// Basic GLnexus command-line interface for use on one compute node

#include <iostream>
#include <exception>
#include <fstream>
#include <getopt.h>
#include <sstream>
#include <cstdlib>
#include <sys/stat.h>
#include "vcf.h"
#include "hfile.h"
#include "service.h"
#include "unifier.h"
#include "BCFKeyValueData.h"
#include "RocksKeyValue.h"
#include "ctpl_stl.h"
#include "spdlog/spdlog.h"
#include "spdlog/sinks/stdout_sinks.h"
#include "cli_utils.h"

// https://gcc.gnu.org/onlinedocs/cpp/Stringizing.html
#define STRINGIFY(x) #x
#define MACRO_TO_STRING(x) STRINGIFY(x)

using namespace std;

auto console = spdlog::stderr_logger_mt("GLnexus");
GLnexus::Status s;
#define H(desc,expr) \
    s = expr; \
    if (s.bad()) { \
        console->error("Failed to {}: {}", desc, s.str()); \
        if (getenv("DX_JOB_ID")) { \
            ofstream joberrorjson("/home/dnanexus/job_error.json"); \
            joberrorjson << "{\"error\": {\"type\": \"AppError\", \"message\": \""; \
            joberrorjson << "Failed to " << desc << ": " << s.str(); \
            joberrorjson << "\"}}"; \
            joberrorjson.close(); \
        } \
        return 1; \
    }

// Perform all the separate GLnexus operations in one go.
// return 0 on success, 1 on failure.
static int all_steps(const vector<string> &vcf_files,
                     const string &bedfilename,
                     const string &dbpath,
                     const string &config_name,
                     bool more_PL, bool squeeze, bool trim_uncalled_alleles,
                     size_t mem_budget, size_t nr_threads,
                     bool debug,
                     bool iter_compare,
                     size_t bucket_size,
                     const string &checkpoint_out_path,
                     const string &checkpoint_in_path,
                     bool checkpoint_only,
                     const string &output_filename) {
    GLnexus::Status s;
    GLnexus::unifier_config unifier_cfg;
    GLnexus::genotyper_config genotyper_cfg;
    string cfg_txt, cfg_crc32c;
    bool is_loaded_from_checkpoint = false;

    if (!checkpoint_in_path.empty()) {
        struct stat info;
        if (stat(dbpath.c_str(), &info) == 0 && (info.st_mode & S_IFDIR)) {
            console->error("Database directory {} already exists. Please remove it before loading from a checkpoint.", dbpath);
            return 1;
        }

        string cmd = "ln -s " + checkpoint_in_path + " " + dbpath;
        console->info("Loading database from checkpoint by creating a symbolic link: {}", cmd);
        if (system(cmd.c_str()) != 0) {
            console->error("Failed to create symbolic link from {} to {}.", checkpoint_in_path, dbpath);
            return 1;
        }
        is_loaded_from_checkpoint = true;
    }

    if (vcf_files.empty() && !is_loaded_from_checkpoint) {
        console->error("No source GVCF files specified");
        return 1;
    }

    H("load unifier/genotyper configuration",
        GLnexus::cli::utils::load_config(console, config_name, unifier_cfg, genotyper_cfg, cfg_txt, cfg_crc32c,
                                         more_PL, squeeze, trim_uncalled_alleles));

    // initilize empty database
    vector<pair<string,size_t> > contigs;
    if (is_loaded_from_checkpoint) {
        H("read the contigs from DB",
          GLnexus::cli::utils::db_get_contigs(console, dbpath, contigs));
    } else {
        H("initialize database", GLnexus::cli::utils::db_init(console, dbpath, vcf_files[0], contigs,
                                                              bucket_size));
        {
            // sanity check, see that we can get the contigs back
            vector<pair<string,size_t> > contigs_dbg;
            H("read the contigs back from DB",
              GLnexus::cli::utils::db_get_contigs(console, dbpath, contigs_dbg));
            if (contigs_dbg != contigs)
                return GLnexus::Status::Invalid("error, contigs read from DB do not match originals");
        }
    }

    if (nr_threads == 0) {
        nr_threads = std::thread::hardware_concurrency();
    }

    // Load the GVCFs into the database
    unique_ptr<GLnexus::KeyValue::DB> db;
    if (!is_loaded_from_checkpoint) {
        // use an empty range filter
        vector<GLnexus::range> ranges;
        H("bulk load into DB",
          GLnexus::cli::utils::db_bulk_load(console, mem_budget, nr_threads, vcf_files, dbpath, ranges, contigs, &db, false));
    } else {
        GLnexus::RocksKeyValue::config cfg;
        cfg.thread_budget = nr_threads;
        H("open database", GLnexus::RocksKeyValue::Open(dbpath, cfg, db));
    }
    assert(db);

    if (iter_compare) {
        H("compare database iteration methods",
          GLnexus::cli::utils::compare_db_itertion_algorithms(console, dbpath, 50));
    }

    // discover alleles
    vector<GLnexus::range> ranges;
    if (bedfilename.empty()) {
        console->warn("Processing full length of {} contigs, as no --bed was provided. Providing a BED file with regions of interest, if applicable, can speed this up.", std::to_string(contigs.size()));
        for (int rid = 0; rid < contigs.size(); ++rid) {
            ranges.push_back(GLnexus::range(rid, 0, contigs[rid].second));
        }
    } else {
        H("parse the bed file", GLnexus::cli::utils::parse_bed_file(console, bedfilename, contigs, ranges));
    }
    GLnexus::discovered_alleles dsals;
    unsigned sample_count = 0;
    auto nr_threads_m2 = nr_threads > 2 ? nr_threads-2 : 1; // reserve threads for DB bg compactions
    H("discover alleles",
      GLnexus::cli::utils::discover_alleles(console, nr_threads_m2, db.get(), ranges, contigs, dsals, sample_count,
                                            unifier_cfg.min_allele_copy_number == 0));
    if (debug) {
        string filename("/tmp/dsals.yml");
        console->info("Writing discovered alleles as YAML to {}", filename);
        H("serialize discovered alleles to a file",
          GLnexus::cli::utils::yaml_write_discovered_alleles_to_file(dsals, contigs, sample_count, filename));
    }

    // partition dsals by contig to reduce peak memory usage in the unifier
    std::vector<GLnexus::discovered_alleles> dsals_by_contig(contigs.size());
    for (auto p = dsals.begin(); p != dsals.end(); dsals.erase(p++)) {
        UNPAIR(*p, al, dai);
        assert(al.pos.rid >= 0 && al.pos.rid < contigs.size());
        dsals_by_contig[al.pos.rid][al] = dai;
    }

    // unify sites (parallel over dsals_by_contig)
    ctpl::thread_pool unify_pool(nr_threads_m2);
    vector<future<GLnexus::Status>> statuses;
    vector<vector<GLnexus::unified_site>> sites_by_contig(contigs.size());
    vector<GLnexus::unifier_stats> stats_by_contig(contigs.size());
    for (size_t i = 0; i < contigs.size(); i++) {
        statuses.push_back(unify_pool.push([&, i](int tid){
            return GLnexus::cli::utils::unify_sites(console, unifier_cfg, contigs, dsals_by_contig[i],
                                                    sample_count, sites_by_contig[i], stats_by_contig[i]);
        }));
    }

    vector<GLnexus::unified_site> sites;
    GLnexus::unifier_stats stats;
    for (size_t i = 0; i < contigs.size(); i++) {
        H("unify sites", statuses[i].get());
        stats += stats_by_contig[i];
        auto& sites_i = sites_by_contig[i];
        sites.insert(sites.end(), make_move_iterator(sites_i.begin()),
                                  make_move_iterator(sites_i.end()));
        sites_i.clear();
    }
    assert(std::is_sorted(sites.begin(), sites.end()));

    console->info("unified to {} sites cleanly with {} ALT alleles. {} ALT alleles were {} and {} were filtered out on quality thresholds.",
                  sites.size(), stats.unified_alleles, stats.lost_alleles,
                  (unifier_cfg.monoallelic_sites_for_lost_alleles ? "additionally included in monoallelic sites" : "lost due to failure to unify"),
                  stats.filtered_alleles);
    if (debug) {
        string filename("/tmp/sites.yml");
        console->info("Writing unified sites as YAML to {}", filename);
        H("write unified sites to file",
          GLnexus::cli::utils::write_unified_sites_to_file(sites, contigs, filename));
    }

    console->info("Finishing database compaction...");
    if (!checkpoint_out_path.empty()) {
        console->info("Creating checkpoint at {}", checkpoint_out_path);
        H("create checkpoint", db->CreateCheckpoint(checkpoint_out_path));
        if (checkpoint_only) {
            console->info("Checkpoint created. Exiting as requested by --checkpoint-only.");
            return 0;
        }
    }
    db.reset();

    // genotype
    genotyper_cfg.output_residuals = debug;
    vector<string> hdr_lines = {
        ("##GLnexusConfigName="+config_name),
        ("##GLnexusConfigCRC32C="+cfg_crc32c),
        ("##GLnexusConfig="+cfg_txt)
    };
    auto DX_JOB_ID = std::getenv("DX_JOB_ID");
    if (DX_JOB_ID) {
        // if running in DNAnexus, record job ID in header
        hdr_lines.push_back(string("##DX_JOB_ID=")+DX_JOB_ID);
    }
    string outfile(output_filename);
    H("genotype",
      GLnexus::cli::utils::genotype(console, mem_budget, nr_threads, dbpath, genotyper_cfg, sites, hdr_lines, outfile));

    return 0;
}


void help(const char* prog) {
    cout << "Usage: " << prog << " [options] /vcf/file/1 .. /vcf/file/N" << endl
         << "Merge and joint-call input gVCF files, emitting multi-sample BCF on standard output." << endl << endl
         << "Options:" << endl
         << "  --dir DIR, -d DIR              scratch directory path (default: ./GLnexus.DB)" << endl
         << "  --config X, -c X               configuration preset name or .yml filename (default: gatk)" << endl
         << "  --bed FILE, -b FILE            three-column BED file with ranges to analyze" << endl << endl
         << "gVCF Input (choose one):" << endl
         << "  --list FILE, -l FILE           text file containing a list of gVCF filenames, one per line" << endl
         << "  (or provide gVCF files directly on the command line)" << endl << endl
         << "Checkpointing:" << endl
         << "  --checkpoint-in PATH           load database from checkpoint at PATH (cannot be used with gVCF input)" << endl
         << "  --checkpoint-out PATH          save database checkpoint to PATH after data loading (requires gVCF input)" << endl
         << "  --checkpoint-only              exit after creating checkpoint (requires --checkpoint-out)" << endl << endl

         << "Other Options:" << endl
         << "  --more-PL, -P                  include PL from reference bands and other cases omitted by default" << endl
         << "  --squeeze, -S                  reduce pVCF size by suppressing detail in cells derived from reference bands" << endl
         << "  --trim-uncalled-alleles, -a    remove alleles with no output GT calls in postprocessing" << endl
         << "  --mem-gbytes X, -m X           memory budget, in gbytes (default: most of system memory)" << endl
         << "  --threads X, -t X              thread budget (default: all hardware threads)" << endl
         << "  --output FILE, -o FILE         BCF output file (default: standard output)" << endl
         << "  --help, -h                     print this help message" << endl
         << endl << "Configuration presets:" << endl;
    cout << GLnexus::cli::utils::describe_config_presets() << endl;
}

// Expected usage:
//    glnexus [vcf files]
//
int main(int argc, char *argv[]) {
    GLnexus::Status s;
    spdlog::set_level(spdlog::level::info);
    spdlog::set_pattern("[%t] %+");
    #ifdef NDEBUG
    #define BUILD_CONFIG "release"
    #else
    #define BUILD_CONFIG "debug"
    #endif
    console->info("glnexus_cli {} {} {}", BUILD_CONFIG, MACRO_TO_STRING(GIT_REVISION), __DATE__);
    GLnexus::cli::utils::detect_jemalloc(console);

    if (argc < 2) {
        help(argv[0]);
        return 1;
    }

    static struct option long_options[] = {
        {"help", no_argument, 0, 'h'},
        {"bed", required_argument, 0, 'b'},
        {"dir", required_argument, 0, 'd'},
        {"config", required_argument, 0, 'c'},
        {"more-PL", no_argument, 0, 'P'},
        {"squeeze", no_argument, 0, 'S'},
        {"trim-uncalled-alleles", no_argument, 0, 'a'},
        {"list", no_argument, 0, 'l'},
        {"mem-gbytes", required_argument, 0, 'm'},
        {"threads", required_argument, 0, 't'},
        {"bucket_size", required_argument, 0, 'x'},
        {"debug", no_argument, 0, 'g'},
        {"iter_compare", no_argument, 0, 'i'},
        {"output", required_argument, 0, 'o'},
        {"checkpoint-out", required_argument, 0, 1},
        {"checkpoint-in", required_argument, 0, 2},
        {"checkpoint-only", no_argument, 0, 3},
        {0, 0, 0, 0}
    };

    int c;
    string dbpath = "GLnexus.DB";
    string config_name = "gatk";
    bool more_PL = false;
    bool squeeze = false;
    bool trim_uncalled_alleles = false;
    bool list_of_files = false;
    bool debug = false;
    bool iter_compare = false;
    string bedfilename;
    size_t mem_budget = 0, nr_threads = 0;
    size_t bucket_size = GLnexus::BCFKeyValueData::default_bucket_size;
    string checkpoint_out_path, checkpoint_in_path;
    bool checkpoint_only = false;
    string output_filename = "-";

    while (-1 != (c = getopt_long(argc, argv, "hPSadil:b:x:m:t:c:o:",
                                  long_options, nullptr))) {
        switch (c) {
            case 1:
                checkpoint_out_path = string(optarg);
                break;

            case 2:
                checkpoint_in_path = string(optarg);
                break;

            case 3:
                checkpoint_only = true;
                break;

            case 'o':
                output_filename = string(optarg);
                break;

            case 'd':
                dbpath = string(optarg);
                break;

            case 'b':
                bedfilename = string(optarg);
                if (bedfilename.size() == 0) {
                    cerr <<  "invalid BED filename" << endl;
                    return 1;
                }
                break;

            case 'l':
                list_of_files = true;
                break;

            case 'c':
                config_name = string(optarg);
                break;

            case 'P':
                more_PL = true;
                break;

            case 'S':
                squeeze = true;
                break;

            case 'a':
                trim_uncalled_alleles = true;
                break;

            case 'g':
                debug = true;
                break;

            case 'h':
            case '?':
                help(argv[0]);
                exit(0);
                break;

            case 'i':
                iter_compare = true;
                break;

            case 'x':
                bucket_size = strtoul(optarg, nullptr, 10);
                if (bucket_size == 0 || bucket_size > 1000000000) {
                    cerr << "bucket size should be in (1,1e9]" << endl;
                    return 1;
                }
                break;

            case 'm':
                mem_budget = strtoull(optarg, nullptr, 10);
                if (mem_budget == 0 || mem_budget > 16*1024) {
                    cerr << "invalid --mem-gbytes" << endl;
                    return 1;
                }
                mem_budget <<= 30;
                break;

            case 't':
                nr_threads = strtoull(optarg, nullptr, 10);
                if (nr_threads == 0 || nr_threads > 1024) {
                    cerr << "invalid --threads" << endl;
                    return 1;
                }
                break;

            default:
                abort ();
        }
    }

    if (optind > argc-1 && checkpoint_in_path.empty()) {
        help(argv[0]);
        return 1;
    }

    bool has_input_files = list_of_files || (optind < argc);

    if (!checkpoint_in_path.empty() && has_input_files) {
        console->error("Cannot provide input files when loading from a checkpoint with --checkpoint-in.");
        return 1;
    }

    if (checkpoint_in_path.empty() && !has_input_files) {
        console->error("Either provide input files or use --checkpoint-in to load from a checkpoint.");
        return 1;
    }

    if (!checkpoint_out_path.empty() && !has_input_files) {
        console->error("--checkpoint-out can only be used when providing input files.");
        return 1;
    }

    if (checkpoint_only && checkpoint_out_path.empty()) {
        console->error("--checkpoint-only can only be used with --checkpoint-out.");
        return 1;
    }

    vector<string> vcf_files, vcf_files_precursor;
    for (int i=optind; i < argc; i++) {
        vcf_files_precursor.push_back(string(argv[i]));
    }

    if (list_of_files) {
        for (const string& fn : vcf_files_precursor) {
            string gvcf;
            ifstream infile(fn);
            while (getline(infile, gvcf)) {
                vcf_files.push_back(gvcf);
            }
            if (infile.bad() || !infile.eof()) {
                H("read input file list", GLnexus::Status::IOError("reading", fn));
            }
        }
    } else {
        vcf_files = vcf_files_precursor;
    }

    return all_steps(vcf_files, bedfilename, dbpath, config_name, more_PL, squeeze, trim_uncalled_alleles,
                     mem_budget, nr_threads, debug, iter_compare, bucket_size, checkpoint_out_path, checkpoint_in_path, checkpoint_only, output_filename);
}
