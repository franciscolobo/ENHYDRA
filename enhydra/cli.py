import os
import sys
import logging
import argparse
import multiprocessing

from .io import read_config_file, read_species_list
from .utils import check_parameters
from .filtering import filter_length, filter_groups
from .alignment import run_mafft, run_trimal
from .tables import make_tables
from .gsea import run_gsea
from .orthofinder import preprocess_orthofinder
from .exceptions import EnhydraConfigError, EnhydraIOError, EnhydraToolError


def _setup_logging(outdir: str):
    """Configure logging to both console and a log file in outdir."""
    log_path = os.path.join(outdir, "enhydra.log")
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        handlers=[
            logging.FileHandler(log_path),
            logging.StreamHandler(sys.stdout),
        ]
    )


def _build_arg_parser():
    parser = argparse.ArgumentParser(
        prog="enhydra",
        description="Gene Set Enrichment Analysis for evolutionary genomics."
    )

    # --- Positional arguments ---
    parser.add_argument("code_config",    help="Path to the code configuration file.")
    parser.add_argument("project_config", help="Path to the project configuration file.")

    # --- Input mode ---
    input_group = parser.add_mutually_exclusive_group()
    input_group.add_argument(
        "--orthofinder-dir",
        help="Path to an OrthoFinder 3 output directory. When provided, "
             "ENHYDRA will preprocess the Orthogroup_Sequences/ directory "
             "into inputdir before running the pipeline."
    )

    # --- Run control ---
    parser.add_argument(
        "--resume",
        action="store_true",
        default=False,
        help="Resume a previously interrupted run. Skips steps whose output "
             "directories already exist. The outdir must already exist."
    )

    # --- Paralog handling ---
    parser.add_argument(
        "--paralogs",
        choices=["all", "remove", "longest"],
        default="all",
        help=(
            "How to handle paralogs (multiple sequences per species per group). "
            "'all' keeps all sequences (default). "
            "'remove' discards any group that contains paralogs. "
            "'longest' keeps only the longest sequence per species, "
            "breaking ties at random."
        )
    )

    # --- Two-list differential mode ---
    diff_group = parser.add_argument_group("two-list differential mode")
    diff_group.add_argument(
        "--list1",
        help="Path to a text file listing species IDs for list 1 "
             "(one species ID per line, '#' lines ignored)."
    )
    diff_group.add_argument(
        "--list2",
        help="Path to a text file listing species IDs for list 2 "
             "(one species ID per line, '#' lines ignored)."
    )
    diff_group.add_argument(
        "--anchor1",
        help="Anchor species ID for list 1. Used for group filtering and "
             "table generation. Does not need to be present in list 2."
    )
    diff_group.add_argument(
        "--anchor2",
        help="Anchor species ID for list 2. Used for group filtering and "
             "table generation. Does not need to be present in list 1."
    )
    diff_group.add_argument(
        "--metric",
        choices=["identity", "rank"],
        default="identity",
        help=(
            "Metric used to compute the differential score between the two lists. "
            "'identity' — difference in mean alignment identity (default). "
            "'rank'     — difference in normalised rank position. "
            "Both metrics produce signed scores: positive values indicate "
            "higher conservation in list 1, negative in list 2."
        )
    )
    gmt_group = parser.add_mutually_exclusive_group(required=True)
    gmt_group.add_argument(
        "--organism",
        help=(
            "g:Profiler organism name (e.g. 'hsapiens', 'athaliana'). "
            "Fetches annotations via the g:Profiler API. "
            "See https://biit.cs.ut.ee/gprofiler for supported organisms."
        )
    )
    gmt_group.add_argument(
        "--gene-sets",
        help=(
            "Path to a local .gmt file for GSEApy prerank. "
            "Use for species not supported by g:Profiler or custom annotations."
        )
    )

    # --- g:Profiler options ---
    parser.add_argument(
        "--sources",
        nargs="+",
        default=["GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC"],
        help="g:Profiler data sources to query. Only used with --organism. "
             "Default: GO:BP GO:MF GO:CC KEGG REAC."
    )

    # --- GSEA options ---
    parser.add_argument(
        "--permutations", type=int, default=1000,
        help="Number of GSEA permutations (default: 1000)."
    )
    parser.add_argument(
        "--min-size", type=int, default=5,
        help="Minimum gene set size to test (default: 5)."
    )
    parser.add_argument(
        "--max-size", type=int, default=500,
        help="Maximum gene set size to test (default: 500)."
    )
    parser.add_argument(
        "--seed", type=int, default=42,
        help="Random seed for reproducibility (default: 42)."
    )

    return parser


def main():
    parser = _build_arg_parser()
    args = parser.parse_args()

    # --- Read and validate config ---
    try:
        with open(args.project_config, "r") as fh_project, \
             open(args.code_config, "r") as fh_code:
            parameters = read_config_file(fh_project, fh_code)
    except OSError as e:
        sys.exit("Could not open configuration file: %s" % e)

    try:
        check_parameters(parameters, args.code_config)
    except (EnhydraConfigError, EnhydraIOError, EnhydraToolError) as e:
        sys.exit("Configuration error: %s" % e)

    # --- Set up output directory and logging ---
    outdir = parameters['outdir']
    if os.path.isdir(outdir) and not args.resume:
        sys.exit(
            "Output directory '%s' already exists. Use --resume to continue "
            "a previous run, or change 'outdir' in your project config." % outdir
        )
    os.makedirs(outdir, exist_ok=True)

    _setup_logging(outdir)
    logger = logging.getLogger(__name__)
    logger.info("Welcome to Enhydra")

    # --- Define subdirectories ---
    length_stats_dir  = os.path.join(outdir, "length_stats")
    length_filter_dir = os.path.join(outdir, "length_filter")
    group_filter_dir  = os.path.join(outdir, "group_filter")
    alignment_dir     = os.path.join(outdir, "alignment")
    ident_dir         = os.path.join(outdir, "ident_alignment")
    tables_dir        = os.path.join(outdir, "tables")
    results_dir       = os.path.join(outdir, "enrichment")

    def _should_skip(step_dir: str, step_name: str) -> bool:
        if args.resume and os.path.isdir(step_dir):
            logger.info(
                "Skipping %s (output already exists: %s)", step_name, step_dir
            )
            return True
        return False

    # --- OrthoFinder preprocessing (optional) ---
    if args.orthofinder_dir:
        logger.info("OrthoFinder mode: preprocessing Orthogroup_Sequences/")
        preprocess_orthofinder(
            orthofinder_dir=args.orthofinder_dir,
            inputdir=parameters['inputdir'],
        )

    # --- Pipeline ---
    inputfiles = os.listdir(parameters['inputdir'])

    logger.info("Step 1: Length filtering")
    if not _should_skip(length_filter_dir, "length filtering"):
        os.makedirs(length_stats_dir, exist_ok=True)
        os.makedirs(length_filter_dir, exist_ok=True)
        args_list = [
            (os.path.join(parameters['inputdir'], f), length_stats_dir, length_filter_dir)
            for f in inputfiles
        ]
        with multiprocessing.Pool(processes=parameters['max_process']) as pool:
            pool.starmap(filter_length, args_list)

    logger.info("Step 2: Group filtering (anchor and min_species)")
    if not _should_skip(group_filter_dir, "group filtering"):
        filter_groups(
            length_filter_dir=length_filter_dir,
            group_filter_dir=group_filter_dir,
            anchor=parameters['anchor'],
            min_species=parameters['min_species'],
            paralog_mode=args.paralogs,
        )

    logger.info("Step 3: Alignment with MAFFT")
    if not _should_skip(alignment_dir, "alignment"):
        run_mafft(
            group_filter_dir=group_filter_dir,
            alignment_dir=alignment_dir,
            mafft_path=parameters['mafft'],
            threads=parameters['max_process']
        )

    logger.info("Step 4: Identity estimation with trimAl")
    if not _should_skip(ident_dir, "identity estimation"):
        run_trimal(
            alignment_dir=alignment_dir,
            ident_dir=ident_dir,
            trimal_path=parameters['trimal']
        )

    logger.info("Step 5: Generating tables")
    if not _should_skip(tables_dir, "table generation"):
        make_tables(
            alignment_dir=alignment_dir,
            ident_dir=ident_dir,
            tables_dir=tables_dir,
            anchor=parameters['anchor']
        )

    logger.info("Step 6: Enrichment analysis")
    if not _should_skip(results_dir, "enrichment analysis"):
        run_gsea(
            anchor2mean_path=os.path.join(tables_dir, "anchor2mean.tsv"),
            results_dir=results_dir,
            gene_sets=args.gene_sets,
            organism=args.organism,
            sources=args.sources,
            permutations=args.permutations,
            min_size=args.min_size,
            max_size=args.max_size,
            seed=args.seed,
        )

    logger.info("Enhydra finished successfully.")
