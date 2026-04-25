import os
import sys
import logging
import argparse
import multiprocessing

from .io import read_config_file
from .utils import check_parameters
from .filtering import filter_length, filter_groups
from .alignment import run_mafft, run_trimal
from .tables import make_tables
from .gsea import download_gmt, run_gsea
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
    parser.add_argument("code_config",    help="Path to the code configuration file.")
    parser.add_argument("project_config", help="Path to the project configuration file.")
    gmt_group = parser.add_mutually_exclusive_group(required=True)
    gmt_group.add_argument(
        "--organism",
        help="g:Profiler organism name (e.g. 'hsapiens', 'athaliana'). "
             "Automatically downloads the corresponding GMT file. "
             "See https://biit.cs.ut.ee/gprofiler for supported organisms."
    )
    gmt_group.add_argument(
        "--gene-sets",
        help="Path to a local .gmt file. Use this for custom annotations "
             "or organisms not supported by g:Profiler."
    )
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
    if os.path.isdir(outdir):
        sys.exit("Output directory '%s' already exists, please change 'outdir'." % outdir)
    os.makedirs(outdir)

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

    os.makedirs(length_stats_dir)
    os.makedirs(length_filter_dir)

    # --- Pipeline ---
    inputfiles = os.listdir(parameters['inputdir'])

    logger.info("Step 1: Length filtering")
    args_list = [
        (os.path.join(parameters['inputdir'], f), length_stats_dir, length_filter_dir)
        for f in inputfiles
    ]
    with multiprocessing.Pool(processes=parameters['max_process']) as pool:
        pool.starmap(filter_length, args_list)

    logger.info("Step 2: Group filtering (anchor and min_species)")
    filter_groups(
        length_filter_dir=length_filter_dir,
        group_filter_dir=group_filter_dir,
        anchor=parameters['anchor'],
        min_species=parameters['min_species']
    )

    logger.info("Step 3: Alignment with MAFFT")
    run_mafft(
        group_filter_dir=group_filter_dir,
        alignment_dir=alignment_dir,
        mafft_path=parameters['mafft'],
        threads=parameters['max_process']
    )

    logger.info("Step 4: Identity estimation with trimAl")
    run_trimal(
        alignment_dir=alignment_dir,
        ident_dir=ident_dir,
        trimal_path=parameters['trimal']
    )

    logger.info("Step 5: Generating tables")
    make_tables(
        alignment_dir=alignment_dir,
        ident_dir=ident_dir,
        tables_dir=tables_dir,
        anchor=parameters['anchor']
    )

    logger.info("Step 6: Running GSEA prerank")
    gsea_dir = os.path.join(outdir, "gsea")

    if args.organism:
        gmt_path = download_gmt(organism=args.organism, outdir=outdir)
    else:
        gmt_path = args.gene_sets

    run_gsea(
        anchor2mean_path=os.path.join(tables_dir, "anchor2mean.tsv"),
        gsea_dir=gsea_dir,
        gene_sets=gmt_path,
        permutations=args.permutations,
        min_size=args.min_size,
        max_size=args.max_size,
        seed=args.seed,
    )

    logger.info("Enhydra finished successfully.")
