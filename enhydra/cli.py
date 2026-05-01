import os
import sys
import logging
import argparse
import multiprocessing

from .io import read_config_file, read_species_list
from .utils import check_parameters
from .filtering import filter_length, filter_groups, subset_groups
from .alignment import run_mafft, run_trimal
from .tables import make_tables
from .gsea import run_gsea
from .orthofinder import preprocess_orthofinder
from .differential import compute_differential, normalise_scores
from .plotting import make_single_list_plots, make_differential_plots
from .report import build_report
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


def _run_single_list(
    inputdir: str,
    listdir: str,
    anchor: str,
    min_species: int,
    mafft_path: str,
    trimal_path: str,
    max_process: int,
    paralog_mode: str,
    require_anchor: bool,
    resume: bool,
    species: list[str] | None = None,
):
    """Run steps 1–5 (filter, align, identity, tables) for one species list.

    Args:
        inputdir:       Directory of input FASTA files (one per orthogroup).
        listdir:        Output subdirectory for this list's intermediate files.
        anchor:         Anchor species ID for table generation.
        min_species:    Minimum species count for group retention.
        mafft_path:     Path to MAFFT executable.
        trimal_path:    Path to trimAl executable.
        max_process:    Number of parallel processes.
        paralog_mode:   Paralog handling strategy ('all', 'remove', 'longest').
        require_anchor: Whether anchor presence is required in each group.
        resume:         Skip steps whose output directories already exist.
        species:        If provided, subset input groups to these species first.
                        Used in two-list differential mode.
    """
    logger = logging.getLogger(__name__)

    def _should_skip(step_dir: str, step_name: str) -> bool:
        if resume and os.path.isdir(step_dir):
            logger.info("Skipping %s (output already exists: %s)", step_name, step_dir)
            return True
        return False

    subset_dir       = os.path.join(listdir, "subset")
    length_stats_dir = os.path.join(listdir, "length_stats")
    length_filter_dir= os.path.join(listdir, "length_filter")
    group_filter_dir = os.path.join(listdir, "group_filter")
    alignment_dir    = os.path.join(listdir, "alignment")
    ident_dir        = os.path.join(listdir, "ident_alignment")
    tables_dir       = os.path.join(listdir, "tables")

    os.makedirs(listdir, exist_ok=True)

    # Step 0 (two-list mode only): subset groups to species list
    if species is not None:
        if not _should_skip(subset_dir, "subsetting"):
            subset_groups(inputdir, subset_dir, species)
        source_dir = subset_dir
    else:
        source_dir = inputdir

    inputfiles = os.listdir(source_dir)

    logger.info("Step 1: Length filtering")
    if not _should_skip(length_filter_dir, "length filtering"):
        os.makedirs(length_stats_dir, exist_ok=True)
        os.makedirs(length_filter_dir, exist_ok=True)
        args_list = [
            (os.path.join(source_dir, f), length_stats_dir, length_filter_dir)
            for f in inputfiles
        ]
        with multiprocessing.Pool(processes=max_process) as pool:
            pool.starmap(filter_length, args_list)

    logger.info("Step 2: Group filtering")
    if not _should_skip(group_filter_dir, "group filtering"):
        filter_groups(
            length_filter_dir=length_filter_dir,
            group_filter_dir=group_filter_dir,
            anchor=anchor,
            min_species=min_species,
            paralog_mode=paralog_mode,
            require_anchor=require_anchor,
        )

    logger.info("Step 3: Alignment with MAFFT")
    if not _should_skip(alignment_dir, "alignment"):
        run_mafft(
            group_filter_dir=group_filter_dir,
            alignment_dir=alignment_dir,
            mafft_path=mafft_path,
            threads=max_process,
        )

    logger.info("Step 4: Identity estimation with trimAl")
    if not _should_skip(ident_dir, "identity estimation"):
        run_trimal(
            alignment_dir=alignment_dir,
            ident_dir=ident_dir,
            trimal_path=trimal_path,
        )

    logger.info("Step 5: Generating tables")
    if not _should_skip(tables_dir, "table generation"):
        make_tables(
            alignment_dir=alignment_dir,
            ident_dir=ident_dir,
            tables_dir=tables_dir,
            anchor=anchor,
        )

    return tables_dir


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

    parser.add_argument(
        "--obo-cache",
        default=None,
        help="Path to the directory containing the cached go-basic.obo file. "
             "Used to add GO term names to the HTML report. "
             "Same directory as used in build_gmt_interproscan.py."
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
        "--metric",
        choices=["zscore", "identity", "rank"],
        default="zscore",
        help=(
            "Metric for the differential score between the two lists. "
            "'zscore'   — z-score normalised identity difference (default). "
            "            Removes the effect of overall divergence level, "
            "            recommended when the two groups have different "
            "            evolutionary depths. "
            "'identity' — raw mean identity difference (identity_1 - identity_2). "
            "            Only appropriate when both groups are evolutionarily "
            "            comparable. "
            "'rank'     — normalised rank difference (rank_1/N - rank_2/N). "
            "            Most robust to non-normality in identity distributions. "
            "Positive scores indicate higher relative conservation in list 1."
        )
    )

    # --- Gene set source (mutually exclusive, one required) ---
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

    # --- Validate two-list arguments ---
    two_list_mode = args.list1 is not None or args.list2 is not None
    if two_list_mode:
        missing = [f for f in ["--list1", "--list2"]
                   if getattr(args, f.lstrip("-").replace("-", "_")) is None]
        if missing:
            parser.error(
                "Two-list mode requires both --list1 and --list2. "
                "Missing: %s" % ", ".join(missing)
            )

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

    # --- OrthoFinder preprocessing (optional) ---
    if args.orthofinder_dir:
        logger.info("OrthoFinder mode: preprocessing Orthogroup_Sequences/")
        preprocess_orthofinder(
            orthofinder_dir=args.orthofinder_dir,
            inputdir=parameters['inputdir'],
        )

    # --- Common kwargs for _run_single_list ---
    common_kwargs = dict(
        inputdir=parameters['inputdir'],
        min_species=parameters['min_species'],
        mafft_path=parameters['mafft'],
        trimal_path=parameters['trimal'],
        max_process=parameters['max_process'],
        paralog_mode=args.paralogs,
        resume=args.resume,
    )

    # --- Single-list mode ---
    if not two_list_mode:
        logger.info("Running in single-list mode.")
        tables_dir = _run_single_list(
            listdir=os.path.join(outdir),
            anchor=parameters['anchor'],
            require_anchor=True,
            species=None,
            **common_kwargs,
        )
        results_dir = os.path.join(outdir, "enrichment")

    # --- Two-list differential mode ---
    else:
        logger.info("Running in two-list differential mode.")

        species1 = read_species_list(args.list1)
        species2 = read_species_list(args.list2)
        logger.info("List 1: %d species", len(species1))
        logger.info("List 2: %d species", len(species2))
        logger.info("Anchor: %s", parameters['anchor'])

        logger.info("--- Processing list 1 ---")
        tables_dir1 = _run_single_list(
            listdir=os.path.join(outdir, "list1"),
            anchor=parameters['anchor'],
            require_anchor=False,
            species=species1,
            **common_kwargs,
        )

        logger.info("--- Processing list 2 ---")
        tables_dir2 = _run_single_list(
            listdir=os.path.join(outdir, "list2"),
            anchor=parameters['anchor'],
            require_anchor=False,
            species=species2,
            **common_kwargs,
        )

        # --- Differential ranking ---
        diff_dir = os.path.join(outdir, "differential")
        results_dir = os.path.join(diff_dir, "enrichment")

        def _should_skip_diff(step_dir, step_name):
            if args.resume and os.path.isdir(step_dir):
                logger.info("Skipping %s (output already exists).", step_name)
                return True
            return False

        logger.info("--- Computing differential scores ---")
        anchor2mean_path = os.path.join(diff_dir, "anchor2mean.tsv")
        if not _should_skip_diff(diff_dir, "differential ranking"):
            compute_differential(
                tables_dir1=tables_dir1,
                tables_dir2=tables_dir2,
                diff_dir=diff_dir,
                metric=args.metric,
            )

        logger.info("Step 6: Enrichment analysis (differential)")
        if not _should_skip_diff(results_dir, "enrichment analysis"):
            run_gsea(
                anchor2mean_path=anchor2mean_path,
                results_dir=results_dir,
                gene_sets=args.gene_sets,
                organism=args.organism,
                sources=args.sources,
                permutations=args.permutations,
                min_size=args.min_size,
                max_size=args.max_size,
                seed=args.seed,
            )

        logger.info("Generating differential plots")
        make_differential_plots(
            tables_dir1=tables_dir1,
            tables_dir2=tables_dir2,
            diff_dir=diff_dir,
            plots_dir=os.path.join(diff_dir, "plots"),
            metric=args.metric,
        )

        logger.info("Building HTML report")
        obo_path = os.path.join(args.obo_cache, "go-basic.obo") \
            if args.obo_cache else None
        build_report(
            results_dir=results_dir,
            plots_dir=os.path.join(diff_dir, "plots"),
            report_path=os.path.join(diff_dir, "report.html"),
            obo_path=obo_path,
            mode="differential",
            metric=args.metric,
        )

        logger.info("Enhydra finished successfully.")
        return

    # --- GSEA (single-list mode) ---
    def _should_skip_gsea():
        if args.resume and os.path.isdir(results_dir):
            logger.info("Skipping enrichment analysis (output already exists).")
            return True
        return False

    logger.info("Step 6: Enrichment analysis")
    raw_anchor2mean = os.path.join(outdir, "tables", "anchor2mean.tsv")

    # Apply metric normalisation if requested
    if args.metric == "identity":
        gsea_anchor2mean = raw_anchor2mean
    else:
        import pandas as _pd
        df = _pd.read_csv(raw_anchor2mean, sep="\t", header=None,
                          names=["gene_id", "score"])
        df["score"] = _pd.to_numeric(df["score"], errors="coerce")
        df = df.dropna(subset=["score"])
        df["score"] = normalise_scores(df.set_index("gene_id")["score"],
                                       args.metric).values
        gsea_anchor2mean = os.path.join(outdir, "tables",
                                        "anchor2mean_%s.tsv" % args.metric)
        df.to_csv(gsea_anchor2mean, sep="\t", index=False, header=False)
        logger.info("Scores normalised using metric: %s", args.metric)

    if not _should_skip_gsea():
        run_gsea(
            anchor2mean_path=gsea_anchor2mean,
            results_dir=results_dir,
            gene_sets=args.gene_sets,
            organism=args.organism,
            sources=args.sources,
            permutations=args.permutations,
            min_size=args.min_size,
            max_size=args.max_size,
            seed=args.seed,
        )

    logger.info("Generating plots")
    make_single_list_plots(
        anchor2mean_path=os.path.join(outdir, "tables", "anchor2mean.tsv"),
        results_dir=results_dir,
        plots_dir=os.path.join(outdir, "plots"),
    )

    logger.info("Building HTML report")
    obo_path = os.path.join(args.obo_cache, "go-basic.obo") \
        if args.obo_cache else None
    build_report(
        results_dir=results_dir,
        plots_dir=os.path.join(outdir, "plots"),
        report_path=os.path.join(outdir, "report.html"),
        obo_path=obo_path,
        mode="single",
    )

    logger.info("Enhydra finished successfully.")
