import os
import sys
import logging
import argparse
import multiprocessing

from tqdm import tqdm

from .io import read_config_file, read_species_list, parse_obo_names
from .utils import check_parameters
from .filtering import filter_length, filter_groups, subset_groups
from .alignment import run_aligner, run_trimal
from .tables import make_tables, tables_complete
from .gsea import run_gsea
from .orthofinder import preprocess_orthofinder
from .differential import compute_differential, normalise_scores
from .plotting import make_single_list_plots, make_differential_plots
from .report import build_report
from .exceptions import EnhydraConfigError, EnhydraIOError, EnhydraToolError


def _setup_logging(outdir: str, quiet: bool = False):
    """Configure logging to both console and a log file in outdir.

    In quiet mode the console handler is raised to WARNING so that only
    warnings and errors appear on stdout; INFO messages are still written
    to enhydra.log regardless.
    """
    log_path = os.path.join(outdir, "enhydra.log")
    fmt = "%(asctime)s [%(levelname)s] %(name)s: %(message)s"
    file_handler = logging.FileHandler(log_path)
    file_handler.setLevel(logging.INFO)
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.ERROR if quiet else logging.INFO)
    logging.basicConfig(level=logging.INFO, format=fmt,
                        handlers=[file_handler, console_handler])


def _resolve(cli_val, config_val, default=None):
    """Return CLI value if explicitly set, otherwise config value, then default."""
    if cli_val is not None:
        return cli_val
    if config_val not in (None, ''):
        return config_val
    return default


def _filter_length_star(args):
    """Unpack args tuple and call filter_length.

    Required at module level so multiprocessing can pickle it when
    pool.imap_unordered is used in place of pool.starmap.
    """
    return filter_length(*args)


def _run_single_list(
    inputdir: str,
    listdir: str,
    anchor: str,
    min_species: int,
    min_sequences: int,
    mafft_path: str,
    trimal_path: str,
    max_process: int,
    paralog_mode: str,
    require_anchor: bool,
    resume: bool,
    sd_multiplier: float = 2.0,
    aligner: str = "mafft",
    mafft_mode: str = "auto",
    parameters: dict = None,
    species: list[str] | None = None,
    show_progress: bool = False,
    label: str = "",
):
    """Run steps 1–5 (filter, align, identity, tables) for one species list."""
    logger = logging.getLogger(__name__)

    def _should_skip(step_dir: str, step_name: str) -> bool:
        if resume and os.path.isdir(step_dir):
            logger.info("Skipping %s (output already exists: %s)", step_name, step_dir)
            return True
        return False

    def _desc(step: str) -> str:
        return ("%s: %s" % (label, step)) if label else step

    subset_dir        = os.path.join(listdir, "subset")
    length_stats_dir  = os.path.join(listdir, "length_stats")
    length_filter_dir = os.path.join(listdir, "length_filter")
    group_filter_dir  = os.path.join(listdir, "group_filter")
    alignment_dir     = os.path.join(listdir, "alignment")
    ident_dir         = os.path.join(listdir, "ident_alignment")
    tables_dir        = os.path.join(listdir, "tables")

    os.makedirs(listdir, exist_ok=True)

    n_steps = 6 if species is not None else 5

    with tqdm(total=n_steps, desc=_desc("starting"),
              unit="step", disable=not show_progress, leave=True) as sbar:

        # Step 0 (two-list only): subset groups to species list
        if species is not None:
            sbar.set_description(_desc("subsetting"))
            if not _should_skip(subset_dir, "subsetting"):
                subset_groups(inputdir, subset_dir, species,
                              show_progress=show_progress)
            source_dir = subset_dir
            sbar.update(1)
        else:
            source_dir = inputdir

        # Step 1: length filtering
        sbar.set_description(_desc("length filter"))
        logger.info("Step 1: Length filtering")
        if not _should_skip(length_filter_dir, "length filtering"):
            os.makedirs(length_stats_dir, exist_ok=True)
            os.makedirs(length_filter_dir, exist_ok=True)
            files = os.listdir(source_dir)
            args_list = [
                (os.path.join(source_dir, f), length_stats_dir,
                 length_filter_dir, sd_multiplier)
                for f in files
            ]
            with multiprocessing.Pool(processes=max_process) as pool:
                try:
                    list(tqdm(
                        pool.imap_unordered(_filter_length_star, args_list),
                        total=len(args_list),
                        desc="  groups",
                        unit="group",
                        leave=False,
                        disable=not show_progress,
                    ))
                except Exception as e:
                    pool.terminate()
                    raise EnhydraToolError(
                        "Length filtering failed: %s" % e
                    ) from e
        sbar.update(1)

        # Step 2: group filtering
        sbar.set_description(_desc("group filter"))
        logger.info("Step 2: Group filtering")
        if not _should_skip(group_filter_dir, "group filtering"):
            filter_groups(
                length_filter_dir=length_filter_dir,
                group_filter_dir=group_filter_dir,
                anchor=anchor,
                min_species=min_species,
                min_sequences=min_sequences,
                paralog_mode=paralog_mode,
                require_anchor=require_anchor,
                show_progress=show_progress,
            )
        sbar.update(1)

        # Step 3: alignment
        sbar.set_description(_desc("alignment"))
        logger.info("Step 3: Alignment with %s", aligner.upper())
        if not _should_skip(alignment_dir, "alignment"):
            run_aligner(
                group_filter_dir=group_filter_dir,
                alignment_dir=alignment_dir,
                aligner=aligner,
                parameters=parameters,
                show_progress=show_progress,
            )
        sbar.update(1)

        # Step 4: identity estimation
        sbar.set_description(_desc("identity"))
        logger.info("Step 4: Identity estimation with trimAl")
        if not _should_skip(ident_dir, "identity estimation"):
            run_trimal(
                alignment_dir=alignment_dir,
                ident_dir=ident_dir,
                trimal_path=trimal_path,
                show_progress=show_progress,
            )
        sbar.update(1)

        # Step 5: tables
        sbar.set_description(_desc("tables"))
        logger.info("Step 5: Generating tables")
        if not (resume and tables_complete(tables_dir)):
            make_tables(
                alignment_dir=alignment_dir,
                ident_dir=ident_dir,
                tables_dir=tables_dir,
                anchor=anchor,
                show_progress=show_progress,
            )
        sbar.update(1)
        sbar.set_description(_desc("done"))

    stats = {
        "n_input":         len(os.listdir(source_dir)),
        "n_length_filter": len(os.listdir(length_filter_dir)) if os.path.isdir(length_filter_dir) else 0,
        "n_group_filter":  len(os.listdir(group_filter_dir))  if os.path.isdir(group_filter_dir)  else 0,
    }
    return tables_dir, stats


def _build_arg_parser():
    parser = argparse.ArgumentParser(
        prog="enhydra",
        description="Gene Set Enrichment Analysis for evolutionary genomics."
    )

    parser.add_argument("code_config",    help="Path to the code configuration file.")
    parser.add_argument("project_config", help="Path to the project configuration file.")

    input_group = parser.add_mutually_exclusive_group()
    input_group.add_argument(
        "--orthofinder-dir",
        help="Path to an OrthoFinder 3 output directory."
    )

    parser.add_argument(
        "--resume", action="store_true", default=False,
        help="Resume a previously interrupted run."
    )
    parser.add_argument(
        "--quiet", action="store_true", default=False,
        help="Suppress INFO messages on the console and show progress bars instead. "
             "All messages are still written to enhydra.log."
    )
    parser.add_argument(
        "--paralogs", choices=["all", "remove", "longest"], default=None,
        help="Paralog handling strategy. Overrides config 'paralogs'."
    )
    parser.add_argument(
        "--min-species", type=int, default=None,
        help="Minimum species per group. Overrides config 'min_species'."
    )

    diff_group = parser.add_argument_group("two-list differential mode")
    diff_group.add_argument(
        "--list1", default=None,
        help="Path to species list file for list 1. Overrides config 'list1'."
    )
    diff_group.add_argument(
        "--list2", default=None,
        help="Path to species list file for list 2. Overrides config 'list2'."
    )
    diff_group.add_argument(
        "--metric", choices=["zscore", "identity", "rank"], default=None,
        help="Differential/ranking metric. Overrides config 'metric'."
    )

    gmt_group = parser.add_mutually_exclusive_group()
    gmt_group.add_argument(
        "--organism", default=None,
        help="g:Profiler organism name. Overrides config 'organism'."
    )
    gmt_group.add_argument(
        "--gene-sets", default=None,
        help="Path to a local .gmt file. Overrides config 'gene_sets'."
    )

    parser.add_argument("--sources", nargs="+", default=None)
    parser.add_argument("--permutations", type=int, default=None)
    parser.add_argument("--min-size", type=int, default=None)
    parser.add_argument("--max-size", type=int, default=None)
    parser.add_argument("--seed", type=int, default=None)
    parser.add_argument("--fdr-threshold", type=float, default=None)
    parser.add_argument("--obo-cache", default=None)

    return parser

def _log_summary(
    stats: dict,
    results_dir: str,
    fdr_threshold: float,
    label: str = "",
):
    """Log a concise end-of-run summary to INFO."""
    logger = logging.getLogger(__name__)
    prefix = ("[%s] " % label) if label else ""
    n_tested = n_sig = 0
    gsea_csv = os.path.join(results_dir, "gseapy.gene_set.prerank.report.csv")
    if os.path.isfile(gsea_csv):
        import pandas as _pd
        df = _pd.read_csv(gsea_csv)
        n_tested = len(df)
        n_sig = int((df["FDR q-val"] < fdr_threshold).sum())
    logger.info("")
    logger.info("=" * 52)
    logger.info("%sENHYDRA RUN SUMMARY", prefix)
    logger.info("=" * 52)
    logger.info("  Input groups:             %d", stats["n_input"])
    logger.info("  After length filter:      %d", stats["n_length_filter"])
    logger.info("  After group filter:       %d", stats["n_group_filter"])
    logger.info("  Gene sets tested:         %d", n_tested)
    logger.info("  Significant (FDR < %.2f): %d", fdr_threshold, n_sig)
    logger.info("=" * 52)


def main():
    parser = _build_arg_parser()
    args = parser.parse_args()

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

    min_species    = _resolve(args.min_species,   parameters['min_species'],   4)
    min_sequences  = parameters['min_sequences']
    paralogs       = _resolve(args.paralogs,      parameters['paralogs'],      'all')
    metric         = _resolve(args.metric,        parameters['metric'],        'zscore')
    gene_sets      = _resolve(args.gene_sets,     parameters['gene_sets'],     None)
    organism       = _resolve(args.organism,      parameters['organism'],      None)
    permutations   = _resolve(args.permutations,  parameters['permutations'],  1000)
    min_size       = _resolve(args.min_size,      parameters['min_size'],      5)
    max_size       = _resolve(args.max_size,      parameters['max_size'],      500)
    seed           = _resolve(args.seed,          parameters['seed'],          42)
    fdr_threshold  = _resolve(args.fdr_threshold, parameters['fdr_threshold'], 0.25)
    list1_path     = _resolve(args.list1,         parameters['list1'],         None)
    list2_path     = _resolve(args.list2,         parameters['list2'],         None)
    sources_raw    = _resolve(args.sources,       parameters['sources'],
                              'GO:BP GO:MF GO:CC KEGG REAC')
    sources        = sources_raw if isinstance(sources_raw, list) \
                     else sources_raw.split()
    sd_multiplier  = parameters['length_filter_sd']
    aligner        = parameters['aligner']
    mafft_mode     = parameters['mafft_mode']
    top_n          = parameters['top_n']
    obo_cache      = _resolve(args.obo_cache, parameters['obo_cache'], None)

    if not gene_sets and not organism:
        parser.error(
            "A gene set source is required. Set 'gene_sets' or 'organism' "
            "in your project config, or use --gene-sets / --organism."
        )

    two_list_mode = bool(list1_path or list2_path)
    if two_list_mode and not (list1_path and list2_path):
        parser.error(
            "Two-list mode requires both list1 and list2 to be specified "
            "(in config or via --list1/--list2)."
        )

    outdir = parameters['outdir']
    if os.path.isdir(outdir) and not args.resume:
        sys.exit(
            "Output directory '%s' already exists. Use --resume to continue "
            "a previous run, or change 'outdir' in your project config." % outdir
        )
    os.makedirs(outdir, exist_ok=True)

    _setup_logging(outdir, quiet=args.quiet)
    logger = logging.getLogger(__name__)
    logger.info("Welcome to Enhydra")
    logger.info(
        "Resolved parameters: metric=%s, paralogs=%s, min_species=%d, "
        "permutations=%d, min_size=%d, max_size=%d, fdr_threshold=%.2f",
        metric, paralogs, min_species, permutations,
        min_size, max_size, fdr_threshold,
    )

    if args.orthofinder_dir:
        logger.info("OrthoFinder mode: preprocessing Orthogroup_Sequences/")
        preprocess_orthofinder(
            orthofinder_dir=args.orthofinder_dir,
            inputdir=parameters['inputdir'],
        )

    common_kwargs = dict(
        inputdir=parameters['inputdir'],
        min_species=min_species,
        min_sequences=min_sequences,
        mafft_path=parameters['mafft'],
        trimal_path=parameters['trimal'],
        max_process=parameters['max_process'],
        paralog_mode=paralogs,
        sd_multiplier=sd_multiplier,
        aligner=aligner,
        mafft_mode=mafft_mode,
        parameters=parameters,
        resume=args.resume,
        show_progress=args.quiet,
    )

    obo_path = os.path.join(obo_cache, "go-basic.obo") if obo_cache else None
    obo_names = parse_obo_names(obo_path) \
        if obo_path and os.path.isfile(obo_path) else {}

    # --- Single-list mode ---
    if not two_list_mode:
        logger.info("Running in single-list mode.")
        tables_dir, stats = _run_single_list(
            listdir=outdir,
            anchor=parameters['anchor'],
            require_anchor=True,
            species=None,
            label="",
            **common_kwargs,
        )
        results_dir = os.path.join(outdir, "enrichment")

        def _should_skip_gsea():
            if args.resume and os.path.isdir(results_dir):
                logger.info("Skipping enrichment analysis (output already exists).")
                return True
            return False

        raw_anchor2mean = os.path.join(outdir, "tables", "anchor2mean.tsv")
        if metric == "identity":
            gsea_anchor2mean = raw_anchor2mean
        else:
            import pandas as _pd
            df = _pd.read_csv(raw_anchor2mean, sep="\t", header=None,
                              names=["gene_id", "score"])
            df["score"] = _pd.to_numeric(df["score"], errors="coerce")
            df = df.dropna(subset=["score"])
            df["score"] = normalise_scores(
                df.set_index("gene_id")["score"], metric
            ).values
            gsea_anchor2mean = os.path.join(
                outdir, "tables", "anchor2mean_%s.tsv" % metric
            )
            df.to_csv(gsea_anchor2mean, sep="\t", index=False, header=False)
            logger.info("Scores normalised using metric: %s", metric)

        logger.info("Step 6: Enrichment analysis")
        if not _should_skip_gsea():
            run_gsea(
                anchor2mean_path=gsea_anchor2mean,
                results_dir=results_dir,
                gene_sets=gene_sets,
                organism=organism,
                sources=sources,
                permutations=permutations,
                min_size=min_size,
                max_size=max_size,
                seed=seed,
                fdr_threshold=fdr_threshold,
            )

        logger.info("Generating plots")
        make_single_list_plots(
            anchor2mean_path=raw_anchor2mean,
            results_dir=results_dir,
            plots_dir=os.path.join(outdir, "plots"),
            obo_names=obo_names,
            fdr_threshold=fdr_threshold,
            top_n=top_n,
        )

        logger.info("Building HTML report")
        build_report(
            results_dir=results_dir,
            plots_dir=os.path.join(outdir, "plots"),
            report_path=os.path.join(outdir, "report.html"),
            obo_path=obo_path,
            mode="single",
            fdr_threshold=fdr_threshold,
        )

    # --- Two-list differential mode ---
    else:
        logger.info("Running in two-list differential mode.")
        species1 = read_species_list(list1_path)
        species2 = read_species_list(list2_path)
        logger.info("List 1: %d species. List 2: %d species. Anchor: %s",
                    len(species1), len(species2), parameters['anchor'])

        logger.info("--- Processing list 1 ---")
        tables_dir1, stats1 = _run_single_list(
            listdir=os.path.join(outdir, "list1"),
            anchor=parameters['anchor'],
            require_anchor=False,
            species=species1,
            label="list 1",
            **common_kwargs,
        )

        logger.info("--- Processing list 2 ---")
        tables_dir2, stats2 = _run_single_list(
            listdir=os.path.join(outdir, "list2"),
            anchor=parameters['anchor'],
            require_anchor=False,
            species=species2,
            label="list 2",
            **common_kwargs,
        )
        _log_summary(stats1, results_dir, fdr_threshold, label="list 1")
        _log_summary(stats2, results_dir, fdr_threshold, label="list 2")

        diff_dir    = os.path.join(outdir, "differential")
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
                metric=metric,
            )

        logger.info("Step 6: Enrichment analysis (differential)")
        if not _should_skip_diff(results_dir, "enrichment analysis"):
            run_gsea(
                anchor2mean_path=anchor2mean_path,
                results_dir=results_dir,
                gene_sets=gene_sets,
                organism=organism,
                sources=sources,
                permutations=permutations,
                min_size=min_size,
                max_size=max_size,
                seed=seed,
                fdr_threshold=fdr_threshold,
            )

        logger.info("Generating differential plots")
        make_differential_plots(
            tables_dir1=tables_dir1,
            tables_dir2=tables_dir2,
            diff_dir=diff_dir,
            plots_dir=os.path.join(diff_dir, "plots"),
            metric=metric,
            obo_names=obo_names,
            fdr_threshold=fdr_threshold,
            top_n=top_n,
        )

        logger.info("Building HTML report")
        build_report(
            results_dir=results_dir,
            plots_dir=os.path.join(diff_dir, "plots"),
            report_path=os.path.join(diff_dir, "report.html"),
            obo_path=obo_path,
            mode="differential",
            metric=metric,
            fdr_threshold=fdr_threshold,
        )

    logger.info("Enhydra finished successfully.")
