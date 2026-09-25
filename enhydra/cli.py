import os
import sys
import logging
import argparse
import multiprocessing

from tqdm import tqdm

from .io import read_config_file, read_species_list, parse_obo_names
from .utils import check_parameters, check_lists, resolve_trim_args, \
    resolve_divergence_filter_sd
from .filtering import filter_length, filter_groups, subset_groups, \
    strip_species_from_alignments, aggregate_length_filter_stats, \
    filter_divergent_sequences
from .alignment import run_aligner, run_trimal, run_trimal_columns
from .tables import make_tables
from .gsea import run_gsea
from .orthofinder import preprocess_orthofinder
from .differential import compute_differential, normalise_scores
from .plotting import make_single_list_plots, make_differential_plots
from .report import build_report, build_multi_metric_report
from .stats import aggregate_pipeline_stats, compute_differential_stats
from .msa_viewer import build_alignment_pages, load_group_anchor
from .exceptions import EnhydraConfigError, EnhydraIOError, EnhydraToolError

ALL_METRICS = ("identity", "zscore", "rank")


# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

def _setup_logging(outdir: str, quiet: bool = False):
    log_path        = os.path.join(outdir, "enhydra.log")
    fmt             = "%(asctime)s [%(levelname)s] %(name)s: %(message)s"
    file_handler    = logging.FileHandler(log_path)
    file_handler.setLevel(logging.INFO)
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.ERROR if quiet else logging.INFO)
    logging.basicConfig(level=logging.INFO, format=fmt,
                        handlers=[file_handler, console_handler])


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------

def _resolve(cli_val, config_val, default=None):
    if cli_val is not None:
        return cli_val
    if config_val not in (None, ''):
        return config_val
    return default


def _step_complete(step_dir: str, sentinel_files: list[str] | None = None) -> bool:
    if not os.path.isdir(step_dir):
        return False
    if sentinel_files:
        return all(
            os.path.isfile(os.path.join(step_dir, f)) and
            os.path.getsize(os.path.join(step_dir, f)) > 0
            for f in sentinel_files
        )
    return len(os.listdir(step_dir)) > 0


def _filter_length_star(args):
    return filter_length(*args)


def _reconstruct_alignment_pages(pages_dir: str) -> dict[str, str]:
    """Rebuild a {group_id: page_path} map from an existing alignments/ dir.

    Used when a --resume/--replot run skips regenerating alignment pages
    because pages_dir already has content from a prior run. Downstream
    steps (report building, in a later commit) still need this mapping
    populated regardless of which branch ran, so it is reconstructed from
    whatever HTML files are already on disk rather than left undefined —
    the same "must stay populated either way" concern already handled for
    tables_dir/stats elsewhere in this module.
    """
    if not os.path.isdir(pages_dir):
        return {}
    return {
        os.path.splitext(f)[0]: os.path.abspath(os.path.join(pages_dir, f))
        for f in os.listdir(pages_dir)
        if f.endswith(".html")
    }


def _normalise_anchor2mean(raw_path: str, metric: str, tables_dir: str) -> str:
    """Normalise anchor2mean.tsv scores for a given metric; return scored path."""
    logger = logging.getLogger(__name__)
    if metric == "identity":
        return raw_path

    import pandas as _pd
    df = _pd.read_csv(raw_path, sep="\t", header=None, names=["gene_id", "score"])
    df["score"] = _pd.to_numeric(df["score"], errors="coerce")
    df = df.dropna(subset=["score"])
    series = df.set_index("gene_id")["score"]

    if metric == "zscore":
        normed = normalise_scores(series, "zscore")
    else:
        n      = len(series)
        normed = series.rank(ascending=True) / n

    df = df.set_index("gene_id")
    df["score"] = normed
    df = df.reset_index()
    out_path = os.path.join(tables_dir, "anchor2mean_%s.tsv" % metric)
    df.to_csv(out_path, sep="\t", index=False, header=False)
    logger.info("Normalised anchor2mean [%s] → %s", metric, out_path)
    return out_path


def _log_summary(
    stats: dict,
    results_dir: str,
    fdr_threshold: float,
    label: str = "",
):
    logger = logging.getLogger(__name__)
    prefix = ("[%s] " % label) if label else ""
    n_tested = n_sig = 0
    gsea_csv = os.path.join(results_dir, "gseapy.gene_set.prerank.report.csv")
    if os.path.isfile(gsea_csv):
        import pandas as _pd
        df       = _pd.read_csv(gsea_csv)
        n_tested = len(df)
        n_sig    = int((df["FDR q-val"] < fdr_threshold).sum())
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


# ---------------------------------------------------------------------------
# Core pipeline: filter → align → (divergence filter) → (trim) → identity → tables
# ---------------------------------------------------------------------------

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
    exclude_from_identity: set[str] | None = None,
    trim_args: list[str] | None = None,
    divergence_filter_sd: float | None = None,
    require_anchor_in_tables: bool | None = None,
) -> tuple[str, dict]:
    logger = logging.getLogger(__name__)

    def _desc(step):
        return ("%s: %s" % (label, step)) if label else step

    def _skip(step_dir, step_name, sentinel_files=None):
        if resume and _step_complete(step_dir, sentinel_files):
            logger.info("Skipping %s (output already exists: %s)", step_name, step_dir)
            return True
        return False

    subset_dir           = os.path.join(listdir, "subset")
    length_stats_dir      = os.path.join(listdir, "length_stats")
    length_filter_dir     = os.path.join(listdir, "length_filter")
    group_filter_dir      = os.path.join(listdir, "group_filter")
    group_stats_dir       = os.path.join(listdir, "group_filter_stats")
    alignment_dir         = os.path.join(listdir, "alignment")
    stripped_dir          = os.path.join(listdir, "alignment_stripped")
    # Divergent-sequence filter (see filtering.filter_divergent_sequences()):
    # a scratch trimAl -sident pass, the filter's own pruned/passthrough
    # output, the de-gapped survivors awaiting realignment, and this step's
    # own drop-reason log. Positioned after anchor-stripping (so an
    # injected anchor is never itself a candidate for removal here) and
    # before column trimming (so trimming subsequently operates on the
    # corrected, post-pruning alignment rather than one still distorted by
    # a divergent sequence's spurious gap placement).
    divergence_sident_dir   = os.path.join(listdir, "divergence_sident")
    divergence_filtered_dir = os.path.join(listdir, "alignment_divergence_filtered")
    divergence_realign_dir  = os.path.join(listdir, "divergence_realign_input")
    divergence_stats_dir    = os.path.join(listdir, "divergence_filter_stats")
    trimmed_dir           = os.path.join(listdir, "alignment_trimmed")
    # Sidecar directory for trimAl's -colnumbering output, capturing which
    # original alignment columns survive trimming. Populated only when
    # trim_args is set (see the trim_args block below). cli.py's main()
    # recomputes this exact same path independently — using the same
    # listdir-based naming convention — when building alignment pages, so
    # the two must be kept in sync if this naming ever changes.
    trim_colnumbering_dir = os.path.join(listdir, "alignment_trimmed_colnumbering")
    ident_dir              = os.path.join(listdir, "ident_alignment")
    tables_dir             = os.path.join(listdir, "tables")

    # These two flags answer genuinely different questions and must not be
    # conflated: 'require_anchor' controls whether filter_groups() *drops*
    # groups lacking the anchor sequence entirely (correctly False for both
    # lists in two-list mode, since most groups won't contain the anchor's
    # own species). 'require_anchor_in_tables' controls whether make_tables()
    # *flags* a missing anchor as an anomaly worth recording in
    # tables/drop_reasons.tsv. For list1 specifically, a group reaching the
    # tables step without an anchor sequence IS worth flagging even though
    # it wasn't dropped upstream — that group will later be silently
    # excluded from anchor2mean.tsv/group2anchor.tsv, and list1's
    # group2anchor.tsv is exactly what compute_differential() reads to map
    # groups to anchor gene IDs. For list2, the same absence is fully
    # expected (list2 never contains the anchor's own species) and should
    # stay silent. Defaulting to require_anchor when not given preserves
    # single-list mode's existing behaviour unchanged.
    #
    # Note this is independent of the divergence filter potentially
    # removing the anchor sequence from the *identity computation*: this
    # step's make_tables() call always reads gene-ID mapping from the
    # original alignment_dir (see below), which the divergence filter
    # never modifies in place, so a group whose anchor happened to be its
    # divergent outlier still resolves to the correct anchor gene ID here
    # — only the identity score itself reflects the outlier's removal.
    if require_anchor_in_tables is None:
        require_anchor_in_tables = require_anchor

    os.makedirs(listdir, exist_ok=True)
    n_steps = (5 + (species is not None) + bool(exclude_from_identity)
               + bool(divergence_filter_sd) + bool(trim_args))

    with tqdm(total=n_steps, desc=_desc("starting"),
              unit="step", disable=not show_progress, leave=True) as sbar:

        if species is not None:
            sbar.set_description(_desc("subsetting"))
            if not _skip(subset_dir, "subsetting"):
                subset_groups(inputdir, subset_dir, species,
                              show_progress=show_progress)
            source_dir = subset_dir
            sbar.update(1)
        else:
            source_dir = inputdir

        sbar.set_description(_desc("length filter"))
        logger.info("Step 1: Length filtering")
        if not _skip(length_filter_dir, "length filtering"):
            os.makedirs(length_stats_dir, exist_ok=True)
            os.makedirs(length_filter_dir, exist_ok=True)
            args_list = [
                (os.path.join(source_dir, f), length_stats_dir,
                 length_filter_dir, sd_multiplier)
                for f in os.listdir(source_dir)
            ]
            pool = multiprocessing.Pool(processes=max_process)
            try:
                list(tqdm(
                    pool.imap_unordered(_filter_length_star, args_list),
                    total=len(args_list), desc="  groups", unit="group",
                    leave=False, disable=not show_progress,
                ))
            except Exception as e:
                raise EnhydraToolError("Length filtering failed: %s" % e) from e
            finally:
                pool.terminate()
                pool.join()
        # Runs unconditionally, whether or not the filtering above was just
        # skipped via --resume. This is deliberate: aggregation is a cheap,
        # idempotent read of whatever per-group '_lengthstats' files already
        # exist in length_stats_dir — it is not tied to whether filtering
        # happened in *this* invocation. Nesting it inside the `if not
        # _skip(...)` block would silently produce empty/stale summary
        # files on any run that resumes past an already-completed length
        # filter step, even though the underlying per-group stats files are
        # present and complete on disk.
        aggregate_length_filter_stats(
            length_stats_dir=length_stats_dir,
            length_filter_stats_dir=os.path.join(listdir, "length_filter_stats"),
        )
        sbar.update(1)

        sbar.set_description(_desc("group filter"))
        logger.info("Step 2: Group filtering")
        if not _skip(group_filter_dir, "group filtering"):
            filter_groups(
                length_filter_dir=length_filter_dir,
                group_filter_dir=group_filter_dir,
                group_stats_dir=group_stats_dir,
                anchor=anchor,
                min_species=min_species,
                min_sequences=min_sequences,
                paralog_mode=paralog_mode,
                require_anchor=require_anchor,
                show_progress=show_progress,
            )
        sbar.update(1)

        sbar.set_description(_desc("alignment"))
        logger.info("Step 3: Alignment with %s", aligner.upper())
        if not _skip(alignment_dir, "alignment"):
            run_aligner(
                group_filter_dir=group_filter_dir,
                alignment_dir=alignment_dir,
                aligner=aligner,
                parameters=parameters,
                show_progress=show_progress,
            )
        sbar.update(1)

        trimal_input_dir = alignment_dir

        if exclude_from_identity:
            sbar.set_description(_desc("stripping anchor"))
            logger.info("Step 3b: Stripping injected species from alignments: %s",
                        exclude_from_identity)
            if not _skip(stripped_dir, "stripping anchor from alignments"):
                strip_species_from_alignments(
                    alignment_dir=alignment_dir,
                    stripped_dir=stripped_dir,
                    exclude=exclude_from_identity,
                    show_progress=show_progress,
                )
            trimal_input_dir = stripped_dir
            sbar.update(1)

        if divergence_filter_sd:
            sbar.set_description(_desc("divergence filter"))
            logger.info(
                "Step 3c: Filtering highly divergent sequences "
                "(sd_multiplier=%s)", divergence_filter_sd,
            )
            # Same two-output resume-guard lesson already applied to column
            # trimming below: check every output this step produces
            # (filtered alignments AND the stats log), not just one, or a
            # resumed run whose divergence_filtered_dir predates this
            # feature could skip the block forever and never produce
            # drop_reasons.tsv at all.
            divergence_ready = (
                _step_complete(divergence_filtered_dir)
                and _step_complete(divergence_stats_dir, ["drop_reasons.tsv"])
            )
            if not (resume and divergence_ready):
                os.makedirs(divergence_sident_dir, exist_ok=True)
                run_trimal(
                    alignment_dir=trimal_input_dir,
                    ident_dir=divergence_sident_dir,
                    trimal_path=trimal_path,
                    n_proc=max_process,
                    show_progress=show_progress,
                )
                filter_divergent_sequences(
                    alignment_dir=trimal_input_dir,
                    sident_dir=divergence_sident_dir,
                    filtered_dir=divergence_filtered_dir,
                    realign_input_dir=divergence_realign_dir,
                    stats_dir=divergence_stats_dir,
                    min_species=min_species,
                    min_sequences=min_sequences,
                    sd_multiplier=divergence_filter_sd,
                    show_progress=show_progress,
                )
                n_to_realign = (
                    len(os.listdir(divergence_realign_dir))
                    if os.path.isdir(divergence_realign_dir) else 0
                )
                if n_to_realign:
                    logger.info(
                        "Realigning %d group(s) after divergent sequence "
                        "removal...", n_to_realign,
                    )
                    run_aligner(
                        group_filter_dir=divergence_realign_dir,
                        alignment_dir=divergence_filtered_dir,
                        aligner=aligner,
                        parameters=parameters,
                        show_progress=show_progress,
                    )
            else:
                logger.info(
                    "Skipping divergence filtering (output already exists: "
                    "%s, %s)", divergence_filtered_dir, divergence_stats_dir,
                )
            trimal_input_dir = divergence_filtered_dir
            sbar.update(1)

        if trim_args:
            sbar.set_description(_desc("trimming columns"))
            logger.info(
                "Step 3d: Trimming alignment columns with trimAl (%s)",
                " ".join(trim_args),
            )
            # Colnumbering capture rides along with the existing -out
            # trimming call (trimAl supports both flags together — see
            # msa_viewer module notes), so it is always captured whenever
            # trimming is enabled at all; there is no separate opt-in flag.
            #
            # The skip-guard here deliberately checks *both* trimmed_dir
            # and trim_colnumbering_dir rather than reusing the plain
            # _skip(trimmed_dir, ...) helper: if trimmed_dir already has
            # contents from a run predating this feature, a resumed run
            # would otherwise skip the whole block forever and never
            # produce colnumbering sidecar files at all. This mirrors the
            # aggregate_length_filter_stats() lesson elsewhere in this
            # function — a resume-guard checking only one of two outputs
            # a step produces can silently strand the other output.
            colnumbering_ready = _step_complete(trim_colnumbering_dir)
            if not (resume and _step_complete(trimmed_dir) and colnumbering_ready):
                run_trimal_columns(
                    alignment_dir=trimal_input_dir,
                    trimmed_dir=trimmed_dir,
                    trimal_path=trimal_path,
                    trim_args=trim_args,
                    n_proc=max_process,
                    show_progress=show_progress,
                    colnumbering_dir=trim_colnumbering_dir,
                )
            else:
                logger.info(
                    "Skipping column trimming (output already exists: %s, %s)",
                    trimmed_dir, trim_colnumbering_dir,
                )
            trimal_input_dir = trimmed_dir
            sbar.update(1)

        sbar.set_description(_desc("identity"))
        logger.info("Step 4: Identity estimation with trimAl")
        if not _skip(ident_dir, "identity estimation"):
            run_trimal(
                alignment_dir=trimal_input_dir,
                ident_dir=ident_dir,
                trimal_path=trimal_path,
                n_proc=max_process,
                show_progress=show_progress,
            )
        sbar.update(1)

        sbar.set_description(_desc("tables"))
        logger.info("Step 5: Generating tables")
        if not _skip(tables_dir, "table generation",
                     sentinel_files=["group2mean.tsv", "anchor2mean.tsv",
                                     "group2anchor.tsv"]):
            make_tables(
                alignment_dir=alignment_dir,
                ident_dir=ident_dir,
                tables_dir=tables_dir,
                anchor=anchor,
                require_anchor=require_anchor_in_tables,
                show_progress=show_progress,
            )
        sbar.update(1)
        sbar.set_description(_desc("done"))

    stats = {
        "n_input":         len(os.listdir(source_dir)),
        "n_length_filter": len(os.listdir(length_filter_dir))
                           if os.path.isdir(length_filter_dir) else 0,
        "n_group_filter":  len(os.listdir(group_filter_dir))
                           if os.path.isdir(group_filter_dir) else 0,
    }
    return tables_dir, stats


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _build_arg_parser():
    parser = argparse.ArgumentParser(
        prog="enhydra",
        description="Gene Set Enrichment Analysis for evolutionary genomics.",
    )
    parser.add_argument("code_config",    help="Path to the code configuration file.")
    parser.add_argument("project_config", help="Path to the project configuration file.")

    input_group = parser.add_mutually_exclusive_group()
    input_group.add_argument("--orthofinder-dir",
                             help="Path to an OrthoFinder 3 output directory.")

    parser.add_argument("--resume",  action="store_true", default=False,
                        help="Resume a previously interrupted run.")
    parser.add_argument("--replot",  action="store_true", default=False,
                        help="Re-run GSEA, plots, and the HTML report without "
                             "repeating alignment. Implies --resume.")
    parser.add_argument("--quiet",   action="store_true", default=False,
                        help="Suppress INFO/WARNING on the console; show progress "
                             "bars instead.")
    parser.add_argument("--paralogs", choices=["all", "remove", "longest"], default=None)
    parser.add_argument("--min-species", type=int, default=None)
    parser.add_argument("--trim", default=None,
                        help="Trim alignment columns with trimAl before "
                             "identity estimation. Accepts a number between "
                             "0 and 1 (used as trimAl's -gt gap threshold), "
                             "or one of: strict, strictplus, automated "
                             "(mapped to trimAl's -strict, -strictplus, "
                             "-automated1 respectively). If unset (default), "
                             "no column trimming is performed.")
    parser.add_argument("--divergence-filter-sd", type=float, default=None,
                        help="Remove sequences whose trimAl -sident identity "
                             "to their closest match falls more than this "
                             "many standard deviations below the group's "
                             "mean, then realign the group without them. "
                             "Runs after alignment (and after anchor "
                             "stripping in two-list mode) but before column "
                             "trimming. Groups with fewer than min_species "
                             "sequences are exempt from this check, since "
                             "too few sequences make the mean/SD unreliable. "
                             "If unset (default), no divergence filtering is "
                             "performed.")
    parser.add_argument("--all-metrics", action="store_true", default=False,
                        help="Run GSEA for all three ranking metrics and produce "
                             "a tabbed HTML report.")

    diff_group = parser.add_argument_group("two-list differential mode")
    diff_group.add_argument("--list1",   default=None)
    diff_group.add_argument("--list2",   default=None)
    diff_group.add_argument("--metric",  choices=["zscore", "identity", "rank"],
                            default=None)

    gmt_group = parser.add_mutually_exclusive_group()
    gmt_group.add_argument("--organism",  default=None)
    gmt_group.add_argument("--gene-sets", default=None)

    parser.add_argument("--sources",       nargs="+", default=None)
    parser.add_argument("--permutations",  type=int,   default=None)
    parser.add_argument("--min-size",      type=int,   default=None)
    parser.add_argument("--max-size",      type=int,   default=None)
    parser.add_argument("--seed",          type=int,   default=None)
    parser.add_argument("--fdr-threshold", type=float, default=None)
    parser.add_argument("--top-n",         type=int,   default=None)
    parser.add_argument("--obo-cache",     default=None)

    return parser


def main():
    parser = _build_arg_parser()
    args   = parser.parse_args()

    try:
        with open(args.project_config) as fh_project, \
             open(args.code_config)    as fh_code:
            parameters = read_config_file(fh_project, fh_code)
    except OSError as e:
        sys.exit("Could not open configuration file: %s" % e)

    try:
        check_parameters(parameters, args.code_config)
    except (EnhydraConfigError, EnhydraIOError, EnhydraToolError) as e:
        sys.exit("Configuration error: %s" % e)

    min_species   = _resolve(args.min_species,   parameters['min_species'],   4)
    min_sequences = parameters['min_sequences']
    paralogs      = _resolve(args.paralogs,      parameters['paralogs'],      'all')
    trim          = _resolve(args.trim,          parameters['trim'],          '')
    divergence_filter_sd_raw = _resolve(
        args.divergence_filter_sd, parameters['divergence_filter_sd'], '')
    metric        = _resolve(args.metric,        parameters['metric'],        'zscore')
    gene_sets     = _resolve(args.gene_sets,     parameters['gene_sets'],     None)
    organism      = _resolve(args.organism,      parameters['organism'],      None)
    permutations  = _resolve(args.permutations,  parameters['permutations'],  1000)
    min_size      = _resolve(args.min_size,      parameters['min_size'],      5)
    max_size      = _resolve(args.max_size,      parameters['max_size'],      500)
    seed          = _resolve(args.seed,          parameters['seed'],          42)
    fdr_threshold = _resolve(args.fdr_threshold, parameters['fdr_threshold'], 0.25)
    top_n         = _resolve(args.top_n,         parameters['top_n'],         20)
    list1_path    = _resolve(args.list1,         parameters['list1'],         None)
    list2_path    = _resolve(args.list2,         parameters['list2'],         None)
    list1_name    = parameters['list1_name']
    list2_name    = parameters['list2_name']
    sources_raw   = _resolve(args.sources,       parameters['sources'],
                             'GO:BP GO:MF GO:CC KEGG REAC')
    sources       = sources_raw if isinstance(sources_raw, list) \
                    else sources_raw.split()
    sd_multiplier = parameters['length_filter_sd']
    aligner       = parameters['aligner']
    mafft_mode    = parameters['mafft_mode']
    obo_cache     = _resolve(args.obo_cache, parameters['obo_cache'], None)
    all_metrics   = args.all_metrics
    replot        = args.replot
    resume        = args.resume or replot

    try:
        trim_args = resolve_trim_args(trim)
    except EnhydraConfigError as e:
        sys.exit("Configuration error: %s" % e)

    try:
        divergence_filter_sd = resolve_divergence_filter_sd(divergence_filter_sd_raw)
    except EnhydraConfigError as e:
        sys.exit("Configuration error: %s" % e)

    if not gene_sets and not organism:
        parser.error(
            "A gene set source is required. Set 'gene_sets' or 'organism' in "
            "your project config, or use --gene-sets / --organism."
        )

    two_list_mode = bool(list1_path or list2_path)
    if two_list_mode and not (list1_path and list2_path):
        parser.error("Two-list mode requires both list1 and list2.")

    outdir = parameters['outdir']
    if os.path.isdir(outdir) and not resume:
        sys.exit(
            "Output directory '%s' already exists. Use --resume to continue "
            "a previous run, --replot to re-run GSEA and regenerate plots, "
            "or change 'outdir' in your project config." % outdir
        )
    os.makedirs(outdir, exist_ok=True)

    _setup_logging(outdir, quiet=args.quiet)
    logger = logging.getLogger(__name__)
    logger.info("Welcome to Enhydra")
    logger.info(
        "Resolved parameters: metric=%s, all_metrics=%s, replot=%s, "
        "paralogs=%s, trim=%s, divergence_filter_sd=%s, min_species=%d, "
        "permutations=%d, fdr_threshold=%.2f",
        metric, all_metrics, replot, paralogs, (trim or "none"),
        (divergence_filter_sd if divergence_filter_sd else "none"),
        min_species, permutations, fdr_threshold,
    )

    if args.orthofinder_dir:
        logger.info("OrthoFinder mode: preprocessing Orthogroup_Sequences/")
        preprocess_orthofinder(
            orthofinder_dir=args.orthofinder_dir,
            inputdir=parameters['inputdir'],
        )

    obo_path  = os.path.join(obo_cache, "go-basic.obo") if obo_cache else None
    obo_names = parse_obo_names(obo_path) \
                if obo_path and os.path.isfile(obo_path) else {}

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
        resume=resume,
        show_progress=args.quiet,
        trim_args=trim_args,
        divergence_filter_sd=divergence_filter_sd,
    )

    # In two-list mode, default to all three metrics for a tabbed comparison.
    # An explicit --metric flag overrides this.
    if two_list_mode and not all_metrics and args.metric is None:
        all_metrics = True
        logger.info(
            "Two-list mode: defaulting to --all-metrics for a tabbed report. "
            "Pass an explicit --metric flag to run a single metric instead."
        )

    metrics_to_run = ALL_METRICS if all_metrics else (metric,)

    gsea_kwargs = dict(
        gene_sets=gene_sets, organism=organism, sources=sources,
        permutations=permutations, min_size=min_size, max_size=max_size,
        seed=seed, fdr_threshold=fdr_threshold,
    )

    # ------------------------------------------------------------------ #
    #  Single-list mode                                                    #
    # ------------------------------------------------------------------ #
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
        aggregate_pipeline_stats(outdir, n_input=stats["n_input"])

        # Alignment pages (bucket 3 = every group in group2anchor.tsv).
        # Rendered once here, not per metric — alignment pages don't
        # depend on the ranking metric at all. Path naming mirrors the
        # convention used inside _run_single_list() for trim_colnumbering_dir
        # (listdir/alignment_trimmed_colnumbering); the two must stay in
        # sync since this is recomputed independently rather than threaded
        # through _run_single_list()'s return value.
        # If the divergence filter ran, its output (post-removal, and
        # re-aligned for affected groups) is the alignment that identity
        # was actually computed from — rendering the raw pre-filter
        # 'alignment' dir here would show a since-removed sequence still
        # present, with a mean-identity figure in the metadata panel that
        # no longer matches what's on screen. This mirrors the existing
        # trim_colnumbering_dir naming duplication note just above: this
        # path must stay in sync with _run_single_list()'s own
        # divergence_filtered_dir naming if that ever changes.
        alignment_dir_single   = os.path.join(
            outdir,
            "alignment_divergence_filtered" if divergence_filter_sd else "alignment",
        )
        alignment_pages_dir    = os.path.join(outdir, "alignments")
        colnumbering_dir_single = (
            os.path.join(outdir, "alignment_trimmed_colnumbering")
            if trim_args else None
        )
        divergence_drop_reasons_path_single = (
            os.path.join(outdir, "divergence_filter_stats", "drop_reasons.tsv")
            if divergence_filter_sd else None
        )
        if not (resume and _step_complete(alignment_pages_dir)):
            logger.info("Rendering alignment pages for report...")
            alignment_pages = build_alignment_pages(
                alignment_dir=alignment_dir_single,
                tables_dir=tables_dir,
                outdir=alignment_pages_dir,
                anchor_species=parameters['anchor'],
                colnumbering_dir=colnumbering_dir_single,
                divergence_drop_reasons_path=divergence_drop_reasons_path_single,
                show_progress=args.quiet,
            )
        else:
            logger.info(
                "Skipping alignment page generation (output already exists: %s)",
                alignment_pages_dir,
            )
            alignment_pages = _reconstruct_alignment_pages(alignment_pages_dir)

        raw_anchor2mean = os.path.join(tables_dir, "anchor2mean.tsv")
        metric_outputs  = {}

        for m in metrics_to_run:
            sfx           = ("_%s" % m) if all_metrics else ""
            results_dir_m = os.path.join(outdir, "enrichment%s" % sfx)
            plots_dir_m   = os.path.join(outdir, "plots%s" % sfx)
            gsea_input    = _normalise_anchor2mean(raw_anchor2mean, m, tables_dir)

            logger.info("Step 6 [%s]: Enrichment analysis", m)
            if replot or not _step_complete(results_dir_m,
                                            ["gseapy.gene_set.prerank.report.csv"]):
                run_gsea(anchor2mean_path=gsea_input,
                         results_dir=results_dir_m, **gsea_kwargs)
            else:
                logger.info("Skipping GSEA for metric '%s' (output exists).", m)

            logger.info("Generating plots [%s]", m)
            make_single_list_plots(
                anchor2mean_path=raw_anchor2mean,
                results_dir=results_dir_m,
                plots_dir=plots_dir_m,
                obo_names=obo_names,
                fdr_threshold=fdr_threshold,
                top_n=top_n,
            )
            metric_outputs[m] = {"results_dir": results_dir_m,
                                  "plots_dir":   plots_dir_m}

        logger.info("Building HTML report")
        if all_metrics:
            build_multi_metric_report(
                metric_data=metric_outputs,
                report_path=os.path.join(outdir, "report.html"),
                obo_path=obo_path,
                fdr_threshold=fdr_threshold,
                mode="single",
                gmt_path=gene_sets,
                tables_dir1=tables_dir,
                pipeline_stats_path=os.path.join(outdir, "pipeline_stats.json"),
                alignment_pages1=alignment_pages,
            )
        else:
            build_report(
                results_dir=metric_outputs[metric]["results_dir"],
                plots_dir=metric_outputs[metric]["plots_dir"],
                report_path=os.path.join(outdir, "report.html"),
                obo_path=obo_path,
                mode="single",
                fdr_threshold=fdr_threshold,
                gmt_path=gene_sets,
                tables_dir1=tables_dir,
                pipeline_stats_path=os.path.join(outdir, "pipeline_stats.json"),
                alignment_pages1=alignment_pages,
            )

        for m in metrics_to_run:
            _log_summary(stats, metric_outputs[m]["results_dir"],
                         fdr_threshold, label=m if all_metrics else "")

    # ------------------------------------------------------------------ #
    #  Two-list differential mode                                          #
    # ------------------------------------------------------------------ #
    else:
        logger.info("Running in two-list differential mode.")
        logger.info("List names: '%s' vs '%s'", list1_name, list2_name)
        anchor   = parameters['anchor']
        species1 = read_species_list(list1_path)
        species2 = read_species_list(list2_path)
        check_lists(species1, species2, anchor)

        if anchor not in species1:
            species1        = list(species1) + [anchor]
            anchor_injected = True
        else:
            anchor_injected = False

        logger.info(
            "%s: %d species. %s: %d species. Anchor: %s",
            list1_name, len(species1), list2_name, len(species2), anchor,
        )

        logger.info("--- Processing %s ---", list1_name)
        tables_dir1, stats1 = _run_single_list(
            listdir=os.path.join(outdir, "list1"),
            anchor=anchor, require_anchor=False,
            species=species1, label=list1_name,
            exclude_from_identity={anchor} if anchor_injected else None,
            require_anchor_in_tables=True,
            **common_kwargs,
        )
        aggregate_pipeline_stats(os.path.join(outdir, "list1"),
                                 n_input=stats1["n_input"])

        logger.info("--- Processing %s ---", list2_name)
        tables_dir2, stats2 = _run_single_list(
            listdir=os.path.join(outdir, "list2"),
            anchor=anchor, require_anchor=False,
            species=species2, label=list2_name,
            **common_kwargs,
        )
        aggregate_pipeline_stats(os.path.join(outdir, "list2"),
                                 n_input=stats2["n_input"])
        compute_differential_stats(
            tables_dir1, tables_dir2,
            list1_name=list1_name, list2_name=list2_name,
            output_path=os.path.join(outdir, "differential_stats.json"),
        )

        metric_outputs = {}

        for m in metrics_to_run:
            sfx           = ("_%s" % m) if all_metrics else ""
            diff_dir_m    = os.path.join(outdir, "differential%s" % sfx)
            results_dir_m = os.path.join(diff_dir_m, "enrichment")
            plots_dir_m   = os.path.join(diff_dir_m, "plots")

            logger.info("--- Computing differential scores [metric=%s] ---", m)
            if replot or not _step_complete(diff_dir_m,
                                            ["anchor2mean.tsv",
                                             "differential_scores.tsv"]):
                compute_differential(
                    tables_dir1=tables_dir1,
                    tables_dir2=tables_dir2,
                    diff_dir=diff_dir_m,
                    metric=m,
                )
            else:
                logger.info("Skipping differential ranking for '%s' (output exists).", m)

            logger.info("Step 6 [%s]: Enrichment analysis (differential)", m)
            if replot or not _step_complete(results_dir_m,
                                            ["gseapy.gene_set.prerank.report.csv"]):
                run_gsea(
                    anchor2mean_path=os.path.join(diff_dir_m, "anchor2mean.tsv"),
                    results_dir=results_dir_m,
                    **gsea_kwargs,
                )
            else:
                logger.info("Skipping GSEA for metric '%s' (output exists).", m)

            logger.info("Generating differential plots [%s]", m)
            make_differential_plots(
                tables_dir1=tables_dir1,
                tables_dir2=tables_dir2,
                diff_dir=diff_dir_m,
                plots_dir=plots_dir_m,
                metric=m,
                obo_names=obo_names,
                fdr_threshold=fdr_threshold,
                top_n=top_n,
                name1=list1_name,
                name2=list2_name,
            )
            metric_outputs[m] = {"results_dir": results_dir_m,
                                  "plots_dir":   plots_dir_m}

        # Alignment pages (bucket 3 in two-list mode = every group_id in
        # differential_scores.tsv, i.e. the intersection compute_differential()
        # already computed). Built once here using metrics_to_run[0]'s
        # differential_scores.tsv, not per metric — the group_id set is
        # expected to be identical across metrics (same intersection of
        # tables_dir1/tables_dir2, just a different score column), a
        # simplifying assumption confirmed acceptable rather than
        # re-validated against every metric's own output.
        first_metric   = metrics_to_run[0]
        first_diff_dir = os.path.join(
            outdir, "differential%s" % (("_%s" % first_metric) if all_metrics else ""),
        )
        diff_scores_path = os.path.join(first_diff_dir, "differential_scores.tsv")

        import pandas as _pd
        diff_group_ids = list(
            _pd.read_csv(diff_scores_path, sep="\t")["group_id"].astype(str)
        )

        list1_alignment_dir = os.path.join(
            outdir, "list1",
            "alignment_divergence_filtered" if divergence_filter_sd else "alignment",
        )
        list1_pages_dir      = os.path.join(outdir, "list1", "alignments")
        list1_colnum_dir     = (
            os.path.join(outdir, "list1", "alignment_trimmed_colnumbering")
            if trim_args else None
        )
        list1_divergence_drop_reasons_path = (
            os.path.join(outdir, "list1", "divergence_filter_stats", "drop_reasons.tsv")
            if divergence_filter_sd else None
        )
        list2_alignment_dir = os.path.join(
            outdir, "list2",
            "alignment_divergence_filtered" if divergence_filter_sd else "alignment",
        )
        list2_pages_dir      = os.path.join(outdir, "list2", "alignments")
        list2_colnum_dir     = (
            os.path.join(outdir, "list2", "alignment_trimmed_colnumbering")
            if trim_args else None
        )
        list2_divergence_drop_reasons_path = (
            os.path.join(outdir, "list2", "divergence_filter_stats", "drop_reasons.tsv")
            if divergence_filter_sd else None
        )

        if not (resume and _step_complete(list1_pages_dir)):
            logger.info("Rendering %s alignment pages for report...", list1_name)
            alignment_pages1 = build_alignment_pages(
                alignment_dir=list1_alignment_dir,
                tables_dir=tables_dir1,
                outdir=list1_pages_dir,
                group_ids=diff_group_ids,
                anchor_species=anchor,
                colnumbering_dir=list1_colnum_dir,
                divergence_drop_reasons_path=list1_divergence_drop_reasons_path,
                show_progress=args.quiet,
            )
        else:
            logger.info(
                "Skipping %s alignment page generation (output already exists: %s)",
                list1_name, list1_pages_dir,
            )
            alignment_pages1 = _reconstruct_alignment_pages(list1_pages_dir)

        # list2 never contains the anchor species (anchor_species=None
        # disables pinning), but list1's own anchor gene mapping is passed
        # through as page context so list2 pages still show which anchor
        # gene this group corresponds to.
        list1_anchor_lookup = load_group_anchor(tables_dir1)
        if not (resume and _step_complete(list2_pages_dir)):
            logger.info("Rendering %s alignment pages for report...", list2_name)
            alignment_pages2 = build_alignment_pages(
                alignment_dir=list2_alignment_dir,
                tables_dir=tables_dir2,
                outdir=list2_pages_dir,
                group_ids=diff_group_ids,
                anchor_species=None,
                anchor_gene_lookup=list1_anchor_lookup,
                colnumbering_dir=list2_colnum_dir,
                divergence_drop_reasons_path=list2_divergence_drop_reasons_path,
                show_progress=args.quiet,
            )
        else:
            logger.info(
                "Skipping %s alignment page generation (output already exists: %s)",
                list2_name, list2_pages_dir,
            )
            alignment_pages2 = _reconstruct_alignment_pages(list2_pages_dir)

        logger.info("Building HTML report")
        if all_metrics:
            build_multi_metric_report(
                metric_data=metric_outputs,
                report_path=os.path.join(outdir, "report.html"),
                obo_path=obo_path,
                fdr_threshold=fdr_threshold,
                mode="differential",
                gmt_path=gene_sets,
                tables_dir1=tables_dir1,
                tables_dir2=tables_dir2,
                label1=list1_name,
                label2=list2_name,
                pipeline_stats_path1=os.path.join(outdir, "list1", "pipeline_stats.json"),
                pipeline_stats_path2=os.path.join(outdir, "list2", "pipeline_stats.json"),
                differential_stats_path=os.path.join(outdir, "differential_stats.json"),
                alignment_pages1=alignment_pages1,
                alignment_pages2=alignment_pages2,
            )
        else:
            build_report(
                results_dir=metric_outputs[metric]["results_dir"],
                plots_dir=metric_outputs[metric]["plots_dir"],
                report_path=os.path.join(outdir, "differential", "report.html"),
                obo_path=obo_path,
                mode="differential",
                metric=metric,
                fdr_threshold=fdr_threshold,
                gmt_path=gene_sets,
                tables_dir1=tables_dir1,
                tables_dir2=tables_dir2,
                label1=list1_name,
                label2=list2_name,
                pipeline_stats_path1=os.path.join(outdir, "list1", "pipeline_stats.json"),
                pipeline_stats_path2=os.path.join(outdir, "list2", "pipeline_stats.json"),
                differential_stats_path=os.path.join(outdir, "differential_stats.json"),
                alignment_pages1=alignment_pages1,
                alignment_pages2=alignment_pages2,
            )

        for m in metrics_to_run:
            lbl1 = ("%s [%s]" % (list1_name, m)) if all_metrics else list1_name
            lbl2 = ("%s [%s]" % (list2_name, m)) if all_metrics else list2_name
            _log_summary(stats1, metric_outputs[m]["results_dir"],
                         fdr_threshold, label=lbl1)
            _log_summary(stats2, metric_outputs[m]["results_dir"],
                         fdr_threshold, label=lbl2)

    logger.info("Enhydra finished successfully.")
