import os
import sys
import json
import shutil
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
from .provenance import collect_run_parameters, write_run_parameters, \
    add_stage_summary
from .stage_markers import check_stage, write_stage_marker, aggregate_stage_markers
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
    """Legacy presence-only completion check.

    Still used for the small number of downstream artefacts that sit
    outside _run_single_list()'s own per-stage marker system (alignment
    page rendering, GSEA results, differential scoring) — these are
    single-shot outputs of one function call each, not a pipeline of
    independently-parameterised stages, so the added bookkeeping of a
    stage marker was judged not to carry its weight there. See
    stage_markers.check_stage() for the parameter-aware version used
    inside _run_single_list().
    """
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


# Every subdirectory _run_single_list() itself creates directly under a
# single list's root directory (i.e. 'outdir' in single-list mode, or
# 'outdir/list1' / 'outdir/list2' in two-list mode). Kept as an explicit
# list, cross-referenced against _run_single_list()'s own local directory
# variables, rather than derived dynamically — --fork-from's seeding step
# (see _seed_fork_source() below) needs to know exactly which
# subdirectories are pipeline *stage* output (eligible for reuse across a
# fork) as opposed to report-facing artefacts that always regenerate
# fresh (alignments/, enrichment*/, plots*/, report.html,
# run_parameters.json, pipeline_stats.json — none of which
# _run_single_list() itself writes; see _seed_fork_source()'s own
# docstring for why those are deliberately excluded).
_FORK_SEEDABLE_STAGE_DIRS = (
    "subset",
    "length_stats", "length_filter", "length_filter_stats",
    "group_filter", "group_filter_stats",
    "alignment", "alignment_stripped",
    "divergence_sident", "alignment_divergence_filtered",
    "divergence_realign_input", "divergence_filter_stats",
    "alignment_trimmed", "alignment_trimmed_colnumbering",
    "ident_alignment", "tables",
)


def _seed_fork_source(fork_from_listdir: str, dest_listdir: str, label: str = "") -> None:
    """Copy every known pipeline-stage subdirectory (with its resume marker)
    from a prior completed run into a fresh destination list directory.

    This is the mechanism behind --fork-from: rather than recomputing an
    entire pipeline from scratch just to change one downstream parameter
    (e.g. trim from 'strict' to 'strictplus'), every stage subdirectory
    present in the fork source is copied wholesale into the new
    (previously nonexistent) destination directory. _run_single_list()'s
    own per-stage marker mechanism (stage_markers.check_stage(), via
    _must_run()) then takes over exactly as it would for an ordinary
    --resume: a copied stage whose marker still matches the current
    invocation's parameters is reused as-is; the first stage whose marker
    does not match is recomputed, which automatically cascades to every
    stage after it regardless of whether *their* own parameters also
    changed — a stage downstream of a changed one is stale via its input
    even if its own configuration did not change.

    Blindly copying every stage — even ones that will turn out to be
    invalidated by the parameter change — is deliberate and safe: a
    stage whose copied-over output later fails its marker check is simply
    recomputed and overwritten in place, exactly as for a same-directory
    --resume that detects a mid-pipeline parameter change. There is no
    need to compute the fork point in advance; the existing cascade does
    it implicitly.

    Deliberately NOT copied here (always regenerated fresh by the normal
    post-_run_single_list() pipeline in main(), regardless of
    --fork-from): rendered alignment pages, GSEA results, plots, and the
    HTML report, since these are report-facing artefacts derived from
    whichever stage output ends up active after _run_single_list()
    returns — potentially including output freshly recomputed because of
    the very parameter change --fork-from exists to apply — rather than
    independently-parameterised pipeline stages of their own.

    Args:
        fork_from_listdir: A single list's root directory from the prior
                           completed run (the old --outdir in single-list
                           mode, or its 'list1'/'list2' subdirectory in
                           two-list mode).
        dest_listdir:      The corresponding directory in the new run.
                           Any individual stage subdirectory within it
                           must not already exist — the caller (main())
                           is responsible for having verified the overall
                           destination outdir was freshly created for
                           this invocation.
        label:             Optional list name for log messages (two-list
                           mode).

    Raises:
        EnhydraIOError: If fork_from_listdir does not exist.
    """
    logger = logging.getLogger(__name__)
    prefix = ("%s: " % label) if label else ""
    if not os.path.isdir(fork_from_listdir):
        raise EnhydraIOError(
            "%sFork source directory not found: %s" % (prefix, fork_from_listdir)
        )
    n_copied = 0
    for stage_name in _FORK_SEEDABLE_STAGE_DIRS:
        src = os.path.join(fork_from_listdir, stage_name)
        dst = os.path.join(dest_listdir, stage_name)
        if os.path.isdir(src):
            shutil.copytree(src, dst)
            n_copied += 1
    logger.info(
        "%sSeeded %d stage director%s from fork source: %s",
        prefix, n_copied, "y" if n_copied == 1 else "ies", fork_from_listdir,
    )


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


def _gsea_weight(two_list_mode: bool, metric: str) -> float:
    """Determine the GSEA prerank 'weight' parameter for a run mode/metric pair.

    weight=0 makes GSEA's running-sum enrichment statistic purely
    rank-based (equivalent to the classic, unweighted Kolmogorov-Smirnov
    statistic) — genes contribute to the running sum based only on their
    position in the ranking, ignoring the magnitude of their score.
    weight=1 (GSEApy/GSEA's standard) instead weights each gene's
    contribution by the magnitude of its own ranking-metric value, so
    genes with larger scores pull the running sum more strongly.

    Current policy: weight=0 only for single-list 'identity' and 'rank'
    metrics; weight=1 for every other combination (single-list 'zscore',
    and all three metrics in two-list/differential mode).

    Args:
        two_list_mode: Whether this is a two-list/differential run.
        metric:        The ranking metric ('identity', 'zscore', 'rank').

    Returns:
        0.0 or 1.0.
    """
    if not two_list_mode and metric in ("identity", "rank"):
        return 0.0
    return 1.0


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

    # ------------------------------------------------------------------ #
    # Per-stage resume-safety: force_recompute cascades forward the       #
    # instant any stage's existing output can't be verified as matching   #
    # this invocation's parameters (stage_markers.check_stage() status    #
    # 'mismatch' or 'missing_marker'), or the instant any stage simply    #
    # has no output yet ('not_run'). Once set, every later stage in this  #
    # list's pipeline runs unconditionally regardless of its own marker,  #
    # since its input has now changed even if its own parameters have    #
    # not. See stage_markers.py's module docstring for the full design.  #
    # ------------------------------------------------------------------ #
    force_recompute = False

    def _must_run(dir_specs, marker_dir, stage_name, stage_params):
        nonlocal force_recompute
        if not resume:
            return True
        if force_recompute:
            return True
        result = check_stage(dir_specs, marker_dir, stage_params)
        if result.status == "ok":
            logger.info(
                "Skipping %s (parameters unchanged, output exists: %s)",
                _desc(stage_name), marker_dir,
            )
            return False
        if result.status == "mismatch":
            changes = "; ".join(
                "%s: %r -> %r" % (k, old, new)
                for k, (old, new) in result.changed_parameters.items()
            )
            logger.warning(
                "%s: parameters changed since the run that produced its "
                "existing output (%s) — recomputing this stage and every "
                "stage downstream of it.", _desc(stage_name), changes,
            )
            force_recompute = True
            return True
        if result.status == "missing_marker":
            logger.warning(
                "%s: output exists but has no parameter marker (likely "
                "produced by a run from before this feature, or a "
                "corrupted marker file) — cannot verify its parameters "
                "match this invocation; recomputing this stage and every "
                "stage downstream of it to be safe.", _desc(stage_name),
            )
            force_recompute = True
            return True
        # 'not_run' — first time this stage has ever produced output here.
        force_recompute = True
        return True

    with tqdm(total=n_steps, desc=_desc("starting"),
              unit="step", disable=not show_progress, leave=True) as sbar:

        if species is not None:
            sbar.set_description(_desc("subsetting"))
            subset_params = {"species": sorted(species)}
            if _must_run([(subset_dir, None)], subset_dir, "subsetting", subset_params):
                subset_groups(inputdir, subset_dir, species,
                              show_progress=show_progress)
                write_stage_marker(subset_dir, "subset", subset_params)
            source_dir = subset_dir
            sbar.update(1)
        else:
            source_dir = inputdir

        sbar.set_description(_desc("length filter"))
        logger.info("Step 1: Length filtering")
        length_filter_params = {"length_filter_sd": sd_multiplier}
        if _must_run([(length_filter_dir, None)], length_filter_dir,
                     "length filter", length_filter_params):
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
            write_stage_marker(length_filter_dir, "length_filter", length_filter_params)
        # Runs unconditionally, whether or not the filtering above was just
        # skipped via --resume. This is deliberate: aggregation is a cheap,
        # idempotent read of whatever per-group '_lengthstats' files already
        # exist in length_stats_dir — it is not tied to whether filtering
        # happened in *this* invocation. Nesting it inside the `if
        # _must_run(...)` block would silently produce empty/stale summary
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
        group_filter_params = {
            "min_species":    min_species,
            "min_sequences":  min_sequences,
            "paralog_mode":   paralog_mode,
            "anchor":         anchor,
            "require_anchor": require_anchor,
        }
        if _must_run([(group_filter_dir, None)], group_filter_dir,
                     "group filter", group_filter_params):
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
            write_stage_marker(group_filter_dir, "group_filter", group_filter_params)
        sbar.update(1)

        sbar.set_description(_desc("alignment"))
        logger.info("Step 3: Alignment with %s", aligner.upper())
        alignment_params = {"aligner": aligner, "mafft_mode": mafft_mode}
        if _must_run([(alignment_dir, None)], alignment_dir,
                     "alignment", alignment_params):
            run_aligner(
                group_filter_dir=group_filter_dir,
                alignment_dir=alignment_dir,
                aligner=aligner,
                parameters=parameters,
                show_progress=show_progress,
            )
            write_stage_marker(alignment_dir, "alignment", alignment_params)
        sbar.update(1)

        trimal_input_dir = alignment_dir

        if exclude_from_identity:
            sbar.set_description(_desc("stripping anchor"))
            logger.info("Step 3b: Stripping injected species from alignments: %s",
                        exclude_from_identity)
            strip_params = {"exclude_from_identity": sorted(exclude_from_identity)}
            if _must_run([(stripped_dir, None)], stripped_dir,
                         "anchor stripping", strip_params):
                strip_species_from_alignments(
                    alignment_dir=alignment_dir,
                    stripped_dir=stripped_dir,
                    exclude=exclude_from_identity,
                    show_progress=show_progress,
                )
                write_stage_marker(stripped_dir, "anchor_stripping", strip_params)
            trimal_input_dir = stripped_dir
            sbar.update(1)

        if divergence_filter_sd:
            sbar.set_description(_desc("divergence filter"))
            logger.info(
                "Step 3c: Filtering highly divergent sequences "
                "(sd_multiplier=%s)", divergence_filter_sd,
            )
            divergence_params = {
                "divergence_filter_sd": divergence_filter_sd,
                "min_species":          min_species,
                "min_sequences":        min_sequences,
            }
            divergence_dir_specs = [
                (divergence_filtered_dir, None),
                (divergence_stats_dir, ["drop_reasons.tsv"]),
            ]
            if _must_run(divergence_dir_specs, divergence_filtered_dir,
                        "divergence filter", divergence_params):
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
                write_stage_marker(divergence_filtered_dir, "divergence_filter",
                                   divergence_params)
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
            trim_params = {"trim": trim_args}
            trim_dir_specs = [
                (trimmed_dir, None),
                (trim_colnumbering_dir, None),
            ]
            if _must_run(trim_dir_specs, trimmed_dir,
                        "column trimming", trim_params):
                run_trimal_columns(
                    alignment_dir=trimal_input_dir,
                    trimmed_dir=trimmed_dir,
                    trimal_path=trimal_path,
                    trim_args=trim_args,
                    n_proc=max_process,
                    show_progress=show_progress,
                    colnumbering_dir=trim_colnumbering_dir,
                )
                write_stage_marker(trimmed_dir, "column_trimming", trim_params)
            trimal_input_dir = trimmed_dir
            sbar.update(1)

        sbar.set_description(_desc("identity"))
        logger.info("Step 4: Identity estimation with trimAl")
        # Identity estimation has no tunable parameters of its own — its
        # correctness depends entirely on which alignment directory feeds
        # it (trimal_input_dir, selected by the stages above). It is still
        # tracked as its own marked stage purely so that force_recompute
        # cascading from any upstream stage correctly forces this one to
        # rerun too, even though its own parameter dict never changes.
        identity_params: dict = {}
        if _must_run([(ident_dir, None)], ident_dir,
                     "identity estimation", identity_params):
            run_trimal(
                alignment_dir=trimal_input_dir,
                ident_dir=ident_dir,
                trimal_path=trimal_path,
                n_proc=max_process,
                show_progress=show_progress,
            )
            write_stage_marker(ident_dir, "identity", identity_params)
        sbar.update(1)

        sbar.set_description(_desc("tables"))
        logger.info("Step 5: Generating tables")
        tables_params = {
            "anchor":                   anchor,
            "require_anchor_in_tables": require_anchor_in_tables,
        }
        tables_dir_specs = [
            (tables_dir, ["group2mean.tsv", "anchor2mean.tsv", "group2anchor.tsv"]),
        ]
        if _must_run(tables_dir_specs, tables_dir, "table generation", tables_params):
            make_tables(
                alignment_dir=alignment_dir,
                ident_dir=ident_dir,
                tables_dir=tables_dir,
                anchor=anchor,
                require_anchor=require_anchor_in_tables,
                show_progress=show_progress,
            )
            write_stage_marker(tables_dir, "tables", tables_params)
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
    input_group.add_argument("--orthofinder-dir", default=None,
                             help="Path to an OrthoFinder 3 output directory. "
                                  "Can also be set via 'orthofinder_dir' in "
                                  "the project config.")

    parser.add_argument("--resume",  action="store_true", default=None,
                        help="Resume a previously interrupted run. Each "
                             "pipeline stage's existing output is reused "
                             "only if a per-stage marker confirms it was "
                             "produced with the same parameters as this "
                             "invocation; otherwise that stage and every "
                             "stage downstream of it are recomputed, with "
                             "a warning. Can also be set via 'resume' in "
                             "the project config.")
    parser.add_argument(
        "--fork-from", default=None, metavar="OUTDIR",
        help="Seed a new output directory from a prior completed run's "
             "pipeline-stage outputs, to change a parameter from some "
             "stage onward without recomputing everything before it "
             "(e.g. changing --trim from 'strict' to 'strictplus' reuses "
             "length/group filtering and alignment, and only recomputes "
             "trimming, identity estimation, and tables). Every stage "
             "whose parameters are unchanged from the fork source is "
             "reused as-is; the first stage whose parameters differ, and "
             "every stage after it, is recomputed. The destination "
             "'outdir' (in the project config) must not already exist. "
             "Implies --resume. Can also be set via 'fork_from' in the "
             "project config."
    )
    parser.add_argument("--replot",  action="store_true", default=None,
                        help="Re-run GSEA, plots, and the HTML report without "
                             "repeating alignment. Implies --resume. Can also "
                             "be set via 'replot' in the project config.")
    parser.add_argument("--quiet",   action="store_true", default=None,
                        help="Suppress INFO/WARNING on the console; show progress "
                             "bars instead. Can also be set via 'quiet' in the "
                             "project config.")
    parser.add_argument("--paralogs", choices=["all", "remove", "longest"], default=None)
    parser.add_argument("--min-species", type=int, default=None)
    parser.add_argument("--length-filter-sd", type=float, default=None,
                        help="Number of standard deviations from a group's "
                             "mean sequence length beyond which a sequence "
                             "is removed by the length filter (default: "
                             "2.0, or the value of 'length_filter_sd' in "
                             "the project config).")
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
    parser.add_argument("--all-metrics", action="store_true", default=None,
                        help="Run GSEA for all three ranking metrics and produce "
                             "a tabbed HTML report. Can also be set via "
                             "'all_metrics' in the project config.")

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
    parser.add_argument(
        "--gene-list-fdr-threshold", type=float, default=None, metavar="FDR",
        help="(Advanced) FDR threshold below which a gene set's full gene "
             "list and leading-edge gene list are embedded in the HTML "
             "report as clickable 'View' links; gene sets at or above "
             "this threshold show a count only, to keep report size "
             "manageable for large (e.g. genome-scale) gene set "
             "collections. Defaults to the same value as --fdr-threshold "
             "(i.e. only significant gene sets get embedded gene lists). "
             "Set to 1.0 or higher to embed gene lists for virtually all "
             "tested gene sets (restores the original, unrestricted "
             "behaviour — can produce very large reports). Can also be "
             "set via 'gene_list_fdr_threshold' in the project config."
    )
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

    # (Advanced) Independent FDR cutoff controlling which gene sets get
    # their full gene / leading-edge gene lists embedded in the report
    # (see report._results_table_html()). Kept separate from
    # fdr_threshold itself (which still governs significance everywhere
    # else) so a person can keep a lenient --fdr-threshold for the results
    # table while still capping report size. '' / unset means "use
    # fdr_threshold" — i.e. only significant gene sets get embedded gene
    # lists, which is the default, restricted behaviour.
    gene_list_fdr_threshold_raw = _resolve(
        args.gene_list_fdr_threshold, parameters['gene_list_fdr_threshold'], '')
    try:
        gene_list_fdr_threshold = (
            float(gene_list_fdr_threshold_raw)
            if gene_list_fdr_threshold_raw not in (None, '') else None
        )
    except (TypeError, ValueError):
        sys.exit(
            "Configuration error: invalid 'gene_list_fdr_threshold' value: "
            "%r. Must be a number." % (gene_list_fdr_threshold_raw,)
        )

    top_n         = _resolve(args.top_n,         parameters['top_n'],         20)
    list1_path    = _resolve(args.list1,         parameters['list1'],         None)
    list2_path    = _resolve(args.list2,         parameters['list2'],         None)
    list1_name    = parameters['list1_name']
    list2_name    = parameters['list2_name']
    sources_raw   = _resolve(args.sources,       parameters['sources'],
                             'GO:BP GO:MF GO:CC KEGG REAC')
    sources       = sources_raw if isinstance(sources_raw, list) \
                    else sources_raw.split()
    sd_multiplier = _resolve(args.length_filter_sd, parameters['length_filter_sd'], 2.0)
    aligner       = parameters['aligner']
    mafft_mode    = parameters['mafft_mode']
    obo_cache     = _resolve(args.obo_cache, parameters['obo_cache'], None)

    # Run-mode flags: previously CLI-only (argparse defaults of False/None
    # with no config fallback). All five now accept a project-config value
    # too — argparse's store_true flags use default=None (not False) so
    # that "flag not passed on the CLI" can be distinguished from
    # "explicitly disabled", letting _resolve() fall through to the config
    # value in the former case. An explicit CLI flag always wins over the
    # config file, consistent with every other parameter in this function.
    orthofinder_dir = _resolve(args.orthofinder_dir, parameters['orthofinder_dir'], None) or None
    fork_from       = _resolve(args.fork_from,       parameters['fork_from'],       None) or None
    resume_flag     = _resolve(args.resume,          parameters['resume'],          False)
    replot          = _resolve(args.replot,          parameters['replot'],          False)
    quiet           = _resolve(args.quiet,           parameters['quiet'],           False)
    all_metrics     = _resolve(args.all_metrics,     parameters['all_metrics'],     False)
    resume          = resume_flag or replot

    try:
        trim_args = resolve_trim_args(trim)
    except EnhydraConfigError as e:
        sys.exit("Configuration error: %s" % e)

    try:
        divergence_filter_sd = resolve_divergence_filter_sd(divergence_filter_sd_raw)
    except EnhydraConfigError as e:
        sys.exit("Configuration error: %s" % e)

    # Snapshot of every effective (post-resolution) parameter value, for the
    # reproducibility record written via provenance.collect_run_parameters()
    # below. Built once, right after all _resolve()/resolve_*() calls above
    # have settled, so it reflects exactly what this run actually used —
    # not the raw config file contents, which a CLI flag may have silently
    # overridden. Keys match collect_run_parameters()'s resolved_params
    # contract exactly.
    #
    # Note: gene_list_fdr_threshold is NOT yet included here, so it is not
    # currently recorded in run_parameters.json / shown in the report's
    # Run Parameters tab — a small follow-up, not yet done.
    resolved_params_for_provenance = {
        "min_species":          min_species,
        "min_sequences":        min_sequences,
        "paralogs":             paralogs,
        "length_filter_sd":     sd_multiplier,
        "divergence_filter_sd": divergence_filter_sd,
        "trim":                 trim,
        "aligner":              aligner,
        "mafft_mode":           mafft_mode,
        "metric":               metric,
        "permutations":         permutations,
        "min_size":             min_size,
        "max_size":             max_size,
        "seed":                 seed,
        "fdr_threshold":        fdr_threshold,
        "top_n":                top_n,
        "gene_sets":            gene_sets,
        "organism":             organism,
        "sources":              sources,
    }

    if not gene_sets and not organism:
        parser.error(
            "A gene set source is required. Set 'gene_sets' or 'organism' in "
            "your project config, or use --gene-sets / --organism."
        )

    two_list_mode = bool(list1_path or list2_path)
    if two_list_mode and not (list1_path and list2_path):
        parser.error("Two-list mode requires both list1 and list2.")

    outdir = parameters['outdir']

    if fork_from:
        if os.path.isdir(outdir):
            sys.exit(
                "--fork-from requires a new output directory: '%s' already "
                "exists. Choose a fresh 'outdir' in your project config, or "
                "omit --fork-from and use --resume instead if you intend to "
                "continue this existing directory in place." % outdir
            )
        if not os.path.isdir(fork_from):
            sys.exit("--fork-from directory not found: %s" % fork_from)
        fork_run_params_path = os.path.join(fork_from, "run_parameters.json")
        if os.path.isfile(fork_run_params_path):
            with open(fork_run_params_path) as fh:
                fork_run_params = json.load(fh)
            fork_inputdir = (fork_run_params.get("input") or {}).get("inputdir")
            this_inputdir = os.path.abspath(parameters['inputdir'])
            if fork_inputdir and os.path.abspath(fork_inputdir) != this_inputdir:
                sys.exit(
                    "--fork-from source was run with a different inputdir "
                    "('%s') than this invocation's ('%s'). Forking across "
                    "different input data is not supported — the "
                    "per-stage markers this feature relies on do not "
                    "themselves track inputdir." % (fork_inputdir, this_inputdir)
                )
        else:
            sys.stderr.write(
                "Warning: --fork-from source has no run_parameters.json "
                "(it predates that feature) — cannot verify it used the "
                "same inputdir as this invocation.\n"
            )
        # The seeded stage directories are only useful if the normal
        # pipeline actually consults their markers instead of blindly
        # recomputing everything — i.e. this must behave as a resumed
        # run regardless of whether the user also passed --resume.
        resume = True
    elif os.path.isdir(outdir) and not resume:
        sys.exit(
            "Output directory '%s' already exists. Use --resume to continue "
            "a previous run, --replot to re-run GSEA and regenerate plots, "
            "--fork-from <path> to seed a new outdir from that run's "
            "outputs while changing selected parameters, or change "
            "'outdir' in your project config." % outdir
        )
    os.makedirs(outdir, exist_ok=True)

    _setup_logging(outdir, quiet=quiet)
    logger = logging.getLogger(__name__)
    logger.info("Welcome to Enhydra")
    logger.info(
        "Resolved parameters: metric=%s, all_metrics=%s, replot=%s, "
        "paralogs=%s, trim=%s, divergence_filter_sd=%s, min_species=%d, "
        "permutations=%d, fdr_threshold=%.2f, gene_list_fdr_threshold=%s",
        metric, all_metrics, replot, paralogs, (trim or "none"),
        (divergence_filter_sd if divergence_filter_sd else "none"),
        min_species, permutations, fdr_threshold,
        (gene_list_fdr_threshold if gene_list_fdr_threshold is not None
         else "same as fdr_threshold"),
    )

    if orthofinder_dir:
        logger.info("OrthoFinder mode: preprocessing Orthogroup_Sequences/")
        preprocess_orthofinder(
            orthofinder_dir=orthofinder_dir,
            inputdir=parameters['inputdir'],
        )

    if fork_from:
        logger.info("Seeding pipeline-stage outputs from fork source: %s", fork_from)
        try:
            if two_list_mode:
                _seed_fork_source(os.path.join(fork_from, "list1"),
                                  os.path.join(outdir, "list1"), label=list1_name)
                _seed_fork_source(os.path.join(fork_from, "list2"),
                                  os.path.join(outdir, "list2"), label=list2_name)
            else:
                _seed_fork_source(fork_from, outdir)
        except EnhydraIOError as e:
            sys.exit("--fork-from error: %s" % e)

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
        show_progress=quiet,
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

    run_parameters_path = os.path.join(outdir, "run_parameters.json")

    # ------------------------------------------------------------------ #
    #  Single-list mode                                                    #
    # ------------------------------------------------------------------ #
    if not two_list_mode:
        logger.info("Running in single-list mode.")

        run_params = collect_run_parameters(
            outdir=outdir,
            code_config_path=args.code_config,
            project_config_path=args.project_config,
            command_line=sys.argv,
            inputdir=parameters['inputdir'],
            anchor=parameters['anchor'],
            resolved_params=resolved_params_for_provenance,
            two_list_mode=False,
            orthofinder_dir=orthofinder_dir,
            all_metrics=all_metrics,
            metrics_run=metrics_to_run,
        )
        write_run_parameters(run_params, outdir)

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
        # Length filtering always runs (unlike the divergence filter), so
        # this path is always populated — no conditional needed. The
        # underlying file is written unconditionally by
        # aggregate_length_filter_stats(), even when nothing was dropped.
        length_filter_drop_reasons_path_single = os.path.join(
            outdir, "length_filter_stats", "drop_reasons.tsv"
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
                length_filter_drop_reasons_path=length_filter_drop_reasons_path_single,
                show_progress=quiet,
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
                         results_dir=results_dir_m,
                         weight=_gsea_weight(two_list_mode, m),
                         **gsea_kwargs)
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
                run_parameters_path=run_parameters_path,
                gene_list_fdr_threshold=gene_list_fdr_threshold,
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
                run_parameters_path=run_parameters_path,
                gene_list_fdr_threshold=gene_list_fdr_threshold,
            )

        # Fold every stage marker produced under outdir into
        # run_parameters.json's 'stages' section, for human/report-tab
        # visibility only — see stage_markers.aggregate_stage_markers()
        # and provenance.add_stage_summary() docstrings. This is purely
        # additive bookkeeping and deliberately happens last, after every
        # stage that could possibly run this invocation already has.
        add_stage_summary(run_parameters_path, aggregate_stage_markers(outdir))

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

        run_params = collect_run_parameters(
            outdir=outdir,
            code_config_path=args.code_config,
            project_config_path=args.project_config,
            command_line=sys.argv,
            inputdir=parameters['inputdir'],
            anchor=anchor,
            resolved_params=resolved_params_for_provenance,
            two_list_mode=True,
            orthofinder_dir=orthofinder_dir,
            list1_path=list1_path,
            list2_path=list2_path,
            list1_name=list1_name,
            list2_name=list2_name,
            species1=species1,
            species2=species2,
            anchor_injected=anchor_injected,
            all_metrics=all_metrics,
            metrics_run=metrics_to_run,
        )
        write_run_parameters(run_params, outdir)

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
                    weight=_gsea_weight(two_list_mode, m),
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
        list1_length_filter_drop_reasons_path = os.path.join(
            outdir, "list1", "length_filter_stats", "drop_reasons.tsv"
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
        list2_length_filter_drop_reasons_path = os.path.join(
            outdir, "list2", "length_filter_stats", "drop_reasons.tsv"
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
                length_filter_drop_reasons_path=list1_length_filter_drop_reasons_path,
                show_progress=quiet,
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
                length_filter_drop_reasons_path=list2_length_filter_drop_reasons_path,
                show_progress=quiet,
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
                run_parameters_path=run_parameters_path,
                gene_list_fdr_threshold=gene_list_fdr_threshold,
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
                run_parameters_path=run_parameters_path,
                gene_list_fdr_threshold=gene_list_fdr_threshold,
            )

        # See the single-list branch's identical call for the full
        # rationale — folds both lists' stage markers into
        # run_parameters.json's 'stages' section, nested per list since
        # list1/list2 run their pipeline stages independently.
        add_stage_summary(run_parameters_path, {
            "list1": aggregate_stage_markers(os.path.join(outdir, "list1")),
            "list2": aggregate_stage_markers(os.path.join(outdir, "list2")),
        })

        for m in metrics_to_run:
            lbl1 = ("%s [%s]" % (list1_name, m)) if all_metrics else list1_name
            lbl2 = ("%s [%s]" % (list2_name, m)) if all_metrics else list2_name
            _log_summary(stats1, metric_outputs[m]["results_dir"],
                         fdr_threshold, label=lbl1)
            _log_summary(stats2, metric_outputs[m]["results_dir"],
                         fdr_threshold, label=lbl2)

    logger.info("Enhydra finished successfully.")
