from __future__ import annotations

import os
import json
import logging

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Small file-reading helpers
# ---------------------------------------------------------------------------

def _count_files(directory: str) -> int:
    if not os.path.isdir(directory):
        return 0
    return len(os.listdir(directory))


def _count_lines(path: str) -> int:
    if not os.path.isfile(path):
        return 0
    with open(path) as fh:
        return sum(1 for line in fh if line.strip())


def _count_reasons_by_column(path: str, reason_col: str = "reason") -> dict[str, int]:
    """Read a drop-reasons-style TSV and return {reason: count}.

    Locates the reason column by header name rather than a fixed index,
    since different drop_reasons.tsv files across stages have different
    column layouts (e.g. length_filter's per-sequence drop_reasons.tsv has
    'reason' as its last column, while group_filter's and tables' have it
    second). Returns an empty dict if the file doesn't exist (e.g. a stage
    that hasn't run yet, or was skipped under --resume).
    """
    if not os.path.isfile(path):
        return {}
    counts: dict[str, int] = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        if reason_col not in header:
            logger.warning(
                "Column '%s' not found in %s — cannot summarise reasons.",
                reason_col, path,
            )
            return {}
        idx = header.index(reason_col)
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) <= idx:
                continue
            reason = fields[idx]
            counts[reason] = counts.get(reason, 0) + 1
    return counts


def _read_group_ids(group2mean_path: str) -> set[str]:
    if not os.path.isfile(group2mean_path):
        return set()
    ids: set[str] = set()
    with open(group2mean_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            ids.add(line.split("\t", 1)[0])
    return ids


# ---------------------------------------------------------------------------
# Per-list aggregation
# ---------------------------------------------------------------------------

def aggregate_pipeline_stats(
    listdir: str,
    n_input: int,
    output_path: str | None = None,
) -> dict:
    """Consolidate every filtering stage's stats for one list into one summary.

    Reads the outputs already written by filter_length()/
    aggregate_length_filter_stats(), filter_groups(),
    filter_divergent_sequences() (if the divergence filter was enabled for
    this run), and make_tables() — all of which live under fixed
    subdirectory names of `listdir` — and produces a single
    JSON-serialisable summary combining funnel counts (how many groups
    survive each stage) with a breakdown of *why* groups or sequences were
    dropped at each stage.

    This function only reads existing output files; it does not touch or
    require any of filtering.py/tables.py internals directly, so it can be
    run independently of (and after) the main pipeline, including on a
    partially-resumed run — stages whose stats files don't exist yet simply
    contribute empty reason breakdowns and are not treated as an error.

    The divergence filter step (filtering.filter_divergent_sequences()) is
    optional and only runs when the pipeline's divergence_filter_sd
    parameter is set — unlike every other stage here, it has no fixed
    on/off state across runs. Its presence is therefore auto-detected from
    whether 'divergence_filter_stats/drop_reasons.tsv' exists: that file is
    written unconditionally by filter_divergent_sequences() whenever it
    runs at all (even if it removes nothing), so its existence is a
    reliable "did this step run" signal, distinct from "did it drop
    anything". When the step did not run, 'n_after_divergence_filter' is
    None and no 'divergence_filter' key appears in drop_summary or
    detail_files at all — callers (e.g. the report's funnel/drop-reason
    rendering) are expected to treat an absent key as "not applicable",
    not as "zero drops", since those mean different things to a reader
    deciding whether the filter was even used for this run.

    Args:
        listdir:     Root directory for this list (e.g. outdir/list1, or
                     outdir itself in single-list mode). Expected to
                     contain 'length_filter/', 'length_filter_stats/',
                     'group_filter/', 'group_filter_stats/', 'tables/',
                     and — only if divergence filtering was enabled —
                     'alignment_divergence_filtered/' and
                     'divergence_filter_stats/', as produced by
                     _run_single_list().
        n_input:     Number of input groups for this list, before any
                     filtering (from _run_single_list()'s returned stats
                     dict, since that's the only stage whose source
                     directory isn't a fixed, re-discoverable path here).
        output_path: Where to write the JSON summary. Defaults to
                     '<listdir>/pipeline_stats.json' if not given.

    Returns:
        The summary dict (also written to output_path as JSON).
    """
    length_filter_dir          = os.path.join(listdir, "length_filter")
    group_filter_dir           = os.path.join(listdir, "group_filter")
    length_filter_stats_dir    = os.path.join(listdir, "length_filter_stats")
    group_filter_stats_dir     = os.path.join(listdir, "group_filter_stats")
    divergence_filtered_dir    = os.path.join(listdir, "alignment_divergence_filtered")
    divergence_filter_stats_dir = os.path.join(listdir, "divergence_filter_stats")
    tables_dir                  = os.path.join(listdir, "tables")

    n_after_length_filter = _count_files(length_filter_dir)
    n_after_group_filter  = _count_files(group_filter_dir)

    divergence_drop_reasons_path = os.path.join(
        divergence_filter_stats_dir, "drop_reasons.tsv"
    )
    divergence_filter_ran = os.path.isfile(divergence_drop_reasons_path)
    n_after_divergence_filter = (
        _count_files(divergence_filtered_dir) if divergence_filter_ran else None
    )

    group2mean_path  = os.path.join(tables_dir, "group2mean.tsv")
    anchor2mean_path = os.path.join(tables_dir, "anchor2mean.tsv")
    n_final_groups   = _count_lines(group2mean_path)
    n_anchor_mapped  = _count_lines(anchor2mean_path)

    length_groups_skipped = _count_reasons_by_column(
        os.path.join(length_filter_stats_dir, "skipped_groups.tsv")
    )
    length_seqs_removed = _count_reasons_by_column(
        os.path.join(length_filter_stats_dir, "drop_reasons.tsv")
    )
    group_filter_dropped = _count_reasons_by_column(
        os.path.join(group_filter_stats_dir, "drop_reasons.tsv")
    )
    tables_dropped = _count_reasons_by_column(
        os.path.join(tables_dir, "drop_reasons.tsv")
    )

    # filter_divergent_sequences() writes both of its reasons
    # ('removed_divergent_sequence' and
    # 'below_min_species_after_divergence_filter') into the same
    # drop_reasons.tsv, one 'reason' column shared with every other
    # stage's file. _count_reasons_by_column() already returns a
    # {reason: count} dict from that single file in one pass; splitting
    # it into the two semantically distinct buckets below (a sequence
    # being removed vs. an entire group being dropped) just means
    # picking each known reason string out of that one dict, the same
    # split already applied to the length filter's own two categories.
    divergence_all_reasons = (
        _count_reasons_by_column(divergence_drop_reasons_path)
        if divergence_filter_ran else {}
    )
    divergence_sequences_removed = {
        k: v for k, v in divergence_all_reasons.items()
        if k == "removed_divergent_sequence"
    }
    divergence_groups_dropped = {
        k: v for k, v in divergence_all_reasons.items()
        if k == "below_min_species_after_divergence_filter"
    }

    def _rel(path: str) -> str | None:
        return os.path.relpath(path, listdir) if os.path.isfile(path) else None

    summary = {
        "n_input":                    n_input,
        "n_after_length_filter":      n_after_length_filter,
        "n_after_group_filter":       n_after_group_filter,
        "n_after_divergence_filter":  n_after_divergence_filter,
        "n_final_groups":             n_final_groups,
        "n_anchor_mapped":            n_anchor_mapped,
        "drop_summary": {
            "length_filter": {
                "groups_skipped":   length_groups_skipped,
                "sequences_removed": length_seqs_removed,
            },
            "group_filter": group_filter_dropped,
            "tables":       tables_dropped,
        },
        "detail_files": {
            "length_filter_skipped_groups": _rel(
                os.path.join(length_filter_stats_dir, "skipped_groups.tsv")),
            "length_filter_drop_reasons": _rel(
                os.path.join(length_filter_stats_dir, "drop_reasons.tsv")),
            "group_filter_drop_reasons": _rel(
                os.path.join(group_filter_stats_dir, "drop_reasons.tsv")),
            "group_filter_species_counts": _rel(
                os.path.join(group_filter_stats_dir, "species_counts.tsv")),
            "tables_drop_reasons": _rel(
                os.path.join(tables_dir, "drop_reasons.tsv")),
        },
    }

    if divergence_filter_ran:
        summary["drop_summary"]["divergence_filter"] = {
            "sequences_removed": divergence_sequences_removed,
            "groups_dropped":    divergence_groups_dropped,
        }
        summary["detail_files"]["divergence_filter_drop_reasons"] = _rel(
            divergence_drop_reasons_path
        )

    if output_path is None:
        output_path = os.path.join(listdir, "pipeline_stats.json")
    with open(output_path, "w") as fh:
        json.dump(summary, fh, indent=2)
    logger.info("Pipeline stats summary written to: %s", output_path)

    return summary


# ---------------------------------------------------------------------------
# Two-list overlap stats
# ---------------------------------------------------------------------------

def compute_differential_stats(
    tables_dir1: str,
    tables_dir2: str,
    list1_name: str = "List 1",
    list2_name: str = "List 2",
    output_path: str | None = None,
) -> dict:
    """Summarise group overlap between two lists ahead of differential scoring.

    This mirrors the intersection compute_differential() performs internally,
    but is computed independently here (by directly reading both lists'
    group2mean.tsv files) so it can be reported once regardless of how many
    ranking metrics are run — compute_differential() itself runs once per
    metric (identity/zscore/rank), and duplicating this computation there
    three times, or threading a "only report stats once" flag through it,
    would entangle reporting concerns with its actual scoring logic.

    Args:
        tables_dir1: Path to list 1's tables/ directory.
        tables_dir2: Path to list 2's tables/ directory.
        list1_name:  Display name for list 1 (from project config).
        list2_name:  Display name for list 2 (from project config).
        output_path: Where to write the JSON summary. If not given, the
                     summary is only returned, not written to disk — since,
                     unlike per-list stats, this doesn't naturally belong
                     inside either list's own directory (it's a property of
                     the pair), the caller is expected to pass an explicit
                     top-level path (e.g. outdir/differential_stats.json).

    Returns:
        The summary dict (also written to output_path as JSON, if given).
    """
    groups1 = _read_group_ids(os.path.join(tables_dir1, "group2mean.tsv"))
    groups2 = _read_group_ids(os.path.join(tables_dir2, "group2mean.tsv"))
    common  = groups1 & groups2

    summary = {
        "list1_name":       list1_name,
        "list2_name":       list2_name,
        "n_groups_list1":   len(groups1),
        "n_groups_list2":   len(groups2),
        "n_common_groups":  len(common),
        "n_list1_only":     len(groups1 - groups2),
        "n_list2_only":     len(groups2 - groups1),
    }

    if output_path:
        with open(output_path, "w") as fh:
            json.dump(summary, fh, indent=2)
        logger.info("Differential overlap stats written to: %s", output_path)

    return summary
