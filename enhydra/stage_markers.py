"""Per-stage resume markers.

ENHYDRA's --resume previously decided whether to skip a pipeline stage
using only filesystem presence (see cli.py's old _step_complete(): "does
this directory exist and have files in it"). That check cannot detect
that the *parameters* used to produce an existing stage's output differ
from the parameters requested in the current invocation — e.g. resuming
with a different --min-species than the run that originally produced
group_filter/ silently reuses the old (now-wrong) output.

This module adds a small JSON marker file written into each stage's own
output directory, recording the parameters that were actually used to
produce that directory's contents. --resume can then verify "does this
stage's marker match what I'd use if I ran it now" before deciding to
skip it, rather than trusting bare directory presence alone.

Design notes:
  - One marker file per stage output directory (hidden file, MARKER_FILENAME),
    so deleting a stage's output directory naturally removes its marker
    too — no separate bookkeeping to keep in sync.
  - This module is the sole authority --resume actually consults at
    decision time. provenance.write_run_parameters()'s aggregated
    "stages" section (see aggregate_stage_markers() below) is a read-only
    summary for human/report-tab inspection and is never itself read back
    to make a skip/rerun decision.
  - A stage whose output exists but has no marker (e.g. produced by a
    run predating this feature) or whose marker file is unreadable/
    corrupted is treated the same as a parameter mismatch: cannot verify
    safety, so recompute rather than silently trust it.
"""
from __future__ import annotations

import os
import json
import logging
import datetime
from dataclasses import dataclass

logger = logging.getLogger(__name__)

MARKER_FILENAME = ".stage_marker.json"


def _enhydra_version() -> str:
    """Best-effort lookup of the installed enhydra package version.

    Duplicated from provenance.py's own private helper of the same
    intent, rather than imported, to keep this module's only coupling to
    the rest of the package at the __version__ constant itself (see
    provenance.py for the identical rationale on falling back to
    'unknown' rather than raising).
    """
    try:
        from . import __version__
        return __version__
    except Exception:
        return "unknown"


@dataclass
class StageCheckResult:
    """Result of check_stage().

    Attributes:
        status: One of:
            'not_run'        — the stage has no output yet at all. Normal
                               state for a stage that hasn't run in this
                               outdir before; the caller should simply run
                               it, no warning needed.
            'ok'             — output exists, a marker exists, and its
                               recorded parameters exactly match the
                               parameters given to check_stage(). Safe to
                               skip.
            'mismatch'       — output exists, a marker exists, but its
                               recorded parameters differ from the
                               parameters given to check_stage(). The
                               caller must recompute this stage (and,
                               because its output will now change,
                               every stage downstream of it too).
            'missing_marker' — output exists (by the same presence check
                               used for 'ok'/'mismatch') but no marker
                               file is present, or the marker file exists
                               but could not be parsed. This cannot be
                               distinguished from a genuine parameter
                               mismatch, so it is treated identically:
                               the caller must recompute this stage and
                               everything downstream of it.
        changed_parameters: For status == 'mismatch' only: a dict mapping
                            each differing parameter key to a
                            (previous_value, current_value) tuple, for
                            constructing a clear warning message. None
                            for every other status.
    """
    status: str
    changed_parameters: dict | None = None


def _dir_has_output(stage_dir: str, sentinel_files: list[str] | None) -> bool:
    """Does stage_dir contain output, ignoring the marker file itself?

    Mirrors cli.py's original _step_complete() exactly, except that the
    marker file (MARKER_FILENAME) is never itself counted as output in
    the generic (sentinel_files=None) case — otherwise a stage directory
    containing nothing but its own marker (e.g. a stage that legitimately
    produces zero output files, however unlikely) would be misreported as
    'has output' purely because the marker was written into it.
    """
    if not os.path.isdir(stage_dir):
        return False
    if sentinel_files:
        return all(
            os.path.isfile(os.path.join(stage_dir, f)) and
            os.path.getsize(os.path.join(stage_dir, f)) > 0
            for f in sentinel_files
        )
    return any(f != MARKER_FILENAME for f in os.listdir(stage_dir))


def list_data_files(directory: str) -> list[str]:
    """List directory contents, excluding this module's own marker file.

    write_stage_marker() writes MARKER_FILENAME directly into a stage's
    own output directory (e.g. alignment/), alongside that stage's real
    output — the only place a fixed, well-known marker filename can live
    such that deleting the stage's output directory naturally removes its
    marker too (see this module's own top docstring for that rationale).
    But that means any *later* stage which reads a prior stage's output
    directory wholesale via plain os.listdir() to process "every file in
    it" will also encounter the marker file, and (having no reason to
    expect it) will typically try to treat it as if it were real pipeline
    data — a divergent-sequence filter or trimAl invocation choking on
    '.stage_marker.json.aln' being a concrete example of exactly this
    happening; a directory-indexing helper deriving a bogus group_id from
    it (e.g. an empty string, if it splits a filename on its first '.')
    being a quieter, more insidious example.

    Every call site anywhere in this package that lists a stage's output
    directory in order to iterate over its *contents as pipeline data*
    (as opposed to merely checking whether the directory has any output
    at all, which _dir_has_output() above already handles safely) should
    use this function instead of calling os.listdir() directly, so that
    adding a new pipeline stage automatically stays safe against this
    class of bug without needing to remember to special-case the marker
    filename at every new call site individually.

    Args:
        directory: Path to list.

    Returns:
        The same list os.listdir(directory) would return, with
        MARKER_FILENAME removed if present. Not sorted — callers that
        need a deterministic order should sort the result themselves,
        exactly as they would have needed to with plain os.listdir().
    """
    return [f for f in os.listdir(directory) if f != MARKER_FILENAME]


def write_stage_marker(marker_dir: str, stage: str, parameters: dict) -> None:
    """Write (or overwrite) a stage's resume marker.

    Args:
        marker_dir: Directory to write the marker file into. Created if
                   it does not already exist.
        stage:      Short stage identifier (e.g. 'length_filter'), stored
                   in the marker purely for human readability when
                   inspecting the file directly or via the aggregated
                   summary in run_parameters.json — not itself checked
                   for consistency by check_stage().
        parameters: The exact parameter values that determine this
                   stage's output. Should be a flat, JSON-serialisable
                   dict of only the parameters that actually affect this
                   stage's output — including irrelevant parameters would
                   cause spurious 'mismatch' results for changes that
                   never actually affected this stage.
    """
    os.makedirs(marker_dir, exist_ok=True)
    marker = {
        "stage":           stage,
        "completed_at":    datetime.datetime.now().isoformat(timespec="seconds"),
        "enhydra_version": _enhydra_version(),
        "parameters":      parameters,
    }
    with open(os.path.join(marker_dir, MARKER_FILENAME), "w") as fh:
        json.dump(marker, fh, indent=2)


def check_stage(
    dir_specs: list[tuple[str, list[str] | None]],
    marker_dir: str,
    parameters: dict,
) -> StageCheckResult:
    """Check whether a stage's existing output can safely be reused.

    Args:
        dir_specs:  List of (directory, sentinel_files) pairs describing
                   every directory this stage's "has output" check must
                   consider. Most stages have a single output directory
                   and pass a one-element list, e.g.
                   [(alignment_dir, None)]. A stage whose completion
                   depends on more than one directory (e.g. the
                   divergence filter, which needs both its filtered
                   alignments dir AND its stats dir to be genuinely
                   complete) lists every one of them; ALL must have
                   output for the stage to be considered as having run
                   at all. sentinel_files, if given for a directory,
                   requires those *specific* named files to exist with
                   nonzero size in that directory (as opposed to "any
                   file at all") — see tables.make_tables()'s three
                   named output files for an example.
        marker_dir: Which one of the directories named in dir_specs
                   actually holds (or should hold) this stage's marker
                   file. Must be one of the directories in dir_specs,
                   conventionally the stage's primary/definitive output
                   directory rather than an auxiliary stats directory.
        parameters: The current invocation's resolved parameter values
                   for this stage, in the same shape write_stage_marker()
                   was given when this stage last actually ran.

    Returns:
        A StageCheckResult — see that class's own docstring for the full
        meaning of each status value.
    """
    if not all(_dir_has_output(d, sf) for d, sf in dir_specs):
        return StageCheckResult("not_run")

    marker_path = os.path.join(marker_dir, MARKER_FILENAME)
    if not os.path.isfile(marker_path):
        return StageCheckResult("missing_marker")

    try:
        with open(marker_path) as fh:
            marker = json.load(fh)
        marker_params = marker["parameters"]
    except (OSError, ValueError, KeyError):
        logger.warning(
            "Stage marker at %s could not be read or is malformed — "
            "treating as missing.", marker_path,
        )
        return StageCheckResult("missing_marker")

    if marker_params == parameters:
        return StageCheckResult("ok")

    keys = set(marker_params) | set(parameters)
    changed = {
        k: (marker_params.get(k, "<unset>"), parameters.get(k, "<unset>"))
        for k in keys
        if marker_params.get(k, "<unset>") != parameters.get(k, "<unset>")
    }
    return StageCheckResult("mismatch", changed)


def aggregate_stage_markers(root_dir: str) -> dict:
    """Collect every stage marker found in root_dir's immediate subdirectories.

    Purely a read-only summary for human/report-tab inspection (see
    provenance.write_run_parameters()'s "stages" section) — never itself
    read back by check_stage() to make a skip/rerun decision; that always
    reads the authoritative marker file directly from its own stage
    directory. Walks only one level deep (root_dir's direct children),
    since every stage directory this module writes to is always a direct
    child of a single list's root output directory. Does not require a
    hardcoded list of stage directory names — any subdirectory containing
    a marker file is picked up automatically, so this stays correct if
    new stages are added later without needing a matching update here.

    Args:
        root_dir: A single list's root output directory (i.e. `outdir`
                  in single-list mode, or `outdir/list1` / `outdir/list2`
                  in two-list mode).

    Returns:
        Dict mapping stage name (the marker's own 'stage' field) to
        {'completed_at':, 'parameters':}. Empty dict if root_dir does not
        exist or contains no stage markers at all. A marker file that
        cannot be parsed is skipped with a logged warning rather than
        raising, since this is purely a reporting aggregation and should
        not be able to fail a run.
    """
    result: dict = {}
    if not os.path.isdir(root_dir):
        return result
    for entry in sorted(os.listdir(root_dir)):
        stage_dir = os.path.join(root_dir, entry)
        if not os.path.isdir(stage_dir):
            continue
        marker_path = os.path.join(stage_dir, MARKER_FILENAME)
        if not os.path.isfile(marker_path):
            continue
        try:
            with open(marker_path) as fh:
                marker = json.load(fh)
            result[marker.get("stage", entry)] = {
                "completed_at": marker.get("completed_at"),
                "parameters":   marker.get("parameters", {}),
            }
        except (OSError, ValueError):
            logger.warning(
                "Could not parse stage marker at %s — omitting from "
                "stage summary.", marker_path,
            )
    return result
