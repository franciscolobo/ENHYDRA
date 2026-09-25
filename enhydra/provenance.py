from __future__ import annotations

import os
import json
import logging
import datetime

logger = logging.getLogger(__name__)

_RUN_PARAMETERS_FILENAME = "run_parameters.json"


def collect_run_parameters(
    outdir: str,
    code_config_path: str,
    project_config_path: str,
    command_line: list[str],
    inputdir: str,
    anchor: str,
    resolved_params: dict,
    two_list_mode: bool = False,
    orthofinder_dir: str | None = None,
    list1_path: str | None = None,
    list2_path: str | None = None,
    list1_name: str | None = None,
    list2_name: str | None = None,
    species1: list[str] | None = None,
    species2: list[str] | None = None,
    anchor_injected: bool | None = None,
    all_metrics: bool = False,
    metrics_run: tuple[str, ...] = (),
) -> dict:
    """Assemble a reproducibility record of the effective parameters for one run.

    This captures parameters *after* CLI-override resolution (i.e. the values
    cli.py actually used, via its _resolve() helper) rather than the raw
    config file contents — a run invoked with e.g. --min-species overriding
    the project config should record the value that was actually used, not
    the config file's original value, since those can differ silently.

    Software/tool version detection (aligner, trimAl, key Python packages)
    is deliberately out of scope for this function for now; only configured
    paths and resolved parameter values are recorded. This may be extended
    later without changing this function's existing keys.

    Args:
        outdir:              This run's output directory.
        code_config_path:    Path to the code configuration file used.
        project_config_path: Path to the project configuration file used.
        command_line:        The full command line invoked, as a list of
                             strings (e.g. sys.argv). Stored as given, not
                             re-parsed or validated.
        inputdir:            Path to the input FASTA directory.
        anchor:              Anchor species ID.
        resolved_params:     Dict of effective (post-resolution) parameter
                             values. Expected keys — any missing key is
                             recorded as null rather than raising, since a
                             caller may reasonably omit a value that isn't
                             relevant to a particular run:
                               min_species, min_sequences, paralogs,
                               length_filter_sd, divergence_filter_sd, trim,
                               aligner, mafft_mode, metric, permutations,
                               min_size, max_size, seed, fdr_threshold,
                               top_n, gene_sets, organism, sources.
        two_list_mode:       Whether this was a two-list differential run.
        orthofinder_dir:     Path to an OrthoFinder output directory, if
                             --orthofinder-dir was used. None otherwise.
        list1_path:          Path to list1's species file (two-list mode).
        list2_path:          Path to list2's species file (two-list mode).
        list1_name:          Display name for list 1 (two-list mode).
        list2_name:          Display name for list 2 (two-list mode).
        species1:            List of species IDs in list 1, after any
                             anchor injection (two-list mode). Only the
                             count is recorded, not the full list, to keep
                             this file small — the original species list
                             files remain the source of truth and are
                             already recorded via list1_path/list2_path.
        species2:            List of species IDs in list 2 (two-list mode).
        anchor_injected:     Whether the anchor species was automatically
                             added to list 1 because it was absent from the
                             user-supplied list (two-list mode only).
        all_metrics:         Whether --all-metrics was set.
        metrics_run:         The actual metric(s) run, e.g. ('identity',
                             'zscore', 'rank') under --all-metrics, or a
                             single-element tuple otherwise.

    Returns:
        A JSON-serialisable dict. Also see write_run_parameters() to persist
        this to disk.
    """
    def _get(key):
        return resolved_params.get(key)

    record: dict = {
        "enhydra_version": _enhydra_version(),
        "timestamp":       datetime.datetime.now().isoformat(timespec="seconds"),
        "command_line":    list(command_line),
        "code_config_path":    os.path.abspath(code_config_path),
        "project_config_path": os.path.abspath(project_config_path),
        "outdir":              os.path.abspath(outdir),
        "input": {
            "inputdir":        os.path.abspath(inputdir) if inputdir else None,
            "orthofinder_dir": os.path.abspath(orthofinder_dir)
                               if orthofinder_dir else None,
            "anchor":          anchor,
            "two_list_mode":   two_list_mode,
        },
        "filtering": {
            "min_species":          _get("min_species"),
            "min_sequences":        _get("min_sequences"),
            "paralogs":             _get("paralogs"),
            "length_filter_sd":     _get("length_filter_sd"),
            "divergence_filter_sd": _get("divergence_filter_sd"),
            "trim":                 _get("trim") or None,
        },
        "alignment": {
            "aligner":    _get("aligner"),
            "mafft_mode": _get("mafft_mode"),
        },
        "ranking": {
            "all_metrics": all_metrics,
            "metrics_run": list(metrics_run),
        },
        "gsea": {
            "permutations":  _get("permutations"),
            "min_size":      _get("min_size"),
            "max_size":      _get("max_size"),
            "seed":          _get("seed"),
            "fdr_threshold": _get("fdr_threshold"),
            "top_n":         _get("top_n"),
        },
        "gene_sets": {
            "gene_sets_path": _get("gene_sets") or None,
            "organism":       _get("organism") or None,
            "sources":        _get("sources"),
        },
    }

    if two_list_mode:
        record["input"]["list1"] = {
            "path":            list1_path,
            "name":            list1_name,
            "n_species":       len(species1) if species1 is not None else None,
            "anchor_injected": anchor_injected,
        }
        record["input"]["list2"] = {
            "path":      list2_path,
            "name":      list2_name,
            "n_species": len(species2) if species2 is not None else None,
        }

    return record


def write_run_parameters(
    params: dict,
    outdir: str,
    output_path: str | None = None,
) -> str:
    """Write a run-parameters record to disk as JSON.

    Mirrors stats.aggregate_pipeline_stats()'s own write pattern: build the
    dict, write it once as pretty-printed JSON, return the path written.

    Args:
        params:      The dict to write, as returned by
                    collect_run_parameters().
        outdir:      This run's output directory. Used to derive the
                    default output_path if none is given.
        output_path: Where to write the JSON file. Defaults to
                    '<outdir>/run_parameters.json'.

    Returns:
        The path the file was written to.
    """
    if output_path is None:
        output_path = os.path.join(outdir, _RUN_PARAMETERS_FILENAME)
    with open(output_path, "w") as fh:
        json.dump(params, fh, indent=2)
    logger.info("Run parameters written to: %s", output_path)
    return output_path


def _enhydra_version() -> str:
    """Best-effort lookup of the installed enhydra package version.

    Falls back to 'unknown' rather than raising — a missing or unresolvable
    version should never prevent a run from completing or its provenance
    record from being written.
    """
    try:
        from . import __version__
        return __version__
    except Exception:
        return "unknown"
