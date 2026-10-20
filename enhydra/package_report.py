"""Package a completed ENHYDRA run's report into one shareable, self-contained folder.

Standalone tool, invoked as:

    python -m enhydra.package_report --outdir outdir_A --package-dir shared_report/

Or used programmatically:

    from enhydra.package_report import package_report
    package_report("outdir_A", "shared_report")

ENHYDRA's report.html inlines all of its own CSS/JS dependencies
(jQuery/DataTables/SheetJS) and all of its own summary plots (base64 PNG or
inline SVG) directly into the HTML at build time (see report.py's
_build_report_impl()/_TEMPLATE). The only files report.html actually reads
from disk at view time are:

  - Per-gene-set enrichment plot PNGs, written by GSEApy into a directory
    literally named 'prerank/' under each metric's own results directory
    (e.g. 'enrichment_identity/prerank/GO_0000001.png'), referenced via a
    relative link and loaded into a modal on click.
  - Rendered per-orthogroup alignment pages, written into a directory
    literally named 'alignments/' (e.g. 'alignments/OG000001.html', or
    'list1/alignments/OG000001.html' in two-list mode). Each such page is
    itself fully self-contained (inline CSS only, no external JS or images
    — see msa_viewer.render_alignment_page()), so copying the file alone
    is sufficient.
  - Standalone "View list" filtering-detail pages written by
    report._tsv_to_details_html() as a sibling of the TSV that produced
    them, inside one of a small set of pipeline stats directories
    ('length_filter_stats/', 'group_filter_stats/',
    'divergence_filter_stats/', 'tables/').
  - run_parameters.json at the root of outdir, linked as a raw-JSON
    download from the Run Parameters tab.

Nothing else in a completed run's output directory — raw FASTA files at any
pipeline stage, alignment files, TSV tables, cached GMT/OBO files, GSEA CSV
results, etc. — is read by report.html at view time.

Because every link inside report.html is computed as a path relative to
report.html's own directory (see report.py's `report_dir` usage throughout),
copying only the four categories above into a new directory while preserving
the exact same relative paths is sufficient: every existing link keeps
working unchanged, with no HTML rewriting required.
"""
from __future__ import annotations

import os
import re
import json
import shutil
import logging
import argparse

from .exceptions import EnhydraIOError

logger = logging.getLogger(__name__)

# Directories whose directly-contained HTML files are the standalone
# "View list" detail pages written by report._tsv_to_details_html() — see
# that function and report._generate_drop_details_pages() for how these are
# produced (always a '<stem>.html' sibling of a '<stem>.tsv' it was built
# from). Matched by directory *name* only, not full path, so this is
# agnostic to nesting (works the same whether the directory is
# 'length_filter_stats/' at the root in single-list mode, or
# 'list1/length_filter_stats/' / 'list2/length_filter_stats/' in two-list
# mode).
_DETAILS_PARENT_DIR_NAMES = (
    "length_filter_stats",
    "group_filter_stats",
    "divergence_filter_stats",
    "tables",
)


def _should_include(rel_path: str) -> bool:
    """Decide whether a file (given its path relative to outdir) belongs in
    the shareable report bundle.

    Matched purely on filename + immediate parent directory name, not full
    path depth, so this works identically across single-list, two-list, and
    multi-metric layouts (e.g. 'enrichment_identity/prerank/x.png',
    'differential_zscore/enrichment/prerank/x.png', and
    'list2/alignments/OG0001.html' are all recognised by the same rules
    that recognise 'enrichment/prerank/x.png' and 'alignments/OG0001.html'
    in a plain single-metric single-list run).

    Args:
        rel_path: A file's path relative to the run's outdir, using the
                 platform's own path separator (as returned by
                 os.path.relpath()).

    Returns:
        True if this file should be copied into the shareable bundle.
    """
    parts    = rel_path.replace(os.sep, "/").split("/")
    filename = parts[-1]
    parent   = parts[-2] if len(parts) >= 2 else ""

    if filename == "report.html" and len(parts) == 1:
        return True
    if filename == "run_parameters.json" and len(parts) == 1:
        return True
    if filename.endswith(".png") and parent == "prerank":
        return True
    if filename.endswith(".html") and parent == "alignments":
        return True
    if filename.endswith(".html") and parent in _DETAILS_PARENT_DIR_NAMES:
        return True
    return False


def package_report(
    outdir: str,
    package_dir: str,
    overwrite: bool = False,
) -> str:
    """Copy only the files needed to view a completed run's report into one
    shareable, self-contained folder.

    See this module's own docstring for the full rationale on exactly which
    files are copied and why that set is sufficient for report.html to
    render and function identically (including the alignment-page links,
    enrichment-plot-on-click modals, and filtering-detail "View list"
    links) when the destination folder is shared on its own, without any of
    the original run's intermediate pipeline output.

    Args:
        outdir:      A completed ENHYDRA run's output directory (must
                    contain report.html at its root).
        package_dir: Destination directory for the shareable bundle.
                    Created if it does not already exist.
        overwrite:   If True and package_dir already exists, it is removed
                    and recreated from scratch. If False (default) and
                    package_dir already exists, raises EnhydraIOError
                    rather than merging into or silently overwriting its
                    existing contents.

    Returns:
        The absolute path to the copied report.html inside package_dir.

    Raises:
        EnhydraIOError: If outdir has no report.html at its root, or if
                       package_dir already exists and overwrite=False.
    """
    outdir = os.path.abspath(outdir)
    report_src = os.path.join(outdir, "report.html")
    if not os.path.isfile(report_src):
        raise EnhydraIOError(
            "No report.html found at the root of outdir: %s\n"
            "package_report() expects a completed ENHYDRA run directory "
            "(i.e. the 'outdir' from your project config, after a full "
            "run has finished)." % outdir
        )

    package_dir = os.path.abspath(package_dir)
    if os.path.exists(package_dir):
        if not overwrite:
            raise EnhydraIOError(
                "Destination directory already exists: %s\n"
                "Pass overwrite=True (or --overwrite on the command line) "
                "to replace it, or choose a new destination path."
                % package_dir
            )
        shutil.rmtree(package_dir)
    os.makedirs(package_dir, exist_ok=True)

    n_copied    = 0
    total_bytes = 0
    for root, _dirs, files in os.walk(outdir):
        for filename in files:
            abs_path = os.path.join(root, filename)
            rel_path = os.path.relpath(abs_path, outdir)
            if not _should_include(rel_path):
                continue
            dest_path = os.path.join(package_dir, rel_path)
            os.makedirs(os.path.dirname(dest_path), exist_ok=True)
            shutil.copy2(abs_path, dest_path)
            n_copied    += 1
            total_bytes += os.path.getsize(abs_path)

    logger.info(
        "Packaged shareable report bundle: %d file(s), %.1f MB, written to: %s",
        n_copied, total_bytes / 1_000_000, package_dir,
    )
    return os.path.join(package_dir, "report.html")


# ---------------------------------------------------------------------------
# compare_runs.py comparison reports
# ---------------------------------------------------------------------------

_ENRICHMENT_PLOTS_MAP_PATTERN = re.compile(
    r"var enrichmentPlotsMap\s*=\s*(\{.*?\});", re.DOTALL
)


def _extract_enrichment_plot_paths(report_html_path: str) -> set[str]:
    """Extract every relative PNG path embedded in a compare_runs.py
    comparison report's own enrichmentPlotsMap JS variable.

    compare_runs.py's comparison report (see
    compare_runs.build_comparison_report()) is otherwise fully
    self-contained — every plot other than per-gene-set enrichment plots is
    rendered as embedded SVG/HTML, and all CSS/JS is inlined — except for
    this one JSON blob, which records nested {.../run_name: {go_id:
    relative_path}} mappings to each run's own on-disk enrichment plot
    PNGs (see compare_runs._build_enrichment_plot_index(), reused from
    report.py), referenced by a path relative to the comparison report's
    own directory, using '/' separators regardless of platform (the same
    convention report.py's own _build_enrichment_plot_index() uses).

    Parsed by regex rather than importing compare_runs.py's own data
    structures, since this only needs to work from an already-written
    report file — the caller does not need to re-supply the original
    run_dirs/run_names used to build it.

    Args:
        report_html_path: Path to a comparison report HTML file, as
                          written by compare_runs.build_comparison_report().

    Returns:
        Set of relative path strings (using '/' separators, as embedded).
        Empty set if the variable could not be found or parsed, or
        contained no paths — treated as "nothing to copy", not an error,
        since a comparison run with no ranking metric shared across at
        least 2 runs produces an empty map by design (see
        build_comparison_report()'s own 'usable_metrics' handling).
    """
    with open(report_html_path, encoding="utf-8") as fh:
        html_content = fh.read()

    match = _ENRICHMENT_PLOTS_MAP_PATTERN.search(html_content)
    if not match:
        logger.warning(
            "Could not find enrichmentPlotsMap in report: %s — no "
            "enrichment plot images will be packaged.", report_html_path,
        )
        return set()

    try:
        plots_map = json.loads(match.group(1))
    except json.JSONDecodeError:
        logger.warning(
            "Could not parse enrichmentPlotsMap JSON in report: %s — no "
            "enrichment plot images will be packaged.", report_html_path,
        )
        return set()

    # Nested dict of arbitrary depth (table_id -> run_name -> {go_id: path}
    # in compare_runs.py's case); traversed generically by recursing into
    # dict values until reaching plain strings, so this stays correct
    # regardless of the exact key names or nesting used to produce it.
    paths: set[str] = set()

    def _collect(node):
        if isinstance(node, dict):
            for v in node.values():
                _collect(v)
        elif isinstance(node, str):
            paths.add(node)

    _collect(plots_map)
    return paths


def package_comparison_report(
    report_path: str,
    package_dir: str,
    overwrite: bool = False,
) -> str:
    """Copy a compare_runs.py comparison report and every enrichment plot
    image it links to into one shareable, self-contained folder.

    compare_runs.py's own comparison report is already self-contained for
    everything except per-gene-set enrichment plot images, which it links
    to by relative path on the filesystem rather than embedding (see that
    module's own build_comparison_report() docstring). This function
    copies exactly that remaining dependency, discovered directly from the
    report's own embedded enrichmentPlotsMap JS variable (see
    _extract_enrichment_plot_paths()) rather than requiring the caller to
    re-supply the original run_dirs/run_names used to build it.

    Args:
        report_path: Path to an existing comparison report HTML file, as
                    written by compare_runs.build_comparison_report().
        package_dir: Destination directory for the shareable bundle.
                    Created if it does not already exist.
        overwrite:   If True and package_dir already exists, it is removed
                    and recreated from scratch. If False (default) and
                    package_dir already exists, raises EnhydraIOError
                    rather than merging into or silently overwriting its
                    existing contents.

    Returns:
        The absolute path to the copied report HTML inside package_dir.

    Raises:
        EnhydraIOError: If report_path does not exist, or if package_dir
                       already exists and overwrite=False.
    """
    report_path = os.path.abspath(report_path)
    if not os.path.isfile(report_path):
        raise EnhydraIOError("Comparison report not found: %s" % report_path)
    report_dir = os.path.dirname(report_path)

    package_dir = os.path.abspath(package_dir)
    if os.path.exists(package_dir):
        if not overwrite:
            raise EnhydraIOError(
                "Destination directory already exists: %s\n"
                "Pass overwrite=True (or --overwrite on the command line) "
                "to replace it, or choose a new destination path."
                % package_dir
            )
        shutil.rmtree(package_dir)
    os.makedirs(package_dir, exist_ok=True)

    report_filename = os.path.basename(report_path)
    shutil.copy2(report_path, os.path.join(package_dir, report_filename))

    rel_paths   = _extract_enrichment_plot_paths(report_path)
    n_copied    = 1   # the report itself
    n_missing   = 0
    total_bytes = os.path.getsize(report_path)
    for rel_path in sorted(rel_paths):
        rel_path_local = rel_path.replace("/", os.sep)
        src  = os.path.join(report_dir, rel_path_local)
        dest = os.path.join(package_dir, rel_path_local)
        if not os.path.isfile(src):
            logger.warning(
                "Enrichment plot referenced by report not found on disk "
                "(skipping): %s", src,
            )
            n_missing += 1
            continue
        os.makedirs(os.path.dirname(dest), exist_ok=True)
        shutil.copy2(src, dest)
        n_copied    += 1
        total_bytes += os.path.getsize(src)

    logger.info(
        "Packaged shareable comparison report bundle: %d file(s) (%d "
        "referenced plot(s) missing on disk), %.1f MB, written to: %s",
        n_copied, n_missing, total_bytes / 1_000_000, package_dir,
    )
    return os.path.join(package_dir, report_filename)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m enhydra.package_report",
        description="Package a completed ENHYDRA run's report into one "
                    "shareable, self-contained folder (no intermediate "
                    "pipeline files included).",
    )
    parser.add_argument(
        "--outdir", required=True, metavar="DIR",
        help="A completed ENHYDRA run's output directory "
             "(must contain report.html at its root).",
    )
    parser.add_argument(
        "--package-dir", required=True, metavar="DIR",
        help="Destination directory for the shareable bundle. "
             "Created if it does not already exist.",
    )
    parser.add_argument(
        "--overwrite", action="store_true", default=False,
        help="Replace --package-dir if it already exists, rather than "
             "raising an error.",
    )
    return parser


def main():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )
    parser = _build_arg_parser()
    args   = parser.parse_args()

    try:
        report_path = package_report(
            outdir=args.outdir,
            package_dir=args.package_dir,
            overwrite=args.overwrite,
        )
    except EnhydraIOError as e:
        raise SystemExit("Error: %s" % e)

    print("Shareable report bundle ready: %s" % report_path)


if __name__ == "__main__":
    main()
