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
