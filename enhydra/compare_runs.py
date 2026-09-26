"""Compare two or more completed ENHYDRA run output directories.

Standalone tool, invoked as:

    python -m enhydra.compare_runs --runs outdir_A outdir_B [outdir_C ...] \\
        --names "min_species=4" "min_species=8" \\
        --output comparison_report.html

Produces a single self-contained HTML report (no companion assets
directory) with four tabs:

    Run parameters          — an N-column diff of each run's
                              run_parameters.json, with differing rows
                              flagged.
    Significance agreement  — for every ranking metric present in ALL
                              selected runs, a majority-agreement table
                              (reusing report.build_significance_agreement_table_html(),
                              with "columns" = runs instead of metrics).
    Rank correlation        — Spearman correlation of GO-term NES between
                              runs, per common metric: a direct scatter
                              for exactly 2 runs, or a correlation matrix
                              plus collapsible pairwise scatters for more.
    Significance transitions — a pairwise Sankey-style diagram per run
                              pair per common metric, showing gene sets
                              gaining, losing, or retaining significance.

Two runs must share at least one ranking metric (identity/zscore/rank) to
be compared at all beyond the Run parameters tab — comparing, say, a
zscore-ranked run against an identity-ranked run would conflate "the runs
differ" with "the ranking metric differs," so this tool deliberately does
not attempt it.

GO term identity (not group/gene identity) is the basis for every
comparison here, since it is the one identifier that remains meaningful
across runs regardless of anchor species, min_species, or other filtering
differences — unlike orthogroup IDs, which are only comparable if two
runs share identical input data.
"""
from __future__ import annotations

import os
import json
import html
import logging
import argparse
from itertools import combinations

import pandas as pd

from .io import parse_obo_names
from .plotting import _rdbu_color
from .report import (
    METRIC_LABELS,
    format_provenance_value,
    get_run_parameter_field_specs,
    build_significance_agreement_table_html,
    _fetch_cached,
    _JQUERY_URL,
    _DATATABLES_JS_URL,
    _DATATABLES_CSS_URL,
    _build_enrichment_plot_index,
    _gmt_term_names,
    _find_gmt_in_dir,
)

logger = logging.getLogger(__name__)

_ALL_METRICS = ("identity", "zscore", "rank")


# ---------------------------------------------------------------------------
# Run discovery
# ---------------------------------------------------------------------------

def _discover_run(outdir: str) -> dict:
    """Load a completed ENHYDRA run's key artifact locations.

    Prefers run_parameters.json (written by provenance.write_run_parameters()
    since the run-parameters-tab feature) to determine the run's mode
    (single-list vs. two-list) and which metrics were run, since that file
    records this unambiguously rather than requiring directory-layout
    guesses. Falls back to sniffing the directory layout (presence of
    list1/, and of <name>_<metric>-suffixed enrichment directories) for
    runs that predate that file, with a logged warning — that fallback
    cannot distinguish *which* metric an un-suffixed 'enrichment/' or
    'differential/' directory was run with (run_parameters.json is the
    only source of that information), so such a run's single metric is
    recorded under the key 'unknown' rather than guessed. Two 'unknown'
    metric runs will still compare against each other (both keyed the
    same), but an 'unknown' run will never match a named-metric run from
    another run's provenance record — a known, accepted limitation of
    comparing against pre-provenance run directories.

    Args:
        outdir: Path to a completed ENHYDRA run's output directory.

    Returns:
        Dict with keys:
            outdir:          Absolute path to outdir.
            run_parameters:  The parsed run_parameters.json dict, or None
                             if not found.
            two_list_mode:   bool.
            all_metrics:     bool.
            metric_paths:    {metric: {"results_dir":, "tables_dir1":,
                             "tables_dir2": or None}} for every metric with
                             a usable GSEA results CSV found on disk.
    """
    outdir = os.path.abspath(outdir)
    run_params_path = os.path.join(outdir, "run_parameters.json")
    run_params = None
    if os.path.isfile(run_params_path):
        with open(run_params_path) as fh:
            run_params = json.load(fh)

    if run_params is not None:
        two_list    = bool((run_params.get("input") or {}).get("two_list_mode"))
        all_metrics = bool((run_params.get("ranking") or {}).get("all_metrics"))
        metrics_run = (run_params.get("ranking") or {}).get("metrics_run") or []
    else:
        logger.warning(
            "run_parameters.json not found in %s — falling back to "
            "directory sniffing. The specific ranking metric of a "
            "non-all-metrics run cannot be determined this way and will "
            "be recorded as 'unknown'.", outdir,
        )
        two_list = os.path.isdir(os.path.join(outdir, "list1"))
        base_prefix = "differential" if two_list else "enrichment"
        suffixed = [
            m for m in _ALL_METRICS
            if os.path.isdir(os.path.join(outdir, "%s_%s" % (base_prefix, m)))
        ]
        all_metrics = bool(suffixed)
        metrics_run = suffixed if suffixed else ["unknown"]

    metric_paths: dict[str, dict] = {}
    for m in metrics_run:
        sfx = ("_%s" % m) if (all_metrics and m != "unknown") else ""
        if two_list:
            diff_dir    = os.path.join(outdir, "differential%s" % sfx)
            results_dir = os.path.join(diff_dir, "enrichment")
            tables_dir1 = os.path.join(outdir, "list1", "tables")
            tables_dir2 = os.path.join(outdir, "list2", "tables")
        else:
            results_dir = os.path.join(outdir, "enrichment%s" % sfx)
            tables_dir1 = os.path.join(outdir, "tables")
            tables_dir2 = None

        csv_path = os.path.join(results_dir, "gseapy.gene_set.prerank.report.csv")
        if not os.path.isfile(csv_path):
            logger.warning(
                "No GSEA results found for metric '%s' in %s (expected: "
                "%s) — this metric will be unavailable for comparison.",
                m, outdir, csv_path,
            )
            continue

        metric_paths[m] = {
            "results_dir": results_dir,
            "tables_dir1": tables_dir1,
            "tables_dir2": tables_dir2,
        }

    return {
        "outdir":         outdir,
        "run_parameters": run_params,
        "two_list_mode":  two_list,
        "all_metrics":    all_metrics,
        "metric_paths":   metric_paths,
    }


def _load_run_metric_df(path_info: dict) -> pd.DataFrame | None:
    """Load a run's GSEA results for one metric as a {Term, NES, FDR q-val} frame."""
    csv_path = os.path.join(
        path_info["results_dir"], "gseapy.gene_set.prerank.report.csv"
    )
    if not os.path.isfile(csv_path):
        return None
    df = pd.read_csv(csv_path)
    required = {"Term", "NES", "FDR q-val"}
    if not required.issubset(df.columns):
        logger.warning(
            "GSEA results at %s are missing expected column(s): %s",
            csv_path, required - set(df.columns),
        )
        return None
    return df[["Term", "NES", "FDR q-val"]].copy()


# ---------------------------------------------------------------------------
# Pairwise comparison math
# ---------------------------------------------------------------------------

def _merge_pair(df_a: pd.DataFrame, df_b: pd.DataFrame) -> pd.DataFrame:
    """Inner-join two runs' {Term, NES, FDR q-val} frames on Term.

    Returns columns: Term, NES_a, FDR_a, NES_b, FDR_b.
    """
    merged = df_a.merge(df_b, on="Term", suffixes=("_a", "_b"), how="inner")
    merged = merged.rename(columns={
        "FDR q-val_a": "FDR_a", "FDR q-val_b": "FDR_b",
    })
    return merged


def _sankey_counts(merged: pd.DataFrame, fdr_threshold: float) -> dict[tuple[bool, bool], int]:
    """Count common gene sets by (significant_in_a, significant_in_b).

    Args:
        merged: Output of _merge_pair() — must have FDR_a/FDR_b columns.
        fdr_threshold: Significance cutoff.

    Returns:
        Dict keyed by (bool, bool) with all four combinations present
        (zero-valued if empty), so callers never need a .get() default.
    """
    sig_a = merged["FDR_a"] < fdr_threshold
    sig_b = merged["FDR_b"] < fdr_threshold
    return {
        (True, True):   int((sig_a & sig_b).sum()),
        (True, False):  int((sig_a & ~sig_b).sum()),
        (False, True):  int((~sig_a & sig_b).sum()),
        (False, False): int((~sig_a & ~sig_b).sum()),
    }


def _spearman_matrix(dfs: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """Pairwise Spearman correlation of NES over each pair's common GO terms.

    Args:
        dfs: {run_name: DataFrame} with 'Term' and 'NES' columns, for one
             fixed ranking metric.

    Returns:
        A square DataFrame indexed and columned by run name. A cell is
        NaN if the corresponding pair of runs has fewer than 2 GO terms
        in common (Spearman correlation is undefined below that).
    """
    names  = list(dfs.keys())
    series = {n: dfs[n].set_index("Term")["NES"] for n in names}
    mat = pd.DataFrame(index=names, columns=names, dtype=float)
    for a in names:
        for b in names:
            common = series[a].index.intersection(series[b].index)
            if len(common) < 2:
                mat.loc[a, b] = float("nan")
                continue
            mat.loc[a, b] = series[a].loc[common].corr(
                series[b].loc[common], method="spearman"
            )
    return mat


# ---------------------------------------------------------------------------
# SVG rendering (hand-rolled, matching plotting.py's existing style —
# element tree building, no external charting dependency)
# ---------------------------------------------------------------------------

def _render_correlation_scatter_svg(
    name_a: str,
    name_b: str,
    merged: pd.DataFrame,
    fdr_threshold: float,
    spearman_r: float,
) -> str:
    """Render an interactive NES-vs-NES scatter for one pair of runs.

    Points are coloured by significance-transition category (stable
    significant / gained / lost / stable not-significant) rather than a
    continuous scale, since that is the information most relevant to a
    robustness check between two runs — the same categorisation used by
    _render_sankey_svg() for the same pair, so the two views agree.

    Args:
        name_a, name_b: Display names for the two runs.
        merged:         Output of _merge_pair() (Term, NES_a, FDR_a,
                        NES_b, FDR_b). Must have at least 2 rows.
        fdr_threshold:  Significance cutoff.
        spearman_r:     Precomputed Spearman correlation, shown in the title.

    Returns:
        A standalone '<svg ...>...</svg>' markup string for direct
        embedding (this tool produces one self-contained HTML file with
        no companion plot files).
    """
    import xml.etree.ElementTree as ET

    ml, mr, mt, mb = 60, 20, 40, 55
    pw, ph = 380, 380
    sw, sh = ml + pw + mr, mt + ph + mb

    nes_a = merged["NES_a"].astype(float)
    nes_b = merged["NES_b"].astype(float)
    data_min = float(min(nes_a.min(), nes_b.min()))
    data_max = float(max(nes_a.max(), nes_b.max()))
    span = data_max - data_min
    pad  = span * 0.06 if span > 0 else 0.5
    lo, hi = data_min - pad, data_max + pad
    rng = hi - lo or 1.0

    def spx(v): return ml + (v - lo) / rng * pw
    def spy(v): return mt + ph - (v - lo) / rng * ph

    svg = ET.Element("svg", {
        "xmlns": "http://www.w3.org/2000/svg",
        "width": str(sw), "height": str(sh),
        "font-family": "Arial, sans-serif",
    })
    ET.SubElement(svg, "text", {
        "x": str(sw // 2), "y": "20", "text-anchor": "middle",
        "font-size": "12", "font-weight": "bold", "fill": "#1a3a5c",
    }).text = "%s vs %s (Spearman r = %.3f)" % (name_a, name_b, spearman_r)

    ax_y = mt + ph
    ET.SubElement(svg, "line", {
        "x1": str(ml), "y1": str(mt), "x2": str(ml), "y2": str(ax_y),
        "stroke": "#555", "stroke-width": "1"})
    ET.SubElement(svg, "line", {
        "x1": str(ml), "y1": str(ax_y), "x2": str(ml + pw), "y2": str(ax_y),
        "stroke": "#555", "stroke-width": "1"})
    ET.SubElement(svg, "line", {
        "x1": str(spx(lo)), "y1": str(spy(lo)),
        "x2": str(spx(hi)), "y2": str(spy(hi)),
        "stroke": "#555", "stroke-width": "0.8",
        "stroke-dasharray": "4,3", "opacity": "0.5"})

    for v in [lo + rng * i / 4 for i in range(5)]:
        ET.SubElement(svg, "line", {
            "x1": str(spx(v)), "y1": str(ax_y),
            "x2": str(spx(v)), "y2": str(ax_y + 4),
            "stroke": "#555", "stroke-width": "1"})
        ET.SubElement(svg, "text", {
            "x": str(spx(v)), "y": str(ax_y + 14),
            "text-anchor": "middle", "font-size": "9", "fill": "#555",
        }).text = "%.2f" % v
        ET.SubElement(svg, "line", {
            "x1": str(ml - 4), "y1": str(spy(v)),
            "x2": str(ml), "y2": str(spy(v)),
            "stroke": "#555", "stroke-width": "1"})
        ET.SubElement(svg, "text", {
            "x": str(ml - 7), "y": str(spy(v) + 4),
            "text-anchor": "end", "font-size": "9", "fill": "#555",
        }).text = "%.2f" % v

    ET.SubElement(svg, "text", {
        "x": str(ml + pw // 2), "y": str(sh - 5),
        "text-anchor": "middle", "font-size": "11", "fill": "#333",
    }).text = "NES (%s)" % name_a
    ET.SubElement(svg, "text", {
        "x": "0", "y": "0",
        "transform": "translate(13,%d) rotate(-90)" % (mt + ph // 2),
        "text-anchor": "middle", "font-size": "11", "fill": "#333",
    }).text = "NES (%s)" % name_b

    for _, row in merged.iterrows():
        sig_a, sig_b = row["FDR_a"] < fdr_threshold, row["FDR_b"] < fdr_threshold
        if sig_a and sig_b:
            color, status = "#2166ac", "significant in both"
        elif sig_a and not sig_b:
            color, status = "#d6604d", "lost significance in %s" % name_b
        elif sig_b and not sig_a:
            color, status = "#2f7d3c", "gained significance in %s" % name_b
        else:
            color, status = "#cccccc", "not significant in either"
        tip = "%s | NES %s: %.3f | NES %s: %.3f | %s" % (
            row["Term"], name_a, row["NES_a"], name_b, row["NES_b"], status,
        )
        ET.SubElement(svg, "circle", {
            "cx": "%.1f" % spx(float(row["NES_a"])),
            "cy": "%.1f" % spy(float(row["NES_b"])),
            "r": "4", "fill": color, "opacity": "0.7",
            "data-tip": html.escape(tip), "style": "cursor:pointer;",
        })

    return ET.tostring(svg, encoding="unicode")


def _render_correlation_matrix_html(mat: pd.DataFrame) -> str:
    """Render a Spearman correlation matrix as a colour-shaded HTML table."""
    names = list(mat.index)
    header = "<th></th>" + "".join(
        "<th>%s</th>" % html.escape(str(n)) for n in names
    )
    rows = []
    for a in names:
        cells = []
        for b in names:
            v = mat.loc[a, b]
            if pd.isna(v):
                cells.append('<td style="text-align:center;">&mdash;</td>')
            else:
                color = _rdbu_color(float(v))
                cells.append(
                    '<td style="background-color:%s;text-align:center;">%.3f</td>'
                    % (color, v)
                )
        rows.append(
            "<tr><td><strong>%s</strong></td>%s</tr>"
            % (html.escape(str(a)), "".join(cells))
        )
    return (
        '<table class="reason-table"><thead><tr>%s</tr></thead>'
        "<tbody>%s</tbody></table>" % (header, "".join(rows))
    )


def _pairwise_correlation_block(
    metric_label: str,
    name_a: str,
    name_b: str,
    df_a: pd.DataFrame,
    df_b: pd.DataFrame,
    fdr_threshold: float,
) -> str:
    """Build one pair's scatter + caption, or an explanatory message if
    there are too few common GO terms to correlate."""
    merged = _merge_pair(df_a, df_b)
    only_a = set(df_a["Term"]) - set(df_b["Term"])
    only_b = set(df_b["Term"]) - set(df_a["Term"])
    if len(merged) < 2:
        return (
            "<p>Fewer than 2 gene sets in common between %s and %s for "
            "%s \u2014 cannot compute a correlation.</p>"
            % (html.escape(name_a), html.escape(name_b), html.escape(metric_label))
        )
    r = float(merged["NES_a"].corr(merged["NES_b"], method="spearman"))
    svg = _render_correlation_scatter_svg(name_a, name_b, merged, fdr_threshold, r)
    caption = (
        "%d gene sets tested in both runs (Spearman r = %.3f). %d tested "
        "only in %s; %d tested only in %s (excluded here)."
        % (len(merged), r, len(only_a), name_a, len(only_b), name_b)
    )
    return '<div class="plot-block">%s<p class="plot-caption">%s</p></div>' % (svg, caption)


def _render_sankey_svg(
    name_a: str,
    name_b: str,
    metric_label: str,
    counts: dict[tuple[bool, bool], int],
) -> str:
    """Render a pairwise significance-transition Sankey diagram.

    Two nodes per side (Significant / Not significant, for run A on the
    left and run B on the right), with four flows between them:
    stable-significant, lost, gained, stable-not-significant.

    Args:
        name_a, name_b: Display names for the two runs.
        metric_label:   Display label for the ranking metric being compared.
        counts:         Output of _sankey_counts().

    Returns:
        A standalone '<svg ...>...</svg>' markup string, or a plain <p>
        message if there is nothing to draw (all four counts are zero).
    """
    import xml.etree.ElementTree as ET

    n_tt = counts[(True, True)]
    n_tf = counts[(True, False)]
    n_ft = counts[(False, True)]
    n_ff = counts[(False, False)]
    grand_total = n_tt + n_tf + n_ft + n_ff
    if grand_total == 0:
        return "<p>No common gene sets to compare between these two runs for this metric.</p>"

    total_a_sig, total_a_notsig = n_tt + n_tf, n_ft + n_ff
    total_b_sig, total_b_notsig = n_tt + n_ft, n_tf + n_ff

    plot_h = 260.0
    gap    = 6.0
    px_per_unit = plot_h / grand_total

    def seg_h(n: int) -> float:
        return max(n * px_per_unit, 1.0) if n > 0 else 0.0

    ml, mr, mt = 150, 150, 50
    node_w  = 18
    right_x = ml + 140
    sw = ml + 140 + mr
    sh = mt + plot_h + gap + 60

    a_sig_h, a_notsig_h = seg_h(total_a_sig), seg_h(total_a_notsig)
    a_sig_y0 = mt
    a_sig_y1 = a_sig_y0 + a_sig_h
    a_notsig_y0 = a_sig_y1 + gap
    a_notsig_y1 = a_notsig_y0 + a_notsig_h

    b_sig_h, b_notsig_h = seg_h(total_b_sig), seg_h(total_b_notsig)
    b_sig_y0 = mt
    b_sig_y1 = b_sig_y0 + b_sig_h
    b_notsig_y0 = b_sig_y1 + gap
    b_notsig_y1 = b_notsig_y0 + b_notsig_h

    svg = ET.Element("svg", {
        "xmlns": "http://www.w3.org/2000/svg",
        "width": str(sw), "height": str(sh),
        "font-family": "Arial, sans-serif",
    })
    ET.SubElement(svg, "text", {
        "x": str(sw // 2), "y": "20", "text-anchor": "middle",
        "font-size": "12", "font-weight": "bold", "fill": "#1a3a5c",
    }).text = "%s: %s \u2192 %s" % (metric_label, name_a, name_b)

    def node_rect(x, y0, y1, label, count, is_left):
        if y1 <= y0:
            return
        ET.SubElement(svg, "rect", {
            "x": str(x), "y": "%.1f" % y0, "width": str(node_w),
            "height": "%.1f" % (y1 - y0), "fill": "#1a3a5c",
        })
        anchor = "end" if is_left else "start"
        tx = x - 6 if is_left else x + node_w + 6
        ET.SubElement(svg, "text", {
            "x": str(tx), "y": "%.1f" % ((y0 + y1) / 2 + 4),
            "text-anchor": anchor, "font-size": "10", "fill": "#333",
        }).text = "%s (%d)" % (label, count)

    node_rect(ml, a_sig_y0, a_sig_y1, "Significant in %s" % name_a, total_a_sig, True)
    node_rect(ml, a_notsig_y0, a_notsig_y1, "Not significant in %s" % name_a, total_a_notsig, True)
    node_rect(right_x, b_sig_y0, b_sig_y1, "Significant in %s" % name_b, total_b_sig, False)
    node_rect(right_x, b_notsig_y0, b_notsig_y1, "Not significant in %s" % name_b, total_b_notsig, False)

    def ribbon(l_y0, l_offset, r_y0, r_offset, n, color, label):
        if n <= 0:
            return l_offset, r_offset
        seg = seg_h(n)
        ly0, ly1 = l_y0 + l_offset, l_y0 + l_offset + seg
        ry0, ry1 = r_y0 + r_offset, r_y0 + r_offset + seg
        x1, x2 = ml + node_w, right_x
        xm = (x1 + x2) / 2
        path_d = (
            "M %.1f %.1f C %.1f %.1f %.1f %.1f %.1f %.1f "
            "L %.1f %.1f C %.1f %.1f %.1f %.1f %.1f %.1f Z"
            % (x1, ly0, xm, ly0, xm, ry0, x2, ry0,
               x2, ry1, xm, ry1, xm, ly1, x1, ly1)
        )
        tip = "%s: %d gene set%s" % (label, n, "" if n == 1 else "s")
        ET.SubElement(svg, "path", {
            "d": path_d, "fill": color, "opacity": "0.55",
            "data-tip": html.escape(tip), "style": "cursor:pointer;",
        })
        return l_offset + seg, r_offset + seg

    a_sig_off = a_notsig_off = b_sig_off = b_notsig_off = 0.0
    a_sig_off, b_sig_off = ribbon(
        a_sig_y0, a_sig_off, b_sig_y0, b_sig_off, n_tt,
        "#2166ac", "Stable: significant in both",
    )
    a_sig_off, b_notsig_off = ribbon(
        a_sig_y0, a_sig_off, b_notsig_y0, b_notsig_off, n_tf,
        "#d6604d", "Lost significance (%s \u2192 %s)" % (name_a, name_b),
    )
    a_notsig_off, b_sig_off = ribbon(
        a_notsig_y0, a_notsig_off, b_sig_y0, b_sig_off, n_ft,
        "#2f7d3c", "Gained significance (%s \u2192 %s)" % (name_a, name_b),
    )
    a_notsig_off, b_notsig_off = ribbon(
        a_notsig_y0, a_notsig_off, b_notsig_y0, b_notsig_off, n_ff,
        "#cccccc", "Stable: not significant in either",
    )

    return ET.tostring(svg, encoding="unicode")


# ---------------------------------------------------------------------------
# Run parameters diff table
# ---------------------------------------------------------------------------

def _diff_rows(
    field_specs: list[tuple[str, str]],
    section_getter,
    records: list[dict | None],
    names: list[str],
) -> str:
    header_cells = "<th>Field</th>" + "".join(
        "<th>%s</th>" % html.escape(n) for n in names
    )
    rows = []
    for key, label in field_specs:
        values = [
            (section_getter(r).get(key) if r is not None else None)
            for r in records
        ]
        differs = len({repr(v) for v in values}) > 1
        row_class = ' class="diff-row"' if differs else ""
        cells = "".join(
            "<td>%s</td>" % format_provenance_value(v) for v in values
        )
        rows.append(
            "<tr%s><td>%s</td>%s</tr>" % (row_class, html.escape(label), cells)
        )
    return (
        '<table class="reason-table"><thead><tr>%s</tr></thead>'
        "<tbody>%s</tbody></table>" % (header_cells, "".join(rows))
    )


def _build_run_params_diff_html(
    records: list[dict | None],
    names: list[str],
) -> str:
    """Build the Run parameters tab: one N-column table per section,
    with rows where values differ across runs flagged (see the
    '.diff-row' CSS class in _COMPARE_TEMPLATE).

    A run whose run_parameters.json could not be found (record is None)
    shows 'not set' in every cell for that run, rather than raising.
    """
    specs = get_run_parameter_field_specs()
    any_two_list = any(
        r and ((r.get("input") or {}).get("two_list_mode")) for r in records
    )

    sections = [
        ("run_info",  "Run info",  lambda r: r or {}),
        ("input",     "Input",     lambda r: (r or {}).get("input") or {}),
        ("filtering", "Filtering", lambda r: (r or {}).get("filtering") or {}),
        ("alignment", "Alignment", lambda r: (r or {}).get("alignment") or {}),
        ("ranking",   "Ranking",   lambda r: (r or {}).get("ranking") or {}),
        ("gsea",      "GSEA",      lambda r: (r or {}).get("gsea") or {}),
        ("gene_sets", "Gene sets", lambda r: (r or {}).get("gene_sets") or {}),
    ]

    parts = []
    for section_key, title, getter in sections:
        parts.append("<h4>%s</h4>" % html.escape(title))
        parts.append(_diff_rows(specs[section_key], getter, records, names))

    if any_two_list:
        parts.append("<h4>Two-list mode \u2014 List 1</h4>")
        parts.append(_diff_rows(
            specs["list"],
            lambda r: ((r or {}).get("input") or {}).get("list1") or {},
            records, names,
        ))
        parts.append("<h4>Two-list mode \u2014 List 2</h4>")
        parts.append(_diff_rows(
            specs["list"],
            lambda r: ((r or {}).get("input") or {}).get("list2") or {},
            records, names,
        ))

    return "".join(parts)


def _build_run_summary_html(runs: dict[str, dict]) -> str:
    rows = []
    for name, info in runs.items():
        mode = "Two-list differential" if info["two_list_mode"] else "Single-list"
        metrics = ", ".join(sorted(info["metric_paths"].keys())) or "(none found)"
        note = ""
        if info["run_parameters"] is None:
            note = (
                ' <span style="color:#d9822b;">(run_parameters.json not '
                "found \u2014 mode/metrics inferred from directory layout)"
                "</span>"
            )
        rows.append(
            "<tr><td>%s</td><td>%s</td><td>%s</td><td>%s%s</td></tr>"
            % (html.escape(name), html.escape(info["outdir"]), mode,
               html.escape(metrics), note)
        )
    return (
        '<table class="reason-table"><thead><tr><th>Run</th>'
        "<th>Directory</th><th>Mode</th><th>Metrics found</th></tr></thead>"
        "<tbody>%s</tbody></table>" % "".join(rows)
    )


# ---------------------------------------------------------------------------
# HTML template
# ---------------------------------------------------------------------------

_COMPARE_TEMPLATE = """\
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8"/>
<meta name="viewport" content="width=device-width, initial-scale=1.0"/>
<title>{title}</title>
<style>
{dt_css}
body {{ font-family: Arial, sans-serif; margin: 0; padding: 0;
        background: #f5f5f5; color: #222; }}
header {{ background: #1a3a5c; color: white; padding: 24px 40px; }}
header h1 {{ margin: 0; font-size: 1.8em; }}
header p  {{ margin: 4px 0 0; font-size: 0.95em; opacity: 0.85; }}
main {{ max-width: 1300px; margin: 32px auto; padding: 0 24px; }}
section {{ background: white; border-radius: 8px;
           box-shadow: 0 1px 4px rgba(0,0,0,0.1);
           padding: 28px 32px; margin-bottom: 28px; }}
h2 {{ margin-top: 0; color: #1a3a5c; border-bottom: 2px solid #e0e0e0;
      padding-bottom: 8px; }}
h4 {{ color: #1a3a5c; margin: 22px 0 8px; }}
.tab-nav {{ display: flex; gap: 0; border-bottom: 3px solid #1a3a5c;
            margin-bottom: 28px; flex-wrap: wrap; }}
.tab-btn {{ padding: 11px 32px; border: none; border-radius: 6px 6px 0 0;
            background: #e2eaf3; cursor: pointer; font-size: 14px;
            font-weight: 600; color: #555; margin-right: 3px; }}
.tab-btn:hover:not(.active) {{ background: #c8d8ec; color: #1a3a5c; }}
.tab-btn.active {{ background: #1a3a5c; color: white; }}
.tab-panel {{ display: none; }}
.tab-panel.active {{ display: block; }}
.metric-desc {{ font-size: 0.9em; color: #444; margin: 0 0 20px;
                padding: 10px 14px; background: #f0f5fa;
                border-left: 3px solid #1a3a5c; border-radius: 0 4px 4px 0; }}
table.reason-table {{ width: 100%; border-collapse: collapse; margin: 6px 0 20px; }}
table.reason-table th, table.reason-table td {{ border: 1px solid #e0e0e0;
    padding: 5px 10px; font-size: 12.5px; text-align: left; }}
table.reason-table th {{ background: #eef2f6; }}
tr.diff-row {{ background-color: #fff4e5 !important; }}
tr.sig-row {{ background-color: #eaf3fb !important; font-weight: bold; }}
a.go-link {{ color: #1a3a5c; text-decoration: underline dotted; cursor: pointer; }}
.col-tip {{ display: inline-block; width: 14px; height: 14px; line-height: 14px;
            font-size: 10px; text-align: center; border-radius: 50%;
            background: #aaa; color: white; cursor: help; margin-left: 3px;
            position: relative; }}
.col-tip .tip-text {{ display: none; position: absolute; bottom: 120%; left: 50%;
                      transform: translateX(-50%); background: #333; color: #fff;
                      padding: 6px 10px; border-radius: 4px; font-size: 11px;
                      white-space: normal; width: 220px; z-index: 999;
                      font-weight: normal; line-height: 1.4; }}
.col-tip:hover .tip-text {{ display: block; }}
thead tr.filter-row th input {{
    width: 100%; box-sizing: border-box; font-size: 11px;
    padding: 3px; border: 1px solid #ccc; border-radius: 3px; }}
thead tr.filter-row th {{ padding: 4px 8px; }}
.svg-tooltip {{ position: fixed; background: #333; color: #fff;
                padding: 8px 12px; border-radius: 4px; font-size: 12px;
                max-width: 360px; line-height: 1.5; pointer-events: none;
                z-index: 2000; display: none; white-space: normal; }}
#modal-overlay {{ display: none; position: fixed; top: 0; left: 0;
                  width: 100%; height: 100%; background: rgba(0,0,0,0.6);
                  z-index: 1000; justify-content: center; align-items: center; }}
#modal-overlay.active {{ display: flex; }}
#modal-box {{ background: white; border-radius: 8px; padding: 24px;
              max-width: 700px; width: 90%; position: relative; }}
#modal-title {{ font-size: 1.1em; font-weight: bold; color: #1a3a5c;
                margin-bottom: 12px; }}
#modal-img {{ width: 100%; border: 1px solid #e0e0e0; border-radius: 4px; }}
#modal-close {{ position: absolute; top: 12px; right: 16px; font-size: 1.4em;
                cursor: pointer; color: #555; background: none; border: none; }}
.plot-block {{ margin: 16px 0; }}
.plot-block svg {{ max-width: 100%; }}
.plot-caption {{ font-size: 0.88em; color: #555; margin: 4px 0 0; }}
details.compare-pair {{ background: #fafbfc; border: 1px solid #e0e0e0;
    border-radius: 6px; margin-bottom: 10px; padding: 8px 14px; }}
details.compare-pair > summary {{ cursor: pointer; font-weight: 600;
    color: #1a3a5c; padding: 4px 0; }}
footer {{ text-align: center; padding: 20px; font-size: 0.85em; color: #888; }}
</style>
</head>
<body>
<header>
  <h1>{title}</h1>
  <p>ENHYDRA &mdash; cross-run comparison report</p>
</header>
<div id="svg-tooltip" class="svg-tooltip"></div>
<div id="modal-overlay">
  <div id="modal-box">
    <button id="modal-close" title="Close">&times;</button>
    <div id="modal-title"></div>
    <img id="modal-img" src="" alt="Enrichment plot"/>
  </div>
</div>
<main>
<section>
  <h2>Runs compared</h2>
  {run_summary_html}
</section>
<section>
  <h2>Comparison</h2>
  <nav class="tab-nav" role="tablist">
{tab_buttons}
  </nav>
{tab_panels}
</section>
</main>
<footer>Generated by ENHYDRA compare_runs</footer>
<script>{jquery_js}</script>
<script>{dt_js}</script>
<script>
var enrichmentPlotsMap = {enrichment_plots_map};
var numericColsMap     = {numeric_cols_map};
var dtInstances        = {{}};
var colFiltersMap      = {{}};

$.fn.dataTable.ext.search.push(function(settings, data) {{
    var tableId = settings.nTable.id;
    var cf = colFiltersMap[tableId] || {{}};
    for (var i in cf) {{
        var f = cf[i];
        if (f.text !== undefined) {{
            if (data[i].toLowerCase().indexOf(f.text) === -1) return false;
        }} else {{
            var v = parseFloat(data[i]);
            if (isNaN(v)) return false;
            if (f.op === '<'  && !(v <  f.num)) return false;
            if (f.op === '<=' && !(v <= f.num)) return false;
            if (f.op === '>'  && !(v >  f.num)) return false;
            if (f.op === '>=' && !(v >= f.num)) return false;
            if ((f.op === '=' || f.op === '==') && v !== f.num) return false;
            if (f.op === '!=' && v === f.num)   return false;
        }}
    }}
    return true;
}});

function initTable(tableId) {{
    if (dtInstances[tableId]) return;
    colFiltersMap[tableId] = {{}};
    var numericCols = numericColsMap[tableId] || [];
    var dt = $('#' + tableId).DataTable({{
        pageLength: 25, orderCellsTop: true, order: [[4, 'asc']],
        columnDefs: [{{ targets: numericCols, type: 'num' }}],
    }});
    $('#' + tableId + ' thead tr.filter-row th').each(function(i) {{
        var isNum = numericCols.indexOf(i) !== -1;
        var inp   = $('<input type="text" placeholder="' +
                      (isNum ? 'e.g. < 0.05' : 'Filter...') + '"/>');
        $(this).html(inp);
        inp.on('keyup change', (function(col) {{
            return function() {{
                var val = $.trim(this.value);
                if (!val) {{
                    delete colFiltersMap[tableId][col];
                }} else if (isNum) {{
                    var m = val.match(/^([<>=!]=?)\\s*([\\d.eE+\\-]+)$/);
                    if (m) colFiltersMap[tableId][col] = {{ op: m[1], num: parseFloat(m[2]) }};
                    else   delete colFiltersMap[tableId][col];
                }} else {{
                    colFiltersMap[tableId][col] = {{ text: val.toLowerCase() }};
                }}
                dt.draw();
            }};
        }})(i));
    }});
    dtInstances[tableId] = dt;
}}

function showImageModal(title, uri) {{
    $('#modal-title').text(title);
    $('#modal-img').attr('src', uri).show();
    $('#modal-overlay').addClass('active');
}}

$(document).ready(function() {{
    var svgTip = document.getElementById('svg-tooltip');
    document.querySelectorAll('[data-tip]').forEach(function(el) {{
        el.addEventListener('mousemove', function(e) {{
            svgTip.innerHTML     = this.getAttribute('data-tip');
            svgTip.style.display = 'block';
            svgTip.style.left    = (e.clientX + 15) + 'px';
            svgTip.style.top     = (e.clientY + 15) + 'px';
        }});
        el.addEventListener('mouseleave', function() {{ svgTip.style.display = 'none'; }});
    }});
    document.querySelectorAll('.tab-btn').forEach(function(btn) {{
        btn.addEventListener('click', function() {{
            var key = this.dataset.tab;
            document.querySelectorAll('.tab-btn').forEach(function(b) {{
                b.classList.remove('active');
            }});
            document.querySelectorAll('.tab-panel').forEach(function(p) {{
                p.classList.remove('active');
            }});
            this.classList.add('active');
            document.getElementById('tab-' + key).classList.add('active');
            document.querySelectorAll('#tab-' + key + ' table[id]').forEach(function(t) {{
                if (numericColsMap.hasOwnProperty(t.id)) initTable(t.id);
            }});
        }});
    }});
    // NES cells emitted by build_significance_agreement_table_html() carry
    // data-goid and data-metric="<run_name>" attributes; the enclosing
    // table's own id disambiguates *which ranking metric's* per-run plot
    // index to look in, since the same run/GO-id pair can point to a
    // different plot file under each metric's own enrichment directory.
    $(document).on('click', '.go-link', function(e) {{
        e.preventDefault();
        var goId    = $(this).data('goid');
        var runName = $(this).data('metric');
        var tableId = $(this).closest('table').attr('id');
        var plots   = (enrichmentPlotsMap[tableId] || {{}})[runName];
        var uri     = plots ? plots[goId] : undefined;
        if (uri) showImageModal(goId + ' (' + runName + ')', uri);
    }});
    $('#modal-close, #modal-overlay').on('click', function(e) {{
        if (e.target === this) $('#modal-overlay').removeClass('active');
    }});
    var firstBtn = document.querySelector('.tab-btn');
    if (firstBtn) firstBtn.click();
}});
</script>
</body>
</html>"""


# ---------------------------------------------------------------------------
# Report assembly
# ---------------------------------------------------------------------------

def build_comparison_report(
    run_dirs: list[str],
    run_names: list[str],
    output_path: str,
    fdr_threshold: float = 0.25,
    obo_path: str | None = None,
) -> str:
    """Build a cross-run comparison HTML report.

    Args:
        run_dirs:      Paths to two or more completed ENHYDRA output
                       directories.
        run_names:     Display names, one per run_dirs entry, in the same
                       order. Must be unique.
        output_path:   Path to write the comparison HTML report to.
        fdr_threshold: Significance cutoff applied throughout.
        obo_path:      Optional path to a cached go-basic.obo file, for
                       resolving GO term names. Falls back to any GMT
                       files found in the runs' own enrichment
                       directories if not given (or for terms the OBO
                       doesn't cover).

    Returns:
        The absolute path the report was written to.

    Note:
        The report links to each run's own per-gene-set enrichment plot
        PNGs by relative path rather than embedding them, so it must
        remain on the same filesystem as the run directories (or be
        moved together with them) for those links to resolve.
    """
    report_dir = os.path.dirname(os.path.abspath(output_path)) or "."
    os.makedirs(report_dir, exist_ok=True)

    logger.info("Discovering %d run(s)...", len(run_dirs))
    runs: dict[str, dict] = {}
    for d, name in zip(run_dirs, run_names):
        runs[name] = _discover_run(d)

    # GO term names: any GMT found in any run's results dirs, then OBO
    # (if given) taking precedence on conflicts — same precedence order
    # used by report.py's own _resolve_term_names().
    term_names: dict[str, str] = {}
    for info in runs.values():
        for mpaths in info["metric_paths"].values():
            gmt = _find_gmt_in_dir(mpaths["results_dir"])
            term_names.update(_gmt_term_names(gmt))
    if obo_path and os.path.isfile(obo_path):
        term_names.update(parse_obo_names(obo_path))

    metric_sets = [set(info["metric_paths"].keys()) for info in runs.values()]
    common_metrics = sorted(set.intersection(*metric_sets)) if metric_sets else []

    dfs: dict[str, dict[str, pd.DataFrame]] = {}
    plot_idx: dict[str, dict[str, dict[str, str]]] = {}
    for metric in common_metrics:
        dfs[metric] = {}
        plot_idx[metric] = {}
        for name, info in runs.items():
            mpaths = info["metric_paths"][metric]
            df = _load_run_metric_df(mpaths)
            if df is not None:
                dfs[metric][name] = df
            plot_idx[metric][name] = _build_enrichment_plot_index(
                mpaths["results_dir"], report_dir,
            )

    usable_metrics = [m for m in common_metrics if len(dfs.get(m, {})) >= 2]
    if not usable_metrics:
        logger.warning(
            "No ranking metric is usable across at least 2 of the selected "
            "runs — Significance agreement, Rank correlation, and "
            "Significance transitions tabs will be empty. This happens "
            "when the selected runs were produced with different ranking "
            "metrics (e.g. one with --metric identity, another with "
            "--metric zscore) — re-run with a shared --metric, or with "
            "--all-metrics everywhere, to enable these comparisons."
        )

    numeric_cols_map: dict[str, list[int]] = {}
    enrichment_plots_map: dict[str, dict[str, dict[str, str]]] = {}

    # --- Tab 1: Run parameters ---
    run_params_html = _build_run_params_diff_html(
        [runs[n]["run_parameters"] for n in run_names], run_names,
    )

    # --- Tab 2: Significance agreement ---
    sig_blocks = []
    for metric in usable_metrics:
        table_id = "sig-table-%s" % metric
        column_labels = {n: n for n in dfs[metric]}
        table_html, num_cols = build_significance_agreement_table_html(
            column_dfs=dfs[metric],
            column_labels=column_labels,
            obo_names=term_names,
            fdr_threshold=fdr_threshold,
            enrichment_plots_map=plot_idx[metric],
            table_id=table_id,
        )
        if num_cols is not None:
            numeric_cols_map[table_id] = num_cols
            enrichment_plots_map[table_id] = plot_idx[metric]
        label = METRIC_LABELS.get(metric, metric.capitalize())
        sig_blocks.append(
            "<h4>%s</h4>%s" % (html.escape(label), table_html)
        )
    sig_html = "".join(sig_blocks) if sig_blocks else (
        "<p>No ranking metric is shared by at least 2 of the selected runs.</p>"
    )

    # --- Tab 3: Rank correlation ---
    corr_blocks = []
    for metric in usable_metrics:
        names_m = list(dfs[metric].keys())
        label = METRIC_LABELS.get(metric, metric.capitalize())
        block = ["<h4>%s</h4>" % html.escape(label)]
        if len(names_m) == 2:
            a, b = names_m
            block.append(_pairwise_correlation_block(
                label, a, b, dfs[metric][a], dfs[metric][b], fdr_threshold,
            ))
        else:
            mat = _spearman_matrix(dfs[metric])
            block.append(_render_correlation_matrix_html(mat))
            for a, b in combinations(names_m, 2):
                pair_html = _pairwise_correlation_block(
                    label, a, b, dfs[metric][a], dfs[metric][b], fdr_threshold,
                )
                block.append(
                    '<details class="compare-pair"><summary>%s vs %s</summary>'
                    "%s</details>" % (html.escape(a), html.escape(b), pair_html)
                )
        corr_blocks.append("".join(block))
    corr_html = "".join(corr_blocks) if corr_blocks else (
        "<p>No ranking metric is shared by at least 2 of the selected runs.</p>"
    )

    # --- Tab 4: Significance transitions (Sankey) ---
    sankey_blocks = []
    for metric in usable_metrics:
        names_m = list(dfs[metric].keys())
        label = METRIC_LABELS.get(metric, metric.capitalize())
        block = ["<h4>%s</h4>" % html.escape(label)]
        for a, b in combinations(names_m, 2):
            merged = _merge_pair(dfs[metric][a], dfs[metric][b])
            only_a = set(dfs[metric][a]["Term"]) - set(dfs[metric][b]["Term"])
            only_b = set(dfs[metric][b]["Term"]) - set(dfs[metric][a]["Term"])
            if merged.empty:
                pair_html = "<p>No common gene sets between these two runs.</p>"
            else:
                counts = _sankey_counts(merged, fdr_threshold)
                svg = _render_sankey_svg(a, b, label, counts)
                caption = (
                    "%d gene sets tested in both runs. %d tested only in "
                    "%s; %d tested only in %s (excluded from this diagram)."
                    % (sum(counts.values()), len(only_a), a, len(only_b), b)
                )
                pair_html = (
                    '<div class="plot-block">%s<p class="plot-caption">%s</p></div>'
                    % (svg, caption)
                )
            block.append(
                '<details class="compare-pair" open><summary>%s vs %s</summary>'
                "%s</details>" % (html.escape(a), html.escape(b), pair_html)
            )
        sankey_blocks.append("".join(block))
    sankey_html = "".join(sankey_blocks) if sankey_blocks else (
        "<p>No ranking metric is shared by at least 2 of the selected runs.</p>"
    )

    # --- Assemble tabs ---
    tab_specs = [
        ("run-parameters", "Run parameters", run_params_html),
        ("significance",   "Significance agreement", sig_html),
        ("correlation",    "Rank correlation", corr_html),
        ("sankey",         "Significance transitions", sankey_html),
    ]
    tab_buttons = "\n".join(
        '    <button class="tab-btn" data-tab="%s" role="tab" '
        'aria-controls="tab-%s">%s</button>' % (key, key, label)
        for key, label, _ in tab_specs
    )
    tab_panels = "\n".join(
        '<div id="tab-%s" class="tab-panel" role="tabpanel">%s</div>'
        % (key, body)
        for key, _, body in tab_specs
    )

    cache_dir = os.path.dirname(obo_path) if obo_path else None
    jquery_js = _fetch_cached(_JQUERY_URL, cache_dir, "jquery.min.js")
    dt_js     = _fetch_cached(_DATATABLES_JS_URL, cache_dir, "datatables.min.js")
    dt_css    = _fetch_cached(_DATATABLES_CSS_URL, cache_dir, "datatables.min.css")

    html_content = _COMPARE_TEMPLATE.format(
        title="ENHYDRA Run Comparison (%d runs)" % len(run_names),
        run_summary_html=_build_run_summary_html(runs),
        tab_buttons=tab_buttons,
        tab_panels=tab_panels,
        jquery_js=jquery_js, dt_js=dt_js, dt_css=dt_css,
        enrichment_plots_map=json.dumps(enrichment_plots_map),
        numeric_cols_map=json.dumps(numeric_cols_map),
    )

    output_path = os.path.abspath(output_path)
    with open(output_path, "w", encoding="utf-8") as fh:
        fh.write(html_content)
    logger.info("Comparison report written to: %s", output_path)
    return output_path


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m enhydra.compare_runs",
        description="Compare two or more completed ENHYDRA run output directories.",
    )
    parser.add_argument(
        "--runs", nargs="+", required=True, metavar="DIR",
        help="Two or more ENHYDRA output directories to compare.",
    )
    parser.add_argument(
        "--names", nargs="+", default=None, metavar="NAME",
        help="Display names for each run, in the same order as --runs. "
             "Defaults to each directory's basename. Must be unique.",
    )
    parser.add_argument(
        "--output", default="comparison_report.html",
        help="Path to write the comparison HTML report (default: "
             "comparison_report.html in the current directory).",
    )
    parser.add_argument(
        "--fdr-threshold", type=float, default=0.25,
        help="FDR significance cutoff applied throughout (default: 0.25).",
    )
    parser.add_argument(
        "--obo-cache", default=None,
        help="Directory containing a cached go-basic.obo file, for "
             "resolving GO term names. Optional — falls back to any GMT "
             "files found in the runs' own enrichment directories.",
    )
    return parser


def main():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )
    parser = _build_arg_parser()
    args = parser.parse_args()

    if len(args.runs) < 2:
        parser.error("At least 2 run directories are required for comparison.")
    for r in args.runs:
        if not os.path.isdir(r):
            parser.error("Not a directory: %s" % r)

    if args.names:
        if len(args.names) != len(args.runs):
            parser.error("--names must list exactly one name per --runs entry.")
        names = args.names
    else:
        names = [os.path.basename(os.path.normpath(r)) for r in args.runs]
    if len(set(names)) != len(names):
        parser.error(
            "Run names must be unique — pass --names explicitly to "
            "disambiguate directories that share a basename."
        )

    obo_path = (
        os.path.join(args.obo_cache, "go-basic.obo") if args.obo_cache else None
    )

    build_comparison_report(
        run_dirs=args.runs,
        run_names=names,
        output_path=args.output,
        fdr_threshold=args.fdr_threshold,
        obo_path=obo_path,
    )


if __name__ == "__main__":
    main()
