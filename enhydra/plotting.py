from __future__ import annotations

import os
import logging
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

logger = logging.getLogger(__name__)

FIGURE_DPI   = 300
FIGURE_STYLE = "seaborn-v0_8-whitegrid"
PALETTE      = {
    "positive": "#2166ac",
    "negative": "#d6604d",
    "neutral":  "#878787",
}


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

def _save(fig: plt.Figure, out_dir: str, name: str):
    for ext in ("png", "pdf"):
        path = os.path.join(out_dir, "%s.%s" % (name, ext))
        fig.savefig(path, dpi=FIGURE_DPI, bbox_inches="tight")
    logger.info("Plot saved: %s (.png/.pdf)", os.path.join(out_dir, name))
    plt.close(fig)


def _load_anchor2mean(anchor2mean_path: str) -> pd.DataFrame:
    return pd.read_csv(
        anchor2mean_path, sep="\t", header=None, names=["gene_id", "score"]
    )


def _load_group2mean(tables_dir: str) -> pd.Series:
    path = os.path.join(tables_dir, "group2mean.tsv")
    df   = pd.read_csv(path, sep="\t", header=None, names=["group_id", "score"])
    df["score"] = pd.to_numeric(df["score"], errors="coerce")
    return df.dropna().set_index("group_id")["score"]


def _load_gsea_results(results_dir: str) -> pd.DataFrame | None:
    path = os.path.join(results_dir, "gseapy.gene_set.prerank.report.csv")
    if not os.path.isfile(path):
        logger.warning("GSEA results not found: %s", path)
        return None
    return pd.read_csv(path)


def _rdbu_color(t: float) -> str:
    """Map t ∈ [-1, 1] to an RdBu hex color (red → white → blue)."""
    t    = max(-1.0, min(1.0, t))
    RED  = (214, 96,  77)
    MID  = (247, 247, 247)
    BLUE = (33,  102, 172)
    src, dst, s = (RED, MID, -t) if t < 0 else (MID, BLUE, t)
    return "rgb(%d,%d,%d)" % tuple(int(src[i] + s * (dst[i] - src[i])) for i in range(3))


# ---------------------------------------------------------------------------
# Single-list plots
# ---------------------------------------------------------------------------

def plot_identity_distribution(
    anchor2mean_path: str,
    plots_dir: str,
    label: str = "Mean alignment identity",
):
    df = _load_anchor2mean(anchor2mean_path)
    df["score"] = pd.to_numeric(df["score"], errors="coerce")
    df = df.dropna(subset=["score"])

    with plt.style.context(FIGURE_STYLE):
        fig, ax = plt.subplots(figsize=(7, 4))
        ax.hist(df["score"], bins=50,
                color=PALETTE["neutral"], edgecolor="white", linewidth=0.5)
        ax.set_xlabel(label, fontsize=12)
        ax.set_ylabel("Number of orthogroups", fontsize=12)
        ax.set_title("Distribution of %s" % label.lower(), fontsize=13)
        ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
        fig.tight_layout()
    _save(fig, plots_dir, "identity_distribution")


def plot_gsea_barplot(
    results_dir: str,
    plots_dir: str,
    obo_names: dict[str, str] | None = None,
    fdr_threshold: float = 0.25,
    top_n: int = 20,
    title: str = "Top enriched gene sets",
):
    path = os.path.join(results_dir, "gseapy.gene_set.prerank.report.csv")
    if not os.path.isfile(path):
        logger.warning("GSEA results not found: %s", path)
        return
    df  = pd.read_csv(path)
    sig = df[df["FDR q-val"] < fdr_threshold].copy()
    if sig.empty:
        logger.warning("No significant gene sets (FDR < %s) to plot.", fdr_threshold)
        return

    pos     = sig[sig["NES"] > 0].nlargest(top_n, "NES")
    neg     = sig[sig["NES"] < 0].nsmallest(top_n, "NES")
    plot_df = pd.concat([pos, neg]).sort_values("NES").reset_index(drop=True)
    colors  = [PALETTE["positive"] if nes > 0 else PALETTE["negative"]
               for nes in plot_df["NES"]]

    with plt.style.context(FIGURE_STYLE):
        fig_height = max(4, len(plot_df) * 0.45)
        fig, ax    = plt.subplots(figsize=(10, fig_height))
        bars       = ax.barh(plot_df["Term"], plot_df["NES"],
                             color=colors, edgecolor="none")
        ax.axvline(0, color="black", linewidth=0.8, linestyle="--")
        ax.set_xlabel("Normalised Enrichment Score (NES)", fontsize=12)
        ax.set_title(title, fontsize=13)

        x_range = max(abs(plot_df["NES"].max()), abs(plot_df["NES"].min()))
        offset  = x_range * 0.03

        for bar, (_, row) in zip(bars, plot_df.iterrows()):
            if row["NES"] > 0:
                x, ha = bar.get_width() + offset, "left"
            else:
                x, ha = bar.get_width() - offset, "right"
            ax.text(x, bar.get_y() + bar.get_height() / 2,
                    "FDR=%.3f" % row["FDR q-val"],
                    va="center", ha=ha, fontsize=7, color="#333")

        ax.set_xlim(ax.get_xlim()[0] - offset * 8,
                    ax.get_xlim()[1] + offset * 8)
        fig.tight_layout()
    _save(fig, plots_dir, "gsea_barplot")
    _save_gsea_barplot_svg(plot_df, colors, offset, title, plots_dir, obo_names or {})


def _save_gsea_barplot_svg(
    plot_df: pd.DataFrame,
    colors: list[str],
    offset: float,
    title: str,
    plots_dir: str,
    obo_names: dict[str, str],
):
    import xml.etree.ElementTree as ET

    n            = len(plot_df)
    bar_h        = 18
    margin_left  = 120
    margin_right = 20
    margin_top   = 50
    margin_bottom = 50
    plot_w       = 400
    svg_w        = margin_left + plot_w + margin_right
    svg_h        = margin_top + n * bar_h + margin_bottom

    nes_vals  = plot_df["NES"].values
    x_abs_max = max(abs(nes_vals.min()), abs(nes_vals.max()))
    x_min     = -x_abs_max * 1.05
    x_max     =  x_abs_max * 1.05

    def to_px(nes):
        return margin_left + (nes - x_min) / (x_max - x_min) * plot_w

    zero_px = to_px(0)

    svg = ET.Element("svg", {
        "xmlns": "http://www.w3.org/2000/svg",
        "width": str(svg_w), "height": str(svg_h),
        "font-family": "Arial, sans-serif",
    })

    ET.SubElement(svg, "text", {
        "x": str(svg_w // 2), "y": "30",
        "text-anchor": "middle", "font-size": "13",
        "font-weight": "bold", "fill": "#1a3a5c",
    }).text = title

    ET.SubElement(svg, "line", {
        "x1": str(zero_px), "y1": str(margin_top),
        "x2": str(zero_px), "y2": str(margin_top + n * bar_h),
        "stroke": "#555", "stroke-width": "1", "stroke-dasharray": "4,3",
    })

    for nes_tick in [-x_abs_max, -x_abs_max / 2, 0, x_abs_max / 2, x_abs_max]:
        tx   = to_px(nes_tick)
        ax_y = margin_top + n * bar_h
        ET.SubElement(svg, "line", {
            "x1": str(tx), "y1": str(ax_y),
            "x2": str(tx), "y2": str(ax_y + 4),
            "stroke": "#555", "stroke-width": "1",
        })
        ET.SubElement(svg, "text", {
            "x": str(tx), "y": str(ax_y + 15),
            "text-anchor": "middle", "font-size": "9", "fill": "#555",
        }).text = "%.2f" % nes_tick

    ET.SubElement(svg, "text", {
        "x": str(margin_left + plot_w // 2),
        "y": str(margin_top + n * bar_h + 35),
        "text-anchor": "middle", "font-size": "11", "fill": "#333",
    }).text = "Normalised Enrichment Score (NES)"

    for idx, (_, row) in enumerate(plot_df.iterrows()):
        y        = margin_top + idx * bar_h
        nes      = row["NES"]
        bar_x    = min(zero_px, to_px(nes))
        bar_w    = max(abs(to_px(nes) - zero_px), 1)
        color    = PALETTE["positive"] if nes > 0 else PALETTE["negative"]
        go_id    = row["Term"]
        defn     = obo_names.get(go_id, "")
        tip_text = "%s | %s | NES=%.3f | FDR=%.3f" % (
            go_id, defn, nes, row["FDR q-val"]
        )

        g = ET.SubElement(svg, "g", {
            "class": "bar-group", "data-tip": tip_text, "style": "cursor:pointer;",
        })
        ET.SubElement(g, "rect", {
            "x": str(bar_x), "y": str(y + 5),
            "width": str(bar_w), "height": str(bar_h - 10),
            "fill": color, "opacity": "0.85", "rx": "2",
        })
        ET.SubElement(g, "text", {
            "x": str(margin_left - 6), "y": str(y + bar_h // 2 + 4),
            "text-anchor": "end", "font-size": "10", "fill": "#1a3a5c",
        }).text = go_id

        fdr_x      = to_px(nes) + 5 if nes > 0 else to_px(nes) - 5
        fdr_anchor = "start" if nes > 0 else "end"
        ET.SubElement(g, "text", {
            "x": str(fdr_x), "y": str(y + bar_h // 2 + 4),
            "text-anchor": fdr_anchor, "font-size": "9", "fill": "#555",
        }).text = "FDR=%.3f" % row["FDR q-val"]

        ET.SubElement(g, "rect", {
            "x": str(margin_left), "y": str(y),
            "width": str(plot_w), "height": str(bar_h), "fill": "transparent",
        })

    svg_path = os.path.join(plots_dir, "gsea_barplot.svg")
    ET.ElementTree(svg).write(svg_path, encoding="unicode", xml_declaration=False)
    logger.info("SVG barplot written to: %s", svg_path)


# ---------------------------------------------------------------------------
# Two-list plots
# ---------------------------------------------------------------------------

def _save_identity_scatter_svg(
    x: pd.Series,
    y: pd.Series,
    c: pd.Series,
    group2anchor: dict[str, str],
    out_path: str,
    label1: str,
    label2: str,
    vmax: float,
):
    """Generate an interactive SVG scatter plot with per-point hover tooltips.

    Each circle represents one orthogroup.  Hovering shows the orthogroup ID,
    anchor gene ID, list 1 score, list 2 score, and differential score.
    Color encodes the differential score on an RdBu scale (blue = more
    conserved in list 1, red = more conserved in list 2).

    Args:
        x:             List 1 scores, indexed by group_id.
        y:             List 2 scores, indexed by group_id.
        c:             Differential scores, indexed by group_id.
        group2anchor:  group_id → anchor gene ID.
        out_path:      Output SVG file path.
        label1:        X-axis label.
        label2:        Y-axis label.
        vmax:          Differential score magnitude used for colour saturation
                       (typically the 95th percentile of abs(c)).
    """
    import xml.etree.ElementTree as ET

    ml, mr, mt, mb = 55, 20, 30, 55
    pw, ph = 370, 370
    sw, sh = ml + pw + mr, mt + ph + mb

    lo  = min(float(x.min()), float(y.min())) - 0.02
    hi  = max(float(x.max()), float(y.max())) + 0.02
    rng = hi - lo or 1.0

    def spx(v):  return ml + (v - lo) / rng * pw
    def spy(v):  return mt + ph - (v - lo) / rng * ph   # y-axis flipped

    svg = ET.Element("svg", {
        "xmlns": "http://www.w3.org/2000/svg",
        "width": str(sw), "height": str(sh),
        "font-family": "Arial, sans-serif",
    })

    # Title
    ET.SubElement(svg, "text", {
        "x": str(sw // 2), "y": "20",
        "text-anchor": "middle", "font-size": "12",
        "font-weight": "bold", "fill": "#1a3a5c",
    }).text = "Identity comparison between lists"

    # Axes
    ax_y = mt + ph
    ET.SubElement(svg, "line", {
        "x1": str(ml), "y1": str(mt), "x2": str(ml), "y2": str(ax_y),
        "stroke": "#555", "stroke-width": "1"})
    ET.SubElement(svg, "line", {
        "x1": str(ml), "y1": str(ax_y), "x2": str(ml + pw), "y2": str(ax_y),
        "stroke": "#555", "stroke-width": "1"})

    # 1:1 diagonal
    ET.SubElement(svg, "line", {
        "x1": str(spx(lo)), "y1": str(spy(lo)),
        "x2": str(spx(hi)), "y2": str(spy(hi)),
        "stroke": "#555", "stroke-width": "0.8",
        "stroke-dasharray": "4,3", "opacity": "0.5"})

    # Ticks
    for v in [lo + rng * i / 4 for i in range(5)]:
        ET.SubElement(svg, "line", {
            "x1": str(spx(v)), "y1": str(ax_y),
            "x2": str(spx(v)), "y2": str(ax_y + 4),
            "stroke": "#555", "stroke-width": "1"})
        ET.SubElement(svg, "text", {
            "x": str(spx(v)), "y": str(ax_y + 14),
            "text-anchor": "middle", "font-size": "9", "fill": "#555"
        }).text = "%.2f" % v
        ET.SubElement(svg, "line", {
            "x1": str(ml - 4), "y1": str(spy(v)),
            "x2": str(ml),     "y2": str(spy(v)),
            "stroke": "#555", "stroke-width": "1"})
        ET.SubElement(svg, "text", {
            "x": str(ml - 7), "y": str(spy(v) + 4),
            "text-anchor": "end", "font-size": "9", "fill": "#555"
        }).text = "%.2f" % v

    # Axis labels
    ET.SubElement(svg, "text", {
        "x": str(ml + pw // 2), "y": str(sh - 5),
        "text-anchor": "middle", "font-size": "11", "fill": "#333"
    }).text = label1
    ET.SubElement(svg, "text", {
        "x": "0", "y": "0",
        "transform": "translate(13,%d) rotate(-90)" % (mt + ph // 2),
        "text-anchor": "middle", "font-size": "11", "fill": "#333"
    }).text = label2

    # Data points — draw below average first so denser regions don't obscure outliers
    for gid in x.index:
        xi   = float(x[gid])
        yi   = float(y[gid])
        ci   = float(c.get(gid, 0.0))
        gene = group2anchor.get(gid, "—")
        col  = _rdbu_color(ci / vmax if vmax else 0.0)
        tip  = (
            "Group: %s | Gene: %s | %s: %.4f | %s: %.4f | diff: %+.4f"
            % (gid, gene, label1, xi, label2, yi, ci)
        )
        ET.SubElement(svg, "circle", {
            "cx": "%.1f" % spx(xi), "cy": "%.1f" % spy(yi),
            "r": "4", "fill": col, "opacity": "0.65",
            "data-tip": tip, "style": "cursor:pointer;",
        })

    ET.ElementTree(svg).write(out_path, encoding="unicode", xml_declaration=False)
    logger.info("Scatter SVG written to: %s", out_path)


def plot_identity_scatter(
    tables_dir1: str,
    tables_dir2: str,
    diff_scores_path: str,
    plots_dir: str,
    label1: str = "List 1 mean identity",
    label2: str = "List 2 mean identity",
    metric: str = "identity",
):
    """Scatter plot of list 1 vs list 2 scores for common orthogroups.

    Produces both a static PNG (for publication) and an SVG with per-point
    hover tooltips showing orthogroup ID, anchor gene ID, and scores.

    Args:
        tables_dir1:      Path to list 1's tables/ directory.
        tables_dir2:      Path to list 2's tables/ directory.
        diff_scores_path: Path to differential_scores.tsv.
        plots_dir:        Output directory for plots.
        label1:           X-axis label.
        label2:           Y-axis label.
        metric:           Metric name (used in hover tooltip labels).
    """
    s1   = _load_group2mean(tables_dir1)
    s2   = _load_group2mean(tables_dir2)
    diff = pd.read_csv(diff_scores_path, sep="\t").set_index("group_id")["score"]

    common = s1.index.intersection(s2.index)
    x = s1.loc[common]
    y = s2.loc[common]
    c = diff.reindex(common).fillna(0)

    mask = x.notna() & y.notna()
    x, y, c = x[mask], y[mask], c[mask]

    if x.empty:
        logger.warning("No common groups with valid scores for scatter plot.")
        return

    # Load group → anchor gene ID mapping for hover labels.
    group2anchor: dict[str, str] = {}
    g2a_path = os.path.join(tables_dir1, "group2anchor.tsv")
    if os.path.isfile(g2a_path):
        with open(g2a_path) as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 2:
                    group2anchor[parts[0]] = parts[1]

    vmax = float(c.abs().quantile(0.95)) or 1.0

    # Static PNG / PDF
    with plt.style.context(FIGURE_STYLE):
        fig, ax = plt.subplots(figsize=(6, 6))
        sc = ax.scatter(x, y, c=c, cmap="RdBu", alpha=0.5, s=8,
                        vmin=-vmax, vmax=vmax)
        lims = [min(x.min(), y.min()) - 0.02, max(x.max(), y.max()) + 0.02]
        ax.plot(lims, lims, "k--", linewidth=0.8, alpha=0.5)
        ax.set_xlim(lims); ax.set_ylim(lims)
        ax.set_xlabel(label1, fontsize=12)
        ax.set_ylabel(label2, fontsize=12)
        ax.set_title("Identity comparison between lists", fontsize=13)
        cbar = fig.colorbar(sc, ax=ax, shrink=0.8)
        cbar.set_label("Differential score", fontsize=10)
        fig.tight_layout()
    _save(fig, plots_dir, "identity_scatter")

    # Interactive SVG for HTML report
    _save_identity_scatter_svg(
        x, y, c, group2anchor,
        out_path=os.path.join(plots_dir, "identity_scatter.svg"),
        label1=label1, label2=label2, vmax=vmax,
    )


def plot_differential_distribution(
    diff_scores_path: str,
    plots_dir: str,
    metric: str = "zscore",
):
    diff   = pd.read_csv(diff_scores_path, sep="\t")
    scores = pd.to_numeric(diff["score"], errors="coerce").dropna()

    label_map = {
        "zscore":   "Z-score difference (list 1 \u2212 list 2)",
        "identity": "Identity difference (list 1 \u2212 list 2)",
        "rank":     "Normalised rank difference (list 1 \u2212 list 2)",
    }
    xlabel = label_map.get(metric, "Differential score")

    with plt.style.context(FIGURE_STYLE):
        fig, ax = plt.subplots(figsize=(7, 4))
        ax.hist(scores[scores >= 0], bins=40, color=PALETTE["positive"],
                edgecolor="white", linewidth=0.5, label="More conserved in list 1")
        ax.hist(scores[scores < 0], bins=40, color=PALETTE["negative"],
                edgecolor="white", linewidth=0.5, label="More conserved in list 2")
        ax.axvline(0, color="black", linewidth=0.8, linestyle="--")
        ax.set_xlabel(xlabel, fontsize=12)
        ax.set_ylabel("Number of orthogroups", fontsize=12)
        ax.set_title("Differential conservation score distribution", fontsize=13)
        ax.legend(fontsize=10)
        fig.tight_layout()
    _save(fig, plots_dir, "differential_distribution")


# ---------------------------------------------------------------------------
# Entry points called from cli.py
# ---------------------------------------------------------------------------

def make_single_list_plots(
    anchor2mean_path: str,
    results_dir: str,
    plots_dir: str,
    obo_names: dict[str, str] | None = None,
    fdr_threshold: float = 0.25,
    top_n: int = 20,
):
    os.makedirs(plots_dir, exist_ok=True)
    plot_identity_distribution(anchor2mean_path, plots_dir)
    plot_gsea_barplot(results_dir, plots_dir,
                      obo_names=obo_names, fdr_threshold=fdr_threshold, top_n=top_n)


def make_differential_plots(
    tables_dir1: str,
    tables_dir2: str,
    diff_dir: str,
    plots_dir: str,
    metric: str,
    obo_names: dict[str, str] | None = None,
    fdr_threshold: float = 0.25,
    top_n: int = 20,
):
    os.makedirs(plots_dir, exist_ok=True)
    diff_scores = os.path.join(diff_dir, "differential_scores.tsv")
    results_dir = os.path.join(diff_dir, "enrichment")

    _mlabels = {"identity": "mean identity", "zscore": "z-score", "rank": "rank"}
    mlabel   = _mlabels.get(metric, metric)

    plot_identity_distribution(
        os.path.join(diff_dir, "anchor2mean.tsv"),
        plots_dir,
        label="Differential conservation score",
    )
    plot_identity_scatter(
        tables_dir1, tables_dir2, diff_scores, plots_dir,
        label1="List 1 %s" % mlabel,
        label2="List 2 %s" % mlabel,
        metric=metric,
    )
    plot_differential_distribution(diff_scores, plots_dir, metric=metric)
    plot_gsea_barplot(
        results_dir, plots_dir,
        obo_names=obo_names, fdr_threshold=fdr_threshold, top_n=top_n,
        title="Top differentially enriched gene sets",
    )
