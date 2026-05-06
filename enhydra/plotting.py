from __future__ import annotations

import os
import logging
import pandas as pd
import matplotlib
matplotlib.use("Agg")  # non-interactive backend, safe for pipeline use
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

logger = logging.getLogger(__name__)

# --- Shared style ---
FIGURE_DPI   = 300
FIGURE_STYLE = "seaborn-v0_8-whitegrid"
PALETTE      = {
    "positive": "#2166ac",  # blue  — more conserved in list 1 / overall
    "negative": "#d6604d",  # red   — faster evolving in list 1 / overall
    "neutral":  "#878787",  # grey  — background
}


def _save(fig: plt.Figure, out_dir: str, name: str):
    """Save a figure as both PNG and PDF."""
    for ext in ("png", "pdf"):
        path = os.path.join(out_dir, "%s.%s" % (name, ext))
        fig.savefig(path, dpi=FIGURE_DPI, bbox_inches="tight")
    logger.info("Plot saved: %s (.png/.pdf)", os.path.join(out_dir, name))
    plt.close(fig)


def _load_anchor2mean(anchor2mean_path: str) -> pd.DataFrame:
    """Load anchor2mean.tsv into a DataFrame."""
    return pd.read_csv(
        anchor2mean_path, sep="\t", header=None, names=["gene_id", "score"]
    )


def _load_group2mean(tables_dir: str) -> pd.Series:
    """Load group2mean.tsv as a Series indexed by group ID."""
    path = os.path.join(tables_dir, "group2mean.tsv")
    df = pd.read_csv(path, sep="\t", header=None, names=["group_id", "score"])
    df["score"] = pd.to_numeric(df["score"], errors="coerce")
    return df.dropna().set_index("group_id")["score"]


def _load_gsea_results(results_dir: str) -> pd.DataFrame | None:
    """Load GSEApy prerank results CSV."""
    path = os.path.join(results_dir, "gseapy.gene_set.prerank.report.csv")
    if not os.path.isfile(path):
        logger.warning("GSEA results not found: %s", path)
        return None
    return pd.read_csv(path)


# --- Single-list plots ---

def plot_identity_distribution(
    anchor2mean_path: str,
    plots_dir: str,
    label: str = "Mean alignment identity",
):
    """Histogram of mean alignment identity scores across all groups.

    Args:
        anchor2mean_path: Path to anchor2mean.tsv.
        plots_dir:        Output directory for plots.
        label:            X-axis label (adjust for differential mode).
    """
    df = _load_anchor2mean(anchor2mean_path)
    df["score"] = pd.to_numeric(df["score"], errors="coerce")
    df = df.dropna(subset=["score"])

    with plt.style.context(FIGURE_STYLE):
        fig, ax = plt.subplots(figsize=(7, 4))
        ax.hist(
            df["score"], bins=50,
            color=PALETTE["neutral"], edgecolor="white", linewidth=0.5
        )
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
    """Horizontal bar plot of NES scores for significant gene sets.

    Generates PNG and PDF for publication and an SVG with GO definition
    tooltips for the HTML report.

    Args:
        results_dir:   Directory containing GSEApy results CSV.
        plots_dir:     Output directory for plots.
        obo_names:     GO ID → term name dict for SVG tooltips (optional).
        fdr_threshold: FDR threshold for significance (default: 0.25).
        top_n:         Maximum number of gene sets to show per direction.
        title:         Plot title.
    """
    path = os.path.join(results_dir, "gseapy.gene_set.prerank.report.csv")
    if not os.path.isfile(path):
        logger.warning("GSEA results not found: %s", path)
        return
    df = pd.read_csv(path)

    sig = df[df["FDR q-val"] < fdr_threshold].copy()
    if sig.empty:
        logger.warning("No significant gene sets (FDR < %s) to plot.", fdr_threshold)
        return

    pos = sig[sig["NES"] > 0].nlargest(top_n, "NES")
    neg = sig[sig["NES"] < 0].nsmallest(top_n, "NES")
    plot_df = pd.concat([pos, neg]).sort_values("NES").reset_index(drop=True)

    colors = [
        PALETTE["positive"] if nes > 0 else PALETTE["negative"]
        for nes in plot_df["NES"]
    ]

    with plt.style.context(FIGURE_STYLE):
        fig_height = max(4, len(plot_df) * 0.45)
        fig, ax = plt.subplots(figsize=(10, fig_height))
        bars = ax.barh(plot_df["Term"], plot_df["NES"], color=colors, edgecolor="none")
        ax.axvline(0, color="black", linewidth=0.8, linestyle="--")
        ax.set_xlabel("Normalised Enrichment Score (NES)", fontsize=12)
        ax.set_title(title, fontsize=13)

        # Fixed offset for FDR labels — always outside the bar
        x_range = max(abs(plot_df["NES"].max()), abs(plot_df["NES"].min()))
        offset = x_range * 0.03

        for bar, (_, row) in zip(bars, plot_df.iterrows()):
            if row["NES"] > 0:
                x = bar.get_width() + offset
                ha = "left"
            else:
                x = bar.get_width() - offset
                ha = "right"
            ax.text(
                x, bar.get_y() + bar.get_height() / 2,
                "FDR=%.3f" % row["FDR q-val"],
                va="center", ha=ha, fontsize=7, color="#333"
            )

        # Add extra margin so labels don't get clipped
        ax.set_xlim(
            ax.get_xlim()[0] - offset * 8,
            ax.get_xlim()[1] + offset * 8
        )
        fig.tight_layout()

    _save(fig, plots_dir, "gsea_barplot")

    # --- SVG version with GO definition tooltips ---
    _save_gsea_barplot_svg(plot_df, colors, offset, title, plots_dir, obo_names or {})


def _save_gsea_barplot_svg(
    plot_df: pd.DataFrame,
    colors: list[str],
    offset: float,
    title: str,
    plots_dir: str,
    obo_names: dict[str, str],
):
    """Save an SVG version of the barplot with browser-native hover tooltips."""
    import xml.etree.ElementTree as ET

    n = len(plot_df)
    bar_h = 22
    margin_left = 260
    margin_right = 130
    margin_top = 50
    margin_bottom = 40
    plot_w = 500
    svg_w = margin_left + plot_w + margin_right
    svg_h = margin_top + n * bar_h + margin_bottom

    x_min = plot_df["NES"].min()
    x_max = plot_df["NES"].max()
    pad = max(abs(x_min), abs(x_max)) * 0.25
    x_min -= pad
    x_max += pad

    def to_px(nes):
        return margin_left + (nes - x_min) / (x_max - x_min) * plot_w

    zero_px = to_px(0)

    svg = ET.Element("svg", {
        "xmlns": "http://www.w3.org/2000/svg",
        "width": str(svg_w), "height": str(svg_h),
        "font-family": "Arial, sans-serif", "font-size": "11",
    })

    # Title
    ET.SubElement(svg, "text", {
        "x": str(svg_w // 2), "y": "30",
        "text-anchor": "middle", "font-size": "13", "font-weight": "bold",
        "fill": "#1a3a5c",
    }).text = title

    # Zero line
    ET.SubElement(svg, "line", {
        "x1": str(zero_px), "y1": str(margin_top),
        "x2": str(zero_px), "y2": str(margin_top + n * bar_h),
        "stroke": "#333", "stroke-width": "1", "stroke-dasharray": "4,3",
    })

    # Bars
    for idx, (_, row) in enumerate(plot_df.iterrows()):
        y = margin_top + idx * bar_h
        nes = row["NES"]
        bar_x = min(zero_px, to_px(nes))
        bar_w = abs(to_px(nes) - zero_px)
        color = PALETTE["positive"] if nes > 0 else PALETTE["negative"]

        go_id = row["Term"]
        definition = obo_names.get(go_id, "")
        tooltip_text = "%s: %s | NES=%.3f FDR=%.3f" % (
            go_id, definition, nes, row["FDR q-val"]
        )

        g = ET.SubElement(svg, "g")
        rect = ET.SubElement(g, "rect", {
            "x": str(bar_x), "y": str(y + 3),
            "width": str(max(bar_w, 1)), "height": str(bar_h - 6),
            "fill": color, "opacity": "0.85",
            "rx": "2",
        })
        ET.SubElement(rect, "title").text = tooltip_text

        # Term label
        label_x = margin_left - 5
        label = ET.SubElement(g, "text", {
            "x": str(label_x), "y": str(y + bar_h // 2 + 4),
            "text-anchor": "end", "fill": "#222", "font-size": "10",
        })
        label.text = go_id
        ET.SubElement(label, "title").text = tooltip_text

        # FDR label
        fdr_x = (to_px(nes) + offset * plot_w * 1.5) if nes > 0 \
                else (to_px(nes) - offset * plot_w * 1.5)
        fdr_anchor = "start" if nes > 0 else "end"
        fdr_el = ET.SubElement(g, "text", {
            "x": str(fdr_x), "y": str(y + bar_h // 2 + 4),
            "text-anchor": fdr_anchor, "fill": "#555", "font-size": "9",
        })
        fdr_el.text = "FDR=%.3f" % row["FDR q-val"]

    # X-axis
    ax_y = margin_top + n * bar_h + 15
    ET.SubElement(svg, "line", {
        "x1": str(margin_left), "y1": str(ax_y),
        "x2": str(margin_left + plot_w), "y2": str(ax_y),
        "stroke": "#333", "stroke-width": "1",
    })
    ET.SubElement(svg, "text", {
        "x": str(margin_left + plot_w // 2),
        "y": str(ax_y + 20),
        "text-anchor": "middle", "fill": "#333",
    }).text = "Normalised Enrichment Score (NES)"

    svg_path = os.path.join(plots_dir, "gsea_barplot.svg")
    ET.ElementTree(svg).write(svg_path, encoding="unicode", xml_declaration=False)
    logger.info("SVG barplot written to: %s", svg_path)


# --- Two-list plots ---

def plot_identity_scatter(
    tables_dir1: str,
    tables_dir2: str,
    diff_scores_path: str,
    plots_dir: str,
    label1: str = "List 1 mean identity",
    label2: str = "List 2 mean identity",
):
    """Scatter plot of list 1 vs list 2 identity for common groups.

    Points are coloured by differential score — blue for genes more conserved
    in list 1, red for genes more conserved in list 2.

    Args:
        tables_dir1:      Path to list 1's tables/ directory.
        tables_dir2:      Path to list 2's tables/ directory.
        diff_scores_path: Path to differential_scores.tsv.
        plots_dir:        Output directory for plots.
        label1:           X-axis label.
        label2:           Y-axis label.
    """
    s1 = _load_group2mean(tables_dir1)
    s2 = _load_group2mean(tables_dir2)
    diff = pd.read_csv(diff_scores_path, sep="\t").set_index("group_id")["score"]

    common = s1.index.intersection(s2.index)
    x = s1.loc[common]
    y = s2.loc[common]
    c = diff.reindex(common).fillna(0)

    mask = x.notna() & y.notna()
    x, y, c = x[mask], y[mask], c[mask]

    if x.empty:
        logger.warning("No common groups with valid identity scores for scatter plot.")
        return

    with plt.style.context(FIGURE_STYLE):
        fig, ax = plt.subplots(figsize=(6, 6))
        vmax = c.abs().quantile(0.95)
        sc = ax.scatter(
            x, y, c=c, cmap="RdBu", alpha=0.5, s=8,
            vmin=-vmax, vmax=vmax
        )
        lims = [min(x.min(), y.min()) - 0.02, max(x.max(), y.max()) + 0.02]
        ax.plot(lims, lims, "k--", linewidth=0.8, alpha=0.5)
        ax.set_xlim(lims)
        ax.set_ylim(lims)
        ax.set_xlabel(label1, fontsize=12)
        ax.set_ylabel(label2, fontsize=12)
        ax.set_title("Identity comparison between lists", fontsize=13)
        cbar = fig.colorbar(sc, ax=ax, shrink=0.8)
        cbar.set_label("Differential score", fontsize=10)
        fig.tight_layout()

    _save(fig, plots_dir, "identity_scatter")


def plot_differential_distribution(
    diff_scores_path: str,
    plots_dir: str,
    metric: str = "zscore",
):
    """Histogram of differential scores between the two lists.

    Args:
        diff_scores_path: Path to differential_scores.tsv.
        plots_dir:        Output directory for plots.
        metric:           Metric used (for axis label).
    """
    diff = pd.read_csv(diff_scores_path, sep="\t")
    scores = pd.to_numeric(diff["score"], errors="coerce").dropna()

    label_map = {
        "zscore":   "Z-score difference (list 1 − list 2)",
        "identity": "Identity difference (list 1 − list 2)",
        "rank":     "Normalised rank difference (list 1 − list 2)",
    }
    xlabel = label_map.get(metric, "Differential score")

    with plt.style.context(FIGURE_STYLE):
        fig, ax = plt.subplots(figsize=(7, 4))
        ax.hist(
            scores[scores >= 0], bins=40,
            color=PALETTE["positive"], edgecolor="white",
            linewidth=0.5, label="More conserved in list 1"
        )
        ax.hist(
            scores[scores < 0], bins=40,
            color=PALETTE["negative"], edgecolor="white",
            linewidth=0.5, label="More conserved in list 2"
        )
        ax.axvline(0, color="black", linewidth=0.8, linestyle="--")
        ax.set_xlabel(xlabel, fontsize=12)
        ax.set_ylabel("Number of orthogroups", fontsize=12)
        ax.set_title("Differential conservation score distribution", fontsize=13)
        ax.legend(fontsize=10)
        fig.tight_layout()

    _save(fig, plots_dir, "differential_distribution")


# --- Entry points called from cli.py ---

def make_single_list_plots(
    anchor2mean_path: str,
    results_dir: str,
    plots_dir: str,
    obo_names: dict[str, str] | None = None,
    fdr_threshold: float = 0.25,
):
    os.makedirs(plots_dir, exist_ok=True)
    plot_identity_distribution(anchor2mean_path, plots_dir)
    plot_gsea_barplot(results_dir, plots_dir,
                      obo_names=obo_names, fdr_threshold=fdr_threshold)


def make_differential_plots(
    tables_dir1: str,
    tables_dir2: str,
    diff_dir: str,
    plots_dir: str,
    metric: str,
    obo_names: dict[str, str] | None = None,
    fdr_threshold: float = 0.25,
):
    os.makedirs(plots_dir, exist_ok=True)
    diff_scores = os.path.join(diff_dir, "differential_scores.tsv")
    results_dir = os.path.join(diff_dir, "enrichment")
    plot_identity_distribution(
        os.path.join(diff_dir, "anchor2mean.tsv"),
        plots_dir,
        label="Differential conservation score",
    )
    plot_identity_scatter(tables_dir1, tables_dir2, diff_scores, plots_dir)
    plot_differential_distribution(diff_scores, plots_dir, metric=metric)
    plot_gsea_barplot(
        results_dir, plots_dir,
        obo_names=obo_names, fdr_threshold=fdr_threshold,
        title="Top differentially enriched gene sets"
    )
