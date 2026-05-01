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
    fdr_threshold: float = 0.25,
    top_n: int = 20,
    title: str = "Top enriched gene sets",
):
    """Horizontal bar plot of NES scores for significant gene sets.

    Positive NES (more conserved) shown in blue, negative (faster evolving)
    in red. Up to top_n gene sets per direction are shown.

    Args:
        results_dir:   Directory containing GSEApy results CSV.
        plots_dir:     Output directory for plots.
        fdr_threshold: FDR threshold for significance (default: 0.25).
        top_n:         Maximum number of gene sets to show per direction.
        title:         Plot title.
    """
    df = _load_gsea_results(results_dir)
    if df is None:
        return

    sig = df[df["FDR q-val"] < fdr_threshold].copy()
    if sig.empty:
        logger.warning("No significant gene sets (FDR < %s) to plot.", fdr_threshold)
        return

    pos = sig[sig["NES"] > 0].nlargest(top_n, "NES")
    neg = sig[sig["NES"] < 0].nsmallest(top_n, "NES")
    plot_df = pd.concat([pos, neg]).sort_values("NES")

    colors = [
        PALETTE["positive"] if nes > 0 else PALETTE["negative"]
        for nes in plot_df["NES"]
    ]

    with plt.style.context(FIGURE_STYLE):
        fig_height = max(4, len(plot_df) * 0.35)
        fig, ax = plt.subplots(figsize=(9, fig_height))
        bars = ax.barh(plot_df["Term"], plot_df["NES"], color=colors, edgecolor="none")
        ax.axvline(0, color="black", linewidth=0.8, linestyle="--")
        ax.set_xlabel("Normalised Enrichment Score (NES)", fontsize=12)
        ax.set_title(title, fontsize=13)
        for bar, (_, row) in zip(bars, plot_df.iterrows()):
            x = row["NES"] + (0.02 if row["NES"] > 0 else -0.02)
            ha = "left" if row["NES"] > 0 else "right"
            ax.text(
                x, bar.get_y() + bar.get_height() / 2,
                "FDR=%.3f" % row["FDR q-val"],
                va="center", ha=ha, fontsize=7, color="black"
            )
        fig.tight_layout()

    _save(fig, plots_dir, "gsea_barplot")


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
):
    """Generate all plots for single-list mode.

    Args:
        anchor2mean_path: Path to anchor2mean.tsv.
        results_dir:      Directory containing GSEApy results.
        plots_dir:        Output directory for plots.
    """
    os.makedirs(plots_dir, exist_ok=True)
    plot_identity_distribution(anchor2mean_path, plots_dir)
    plot_gsea_barplot(results_dir, plots_dir)


def make_differential_plots(
    tables_dir1: str,
    tables_dir2: str,
    diff_dir: str,
    plots_dir: str,
    metric: str,
):
    """Generate all plots for two-list differential mode.

    Args:
        tables_dir1: Path to list 1's tables/ directory.
        tables_dir2: Path to list 2's tables/ directory.
        diff_dir:    Path to the differential/ output directory.
        plots_dir:   Output directory for plots.
        metric:      Differential metric used ('zscore', 'identity', 'rank').
    """
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
        title="Top differentially enriched gene sets"
    )
