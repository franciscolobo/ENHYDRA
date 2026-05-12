from __future__ import annotations

import os
import logging
import pandas as pd

from .exceptions import EnhydraIOError

logger = logging.getLogger(__name__)


def normalise_scores(
    scores: pd.Series,
    metric: str,
) -> pd.Series:
    """Normalise a Series of identity scores according to the chosen metric.

    Args:
        scores: Series of raw mean alignment identity values.
        metric: Normalisation method:
                'identity' — no transformation, raw identity values.
                'zscore'   — z-score normalisation: (x - mean) / sd.
                'rank'     — normalised rank: rank / N (1 = most conserved).

    Returns:
        Normalised Series with the same index.
    """
    if metric == "identity":
        return scores
    elif metric == "zscore":
        return (scores - scores.mean()) / scores.std()
    elif metric == "rank":
        n = len(scores)
        return scores.rank(ascending=False) / n
    else:
        raise ValueError(
            "Unknown metric '%s'. Choose 'identity', 'zscore', or 'rank'." % metric
        )


def _load_group2mean(tables_dir: str) -> pd.Series:
    """Load group2mean.tsv as a Series indexed by group ID.

    Args:
        tables_dir: Path to a list's tables/ directory.

    Returns:
        Series: group_id → mean_identity (float).
    """
    path = os.path.join(tables_dir, "group2mean.tsv")
    if not os.path.isfile(path):
        raise EnhydraIOError("group2mean.tsv not found: %s" % path)
    df = pd.read_csv(path, sep="\t", header=None, names=["group_id", "identity"])
    df["identity"] = pd.to_numeric(df["identity"], errors="coerce")
    df = df.dropna(subset=["identity"])
    return df.set_index("group_id")["identity"]


def _load_group2anchor(tables_dir: str) -> pd.Series:
    """Load group2anchor.tsv as a Series indexed by group ID.

    Args:
        tables_dir: Path to a list's tables/ directory.

    Returns:
        Series: group_id → anchor_gene_id (str).
    """
    path = os.path.join(tables_dir, "group2anchor.tsv")
    if not os.path.isfile(path):
        raise EnhydraIOError("group2anchor.tsv not found: %s" % path)
    df = pd.read_csv(path, sep="\t", header=None, names=["group_id", "gene_id"])
    return df.set_index("group_id")["gene_id"]


def compute_differential(
    tables_dir1: str,
    tables_dir2: str,
    diff_dir: str,
    metric: str = "identity",
) -> str:
    """Compute differential conservation scores between two ranked lists.

    Subsets both lists to common orthogroups, computes a signed differential
    score per group, and maps groups to anchor gene IDs using list 1's
    group2anchor.tsv (so that gene IDs match the GMT file built from anchor 1).

    Positive scores indicate higher conservation in list 1 (relative to list 2).
    Negative scores indicate higher conservation in list 2.

    Args:
        tables_dir1: Path to list 1's tables/ directory.
        tables_dir2: Path to list 2's tables/ directory.
        diff_dir:    Directory where differential output files are written.
        metric:      Scoring metric:
                     'identity' — mean identity difference (identity_1 - identity_2).
                     'zscore'   — z-score normalised identity difference.
                     'rank'     — normalised rank difference (rank_1/N - rank_2/N),
                                  where rank 1 = most conserved and N = list length.

    Returns:
        Path to the differential anchor2mean.tsv file for use with run_gsea().

    Raises:
        EnhydraIOError: If required input files are missing.
        ValueError:     If no common orthogroups are found between the two lists.
    """
    os.makedirs(diff_dir, exist_ok=True)

    logger.info("Loading group identity scores...")
    mean1 = _load_group2mean(tables_dir1)
    mean2 = _load_group2mean(tables_dir2)

    common_groups = mean1.index.intersection(mean2.index)
    n_common = len(common_groups)
    logger.info(
        "List 1: %d groups. List 2: %d groups. Common: %d groups.",
        len(mean1), len(mean2), n_common,
    )
    if n_common == 0:
        raise ValueError(
            "No common orthogroups found between the two lists. "
            "Check that both lists were derived from the same OrthoFinder run."
        )

    mean1 = mean1.loc[common_groups]
    mean2 = mean2.loc[common_groups]

    if metric == "identity":
        scores = mean1 - mean2
        logger.info("Metric: raw identity difference (identity_1 - identity_2).")
    elif metric == "zscore":
        scores = normalise_scores(mean1, "zscore") - normalise_scores(mean2, "zscore")
        logger.info("Metric: z-score normalised identity difference.")
    elif metric == "rank":
        scores = normalise_scores(mean1, "rank") - normalise_scores(mean2, "rank")
        logger.info("Metric: normalised rank difference (rank_1/N - rank_2/N).")
    else:
        raise ValueError(
            "Unknown metric '%s'. Choose 'identity', 'zscore', or 'rank'." % metric
        )

    scores.name = "differential_score"

    logger.info("Mapping orthogroups to anchor gene IDs (from list 1)...")
    group2anchor = _load_group2anchor(tables_dir1)
    mapped = scores.to_frame().join(group2anchor, how="left")
    n_unmapped = mapped["gene_id"].isna().sum()
    if n_unmapped:
        logger.warning(
            "%d common groups could not be mapped to an anchor gene ID "
            "and will be excluded from GSEA.", n_unmapped
        )
    mapped = mapped.dropna(subset=["gene_id"])

    diff_scores_path = os.path.join(diff_dir, "differential_scores.tsv")
    scores.reset_index().rename(
        columns={"group_id": "group_id", "differential_score": "score"}
    ).to_csv(diff_scores_path, sep="\t", index=False)
    logger.info("Differential scores written to: %s", diff_scores_path)

    anchor2mean_path = os.path.join(diff_dir, "anchor2mean.tsv")
    mapped[["gene_id", "differential_score"]].to_csv(
        anchor2mean_path, sep="\t", index=False, header=False
    )
    logger.info(
        "Differential anchor2mean.tsv written: %d genes. Path: %s",
        len(mapped), anchor2mean_path,
    )

    return anchor2mean_path
