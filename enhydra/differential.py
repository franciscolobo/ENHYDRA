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
