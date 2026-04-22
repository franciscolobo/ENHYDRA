import os
import logging
import pandas as pd
import gseapy as gp

from .exceptions import EnhydraIOError

logger = logging.getLogger(__name__)


def run_gsea(
    anchor2mean_path: str,
    gsea_dir: str,
    gene_sets: str,
    permutations: int = 1000,
    min_size: int = 5,
    max_size: int = 500,
    seed: int = 42,
) -> gp.PreRank:
    """Run GSEApy prerank on the anchor gene identity ranking.

    Args:
        anchor2mean_path: Path to anchor2mean.tsv (gene_id, mean_identity).
        gsea_dir:         Directory where GSEA results are written.
        gene_sets:        Gene set database name (e.g. 'GO_Biological_Process_2023')
                          or path to a local .gmt file.
        permutations:     Number of permutations for significance estimation.
        min_size:         Minimum gene set size to test.
        max_size:         Maximum gene set size to test.
        seed:             Random seed for reproducibility.

    Returns:
        A gseapy PreRank result object. Results are also written to gsea_dir.
    """
    if not os.path.isfile(anchor2mean_path):
        raise EnhydraIOError("anchor2mean file not found: %s" % anchor2mean_path)

    logger.info("Loading ranked gene list from %s", anchor2mean_path)
    ranked = pd.read_csv(
        anchor2mean_path,
        sep="\t",
        header=None,
        names=["gene_id", "identity"],
    )
    ranked["identity"] = pd.to_numeric(ranked["identity"], errors="coerce")
    ranked = ranked.dropna(subset=["identity"])
    ranked = ranked.sort_values("identity", ascending=False)

    logger.info(
        "Running GSEA prerank: %d genes, gene sets: %s, permutations: %d",
        len(ranked), gene_sets, permutations,
    )

    results = gp.prerank(
        rnk=ranked,
        gene_sets=gene_sets,
        outdir=gsea_dir,
        permutation_num=permutations,
        min_size=min_size,
        max_size=max_size,
        seed=seed,
        verbose=False,
    )

    n_sig = (results.res2d["FDR q-val"] < 0.25).sum()
    logger.info("GSEA complete. %d significant gene sets (FDR < 0.25).", n_sig)

    return results
