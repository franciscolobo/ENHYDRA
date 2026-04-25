from __future__ import annotations

import os
import logging
import pandas as pd
import gseapy as gp
from gprofiler import GProfiler

from .exceptions import EnhydraIOError

logger = logging.getLogger(__name__)


def _load_ranked(anchor2mean_path: str) -> pd.DataFrame:
    """Load and sort the anchor2mean table by identity (descending)."""
    if not os.path.isfile(anchor2mean_path):
        raise EnhydraIOError("anchor2mean file not found: %s" % anchor2mean_path)
    ranked = pd.read_csv(
        anchor2mean_path,
        sep="\t",
        header=None,
        names=["gene_id", "identity"],
    )
    ranked["identity"] = pd.to_numeric(ranked["identity"], errors="coerce")
    ranked = ranked.dropna(subset=["identity"])
    return ranked.sort_values("identity", ascending=False)


def build_gmt_from_gprofiler(
    gene_ids: list[str],
    organism: str,
    gmt_path: str,
    sources: list[str] | None = None,
) -> str:
    """Fetch annotations from g:Profiler and write a GMT file.

    Queries the g:Profiler API with the full list of anchor gene IDs to
    retrieve all annotated terms and their gene members. The result is
    written as a GMT file for use with GSEApy prerank.

    Args:
        gene_ids: List of anchor gene IDs (Ensembl format).
        organism: g:Profiler organism name (e.g. 'hsapiens', 'athaliana').
        gmt_path: Path where the GMT file will be written.
        sources:  Data sources to query. Defaults to GO:BP, GO:MF, GO:CC,
                  KEGG, and Reactome.

    Returns:
        Path to the written GMT file.
    """
    if sources is None:
        sources = ["GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC"]

    logger.info(
        "Fetching annotations from g:Profiler for %d genes, organism: %s, sources: %s",
        len(gene_ids), organism, sources,
    )

    gprofiler = GProfiler(user_agent="ENHYDRA", return_dataframe=True)
    results = gprofiler.profile(
        organism=organism,
        query=gene_ids,
        sources=sources,
        significant=False,      # return all terms, not just significant ones
        no_evidences=False,     # include intersections column (gene members)
    )

    if results.empty:
        raise EnhydraIOError(
            "g:Profiler returned no annotations for organism '%s'. "
            "Please check the organism name at https://biit.cs.ut.ee/gprofiler."
            % organism
        )

    # Build GMT: each line is term_id <TAB> term_name <TAB> gene1 <TAB> gene2 ...
    n_terms = 0
    with open(gmt_path, "w") as fh:
        for _, row in results.iterrows():
            genes = row.get("intersections")
            if not genes:
                continue
            gene_str = "\t".join(genes)
            fh.write("%s\t%s\t%s\n" % (row["native"], row["name"], gene_str))
            n_terms += 1

    logger.info("GMT file written: %d terms, path: %s", n_terms, gmt_path)
    return gmt_path


def run_gsea(
    anchor2mean_path: str,
    results_dir: str,
    gene_sets: str,
    organism: str | None = None,
    sources: list[str] | None = None,
    permutations: int = 1000,
    min_size: int = 5,
    max_size: int = 500,
    seed: int = 42,
) -> gp.Prerank:
    """Run GSEApy prerank, optionally fetching annotations from g:Profiler.

    If organism is provided, annotations are fetched automatically from the
    g:Profiler API and a GMT file is generated in results_dir — no manual
    download required. If gene_sets is provided instead, that local GMT file
    is used directly.

    Args:
        anchor2mean_path: Path to anchor2mean.tsv (gene_id, mean_identity).
        results_dir:      Directory where results and GMT are written.
        gene_sets:        Path to a local GMT file. Used when organism is None.
        organism:         g:Profiler organism name. When provided, annotations
                          are fetched via the API and gene_sets is ignored.
        sources:          Data sources for g:Profiler (only used with organism).
        permutations:     Number of GSEA permutations.
        min_size:         Minimum gene set size to test.
        max_size:         Maximum gene set size to test.
        seed:             Random seed for reproducibility.

    Returns:
        A gseapy Prerank result object.
    """
    os.makedirs(results_dir, exist_ok=True)
    ranked = _load_ranked(anchor2mean_path)

    if organism is not None:
        gmt_path = os.path.join(results_dir, "gprofiler_%s.gmt" % organism)
        if os.path.isfile(gmt_path):
            logger.info("GMT file already exists, skipping API call: %s", gmt_path)
        else:
            build_gmt_from_gprofiler(
                gene_ids=ranked["gene_id"].tolist(),
                organism=organism,
                gmt_path=gmt_path,
                sources=sources,
            )
        gene_sets = gmt_path

    logger.info(
        "Running GSEApy prerank: %d genes, gene sets: %s, permutations: %d",
        len(ranked), gene_sets, permutations,
    )

    results = gp.prerank(
        rnk=ranked,
        gene_sets=gene_sets,
        outdir=results_dir,
        permutation_num=permutations,
        min_size=min_size,
        max_size=max_size,
        seed=seed,
        verbose=False,
    )

    n_sig = (results.res2d["FDR q-val"] < 0.25).sum()
    logger.info("GSEA complete. %d significant gene sets (FDR < 0.25).", n_sig)

    return results
