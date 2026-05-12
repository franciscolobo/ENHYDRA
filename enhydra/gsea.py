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


def _gmt_gene_set(gmt_path: str) -> set[str]:
    """Return the set of all gene IDs present in a GMT file."""
    genes: set[str] = set()
    with open(gmt_path) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            # GMT: term_id <TAB> description <TAB> gene1 <TAB> gene2 ...
            genes.update(fields[2:])
    return genes


def _check_overlap(ranked: pd.DataFrame, gmt_path: str) -> None:
    """Raise a clear error if the ranked list barely overlaps with the GMT.

    Args:
        ranked:   Ranked gene DataFrame (gene_id, identity).
        gmt_path: Path to the GMT file.

    Raises:
        EnhydraIOError: If fewer than 2 genes overlap, with a diagnostic
                        message showing example IDs from both sides.
    """
    ranked_ids = set(ranked["gene_id"].astype(str))
    gmt_ids    = _gmt_gene_set(gmt_path)
    overlap    = ranked_ids & gmt_ids
    n_overlap  = len(overlap)

    logger.info(
        "Preflight: %d ranked genes, %d genes in GMT, %d in common.",
        len(ranked_ids), len(gmt_ids), n_overlap,
    )

    if n_overlap < 2:
        ranked_ex = ", ".join(list(ranked_ids)[:5])
        gmt_ex    = ", ".join(list(gmt_ids)[:5])
        raise EnhydraIOError(
            "Too few genes overlap between the ranked list and the GMT file "
            "(%d in common — GSEApy requires at least 2).\n\n"
            "This usually means the gene IDs in anchor2mean.tsv do not match "
            "those in the GMT file.\n\n"
            "Example ranked IDs : %s\n"
            "Example GMT IDs    : %s\n\n"
            "Check that the GMT was built from the same anchor proteome used "
            "to run the pipeline."
            % (n_overlap, ranked_ex, gmt_ex)
        )


def _convert_gsea_plots_to_png(results_dir: str):
    """Convert GSEApy-generated PDFs to PNGs for HTML embedding."""
    prerank_dir = os.path.join(results_dir, "prerank")
    if not os.path.isdir(prerank_dir):
        return
    pdfs = [f for f in os.listdir(prerank_dir) if f.endswith(".pdf")]
    if not pdfs:
        return
    try:
        from pdf2image import convert_from_path
    except ImportError:
        logger.warning(
            "pdf2image not installed — enrichment plots will not be embedded "
            "in the HTML report. Install with: pip install 'enhydra[plots]'\n"
            "Also requires poppler: brew install poppler (macOS) or "
            "apt install poppler-utils (Linux)."
        )
        return
    logger.info("Converting %d GSEApy PDF plots to PNG...", len(pdfs))
    for pdf_file in pdfs:
        pdf_path = os.path.join(prerank_dir, pdf_file)
        png_path = pdf_path.replace(".pdf", ".png")
        if os.path.isfile(png_path):
            continue
        try:
            images = convert_from_path(pdf_path, dpi=150)
            if images:
                images[0].save(png_path, "PNG")
        except Exception as e:
            logger.warning("Failed to convert %s: %s", pdf_file, e)


def build_gmt_from_gprofiler(
    gene_ids: list[str],
    organism: str,
    gmt_path: str,
    sources: list[str] | None = None,
) -> str:
    """Fetch annotations from g:Profiler and write a GMT file."""
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
        user_threshold=1.0,
        no_evidences=False,
    )
    if results.empty:
        raise EnhydraIOError(
            "g:Profiler returned no annotations for organism '%s'. "
            "Please check the organism name at https://biit.cs.ut.ee/gprofiler."
            % organism
        )
    n_terms = 0
    with open(gmt_path, "w") as fh:
        for _, row in results.iterrows():
            genes = row.get("intersections")
            if not genes:
                continue
            fh.write("%s\t%s\t%s\n" % (row["native"], row["name"], "\t".join(genes)))
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
    fdr_threshold: float = 0.25,
) -> gp.Prerank:
    """Run GSEApy prerank, optionally fetching annotations from g:Profiler."""
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

    # Preflight: fail early with a clear message if IDs don't match
    _check_overlap(ranked, gene_sets)

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
        graph_num=0,
        verbose=False,
    )

    n_sig = (results.res2d["FDR q-val"] < fdr_threshold).sum()
    logger.info("GSEA complete. %d significant gene sets (FDR < %s).", n_sig, fdr_threshold)

    if n_sig > 0:
        logger.info("Generating enrichment plots for %d significant gene sets...", n_sig)
        gp.prerank(
            rnk=ranked,
            gene_sets=gene_sets,
            outdir=results_dir,
            permutation_num=permutations,
            min_size=min_size,
            max_size=max_size,
            seed=seed,
            graph_num=int(n_sig),
            verbose=False,
        )
        _convert_gsea_plots_to_png(results_dir)

    return results
