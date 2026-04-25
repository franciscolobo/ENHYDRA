import os
import logging
import requests
import pandas as pd
import gseapy as gp

from .exceptions import EnhydraIOError

logger = logging.getLogger(__name__)

_GPROFILER_GMT_URL = "https://biit.cs.ut.ee/gprofiler/static/gprofiler_full_{organism}.ENSG.gmt"


def download_gmt(organism: str, outdir: str) -> str:
    """Download a GMT file for a given organism from g:Profiler.

    GMT files use Ensembl gene IDs, consistent with the ID format expected
    in ENHYDRA FASTA headers.

    Args:
        organism: g:Profiler organism name (e.g. 'hsapiens', 'mmusculus',
                  'athaliana', 'ecoli'). See https://biit.cs.ut.ee/gprofiler
                  for the full list of supported organisms.
        outdir:   Directory where the GMT file will be saved.

    Returns:
        Path to the downloaded GMT file.

    Raises:
        EnhydraIOError: If the download fails or the organism is not found.
    """
    url = _GPROFILER_GMT_URL.format(organism=organism)
    gmt_path = os.path.join(outdir, "gprofiler_%s.ENSG.gmt" % organism)

    if os.path.isfile(gmt_path):
        logger.info("GMT file already exists, skipping download: %s", gmt_path)
        return gmt_path

    logger.info("Downloading GMT file for '%s' from g:Profiler...", organism)
    try:
        response = requests.get(url, timeout=60)
        response.raise_for_status()
    except requests.exceptions.HTTPError as e:
        raise EnhydraIOError(
            "Could not download GMT for organism '%s'. "
            "Please check the organism name at https://biit.cs.ut.ee/gprofiler.\n"
            "HTTP error: %s" % (organism, e)
        )
    except requests.exceptions.RequestException as e:
        raise EnhydraIOError("GMT download failed: %s" % e)

    with open(gmt_path, "w") as fh:
        fh.write(response.text)

    logger.info("GMT file saved to %s", gmt_path)
    return gmt_path


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
        gene_sets:        Path to a local .gmt file (e.g. downloaded via
                          download_gmt()) or an Enrichr database name.
        permutations:     Number of permutations for significance estimation.
        min_size:         Minimum gene set size to test.
        max_size:         Maximum gene set size to test.
        seed:             Random seed for reproducibility.

    Returns:
        A gseapy PreRank result object. Results are also written to gsea_dir.

    Raises:
        EnhydraIOError: If the anchor2mean file is not found.
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
