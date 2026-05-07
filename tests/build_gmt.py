"""
Build a GMT file for an anchor strain using NCBI gene2go and gene2refseq.

Workflow:
    anchor genome accession (e.g. NC_004431)
        → NCBI Entrez: resolve taxonomy ID
        → NCBI FTP: download and cache gene2go.gz + gene2refseq.gz
        → filter both files by taxid
        → join on GeneID: protein_accession → GO terms
        → write one GMT file per requested namespace

GMT format (tab-separated):
    GO_ID <TAB> GO_term <TAB> protein1 <TAB> protein2 ...

Usage:
    python tests/build_gmt.py \\
        --anchor    NC_004431 \\
        --email     your@email.com \\
        --outdir    test_data/gmt/ \\
        --cache     test_data/ncbi_cache/ \\
        --namespaces GO_BP

Note:
    gene2refseq.gz is ~1.2 GB compressed. It is downloaded once and cached.
    Subsequent runs reuse the cached file.

Dependencies: Biopython
"""

from __future__ import annotations

import os
import sys
import gzip
import logging
import argparse
import urllib.request
from collections import defaultdict
from Bio import Entrez

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)]
)
logger = logging.getLogger(__name__)

# Mapping from namespace name to the Category string used in gene2go
NAMESPACE_CATEGORIES: dict[str, str] = {
    "GO_BP": "Process",
    "GO_MF": "Function",
    "GO_CC": "Component",
}

NCBI_FTP_BASE = "https://ftp.ncbi.nlm.nih.gov/gene/DATA"


# --- Step 1: Taxonomy ID resolution ---

def get_taxid(genome_accession: str, email: str) -> str:
    """Resolve a genome accession to an NCBI taxonomy ID.

    Args:
        genome_accession: Genome sequence accession (e.g. NC_004431).
        email:            Email address for NCBI Entrez (required by NCBI).

    Returns:
        Taxonomy ID as a string.

    Raises:
        ValueError: If the accession cannot be found in NCBI nuccore.
    """
    Entrez.email = email
    logger.info("Resolving taxid for %s...", genome_accession)

    handle = Entrez.esearch(db="nuccore", term=genome_accession, retmax=1)
    record = Entrez.read(handle)
    handle.close()

    if not record["IdList"]:
        raise ValueError(
            "Could not find '%s' in NCBI nuccore. "
            "Please check the accession." % genome_accession
        )
    uid = record["IdList"][0]

    handle = Entrez.esummary(db="nuccore", id=uid)
    summary = Entrez.read(handle)
    handle.close()

    taxid = str(int(summary[0]["TaxId"]))
    logger.info("Taxid for %s: %s", genome_accession, taxid)
    return taxid


# --- Step 2: NCBI FTP downloads ---

def _download_with_progress(url: str, local_path: str):
    """Download a file from url to local_path with progress reporting."""
    def _progress(block_count, block_size, total_size):
        if total_size > 0 and block_count % 1000 == 0:
            mb_done = block_count * block_size / 1_000_000
            mb_total = total_size / 1_000_000
            logger.info("  %.0f / %.0f MB", mb_done, mb_total)

    urllib.request.urlretrieve(url, local_path, reporthook=_progress)


def download_ncbi_file(filename: str, cache_dir: str) -> str:
    """Download a file from the NCBI FTP gene/DATA directory if not cached.

    Args:
        filename:  Filename on the FTP server (e.g. 'gene2go.gz').
        cache_dir: Local directory to cache the file.

    Returns:
        Path to the local (possibly cached) file.
    """
    local_path = os.path.join(cache_dir, filename)
    if os.path.isfile(local_path):
        logger.info("Using cached file: %s", local_path)
        return local_path

    url = "%s/%s" % (NCBI_FTP_BASE, filename)
    logger.info("Downloading %s (this may take a while)...", filename)
    _download_with_progress(url, local_path)
    logger.info("Saved to: %s", local_path)
    return local_path


# --- Step 3: Load and filter NCBI data files ---

def load_gene2go(
    gene2go_path: str,
    taxid: str,
    namespaces: list[str],
) -> dict[str, dict[str, list[tuple[str, str]]]]:
    """Load gene2go.gz filtered by taxid and namespaces.

    Args:
        gene2go_path: Path to gene2go.gz.
        taxid:        Taxonomy ID to filter on.
        namespaces:   List of namespace keys (e.g. ['GO_BP']).

    Returns:
        Nested dict: gene_id → {namespace → [(go_id, go_term)]}.
    """
    target_categories = {NAMESPACE_CATEGORIES[ns]: ns for ns in namespaces}
    gene_to_go: dict[str, dict[str, list]] = defaultdict(lambda: defaultdict(list))

    logger.info("Loading gene2go for taxid %s, namespaces: %s...", taxid, namespaces)
    with gzip.open(gene2go_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            if fields[0] != taxid:
                continue
            category = fields[7]
            if category not in target_categories:
                continue
            ns = target_categories[category]
            gene_id, go_id, go_term = fields[1], fields[2], fields[5]
            gene_to_go[gene_id][ns].append((go_id, go_term))

    logger.info(
        "Loaded GO annotations for %d genes (taxid %s).", len(gene_to_go), taxid
    )
    return dict(gene_to_go)


def load_gene2refseq(gene2refseq_path: str, taxid: str) -> dict[str, list[str]]:
    """Stream gene2refseq.gz and extract gene_id → protein accessions for taxid.

    Both versioned (WP_000001.1) and unversioned (WP_000001) accessions
    are retained. Entries without a protein accession ('-') are skipped.

    Args:
        gene2refseq_path: Path to gene2refseq.gz.
        taxid:            Taxonomy ID to filter on.

    Returns:
        Dict: gene_id → list of protein accessions.
    """
    gene_to_proteins: dict[str, list[str]] = defaultdict(list)

    logger.info(
        "Streaming gene2refseq for taxid %s (file is large, please wait)...", taxid
    )
    with gzip.open(gene2refseq_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                continue
            if fields[0] != taxid:
                continue
            gene_id = fields[1]
            protein_acc = fields[5]  # protein_accession.version column
            if protein_acc != "-":
                gene_to_proteins[gene_id].append(protein_acc)

    logger.info(
        "Loaded protein accessions for %d genes (taxid %s).", len(gene_to_proteins), taxid
    )
    return dict(gene_to_proteins)


# --- Step 4: Build and write GMT ---

def build_gmt(
    gene_to_go: dict[str, dict[str, list[tuple[str, str]]]],
    gene_to_proteins: dict[str, list[str]],
    namespace: str,
    out_path: str,
):
    """Build and write a GMT file for a single namespace.

    Args:
        gene_to_go:       gene_id → {namespace → [(go_id, go_term)]}.
        gene_to_proteins: gene_id → [protein_accessions].
        namespace:        Namespace key to write (e.g. 'GO_BP').
        out_path:         Output GMT file path.
    """
    go_to_proteins: dict[str, tuple[str, set[str]]] = {}

    n_genes_mapped = 0
    for gene_id, ns_dict in gene_to_go.items():
        go_terms = ns_dict.get(namespace, [])
        proteins = gene_to_proteins.get(gene_id, [])
        if not go_terms or not proteins:
            continue
        n_genes_mapped += 1
        for go_id, go_term in go_terms:
            if go_id not in go_to_proteins:
                go_to_proteins[go_id] = (go_term, set())
            go_to_proteins[go_id][1].update(proteins)

    logger.info(
        "Namespace %s: %d GO terms, %d genes with protein annotations.",
        namespace, len(go_to_proteins), n_genes_mapped,
    )

    with open(out_path, "w") as fh:
        for go_id, (go_term, proteins) in sorted(go_to_proteins.items()):
            fh.write("%s\t%s\t%s\n" % (go_id, go_term, "\t".join(sorted(proteins))))

    logger.info("GMT written to: %s", out_path)


# --- Main ---

def main():
    parser = argparse.ArgumentParser(
        description="Build a GMT file for an anchor strain using NCBI gene2go and gene2refseq."
    )
    parser.add_argument(
        "--anchor", required=True,
        help="Anchor genome accession (e.g. NC_004431)."
    )
    parser.add_argument(
        "--email", required=True,
        help="Email address for NCBI Entrez (required by NCBI policy)."
    )
    parser.add_argument(
        "--outdir", required=True,
        help="Directory where GMT files will be written."
    )
    parser.add_argument(
        "--cache", required=True,
        help="Directory for caching downloaded NCBI data files."
    )
    parser.add_argument(
        "--namespaces",
        nargs="+",
        default=["GO_BP"],
        choices=list(NAMESPACE_CATEGORIES.keys()),
        help="GO namespaces to include. Default: GO_BP. "
             "Options: GO_BP GO_MF GO_CC."
    )
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    os.makedirs(args.cache, exist_ok=True)

    # Step 1: Resolve taxid
    taxid = get_taxid(args.anchor, args.email)

    # Step 2: Download NCBI data files (cached)
    gene2go_path     = download_ncbi_file("gene2go.gz",     args.cache)
    gene2refseq_path = download_ncbi_file("gene2refseq.gz", args.cache)

    # Step 3: Load and filter
    gene_to_go       = load_gene2go(gene2go_path, taxid, args.namespaces)
    gene_to_proteins = load_gene2refseq(gene2refseq_path, taxid)

    # Step 4: Build one GMT per namespace
    for namespace in args.namespaces:
        out_path = os.path.join(
            args.outdir, "%s_%s.gmt" % (args.anchor, namespace)
        )
        build_gmt(gene_to_go, gene_to_proteins, namespace, out_path)

    logger.info("Done.")


if __name__ == "__main__":
    main()
