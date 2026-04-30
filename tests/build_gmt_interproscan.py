"""
Build GMT files from InterProScan TSV output and the GO OBO file.

Workflow:
    InterProScan TSV (for anchor proteome)
        → parse protein → GO term mappings
        → download GO OBO file (cached)
        → resolve GO term names and namespaces
        → write one GMT file per requested namespace

GMT format (tab-separated):
    GO_ID <TAB> GO_term_name <TAB> protein1 <TAB> protein2 ...

Protein IDs are extracted from the InterProScan TSV sequence column.
If your FASTA headers follow the ENHYDRA format (speciesID|proteinID),
only the proteinID portion is used.

Usage:
    python tests/build_gmt_interproscan.py \\
        --interproscan  results/NC_004431.tsv \\
        --anchor        NC_004431 \\
        --outdir        test_data/gmt/ \\
        --cache         test_data/obo_cache/ \\
        --namespaces    GO_BP

Dependencies: none beyond the standard library
"""

from __future__ import annotations

import os
import ssl
import sys
import logging
import argparse
import urllib.request
from collections import defaultdict

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)]
)
logger = logging.getLogger(__name__)

GO_OBO_URL = "https://purl.obolibrary.org/obo/go/go-basic.obo"

# Maps CLI namespace argument to GO OBO namespace string
NAMESPACE_MAP: dict[str, str] = {
    "GO_BP": "biological_process",
    "GO_MF": "molecular_function",
    "GO_CC": "cellular_component",
}


# --- Step 1: Download and parse GO OBO ---

def download_obo(cache_dir: str) -> str:
    """Download go-basic.obo if not already cached.

    Args:
        cache_dir: Directory to cache the OBO file.

    Returns:
        Path to the local OBO file.
    """
    local_path = os.path.join(cache_dir, "go-basic.obo")
    if os.path.isfile(local_path):
        logger.info("Using cached OBO file: %s", local_path)
        return local_path
    logger.info("Downloading GO OBO file from %s...", GO_OBO_URL)
    ctx = ssl.create_default_context()
    ctx.check_hostname = False
    ctx.verify_mode = ssl.CERT_NONE
    with urllib.request.urlopen(GO_OBO_URL, context=ctx) as response, \
         open(local_path, "wb") as fh:
        fh.write(response.read())
    logger.info("Saved to: %s", local_path)
    return local_path


def parse_obo(obo_path: str) -> dict[str, dict[str, str]]:
    """Parse go-basic.obo into a lookup dict.

    Args:
        obo_path: Path to the OBO file.

    Returns:
        Dict mapping GO ID → {'name': str, 'namespace': str}.
        Only includes Terms that are not obsolete.
    """
    logger.info("Parsing OBO file: %s", obo_path)
    go_terms: dict[str, dict[str, str]] = {}
    current_id = None
    current_name = None
    current_namespace = None
    is_obsolete = False

    with open(obo_path) as fh:
        for line in fh:
            line = line.rstrip()
            if line == "[Term]":
                # Save previous term
                if current_id and not is_obsolete:
                    go_terms[current_id] = {
                        "name":      current_name or "",
                        "namespace": current_namespace or "",
                    }
                current_id = None
                current_name = None
                current_namespace = None
                is_obsolete = False
            elif line.startswith("id: GO:"):
                current_id = line[4:]
            elif line.startswith("name: "):
                current_name = line[6:]
            elif line.startswith("namespace: "):
                current_namespace = line[11:]
            elif line == "is_obsolete: true":
                is_obsolete = True

    # Save the last term
    if current_id and not is_obsolete:
        go_terms[current_id] = {
            "name":      current_name or "",
            "namespace": current_namespace or "",
        }

    logger.info("Parsed %d non-obsolete GO terms.", len(go_terms))
    return go_terms


# --- Step 2: Parse InterProScan TSV ---

def parse_interproscan_tsv(
    tsv_path: str,
) -> dict[str, set[str]]:
    """Parse an InterProScan TSV file into a protein → GO terms mapping.

    InterProScan TSV columns (0-indexed):
        0:  Protein accession
        1:  Sequence MD5
        2:  Sequence length
        3:  Analysis
        4:  Signature accession
        5:  Signature description
        6:  Start position
        7:  Stop position
        8:  Score
        9:  Status
        10: Date
        11: InterPro accession (optional)
        12: InterPro description (optional)
        13: GO annotations (optional, pipe-separated)
        14: Pathways (optional)

    Protein IDs are taken from column 0. If the ID follows ENHYDRA's
    'speciesID|proteinID' format, only the proteinID portion is used.

    Args:
        tsv_path: Path to the InterProScan TSV file.

    Returns:
        Dict mapping protein_id → set of GO IDs.
    """
    logger.info("Parsing InterProScan TSV: %s", tsv_path)
    protein_to_go: dict[str, set[str]] = defaultdict(set)
    n_rows = 0
    n_with_go = 0

    with open(tsv_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            n_rows += 1

            # Extract protein ID from lcl-style ID:
            # 'lcl|NC_004431.1_prot_WP_001366457.1_448' → 'WP_001366457.1'
            # Falls back to the raw ID if the pattern is not matched
            raw_id = fields[0]
            if "|" in raw_id and "_prot_" in raw_id:
                # lcl|<accession>_prot_<protein_id>_<index>
                prot_part = raw_id.split("_prot_", 1)[1]
                # Strip the trailing _<index> (last underscore-delimited field)
                protein_id = "_".join(prot_part.split("_")[:-1])
            elif "|" in raw_id:
                protein_id = raw_id.split("|")[-1]
            else:
                protein_id = raw_id

            # GO annotations are in column 13 (pipe-separated), if present
            if len(fields) < 14 or not fields[13].strip():
                continue
            go_ids = [g.strip().split("(")[0] for g in fields[13].split("|") if g.strip()]
            if go_ids:
                protein_to_go[protein_id].update(go_ids)
                n_with_go += 1

    logger.info(
        "Parsed %d rows. %d unique proteins with GO annotations.",
        n_rows, len(protein_to_go)
    )
    return dict(protein_to_go)


# --- Step 3: Build and write GMT ---

def build_and_write_gmt(
    protein_to_go: dict[str, set[str]],
    go_terms: dict[str, dict[str, str]],
    namespace: str,
    out_path: str,
):
    """Build and write a GMT file for a single GO namespace.

    Args:
        protein_to_go: protein_id → set of GO IDs.
        go_terms:      GO ID → {'name': str, 'namespace': str}.
        namespace:     Namespace key to write (e.g. 'GO_BP').
        out_path:      Output GMT file path.
    """
    target_namespace = NAMESPACE_MAP[namespace]
    go_to_proteins: dict[str, tuple[str, set[str]]] = defaultdict(
        lambda: ("", set())
    )

    n_skipped = 0
    for protein_id, go_ids in protein_to_go.items():
        for go_id in go_ids:
            term = go_terms.get(go_id)
            if term is None:
                n_skipped += 1
                continue
            if term["namespace"] != target_namespace:
                continue
            name, proteins = go_to_proteins[go_id]
            if not name:
                go_to_proteins[go_id] = (term["name"], proteins | {protein_id})
            else:
                go_to_proteins[go_id] = (name, proteins | {protein_id})

    if n_skipped:
        logger.warning(
            "%d GO ID(s) not found in OBO file (obsolete or malformed).", n_skipped
        )

    logger.info(
        "Namespace %s: %d GO terms covering %d unique proteins.",
        namespace,
        len(go_to_proteins),
        len({p for _, proteins in go_to_proteins.values() for p in proteins}),
    )

    with open(out_path, "w") as fh:
        for go_id, (go_name, proteins) in sorted(go_to_proteins.items()):
            fh.write("%s\t%s\t%s\n" % (go_id, go_name, "\t".join(sorted(proteins))))

    logger.info("GMT written to: %s", out_path)


# --- Main ---

def main():
    parser = argparse.ArgumentParser(
        description="Build GMT files from InterProScan TSV and GO OBO."
    )
    parser.add_argument(
        "--interproscan", required=True,
        help="Path to the InterProScan TSV output file."
    )
    parser.add_argument(
        "--anchor", required=True,
        help="Anchor genome ID (used to name the output GMT files, e.g. NC_004431)."
    )
    parser.add_argument(
        "--outdir", required=True,
        help="Directory where GMT files will be written."
    )
    parser.add_argument(
        "--cache", required=True,
        help="Directory for caching the GO OBO file."
    )
    parser.add_argument(
        "--namespaces",
        nargs="+",
        default=["GO_BP"],
        choices=list(NAMESPACE_MAP.keys()),
        help="GO namespaces to build GMTs for. "
             "Default: GO_BP. Options: GO_BP GO_MF GO_CC."
    )
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    os.makedirs(args.cache, exist_ok=True)

    # Step 1: Download and parse OBO
    obo_path = download_obo(args.cache)
    go_terms = parse_obo(obo_path)

    # Step 2: Parse InterProScan TSV
    protein_to_go = parse_interproscan_tsv(args.interproscan)

    # Step 3: Build one GMT per namespace
    for namespace in args.namespaces:
        out_path = os.path.join(
            args.outdir, "%s_%s.gmt" % (args.anchor, namespace)
        )
        build_and_write_gmt(protein_to_go, go_terms, namespace, out_path)

    logger.info("Done.")


if __name__ == "__main__":
    main()
