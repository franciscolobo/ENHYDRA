"""
Build a GMT file from a GAF (Gene Association File) annotation file.

The script reads GO annotations from a GAF 2.x file, optionally propagates
them up the GO hierarchy (true-path rule), resolves GO term names and
namespaces from the GO OBO file, and writes one GMT file per requested
namespace.

GAF column layout (1-indexed, tab-separated):
    1   DB
    2   DB Object ID    ← gene ID used by ENHYDRA
    3   DB Object Symbol
    4   Qualifier
    5   GO ID
    6   DB:Reference
    7   Evidence Code
    8   With/From
    9   Aspect           P = biological_process
                         F = molecular_function
                         C = cellular_component
    ...

GMT format (tab-separated):
    GO_ID <TAB> GO_term_name <TAB> gene1 <TAB> gene2 ...

Usage:
    python tests/build_gmt_gaf.py \\
        --gaf        /path/to/Amel.gaf \\
        --anchor     Apis_mellifera \\
        --outdir     test_data/gmt/ \\
        --cache      test_data/obo_cache/ \\
        --namespaces GO_BP

    # To exclude IEA (inferred from electronic annotation) evidence:
    python tests/build_gmt_gaf.py \\
        --gaf        /path/to/Amel.gaf \\
        --anchor     Apis_mellifera \\
        --outdir     test_data/gmt/ \\
        --cache      test_data/obo_cache/ \\
        --namespaces GO_BP GO_MF GO_CC \\
        --exclude-evidence IEA

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
    handlers=[logging.StreamHandler(sys.stdout)],
)
logger = logging.getLogger(__name__)

GO_OBO_URL = "https://purl.obolibrary.org/obo/go/go-basic.obo"

NAMESPACE_MAP: dict[str, str] = {
    "GO_BP": "biological_process",
    "GO_MF": "molecular_function",
    "GO_CC": "cellular_component",
}

# GAF column indices (0-based)
_COL_GENE_ID       = 1
_COL_QUALIFIER     = 3
_COL_GO_ID         = 4
_COL_EVIDENCE_CODE = 6
_COL_ASPECT        = 8

# Aspect letter → OBO namespace string
_ASPECT_TO_NAMESPACE = {
    "P": "biological_process",
    "F": "molecular_function",
    "C": "cellular_component",
}


# ---------------------------------------------------------------------------
# GO OBO parsing
# ---------------------------------------------------------------------------

def download_obo(cache_dir: str) -> str:
    """Download go-basic.obo if not already cached."""
    local_path = os.path.join(cache_dir, "go-basic.obo")
    if os.path.isfile(local_path):
        logger.info("Using cached OBO file: %s", local_path)
        return local_path
    logger.info("Downloading GO OBO from %s ...", GO_OBO_URL)
    ctx = ssl.create_default_context()
    ctx.check_hostname = False
    ctx.verify_mode    = ssl.CERT_NONE
    with urllib.request.urlopen(GO_OBO_URL, context=ctx) as r, \
         open(local_path, "wb") as fh:
        fh.write(r.read())
    logger.info("Saved to: %s", local_path)
    return local_path


def parse_obo(obo_path: str) -> tuple[
    dict[str, dict[str, str]],   # go_id → {"name": str, "namespace": str}
    dict[str, set[str]],          # go_id → set of direct parent go_ids
]:
    """Parse go-basic.obo into a term dict and a parent-relationship dict."""
    logger.info("Parsing OBO file: %s", obo_path)
    go_terms: dict[str, dict[str, str]] = {}
    parents:  dict[str, set[str]]       = {}

    cur_id, cur_name, cur_ns, cur_parents, obsolete = None, None, None, set(), False

    def _flush():
        if cur_id and not obsolete:
            go_terms[cur_id] = {"name": cur_name or "", "namespace": cur_ns or ""}
            parents[cur_id]  = cur_parents

    with open(obo_path) as fh:
        for line in fh:
            line = line.rstrip()
            if line == "[Term]":
                _flush()
                cur_id, cur_name, cur_ns, cur_parents, obsolete = \
                    None, None, None, set(), False
            elif line.startswith("id: GO:"):
                cur_id = line[4:]
            elif line.startswith("name: "):
                cur_name = line[6:]
            elif line.startswith("namespace: "):
                cur_ns = line[11:]
            elif line.startswith("is_a:") or line.startswith("part_of:"):
                parts = line.split()
                if len(parts) >= 2 and parts[1].startswith("GO:"):
                    cur_parents.add(parts[1])
            elif line == "is_obsolete: true":
                obsolete = True
    _flush()

    logger.info("Parsed %d non-obsolete GO terms.", len(go_terms))
    return go_terms, parents


def get_ancestors(go_id: str, parents: dict[str, set[str]]) -> set[str]:
    """Return all ancestor GO IDs for a given term (transitive closure)."""
    ancestors: set[str] = set()
    queue = list(parents.get(go_id, []))
    while queue:
        p = queue.pop()
        if p not in ancestors:
            ancestors.add(p)
            queue.extend(parents.get(p, []))
    return ancestors


# ---------------------------------------------------------------------------
# GAF parsing
# ---------------------------------------------------------------------------

def parse_gaf(
    gaf_path: str,
    exclude_evidence: set[str] | None = None,
) -> dict[str, set[str]]:
    """Parse a GAF file and return {gene_id: set of GO IDs}.

    Header lines (starting with '!') are skipped.
    Negated annotations (qualifier contains 'NOT') are excluded.

    Args:
        gaf_path:         Path to the GAF annotation file.
        exclude_evidence: Set of evidence codes to exclude (e.g. {'IEA'}).

    Returns:
        Dict mapping gene_id (field 2) to a set of directly annotated GO IDs.
    """
    gene_to_go: dict[str, set[str]] = defaultdict(set)
    exclude_evidence = exclude_evidence or set()
    n_rows = n_skipped_neg = n_skipped_ev = 0

    with open(gaf_path) as fh:
        for line in fh:
            if line.startswith("!"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            n_rows += 1

            qualifier = fields[_COL_QUALIFIER].upper()
            if "NOT" in qualifier:
                n_skipped_neg += 1
                continue

            ev_code = fields[_COL_EVIDENCE_CODE].strip()
            if ev_code in exclude_evidence:
                n_skipped_ev += 1
                continue

            gene_id = fields[_COL_GENE_ID].strip()
            go_id   = fields[_COL_GO_ID].strip()
            if gene_id and go_id.startswith("GO:"):
                gene_to_go[gene_id].add(go_id)

    logger.info(
        "GAF: %d annotation rows; negated excluded: %d; "
        "evidence-filtered: %d; unique genes: %d.",
        n_rows, n_skipped_neg, n_skipped_ev, len(gene_to_go),
    )
    return dict(gene_to_go)


# ---------------------------------------------------------------------------
# GMT building
# ---------------------------------------------------------------------------

def build_and_write_gmt(
    gene_to_go: dict[str, set[str]],
    go_terms: dict[str, dict[str, str]],
    parents: dict[str, set[str]],
    namespace: str,
    out_path: str,
    propagate: bool = True,
) -> None:
    """Build and write a GMT file for one GO namespace.

    Args:
        gene_to_go: gene_id → set of directly annotated GO IDs.
        go_terms:   GO ID → {'name', 'namespace'}.
        parents:    GO ID → set of direct parent GO IDs.
        namespace:  Namespace key (e.g. 'GO_BP').
        out_path:   Output GMT file path.
        propagate:  If True, propagate annotations up the GO hierarchy
                    (true-path rule) so parent terms also receive the gene.
    """
    target_ns = NAMESPACE_MAP[namespace]
    go_to_genes: dict[str, tuple[str, set[str]]] = {}

    n_direct = n_propagated = n_skipped = 0
    for gene_id, go_ids in gene_to_go.items():
        expanded = set(go_ids)
        if propagate:
            for go_id in go_ids:
                if go_id in go_terms:
                    expanded.update(get_ancestors(go_id, parents))

        for go_id in expanded:
            term = go_terms.get(go_id)
            if term is None:
                n_skipped += 1
                continue
            if term["namespace"] != target_ns:
                continue
            if go_id not in go_to_genes:
                go_to_genes[go_id] = (term["name"], set())
            go_to_genes[go_id][1].add(gene_id)

        n_direct     += len(go_ids)
        n_propagated += len(expanded) - len(go_ids)

    if n_skipped:
        logger.warning(
            "%d GO ID(s) not found in OBO (obsolete or malformed).", n_skipped
        )
    logger.info(
        "Namespace %s: %d GO terms, %d direct annotations, "
        "%d propagated, %d unique genes.",
        namespace, len(go_to_genes), n_direct, n_propagated,
        len({g for _, gs in go_to_genes.values() for g in gs}),
    )

    with open(out_path, "w") as fh:
        for go_id, (go_name, genes) in sorted(go_to_genes.items()):
            fh.write("%s\t%s\t%s\n" % (go_id, go_name, "\t".join(sorted(genes))))
    logger.info("GMT written to: %s", out_path)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build a GMT file from a GAF annotation file."
    )
    parser.add_argument(
        "--gaf", required=True,
        help="Path to the GAF annotation file.",
    )
    parser.add_argument(
        "--anchor", required=True,
        help="Anchor species ID used to name the output GMT files "
             "(e.g. Apis_mellifera).",
    )
    parser.add_argument(
        "--outdir", required=True,
        help="Directory where GMT files will be written.",
    )
    parser.add_argument(
        "--cache", required=True,
        help="Directory for caching the GO OBO file.",
    )
    parser.add_argument(
        "--namespaces", nargs="+", default=["GO_BP"],
        choices=list(NAMESPACE_MAP.keys()),
        help="GO namespaces to build GMTs for. "
             "Default: GO_BP. Options: GO_BP GO_MF GO_CC.",
    )
    parser.add_argument(
        "--exclude-evidence", nargs="*", default=[],
        metavar="CODE",
        help="Evidence codes to exclude (e.g. IEA ND). "
             "Default: none excluded.",
    )
    parser.add_argument(
        "--no-propagation", action="store_true", default=False,
        help="Disable annotation propagation up the GO hierarchy. "
             "By default, annotations are propagated (true-path rule).",
    )
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    os.makedirs(args.cache,  exist_ok=True)

    exclude_ev = set(args.exclude_evidence)
    if exclude_ev:
        logger.info("Excluding evidence codes: %s", sorted(exclude_ev))

    # Parse GAF
    gene_to_go = parse_gaf(args.gaf, exclude_evidence=exclude_ev)

    # Load OBO
    obo_path          = download_obo(args.cache)
    go_terms, parents = parse_obo(obo_path)

    # Build one GMT per namespace
    propagate = not args.no_propagation
    if propagate:
        logger.info("Annotation propagation: ON (true-path rule).")
    else:
        logger.info("Annotation propagation: OFF.")

    for namespace in args.namespaces:
        out_path = os.path.join(
            args.outdir, "%s_%s.gmt" % (args.anchor, namespace)
        )
        build_and_write_gmt(
            gene_to_go, go_terms, parents,
            namespace=namespace,
            out_path=out_path,
            propagate=propagate,
        )

    logger.info("Done.")


if __name__ == "__main__":
    main()
