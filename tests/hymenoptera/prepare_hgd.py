"""
Prepare ENHYDRA input from Hymenoptera Genome Database (HGD) data.

Converts HGD proteomes and homolog group table into a directory of per-group
FASTA files compatible with ENHYDRA's expected input format:

    >speciesID|geneID
    SEQUENCE

where speciesID is derived from the proteome filename as "Genus_species"
(e.g. "Apis_mellifera") and geneID is the gene identifier from the homolog
table (the third tab-separated field).

Homolog table format (tab-separated, no header):
    internal_id  group_id  gene_id  gene_id  species_name  taxon_id  root_taxon  -  -

Proteome FASTA header format (space-separated fields after '>'):
    protein_id  gene_id  transcript_id
    e.g. >NP_001010975.1 406075 NM_001010975.1

Usage:
    python tests/prepare_hgd.py \\
        --proteomes  /path/to/hgd_proteomes/ \\
        --homologs   /path/to/hgd_homologs.txt \\
        --outdir     /path/to/enhydra_input/

Dependencies: Biopython
"""

from __future__ import annotations

import os
import sys
import logging
import argparse
from collections import defaultdict

from Bio import SeqIO

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Step 1: Discover proteome files and derive species IDs
# ---------------------------------------------------------------------------

def _species_id_from_filename(filename: str) -> str:
    """Derive a Genus_species ID from an HGD proteome filename.

    Examples:
        Apis_mellifera_GCF_003254395.2_Amel_HAv3.1_AR104_RefSeq_protein.fa
        → 'Apis_mellifera'
    """
    stem = filename.replace(".fa", "").replace(".fasta", "").replace(".faa", "")
    parts = stem.split("_")
    return "%s_%s" % (parts[0], parts[1])


def _species_name_to_id(species_name: str) -> str:
    """Convert a species name from the homolog table to a species ID.

    'Apis mellifera' → 'Apis_mellifera'
    """
    return species_name.strip().replace(" ", "_")


def discover_proteomes(proteomes_dir: str) -> dict[str, str]:
    """Scan proteomes_dir and return {species_id: filepath}.

    Only files ending in .fa, .fasta, or .faa are considered.

    Raises:
        SystemExit: If no proteome files are found.
    """
    exts = (".fa", ".fasta", ".faa")
    mapping: dict[str, str] = {}
    for fname in sorted(os.listdir(proteomes_dir)):
        if not any(fname.endswith(ext) for ext in exts):
            continue
        sid  = _species_id_from_filename(fname)
        path = os.path.join(proteomes_dir, fname)
        if sid in mapping:
            logger.warning(
                "Duplicate species ID '%s' — skipping '%s' (keeping '%s').",
                sid, path, mapping[sid],
            )
            continue
        mapping[sid] = path
    if not mapping:
        sys.exit("No proteome FASTA files found in: %s" % proteomes_dir)
    logger.info("Discovered %d proteome files.", len(mapping))
    return mapping


# ---------------------------------------------------------------------------
# Step 2: Index all proteomes  (gene_id → sequence)
# ---------------------------------------------------------------------------

def index_proteomes(
    species_map: dict[str, str],
) -> dict[str, dict[str, str]]:
    """Parse every proteome FASTA and return {species_id: {gene_id: sequence}}.

    The gene ID is taken from the second whitespace-separated field of each
    FASTA header (index 1 after splitting on whitespace).  When multiple
    sequences share a gene ID, the longest is retained and a warning is
    emitted.
    """
    index: dict[str, dict[str, str]] = {}
    for sid, path in species_map.items():
        gene_map: dict[str, str] = {}
        n_dup = 0
        for record in SeqIO.parse(path, "fasta"):
            fields = record.description.split()
            if len(fields) < 2:
                logger.warning(
                    "Skipping malformed header in %s: '%s'",
                    os.path.basename(path), record.description,
                )
                continue
            gene_id = fields[1]
            seq     = str(record.seq)
            if gene_id in gene_map:
                n_dup += 1
                if len(seq) > len(gene_map[gene_id]):
                    gene_map[gene_id] = seq   # keep longest
            else:
                gene_map[gene_id] = seq
        if n_dup:
            logger.warning(
                "%s: %d duplicate gene ID(s) — longest sequence retained.",
                sid, n_dup,
            )
        logger.info("%s: indexed %d sequences.", sid, len(gene_map))
        index[sid] = gene_map
    return index


# ---------------------------------------------------------------------------
# Step 3: Parse the homolog group table
# ---------------------------------------------------------------------------

def parse_homolog_table(
    table_path: str,
    known_species: set[str],
) -> dict[str, list[tuple[str, str]]]:
    """Parse the HGD homolog table into {group_id: [(species_id, gene_id)]}.

    Tab-separated columns (0-indexed):
        0  internal_id
        1  group_id
        2  gene_id
        3  gene_id (duplicate)
        4  species_name
        5  taxon_id
        6  root_taxon
        7  unknown
        8  unknown

    Only rows whose species_id matches a discovered proteome are included.
    Rows with fewer than 5 fields or a '-' group_id are skipped.

    Args:
        table_path:    Path to the HGD homolog table.
        known_species: Set of species IDs from discover_proteomes().

    Returns:
        Dict mapping group_id → list of (species_id, gene_id) tuples.
    """
    groups: dict[str, list[tuple[str, str]]] = defaultdict(list)
    n_skipped_species = 0
    n_skipped_format  = 0

    with open(table_path) as fh:
        for lineno, line in enumerate(fh, 1):
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 5:
                n_skipped_format += 1
                if n_skipped_format <= 5:
                    logger.warning("Line %d: too few fields — skipping.", lineno)
                continue
            group_id    = fields[1].strip()
            gene_id     = fields[2].strip()
            species_name = fields[4].strip()
            if not group_id or group_id == "-":
                n_skipped_format += 1
                continue
            sid = _species_name_to_id(species_name)
            if sid not in known_species:
                n_skipped_species += 1
                continue
            groups[group_id].append((sid, gene_id))

    logger.info(
        "Parsed homolog table: %d groups, %d rows skipped "
        "(unmatched species: %d, malformed: %d).",
        len(groups), n_skipped_species + n_skipped_format,
        n_skipped_species, n_skipped_format,
    )
    return dict(groups)


# ---------------------------------------------------------------------------
# Step 4: Write per-group FASTA files
# ---------------------------------------------------------------------------

def write_group_fastas(
    groups: dict[str, list[tuple[str, str]]],
    index: dict[str, dict[str, str]],
    outdir: str,
) -> None:
    """Write one FASTA per homolog group to outdir.

    Each sequence header follows the ENHYDRA format:
        >speciesID|geneID

    Groups where no sequence can be resolved are skipped.

    Args:
        groups: {group_id: [(species_id, gene_id)]}
        index:  {species_id: {gene_id: sequence}}
        outdir: Output directory.
    """
    os.makedirs(outdir, exist_ok=True)

    n_written      = 0
    n_skipped      = 0
    n_missing_seqs = 0

    for group_id, members in groups.items():
        entries: list[tuple[str, str, str]] = []   # (species_id, gene_id, seq)
        for sid, gene_id in members:
            seq = index.get(sid, {}).get(gene_id)
            if seq is None:
                n_missing_seqs += 1
                logger.debug(
                    "Group %s: gene %s not found in proteome of %s — skipping sequence.",
                    group_id, gene_id, sid,
                )
            else:
                entries.append((sid, gene_id, seq))

        if not entries:
            n_skipped += 1
            logger.warning(
                "Group %s: no sequences resolved — group skipped.", group_id,
            )
            continue

        out_path = os.path.join(outdir, group_id)
        with open(out_path, "w") as fh:
            for sid, gene_id, seq in entries:
                fh.write(">%s|%s\n%s\n" % (sid, gene_id, seq))
        n_written += 1

    logger.info(
        "Written: %d groups. Skipped (no sequences): %d. "
        "Missing individual sequences: %d.",
        n_written, n_skipped, n_missing_seqs,
    )


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

def _print_summary(
    species_map: dict[str, str],
    groups: dict[str, list[tuple[str, str]]],
    outdir: str,
) -> None:
    """Log a concise summary and a few example species IDs for the user."""
    n_files = len(os.listdir(outdir)) if os.path.isdir(outdir) else 0
    example_ids = sorted(species_map.keys())[:5]
    logger.info("")
    logger.info("=" * 56)
    logger.info("SUMMARY")
    logger.info("=" * 56)
    logger.info("  Species (proteomes):   %d", len(species_map))
    logger.info("  Homolog groups parsed: %d", len(groups))
    logger.info("  Group FASTAs written:  %d  →  %s", n_files, outdir)
    logger.info("")
    logger.info("  Example species IDs (use one as 'anchor' in ENHYDRA config):")
    for sid in example_ids:
        logger.info("    %s", sid)
    logger.info("=" * 56)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Prepare ENHYDRA input from Hymenoptera Genome Database data."
    )
    parser.add_argument(
        "--proteomes", required=True,
        help="Directory containing HGD proteome FASTA files.",
    )
    parser.add_argument(
        "--homologs", required=True,
        help="Path to the HGD homolog group table (tab-separated).",
    )
    parser.add_argument(
        "--outdir", required=True,
        help="Output directory for ENHYDRA-compatible group FASTAs.",
    )
    args = parser.parse_args()

    # Step 1: discover proteomes
    species_map = discover_proteomes(args.proteomes)

    # Step 2: index all proteomes
    index = index_proteomes(species_map)

    # Step 3: parse homolog table
    groups = parse_homolog_table(args.homologs, known_species=set(species_map.keys()))

    # Step 4: write group FASTAs
    write_group_fastas(groups, index, args.outdir)

    _print_summary(species_map, groups, args.outdir)
    logger.info("Done.")


if __name__ == "__main__":
    main()
