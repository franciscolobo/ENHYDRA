"""
Reconstruct orthogroup FASTA files from OrthoFinder output and E. coli proteomes.

Produces three sets of per-orthogroup FASTA files:
  all/             one FASTA per OG, all species
  pathogenic/      one FASTA per OG, pathogenic species only
  non_pathogenic/  one FASTA per OG, non-pathogenic species only

Sequence headers follow ENHYDRA's expected format:
  >genomeID|proteinID

Usage:
    python tests/reconstruct_orthogroups.py \\
        --table    tests/E_coli_genomes.tsv \\
        --proteomes test_data/ecoli_proteomes/ \\
        --orthogroups <orthofinder_dir>/Orthogroups/Orthogroups.txt \\
        --outdir   test_data/ecoli_orthogroups/

Dependencies: Biopython
"""

from __future__ import annotations

import os
import re
import sys
import logging
import argparse
from Bio import SeqIO

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)]
)
logger = logging.getLogger(__name__)


# --- Parsing ---

def parse_table(table_path: str) -> dict[str, str]:
    """Parse the strain table.

    Returns:
        pathogenicity: maps genome_id → 'pathogenic' | 'non-pathogenic'
    """
    pathogenicity = {}
    with open(table_path) as fh:
        header = fh.readline().rstrip().split("\t")
        for line in fh:
            fields = line.rstrip().split("\t")
            record = dict(zip(header, fields))
            pathogenicity[record["genomeID"]] = record["Pathogenicity"]
    return pathogenicity


def parse_lcl_id(lcl_id: str) -> str | None:
    """Extract the sequence accession from a lcl-style sequence ID.

    'lcl|NZ_CP012380.1_prot_WP_000870373.1_847' → 'NZ_CP012380.1'

    Returns None if the format is not recognised.
    """
    match = re.match(r"lcl\|(.+?)_prot_", lcl_id)
    if match:
        return match.group(1)
    return None


def parse_protein_id(full_header: str, lcl_id: str) -> str:
    """Extract protein_id from a full efetch FASTA header.

    First tries to extract from the [protein_id=...] field. If not found,
    falls back to parsing the protein ID from the lcl-style sequence ID:
        lcl|NC_004431.1_prot_WP_000361612.1_1 → WP_000361612.1

    Args:
        full_header: Full sequence description string.
        lcl_id:      Raw sequence ID (lcl|... style).

    Returns:
        Protein ID string.
    """
    # Try [protein_id=...] field first
    match = re.search(r"\[protein_id=([^\]]+)\]", full_header)
    if match:
        return match.group(1)

    # Fall back to parsing from lcl ID:
    # lcl|<accession>_prot_<protein_id>_<index> → <protein_id>
    if "_prot_" in lcl_id:
        prot_part = lcl_id.split("_prot_", 1)[1]
        return "_".join(prot_part.split("_")[:-1])

    # Last resort: return the lcl_id as-is
    return lcl_id


def index_proteomes(proteomes_dir: str) -> dict[str, tuple[str, str, str]]:
    """Build a lookup index from proteome FASTAs.

    Indexes sequences by their lcl ID (as used in Orthogroups.txt).
    The genome ID is determined first by the proteome filename (authoritative),
    then by the accession2genome table as a fallback.

    Args:
        proteomes_dir:    Directory of per-strain proteome FASTAs.
        accession2genome: Maps any sequence accession → genome_id.

    Returns:
        Dictionary mapping lcl_id → (genome_id, protein_id, sequence).
    """
    index = {}
    files = [f for f in os.listdir(proteomes_dir) if f.endswith(".fa")]
    logger.info("Indexing %d proteome files...", len(files))

    for filename in files:
        # Genome ID is always the filename stem — authoritative regardless
        # of what accessions appear inside the lcl headers
        genome_id = filename.replace(".fa", "")
        path = os.path.join(proteomes_dir, filename)

        for record in SeqIO.parse(path, "fasta"):
            lcl_id = record.id
            protein_id = parse_protein_id(record.description, lcl_id)
            index[lcl_id] = (genome_id, protein_id, str(record.seq))

    logger.info("Indexed %d sequences.", len(index))
    return index


def parse_orthogroups(orthogroups_path: str) -> dict[str, list[str]]:
    """Parse Orthogroups.txt into a dict mapping OG_id → list of lcl_ids."""
    orthogroups = {}
    with open(orthogroups_path) as fh:
        for line in fh:
            line = line.rstrip()
            if not line:
                continue
            og_id, _, members = line.partition(": ")
            orthogroups[og_id] = members.split()
    logger.info("Parsed %d orthogroups.", len(orthogroups))
    return orthogroups


# --- Writing ---

def write_orthogroup_fastas(
    orthogroups: dict[str, list[str]],
    seq_index: dict[str, tuple[str, str, str]],
    pathogenicity: dict[str, str],
    outdir: str,
):
    """Write three sets of per-orthogroup FASTA files.

    Args:
        orthogroups:   OG_id → list of lcl_ids.
        seq_index:     lcl_id → (genome_id, protein_id, sequence).
        pathogenicity: genome_id → 'pathogenic' | 'non-pathogenic'.
        outdir:        Root output directory.
    """
    dirs = {
        "all":            os.path.join(outdir, "all"),
        "pathogenic":     os.path.join(outdir, "pathogenic"),
        "non_pathogenic": os.path.join(outdir, "non_pathogenic"),
    }
    for d in dirs.values():
        os.makedirs(d, exist_ok=True)

    n_written = 0
    n_skipped = 0

    for og_id, lcl_ids in orthogroups.items():
        entries = []
        for lcl_id in lcl_ids:
            if lcl_id not in seq_index:
                logger.debug("lcl_id not in index: %s", lcl_id)
                n_skipped += 1
                continue
            genome_id, protein_id, seq = seq_index[lcl_id]
            entries.append((genome_id, protein_id, seq))

        if not entries:
            continue

        # Write all
        _write_fasta(entries, os.path.join(dirs["all"], og_id))

        # Write pathogenic subset
        pathogenic_entries = [
            e for e in entries if pathogenicity.get(e[0]) == "pathogenic"
        ]
        if pathogenic_entries:
            _write_fasta(pathogenic_entries, os.path.join(dirs["pathogenic"], og_id))

        # Write non-pathogenic subset
        non_pathogenic_entries = [
            e for e in entries if pathogenicity.get(e[0]) == "non-pathogenic"
        ]
        if non_pathogenic_entries:
            _write_fasta(non_pathogenic_entries, os.path.join(dirs["non_pathogenic"], og_id))

        n_written += 1

    logger.info(
        "Written %d orthogroups. Sequences not found in index: %d.",
        n_written, n_skipped
    )


def _write_fasta(entries: list[tuple[str, str, str]], out_path: str):
    """Write a list of (genome_id, protein_id, sequence) to a FASTA file."""
    with open(out_path, "w") as fh:
        for genome_id, protein_id, seq in entries:
            fh.write(">%s|%s\n%s\n" % (genome_id, protein_id, seq))


# --- Main ---

def main():
    parser = argparse.ArgumentParser(
        description="Reconstruct orthogroup FASTAs from OrthoFinder output."
    )
    parser.add_argument("--table",       required=True, help="Strain TSV table.")
    parser.add_argument("--proteomes",   required=True, help="Directory of proteome FASTAs.")
    parser.add_argument("--orthogroups", required=True, help="Path to Orthogroups.txt.")
    parser.add_argument("--outdir",      required=True, help="Output directory.")
    args = parser.parse_args()

    logger.info("Parsing strain table...")
    pathogenicity = parse_table(args.table)
    logger.info("Loaded %d strains.", len(pathogenicity))

    logger.info("Indexing proteomes...")
    seq_index = index_proteomes(args.proteomes)

    logger.info("Parsing Orthogroups.txt...")
    orthogroups = parse_orthogroups(args.orthogroups)

    logger.info("Writing orthogroup FASTAs...")
    write_orthogroup_fastas(orthogroups, seq_index, pathogenicity, args.outdir)

    logger.info("Done. Output written to: %s", args.outdir)
    logger.info("  all/            → %s", os.path.join(args.outdir, "all"))
    logger.info("  pathogenic/     → %s", os.path.join(args.outdir, "pathogenic"))
    logger.info("  non_pathogenic/ → %s", os.path.join(args.outdir, "non_pathogenic"))


if __name__ == "__main__":
    main()
