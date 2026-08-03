"""
Subset ENHYDRA homolog group FASTAs by hymenopteran clade.

Reads a table that maps species IDs to clade names (Ants, Bees, Sawflies,
Wasps) and produces one output directory per clade, each containing the
subset of group FASTAs restricted to species from that clade.

Output layout:
    outdir/
    ├── Ants/
    │   ├── HGDOG00001at7399
    │   ├── HGDOG00002at7399
    │   └── ...
    ├── Bees/
    ├── Sawflies/
    └── Wasps/

Each group file within a clade directory contains only sequences whose
species ID belongs to that clade.  Groups with fewer than --min-species
distinct species after subsetting are omitted.

Clade table format (tab-separated, optional header):
    species_id    clade
    Apis_mellifera    Bees
    Bombus_terrestris    Bees
    Acromyrmex_echinatior    Ants
    ...

The species_id must match the prefix used in the FASTA sequence headers,
i.e. the part before '|' (as written by prepare_hgd.py).

Usage:
    python tests/subset_hgd_groups.py \\
        --inputdir  /path/to/enhydra_input/ \\
        --table     /path/to/hymenoptera_clades.tsv \\
        --outdir    /path/to/enhydra_subsets/ \\
        --min-species 4

    # If the table has no header row:
    python tests/subset_hgd_groups.py \\
        --inputdir  /path/to/enhydra_input/ \\
        --table     /path/to/hymenoptera_clades.tsv \\
        --outdir    /path/to/enhydra_subsets/ \\
        --no-header

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
# Step 1: Parse the clade table
# ---------------------------------------------------------------------------

def parse_clade_table(
    table_path: str,
    has_header: bool = True,
    species_col: int = 0,
    clade_col: int = 1,
) -> tuple[dict[str, str], dict[str, set[str]]]:
    """Parse a tab-separated clade table.

    Args:
        table_path:  Path to the table file.
        has_header:  Whether the first non-comment line is a header.
        species_col: 0-based column index for the species ID.
        clade_col:   0-based column index for the clade name.

    Returns:
        Tuple of:
            species_to_clade: {species_id: clade_name}
            clade_to_species: {clade_name: set of species_ids}
    """
    species_to_clade: dict[str, str]       = {}
    clade_to_species: dict[str, set[str]]  = defaultdict(set)
    first_data_line = True

    with open(table_path) as fh:
        for lineno, line in enumerate(fh, 1):
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            if has_header and first_data_line:
                first_data_line = False
                logger.info("Skipping header: %s", line)
                continue
            first_data_line = False

            fields = line.split("\t")
            if len(fields) <= max(species_col, clade_col):
                logger.warning("Line %d: too few fields — skipping.", lineno)
                continue

            sid   = fields[species_col].strip()
            clade = fields[clade_col].strip()
            if not sid or not clade:
                continue

            if sid in species_to_clade and species_to_clade[sid] != clade:
                logger.warning(
                    "Species '%s' assigned to multiple clades "
                    "('%s' and '%s') — keeping first assignment.",
                    sid, species_to_clade[sid], clade,
                )
                continue

            species_to_clade[sid] = clade
            clade_to_species[clade].add(sid)

    if not species_to_clade:
        sys.exit("No species–clade assignments parsed from: %s" % table_path)

    clades = sorted(clade_to_species)
    logger.info(
        "Parsed %d species across %d clade(s): %s",
        len(species_to_clade), len(clades), ", ".join(clades),
    )
    for clade in clades:
        logger.info("  %s: %d species", clade, len(clade_to_species[clade]))

    return dict(species_to_clade), dict(clade_to_species)


# ---------------------------------------------------------------------------
# Step 2: Subset and write group FASTAs
# ---------------------------------------------------------------------------

def subset_groups(
    inputdir: str,
    clade_to_species: dict[str, set[str]],
    outdir: str,
    min_species: int = 4,
) -> None:
    """Read each group FASTA once and write clade-specific subsets.

    Each group file is read a single time.  For every clade, sequences
    whose species ID (the part before '|' in the FASTA header) belongs to
    that clade are collected.  Groups with fewer than min_species distinct
    species after subsetting are omitted.

    Args:
        inputdir:         Directory of ENHYDRA-format group FASTAs.
        clade_to_species: {clade_name: set of species_ids}.
        outdir:           Root output directory; one subdirectory per clade.
        min_species:      Minimum number of distinct species for a group to
                          be written to a clade directory.
    """
    # Create one output directory per clade.
    clade_dirs = {}
    for clade in clade_to_species:
        d = os.path.join(outdir, clade)
        os.makedirs(d, exist_ok=True)
        clade_dirs[clade] = d

    group_files = sorted(os.listdir(inputdir))
    logger.info("Processing %d group files ...", len(group_files))

    # Counters: n_written[clade], n_skipped[clade]
    n_written  = defaultdict(int)
    n_skipped  = defaultdict(int)

    for group_file in group_files:
        in_path = os.path.join(inputdir, group_file)
        if not os.path.isfile(in_path):
            continue

        # Parse once, bucket sequences by clade.
        clade_records: dict[str, list] = defaultdict(list)
        for record in SeqIO.parse(in_path, "fasta"):
            species_id = record.id.split("|")[0]
            for clade, species_set in clade_to_species.items():
                if species_id in species_set:
                    clade_records[clade].append(record)
                    break   # a species belongs to exactly one clade

        # Write subset for each clade.
        for clade, records in clade_records.items():
            n_distinct = len({r.id.split("|")[0] for r in records})
            if n_distinct < min_species:
                n_skipped[clade] += 1
                continue
            out_path = os.path.join(clade_dirs[clade], group_file)
            with open(out_path, "w") as fh:
                for r in records:
                    fh.write(">%s\n%s\n" % (r.id, r.seq))
            n_written[clade] += 1

    logger.info("")
    logger.info("=" * 56)
    logger.info("SUBSET RESULTS (min_species=%d)", min_species)
    logger.info("=" * 56)
    for clade in sorted(clade_to_species):
        logger.info(
            "  %-16s written: %5d   skipped: %5d",
            clade, n_written[clade], n_skipped[clade],
        )
    logger.info("=" * 56)
    logger.info("Output: %s", outdir)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Subset ENHYDRA group FASTAs by hymenopteran clade."
    )
    parser.add_argument(
        "--inputdir", required=True,
        help="Directory of ENHYDRA-format group FASTAs (output of prepare_hgd.py).",
    )
    parser.add_argument(
        "--table", required=True,
        help="Tab-separated table mapping species_id to clade name.",
    )
    parser.add_argument(
        "--outdir", required=True,
        help="Root output directory; one subdirectory is created per clade.",
    )
    parser.add_argument(
        "--min-species", type=int, default=4,
        help="Minimum number of distinct species required to write a group "
             "to a clade directory (default: 4).",
    )
    parser.add_argument(
        "--no-header", action="store_true", default=False,
        help="Treat the first line of the table as data (no header row).",
    )
    parser.add_argument(
        "--species-col", type=int, default=0,
        help="0-based column index for species ID (default: 0).",
    )
    parser.add_argument(
        "--clade-col", type=int, default=1,
        help="0-based column index for clade name (default: 1).",
    )
    args = parser.parse_args()

    species_to_clade, clade_to_species = parse_clade_table(
        table_path=args.table,
        has_header=not args.no_header,
        species_col=args.species_col,
        clade_col=args.clade_col,
    )

    subset_groups(
        inputdir=args.inputdir,
        clade_to_species=clade_to_species,
        outdir=args.outdir,
        min_species=args.min_species,
    )
    logger.info("Done.")


if __name__ == "__main__":
    main()
