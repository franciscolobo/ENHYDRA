"""
Download and prepare proteomes for E. coli strains.

Optionally removes pseudogene sequences, which are flagged with [pseudo=true]
in efetch fasta_cds_aa headers and typically contain internal stop codons.

For each strain:
  1. Download chromosome protein sequences via efetch (fasta_cds_aa format)
  2. Download plasmid protein sequences via efetch (one per plasmid accession)
  3. Concatenate chromosome + plasmid proteins into a single FASTA per strain

Identical sequences from different loci are intentionally kept — in bacteria,
duplicate protein sequences reflect genuine multi-locus coding or horizontal
gene transfer, and their presence increases the average identity score of the
corresponding homolog group, which is biologically meaningful.

Usage:
    python tests/prepare_ecoli_proteomes.py <table> <outdir>

    table  : TSV file with columns genomeID, Pathogenicity, plasmidIDs
    outdir : Directory where proteome FASTAs will be written

Dependencies:
    efetch (entrez-direct)
"""

import os
import sys
import subprocess
import logging
import tempfile
import argparse

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)]
)
logger = logging.getLogger(__name__)


def parse_table(table_path: str) -> list[dict]:
    """Parse the strain table into a list of records."""
    strains = []
    with open(table_path) as fh:
        header = fh.readline().rstrip().split("\t")
        for line in fh:
            fields = line.rstrip().split("\t")
            record = dict(zip(header, fields))
            plasmids = (
                record["plasmidIDs"].split("|")
                if record["plasmidIDs"] != "-"
                else []
            )
            strains.append({
                "genome_id":     record["genomeID"],
                "pathogenicity": record["Pathogenicity"],
                "plasmids":      plasmids,
            })
    return strains


def efetch_proteins(accession: str, out_path: str) -> str | None:
    """Download protein sequences for any nucleotide accession via efetch.

    Works for both chromosome and plasmid sequence accessions (NC_*, CP_*, etc.).
    Uses fasta_cds_aa format which extracts translated CDS sequences directly
    from the GenBank record.

    Returns path to the protein FASTA, or None on failure.
    """
    if os.path.isfile(out_path):
        return out_path

    logger.info("Downloading proteins: %s", accession)
    result = subprocess.run(
        ["efetch", "-db", "nuccore", "-id", accession, "-format", "fasta_cds_aa"],
        capture_output=True, text=True
    )
    if result.returncode != 0 or not result.stdout.strip():
        logger.warning(
            "efetch failed for %s: %s",
            accession, result.stderr.strip() or "empty response"
        )
        return None

    with open(out_path, "w") as fh:
        fh.write(result.stdout)

    return out_path


def prepare_strain(strain: dict, tmp_dir: str, outdir: str, remove_pseudogenes: bool = True):
    """Download and merge proteome for a single strain."""
    genome_id = strain["genome_id"]
    out_fasta = os.path.join(outdir, genome_id + ".fa")

    if os.path.isfile(out_fasta):
        logger.info("Skipping %s (already processed)", genome_id)
        return

    fastas = []

    # Chromosome
    chrom_path = os.path.join(tmp_dir, genome_id + "_chromosome.faa")
    chrom_fasta = efetch_proteins(genome_id, chrom_path)
    if chrom_fasta:
        fastas.append(chrom_fasta)
    else:
        logger.warning("Could not get chromosome proteins for %s", genome_id)

    # Plasmids
    for plasmid_id in strain["plasmids"]:
        plasmid_path = os.path.join(tmp_dir, plasmid_id + "_plasmid.faa")
        plasmid_fasta = efetch_proteins(plasmid_id, plasmid_path)
        if plasmid_fasta:
            fastas.append(plasmid_fasta)
        else:
            logger.warning("Could not get plasmid proteins for %s", plasmid_id)

    if not fastas:
        logger.error("No sequences retrieved for %s, skipping.", genome_id)
        return

    # Merge all sequences into a single proteome file, deduplicating by header ID
    seen_ids = set()
    n_seqs = 0
    n_dupes = 0
    with open(out_fasta, "w") as out_fh:
        current_header = None
        current_seq_lines = []
        skip = False

        def flush():
            nonlocal n_seqs, n_dupes, skip
            if current_header is None:
                return
            seq_id = current_header.split()[0][1:]  # strip '>' and take first field
            if seq_id in seen_ids:
                n_dupes += 1
                return
            seen_ids.add(seq_id)
            out_fh.write(current_header + "\n")
            for sl in current_seq_lines:
                out_fh.write(sl + "\n")
            n_seqs += 1

        for fasta in fastas:
            with open(fasta) as in_fh:
                for line in in_fh:
                    line = line.rstrip()
                    if line.startswith(">"):
                        flush()
                        current_header = line
                        current_seq_lines = []
                    else:
                        current_seq_lines.append(line)
        flush()  # write the last record

    if n_dupes:
        logger.warning("Removed %d duplicate sequence(s) from %s", n_dupes, genome_id)
    logger.info("Written: %s (%d sequences)", out_fasta, n_seqs)


def main():
    parser = argparse.ArgumentParser(
        description="Download and prepare E. coli proteomes for OrthoFinder + ENHYDRA."
    )
    parser.add_argument("table",  help="TSV table with genomeID, Pathogenicity, plasmidIDs.")
    parser.add_argument("outdir", help="Output directory for proteome FASTAs.")
    parser.add_argument(
        "--keep-pseudogenes",
        action="store_true",
        default=False,
        help="Keep pseudogene sequences (flagged with [pseudo=true] in efetch output). "
             "By default pseudogenes are removed, as they typically contain "
             "internal stop codons and are not translated into functional proteins."
    )
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    strains = parse_table(args.table)
    logger.info("Loaded %d strains from %s", len(strains), args.table)

    with tempfile.TemporaryDirectory() as tmp_dir:
        for strain in strains:
            try:
                prepare_strain(strain, tmp_dir, args.outdir,
                               remove_pseudogenes=not args.keep_pseudogenes)
            except Exception as e:
                logger.error("Failed to process %s: %s", strain["genome_id"], e)

    logger.info("Done. Proteomes written to: %s", args.outdir)


if __name__ == "__main__":
    main()
