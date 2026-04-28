import os
import logging
from Bio import SeqIO

logger = logging.getLogger(__name__)


def make_tables(alignment_dir: str, ident_dir: str, tables_dir: str, anchor: str):
    """Generate ranked output tables from alignments and identity reports.

    Args:
        alignment_dir: Directory of MAFFT alignment files.
        ident_dir:     Directory of trimAl identity report files.
        tables_dir:    Directory where output tables are written.
        anchor:        Anchor species ID used to map group → gene ID.
    """
    os.makedirs(tables_dir, exist_ok=True)
    ortho_mean = {}
    with open(os.path.join(tables_dir, "group2mean.tsv"), "a") as group2mean, \
         open(os.path.join(tables_dir, "anchor2mean.tsv"), "a") as anchor2mean, \
         open(os.path.join(tables_dir, "group2anchor.tsv"), "a") as group2anchor:
        for file in os.listdir(ident_dir):
            group_name = file.split(".")[0]
            ident_path = os.path.join(ident_dir, file)
            with open(ident_path, "r") as ident_file:
                for line in ident_file:
                    line = line.rstrip()
                    if line.startswith("## AverageIdentity"):
                        mean_percent = line.split()[-1]
                        ortho_mean[group_name] = mean_percent
                        group2mean.write("%s\t%s\n" % (group_name, mean_percent))
                        break
                else:
                    logger.warning("No AverageIdentity found in %s", ident_path)
        for file in os.listdir(alignment_dir):
            group_name = file.split(".")[0]
            seq_file = os.path.join(alignment_dir, file)
            if group_name not in ortho_mean:
                logger.warning("No identity score found for group %s, skipping.", group_name)
                continue
            for seq_record in SeqIO.parse(seq_file, "fasta"):
                ids_fields = seq_record.id.split("|")
                species, gene_id = ids_fields[0], ids_fields[1]
                if species == anchor:
                    anchor2mean.write("%s\t%s\n" % (gene_id, ortho_mean[group_name]))
                    group2anchor.write("%s\t%s\n" % (group_name, gene_id))
