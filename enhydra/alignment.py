import os
import subprocess
from .exceptions import EnhydraToolError


def run_mafft(group_filter_dir: str, alignment_dir: str, mafft_path: str, threads: int = 1):
    """Align each group FASTA file using MAFFT.

    Args:
        group_filter_dir: Directory of group-filtered FASTA files.
        alignment_dir:    Directory where alignments are written.
        mafft_path:       Path to the MAFFT executable.
        threads:          Number of threads passed to MAFFT.
    """
    os.makedirs(alignment_dir, exist_ok=True)
    for file in os.listdir(group_filter_dir):
        input_seq = os.path.join(group_filter_dir, file)
        output_seq = os.path.join(alignment_dir, file + ".aln")
        try:
            with open(output_seq, 'w') as outfile:
                subprocess.run(
                    [mafft_path, '--auto', '--quiet', '--anysymbol', '--thread', str(threads), input_seq],
                    stdout=outfile,
                    stderr=subprocess.PIPE,
                    check=True
                )
        except subprocess.CalledProcessError as e:
            raise EnhydraToolError("mafft failed on %s:\n%s" % (file, e.stderr.decode()))


def run_trimal(alignment_dir: str, ident_dir: str, trimal_path: str):
    """Compute alignment identity for each alignment using trimAl.

    Args:
        alignment_dir: Directory of MAFFT alignment files.
        ident_dir:     Directory where identity reports are written.
        trimal_path:   Path to the trimAl executable.
    """
    os.makedirs(ident_dir, exist_ok=True)
    for file in os.listdir(alignment_dir):
        input_seq = os.path.join(alignment_dir, file)
        outident = os.path.join(ident_dir, file + ".ident")
        try:
            with open(outident, 'w') as outfile:
                subprocess.run(
                    [trimal_path, '-sident', '-in', input_seq],
                    stdout=outfile,
                    stderr=subprocess.PIPE,
                    check=True
                )
        except subprocess.CalledProcessError as e:
            raise EnhydraToolError("trimal failed on %s:\n%s" % (file, e.stderr.decode()))
