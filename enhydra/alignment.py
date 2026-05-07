from __future__ import annotations

import os
import subprocess
import logging
from .exceptions import EnhydraToolError

logger = logging.getLogger(__name__)

ALIGNERS = ("mafft", "muscle", "prank")


def run_mafft(
    group_filter_dir: str,
    alignment_dir: str,
    mafft_path: str,
    threads: int = 1,
    mode: str = "auto",
):
    """Align each group FASTA file using MAFFT.

    Args:
        group_filter_dir: Directory of group-filtered FASTA files.
        alignment_dir:    Directory where alignments are written.
        mafft_path:       Path to the MAFFT executable.
        threads:          Number of threads.
        mode:             MAFFT alignment mode (auto | linsi | ginsi | einsi |
                          fftnsi | fftns | nwnsi | nwns). Default: auto.
    """
    os.makedirs(alignment_dir, exist_ok=True)
    for file in os.listdir(group_filter_dir):
        input_seq  = os.path.join(group_filter_dir, file)
        output_seq = os.path.join(alignment_dir, file + ".aln")

        if mode == "auto":
            cmd = [mafft_path, "--auto", "--quiet",
                   "--thread", str(threads), input_seq]
        else:
            # Named modes are standalone executables or flags in newer MAFFT
            cmd = [mafft_path, "--%s" % mode, "--quiet",
                   "--thread", str(threads), input_seq]

        try:
            with open(output_seq, "w") as outfile:
                subprocess.run(cmd, stdout=outfile, stderr=subprocess.PIPE, check=True)
        except subprocess.CalledProcessError as e:
            raise EnhydraToolError("mafft failed on %s:\n%s" % (file, e.stderr.decode()))


def run_muscle(
    group_filter_dir: str,
    alignment_dir: str,
    muscle_path: str,
    threads: int = 1,
):
    """Align each group FASTA file using MUSCLE.

    Args:
        group_filter_dir: Directory of group-filtered FASTA files.
        alignment_dir:    Directory where alignments are written.
        muscle_path:      Path to the MUSCLE executable.
        threads:          Number of threads.
    """
    os.makedirs(alignment_dir, exist_ok=True)
    for file in os.listdir(group_filter_dir):
        input_seq  = os.path.join(group_filter_dir, file)
        output_seq = os.path.join(alignment_dir, file + ".aln")
        try:
            subprocess.run(
                [muscle_path, "-align", input_seq, "-output", output_seq,
                 "-threads", str(threads)],
                stderr=subprocess.PIPE, check=True
            )
        except subprocess.CalledProcessError as e:
            raise EnhydraToolError("muscle failed on %s:\n%s" % (file, e.stderr.decode()))


def run_prank(
    group_filter_dir: str,
    alignment_dir: str,
    prank_path: str,
):
    """Align each group FASTA file using PRANK.

    Args:
        group_filter_dir: Directory of group-filtered FASTA files.
        alignment_dir:    Directory where alignments are written.
        prank_path:       Path to the PRANK executable.
    """
    os.makedirs(alignment_dir, exist_ok=True)
    for file in os.listdir(group_filter_dir):
        input_seq  = os.path.join(group_filter_dir, file)
        output_stem = os.path.join(alignment_dir, file)
        try:
            subprocess.run(
                [prank_path, "-d=%s" % input_seq, "-o=%s" % output_stem,
                 "-protein", "-quiet"],
                stderr=subprocess.PIPE, check=True
            )
            # PRANK appends .best.fas to the output stem
            prank_out = output_stem + ".best.fas"
            final_out = output_stem + ".aln"
            if os.path.isfile(prank_out):
                os.rename(prank_out, final_out)
        except subprocess.CalledProcessError as e:
            raise EnhydraToolError("prank failed on %s:\n%s" % (file, e.stderr.decode()))


def run_aligner(
    group_filter_dir: str,
    alignment_dir: str,
    aligner: str,
    parameters: dict,
):
    """Dispatch alignment to the configured aligner.

    Args:
        group_filter_dir: Directory of group-filtered FASTA files.
        alignment_dir:    Directory where alignments are written.
        aligner:          Aligner name ('mafft', 'muscle', 'prank').
        parameters:       Full parameters dict from config.
    """
    if aligner not in ALIGNERS:
        raise ValueError(
            "Unknown aligner '%s'. Choose from: %s" % (aligner, ALIGNERS)
        )

    if aligner == "mafft":
        if not parameters.get('mafft'):
            raise EnhydraToolError("No path to mafft specified in code config.")
        run_mafft(
            group_filter_dir=group_filter_dir,
            alignment_dir=alignment_dir,
            mafft_path=parameters['mafft'],
            threads=parameters['max_process'],
            mode=parameters.get('mafft_mode', 'auto'),
        )
    elif aligner == "muscle":
        if not parameters.get('muscle'):
            raise EnhydraToolError("No path to muscle specified in code config.")
        run_muscle(
            group_filter_dir=group_filter_dir,
            alignment_dir=alignment_dir,
            muscle_path=parameters['muscle'],
            threads=parameters['max_process'],
        )
    elif aligner == "prank":
        if not parameters.get('prank'):
            raise EnhydraToolError("No path to prank specified in code config.")
        run_prank(
            group_filter_dir=group_filter_dir,
            alignment_dir=alignment_dir,
            prank_path=parameters['prank'],
        )

    logger.info("Alignment complete using %s.", aligner)


def run_trimal(alignment_dir: str, ident_dir: str, trimal_path: str):
    """Compute alignment identity for each alignment using trimAl."""
    os.makedirs(ident_dir, exist_ok=True)
    for file in os.listdir(alignment_dir):
        input_seq = os.path.join(alignment_dir, file)
        outident  = os.path.join(ident_dir, file + ".ident")
        try:
            with open(outident, "w") as outfile:
                subprocess.run(
                    [trimal_path, "-sident", "-in", input_seq],
                    stdout=outfile, stderr=subprocess.PIPE, check=True
                )
        except subprocess.CalledProcessError as e:
            raise EnhydraToolError("trimal failed on %s:\n%s" % (file, e.stderr.decode()))
