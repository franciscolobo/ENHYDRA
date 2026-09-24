from __future__ import annotations

import os
import subprocess
import logging
import multiprocessing
from tqdm import tqdm
from .exceptions import EnhydraToolError

logger = logging.getLogger(__name__)

ALIGNERS = ("mafft", "muscle", "prank")


# ---------------------------------------------------------------------------
# Module-level worker functions (must be picklable for multiprocessing)
# ---------------------------------------------------------------------------

def _mafft_worker(args: tuple) -> None:
    input_seq, output_seq, mafft_path, mode = args
    if mode == "auto":
        cmd = [mafft_path, "--auto", "--quiet", "--anysymbol",
               "--thread", "1", input_seq]
    else:
        cmd = [mafft_path, "--%s" % mode, "--quiet", "--anysymbol",
               "--thread", "1", input_seq]
    try:
        with open(output_seq, "w") as fh:
            subprocess.run(cmd, stdout=fh, stderr=subprocess.PIPE, check=True)
    except subprocess.CalledProcessError as e:
        raise EnhydraToolError(
            "mafft failed on %s:\n%s" % (input_seq, e.stderr.decode())
        )


def _muscle_worker(args: tuple) -> None:
    input_seq, output_seq, muscle_path = args
    try:
        subprocess.run(
            [muscle_path, "-align", input_seq, "-output", output_seq,
             "-threads", "1"],
            stderr=subprocess.PIPE, check=True,
        )
    except subprocess.CalledProcessError as e:
        raise EnhydraToolError(
            "muscle failed on %s:\n%s" % (input_seq, e.stderr.decode())
        )


def _prank_worker(args: tuple) -> None:
    input_seq, output_stem, prank_path = args
    try:
        subprocess.run(
            [prank_path, "-d=%s" % input_seq, "-o=%s" % output_stem,
             "-protein", "-quiet"],
            stderr=subprocess.PIPE, check=True,
        )
        prank_out = output_stem + ".best.fas"
        final_out = output_stem + ".aln"
        if os.path.isfile(prank_out):
            os.rename(prank_out, final_out)
    except subprocess.CalledProcessError as e:
        raise EnhydraToolError(
            "prank failed on %s:\n%s" % (input_seq, e.stderr.decode())
        )


def _trimal_worker(args: tuple) -> None:
    input_seq, outident, trimal_path = args
    try:
        with open(outident, "w") as fh:
            subprocess.run(
                [trimal_path, "-sident", "-in", input_seq],
                stdout=fh, stderr=subprocess.PIPE, check=True,
            )
    except subprocess.CalledProcessError as e:
        raise EnhydraToolError(
            "trimal failed on %s:\n%s" % (input_seq, e.stderr.decode())
        )

def _extract_colnumbering_line(stdout_text: str) -> str:
    """Return the '#ColumnsMap' line from trimAl -colnumbering stdout.

    trimAl prints the '#ColumnsMap' line (retained original-column indices)
    to stdout even when -out is also given and the trimmed alignment itself
    is written to file. Isolating just this line for the sidecar file keeps
    it to a single line instead of duplicating alignment data, and avoids
    the column-index parser ever seeing unrelated digits.

    Falls back to the last non-empty line if no '#ColumnsMap'-prefixed line
    is found, in case a different trimAl version omits the label.
    """
    lines = [l for l in stdout_text.splitlines() if l.strip()]
    for line in lines:
        if line.startswith("#ColumnsMap"):
            return line
    return lines[-1] if lines else ""


def _trimal_columns_worker(args: tuple) -> None:
    input_seq, output_seq, trimal_path, trim_args, colnumbering_path = args
    cmd = [trimal_path, "-in", input_seq, "-out", output_seq] + list(trim_args)
    if colnumbering_path is not None:
        cmd = cmd + ["-colnumbering"]
    try:
        result = subprocess.run(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True
        )
    except subprocess.CalledProcessError as e:
        raise EnhydraToolError(
            "trimal column trimming failed on %s:\n%s"
            % (input_seq, e.stderr.decode())
        )
    if colnumbering_path is not None:
        colmap_line = _extract_colnumbering_line(result.stdout.decode())
        with open(colnumbering_path, "w") as fh:
            fh.write(colmap_line + "\n")


def run_trimal_columns(
    alignment_dir: str,
    trimmed_dir: str,
    trimal_path: str,
    trim_args: list[str],
    n_proc: int = 1,
    show_progress: bool = False,
    colnumbering_dir: str | None = None,
):
    """Trim alignment columns with trimAl before identity estimation.

    Args:
        alignment_dir:    Directory of (untrimmed) alignment files.
        trimmed_dir:       Directory where column-trimmed alignments are written.
        trimal_path:       Path to the trimAl executable.
        trim_args:         trimAl CLI flags controlling the trimming mode, as
                           returned by utils.resolve_trim_args() — e.g.
                           ['-gt', '0.5'], ['-strict'], ['-strictplus'], or
                           ['-automated1'].
        n_proc:            Number of parallel worker processes.
        show_progress:     Show a tqdm progress bar.
        colnumbering_dir:  If given, also capture -colnumbering output per
                           group as a sidecar file
                           '<colnumbering_dir>/<file>.colnumbering', recording
                           which original-alignment columns survived
                           trimming (for later use masking trimmed columns
                           when rendering an alignment). This runs alongside
                           the existing -out trimming call rather than as a
                           second trimAl invocation, since trimAl supports
                           both flags together. Defaults to None (no
                           colnumbering capture), preserving prior behaviour
                           for callers that don't need it.
    """
    os.makedirs(trimmed_dir, exist_ok=True)
    if colnumbering_dir is not None:
        os.makedirs(colnumbering_dir, exist_ok=True)

    args_list = [
        (os.path.join(alignment_dir, f),
         os.path.join(trimmed_dir, f),
         trimal_path, trim_args,
         os.path.join(colnumbering_dir, f + ".colnumbering")
         if colnumbering_dir is not None else None)
        for f in os.listdir(alignment_dir)
    ]
    _run_pool(_trimal_columns_worker, args_list, n_proc, show_progress)

# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def _run_pool(worker, args_list: list, n_proc: int, show_progress: bool) -> None:
    """Run worker over args_list in parallel with optional tqdm progress.

    Uses an explicit try/finally rather than the Pool context manager so that
    pool.terminate() is guaranteed to run — and in-flight worker processes are
    killed — whether the tqdm iteration completes normally, raises an
    exception, or is interrupted externally.
    """
    pool = multiprocessing.Pool(processes=n_proc)
    try:
        list(tqdm(
            pool.imap_unordered(worker, args_list),
            total=len(args_list),
            desc="  groups", unit="group",
            leave=False, disable=not show_progress,
        ))
    except Exception as e:
        raise EnhydraToolError(str(e)) from e
    finally:
        pool.terminate()
        pool.join()


def run_mafft(
    group_filter_dir: str,
    alignment_dir: str,
    mafft_path: str,
    n_proc: int = 1,
    mode: str = "auto",
    show_progress: bool = False,
):
    """Align each group FASTA in parallel using MAFFT.

    Args:
        group_filter_dir: Directory of group-filtered FASTA files.
        alignment_dir:    Directory where alignments are written.
        mafft_path:       Path to the MAFFT executable.
        n_proc:           Number of parallel worker processes.
        mode:             MAFFT alignment mode (auto | linsi | ginsi | einsi |
                          fftnsi | fftns | nwnsi | nwns). Default: auto.
        show_progress:    Show a tqdm progress bar.
    """
    os.makedirs(alignment_dir, exist_ok=True)
    args_list = [
        (os.path.join(group_filter_dir, f),
         os.path.join(alignment_dir, f + ".aln"),
         mafft_path, mode)
        for f in os.listdir(group_filter_dir)
    ]
    _run_pool(_mafft_worker, args_list, n_proc, show_progress)


def run_muscle(
    group_filter_dir: str,
    alignment_dir: str,
    muscle_path: str,
    n_proc: int = 1,
    show_progress: bool = False,
):
    """Align each group FASTA in parallel using MUSCLE5.

    Args:
        group_filter_dir: Directory of group-filtered FASTA files.
        alignment_dir:    Directory where alignments are written.
        muscle_path:      Path to the MUSCLE executable.
        n_proc:           Number of parallel worker processes.
        show_progress:    Show a tqdm progress bar.
    """
    os.makedirs(alignment_dir, exist_ok=True)
    args_list = [
        (os.path.join(group_filter_dir, f),
         os.path.join(alignment_dir, f + ".aln"),
         muscle_path)
        for f in os.listdir(group_filter_dir)
    ]
    _run_pool(_muscle_worker, args_list, n_proc, show_progress)


def run_prank(
    group_filter_dir: str,
    alignment_dir: str,
    prank_path: str,
    n_proc: int = 1,
    show_progress: bool = False,
):
    """Align each group FASTA in parallel using PRANK.

    Args:
        group_filter_dir: Directory of group-filtered FASTA files.
        alignment_dir:    Directory where alignments are written.
        prank_path:       Path to the PRANK executable.
        n_proc:           Number of parallel worker processes.
        show_progress:    Show a tqdm progress bar.
    """
    os.makedirs(alignment_dir, exist_ok=True)
    args_list = [
        (os.path.join(group_filter_dir, f),
         os.path.join(alignment_dir, f),
         prank_path)
        for f in os.listdir(group_filter_dir)
    ]
    _run_pool(_prank_worker, args_list, n_proc, show_progress)


def run_trimal(
    alignment_dir: str,
    ident_dir: str,
    trimal_path: str,
    n_proc: int = 1,
    show_progress: bool = False,
):
    """Compute alignment identity for each alignment in parallel using trimAl.

    Args:
        alignment_dir: Directory of alignment files.
        ident_dir:     Directory where identity reports are written.
        trimal_path:   Path to the trimAl executable.
        n_proc:        Number of parallel worker processes.
        show_progress: Show a tqdm progress bar.
    """
    os.makedirs(ident_dir, exist_ok=True)
    args_list = [
        (os.path.join(alignment_dir, f),
         os.path.join(ident_dir, f + ".ident"),
         trimal_path)
        for f in os.listdir(alignment_dir)
    ]
    _run_pool(_trimal_worker, args_list, n_proc, show_progress)


def run_aligner(
    group_filter_dir: str,
    alignment_dir: str,
    aligner: str,
    parameters: dict,
    show_progress: bool = False,
):
    """Dispatch alignment to the configured aligner.

    Args:
        group_filter_dir: Directory of group-filtered FASTA files.
        alignment_dir:    Directory where alignments are written.
        aligner:          Aligner name ('mafft', 'muscle', 'prank').
        parameters:       Full parameters dict from config.
        show_progress:    Show a tqdm progress bar.
    """
    if aligner not in ALIGNERS:
        raise ValueError(
            "Unknown aligner '%s'. Choose from: %s" % (aligner, ALIGNERS)
        )

    n_proc = parameters['max_process']

    if aligner == "mafft":
        if not parameters.get('mafft'):
            raise EnhydraToolError("No path to mafft specified in code config.")
        run_mafft(
            group_filter_dir=group_filter_dir,
            alignment_dir=alignment_dir,
            mafft_path=parameters['mafft'],
            n_proc=n_proc,
            mode=parameters.get('mafft_mode', 'auto'),
            show_progress=show_progress,
        )
    elif aligner == "muscle":
        if not parameters.get('muscle'):
            raise EnhydraToolError("No path to muscle specified in code config.")
        run_muscle(
            group_filter_dir=group_filter_dir,
            alignment_dir=alignment_dir,
            muscle_path=parameters['muscle'],
            n_proc=n_proc,
            show_progress=show_progress,
        )
    elif aligner == "prank":
        if not parameters.get('prank'):
            raise EnhydraToolError("No path to prank specified in code config.")
        run_prank(
            group_filter_dir=group_filter_dir,
            alignment_dir=alignment_dir,
            prank_path=parameters['prank'],
            n_proc=n_proc,
            show_progress=show_progress,
        )

    logger.info("Alignment complete using %s.", aligner)
