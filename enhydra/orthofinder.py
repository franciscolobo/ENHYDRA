from __future__ import annotations

import os
import logging
from Bio import SeqIO

from .exceptions import EnhydraIOError

logger = logging.getLogger(__name__)

_OG_SEQUENCES_SUBDIR = "Orthogroup_Sequences"


def _find_og_sequences_dir(orthofinder_dir: str) -> str:
    """Locate the Orthogroup_Sequences directory in an OrthoFinder 3 output tree.

    OrthoFinder 3 nests results under a dated run directory:
        <orthofinder_dir>/Results_<date>/Orthogroup_Sequences/

    If the user points directly at the Results directory, we also handle that.

    Args:
        orthofinder_dir: Path to the OrthoFinder output directory.

    Returns:
        Path to the Orthogroup_Sequences directory.

    Raises:
        EnhydraIOError: If the directory cannot be found.
    """
    # Case 1: user pointed directly at Results_<date>/
    direct = os.path.join(orthofinder_dir, _OG_SEQUENCES_SUBDIR)
    if os.path.isdir(direct):
        return direct

    # Case 2: user pointed at the top-level OrthoFinder output dir,
    # which contains one or more Results_<date>/ subdirectories
    candidates = []
    for entry in os.listdir(orthofinder_dir):
        if entry.startswith("Results_"):
            candidate = os.path.join(orthofinder_dir, entry, _OG_SEQUENCES_SUBDIR)
            if os.path.isdir(candidate):
                candidates.append(candidate)

    if not candidates:
        raise EnhydraIOError(
            "Could not find '%s' in OrthoFinder output at: %s\n"
            "Please point --orthofinder-dir at either the top-level OrthoFinder "
            "output directory or the Results_<date> subdirectory."
            % (_OG_SEQUENCES_SUBDIR, orthofinder_dir)
        )
    if len(candidates) > 1:
        # Use the most recent run
        candidates.sort()
        logger.warning(
            "Multiple OrthoFinder runs found, using the most recent: %s",
            candidates[-1]
        )
    return candidates[-1]


def _reformat_header(header: str) -> tuple[str, str]:
    """Parse an OrthoFinder 3 sequence header into (species_id, gene_id).

    OrthoFinder 3 prefixes gene IDs with the proteome filename stem:
        >Hsapiens_ENSG00000141510  →  species='Hsapiens', gene='ENSG00000141510'

    The split is on the first underscore only, so gene IDs containing
    underscores are preserved intact.

    Args:
        header: Sequence ID string (without the leading '>').

    Returns:
        Tuple of (species_id, gene_id).

    Raises:
        EnhydraIOError: If the header does not match the expected format.
    """
    parts = header.split("_", 1)
    if len(parts) != 2:
        raise EnhydraIOError(
            "Unexpected sequence header format in OrthoFinder output: '%s'. "
            "Expected 'speciesname_geneID'." % header
        )
    return parts[0], parts[1]


def preprocess_orthofinder(orthofinder_dir: str, inputdir: str) -> int:
    """Convert OrthoFinder 3 output into ENHYDRA-compatible input files.

    Reads protein FASTA files from OrthoFinder's Orthogroup_Sequences/
    directory, reformats sequence headers from 'speciesname_geneID' to
    'speciesID|geneID', and writes one FASTA per orthogroup to inputdir.

    Args:
        orthofinder_dir: Path to the OrthoFinder output directory (either
                         the top-level directory or a Results_<date> subdir).
        inputdir:        Path to the directory where reformatted FASTAs will
                         be written. Created if it does not exist.

    Returns:
        Number of orthogroup files written.

    Raises:
        EnhydraIOError: If the OrthoFinder output structure is not recognised.
    """
    og_dir = _find_og_sequences_dir(orthofinder_dir)
    os.makedirs(inputdir, exist_ok=True)

    files = [f for f in os.listdir(og_dir) if f.endswith(".fa")]
    if not files:
        raise EnhydraIOError(
            "No .fa files found in Orthogroup_Sequences directory: %s" % og_dir
        )

    n_written = 0
    for filename in files:
        og_id = filename.replace(".fa", "")
        in_path = os.path.join(og_dir, filename)
        out_path = os.path.join(inputdir, og_id)
        with open(out_path, "w") as out_fh:
            for seq_record in SeqIO.parse(in_path, "fasta"):
                try:
                    species_id, gene_id = _reformat_header(seq_record.id)
                except EnhydraIOError as e:
                    logger.warning("Skipping malformed header in %s: %s", filename, e)
                    continue
                out_fh.write(">%s|%s\n%s\n" % (species_id, gene_id, seq_record.seq))
        n_written += 1

    logger.info(
        "OrthoFinder preprocessing complete: %d orthogroups written to %s",
        n_written, inputdir
    )
    return n_written
