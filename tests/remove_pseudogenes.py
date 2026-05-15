"""Remove pseudogene sequences from existing proteome FASTAs in place.

Filters out any sequence whose header contains [pseudo=true].
Overwrites each .fa file in the given directory.

Usage:
    python tests/remove_pseudogenes.py <proteomes_dir>
"""

import os
import sys


def _filter_pseudogenes(path: str) -> tuple[list[str], int]:
    """Read a FASTA file and return (lines_to_keep, n_removed).

    Uses a sentinel-based loop to avoid nested functions and closures
    entirely — all state is local to this function.

    Args:
        path: Path to the FASTA file.

    Returns:
        Tuple of (list of lines to write back, number of sequences removed).
    """
    kept: list[str] = []
    n_removed       = 0
    current_header: str | None = None
    current_seq:    list[str]  = []

    with open(path) as fh:
        raw_lines = [line.rstrip() for line in fh]

    # Append a sentinel header so the final record is always flushed in-loop,
    # removing the need for a separate post-loop flush call.
    for line in raw_lines + [">__sentinel__"]:
        if line.startswith(">"):
            if current_header is not None:
                if "[pseudo=true]" in current_header:
                    n_removed += 1
                else:
                    kept.append(current_header)
                    kept.extend(current_seq)
            current_header = line
            current_seq    = []
        else:
            current_seq.append(line)

    return kept, n_removed


def main(proteomes_dir: str) -> None:
    for filename in sorted(os.listdir(proteomes_dir)):
        if not filename.endswith(".fa"):
            continue
        path = os.path.join(proteomes_dir, filename)
        kept, n_removed = _filter_pseudogenes(path)

        with open(path, "w") as fh:
            fh.write("\n".join(kept) + "\n")

        if n_removed:
            n_kept = sum(1 for line in kept if line.startswith(">"))
            print("%s: removed %d pseudogene(s), kept %d sequences"
                  % (filename, n_removed, n_kept))


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("Usage: remove_pseudogenes.py <proteomes_dir>")
    main(sys.argv[1])
