from __future__ import annotations

import re
import math
from collections import Counter

__all__ = ["parse_trimal_colnumbering", "compute_column_stats"]

_GAP_CHARS = frozenset("-.")
# Reference max entropy: a uniform distribution over the 20 standard amino
# acids. Used as a fixed denominator so normalised entropy values stay
# comparable across different columns and alignments, rather than being
# rescaled relative to however many distinct residues happen to appear in
# one particular column.
_MAX_ENTROPY_BITS = math.log2(20)


def parse_trimal_colnumbering(path: str) -> set[int]:
    """Parse a trimAl -colnumbering sidecar file into retained column indices.

    The sidecar file is expected to contain a single line as isolated by
    alignment._extract_colnumbering_line(): trimAl's '#ColumnsMap' label
    (if present) followed by a comma-separated list of 0-based column
    indices from the *original* (pre-trimming) alignment that survive
    trimming, e.g.:

        #ColumnsMap\t0, 1, 2, 3, 4, 5, 6, 7, ...

    Older or different trimAl builds may omit the '#ColumnsMap' label
    entirely and just list the numbers, or use a different label —
    parsing is therefore defensive: any '#ColumnsMap' prefix is stripped
    if present, then every integer remaining in the line is extracted.
    Because the sidecar is expected to contain only this one isolated
    line (not the full alignment), there is no risk of this regex
    picking up unrelated digits such as sequence header gene IDs.

    Args:
        path: Path to a sidecar file written by
              alignment.run_trimal_columns() with colnumbering_dir set.

    Returns:
        Set of 0-based column indices (relative to the original,
        untrimmed alignment) that were retained after trimming. A file
        with no parseable integers (e.g. trimAl produced no output for a
        degenerate alignment) returns an empty set.

    Raises:
        FileNotFoundError: If path does not exist.
    """
    with open(path) as fh:
        content = fh.read()

    match = re.search(r"#ColumnsMap\s*(.*)", content, re.DOTALL)
    if match:
        content = match.group(1)

    return {int(tok) for tok in re.findall(r"\d+", content)}


def compute_column_stats(seqs: list[str]) -> tuple[list[float], list[float]]:
    """Compute per-column identity fraction and normalised Shannon entropy.

    For each alignment column:
      - Identity fraction: count of the most common non-gap residue,
        divided by the total number of non-gap residues in that column.
        A tie for most common still just uses that shared count (which
        residue "wins" the tie doesn't affect the fraction).
      - Normalised entropy: Shannon entropy (base 2) of the non-gap
        residue frequency distribution in that column, divided by
        log2(20) so values fall in [0, 1] and are comparable across
        columns/alignments regardless of how many distinct residues
        happen to appear in any single column. Clipped to 1.0 as a
        defensive measure in the unusual case where more than 20
        distinct symbols appear in one column (e.g. heavy use of
        ambiguity codes), which would otherwise push the raw ratio
        above 1.0.

    Both '-' and '.' are treated as gap characters and excluded from
    both calculations. All other characters (including ambiguity codes
    such as 'X') are treated as ordinary residues — this is a
    deliberate simplification, not an oversight. Residues are
    case-normalised (uppercased) before counting.

    A column with no non-gap residues at all (all sequences have a gap
    there) has no conservation signal to measure. Rather than return
    NaN or None — which would force every caller to special-case
    missing values (e.g. when rendering a bounded color scale) — both
    stats are defined as 0.0 for an all-gap column.

    Args:
        seqs: Aligned sequences (equal length, gap-padded). An empty
              list returns ([], []).

    Returns:
        Tuple of (identity_fractions, entropies), each a list with one
        entry per alignment column, in column order.

    Raises:
        ValueError: If the sequences are not all the same length (i.e.
                    not a valid alignment).
    """
    if not seqs:
        return [], []

    width = len(seqs[0])
    lengths = {len(s) for s in seqs}
    if len(lengths) > 1:
        raise ValueError(
            "All sequences must be the same length to compute column "
            "statistics (this must be an alignment, not raw unaligned "
            "sequences). Got lengths: %s" % sorted(lengths)
        )

    identities: list[float] = []
    entropies: list[float] = []

    for col_idx in range(width):
        residues = [
            s[col_idx].upper() for s in seqs if s[col_idx] not in _GAP_CHARS
        ]
        if not residues:
            identities.append(0.0)
            entropies.append(0.0)
            continue

        counts = Counter(residues)
        n = len(residues)

        identities.append(max(counts.values()) / n)

        entropy_bits = -sum(
            (count / n) * math.log2(count / n) for count in counts.values()
        )
        entropies.append(min(entropy_bits / _MAX_ENTROPY_BITS, 1.0))

    return identities, entropies
