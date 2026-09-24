from __future__ import annotations

import re

__all__ = ["parse_trimal_colnumbering"]


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
              alignment.run_trimal_colnumbering() (or the equivalent
              step folded into run_trimal_columns()).

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
