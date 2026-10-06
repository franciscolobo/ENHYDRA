from __future__ import annotations

import re
from io import StringIO

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def _parse_config(fh) -> dict:
    """Parse a config file handle into a dict, ignoring comments and blank lines."""
    pattern = re.compile(r'^\s*(\w+)\s*=\s*([^#]*)')
    config = {}
    for line in fh:
        m = pattern.match(line)
        if m:
            config[m.group(1)] = m.group(2).strip()
    return config


def _parse_bool(value: str, default: bool = False) -> bool:
    """Parse a config-file string as a boolean.

    Accepts (case-insensitive): 'true', 'yes', 'on', '1' as True, and
    'false', 'no', 'off', '0' as False. An empty or unset value returns
    default rather than raising, since these are all optional flags that
    default to off.
    """
    if value is None or value == '':
        return default
    v = str(value).strip().lower()
    if v in ('true', 'yes', 'on', '1'):
        return True
    if v in ('false', 'no', 'off', '0'):
        return False
    return default


def read_config_file(fh_project, fh_code) -> dict:
    """Read project and code config files into a parameters dict."""
    project = _parse_config(fh_project)
    code    = _parse_config(fh_code)

    parameters = {
        # Required
        'inputdir':         project.get('inputdir', ''),
        'outdir':           project.get('outdir', ''),
        'anchor':           project.get('anchor', ''),
        'max_process':      int(project.get('max_process', 1)),
        'mafft':            code.get('mafft', ''),
        'trimal':           code.get('trimal', ''),
        # Optional tool paths
        'muscle':           code.get('muscle', ''),
        'prank':            code.get('prank', ''),
        # Filtering
        'min_species':      int(project.get('min_species', 4)),
        'min_sequences':    int(project.get('min_sequences', 2)),
        'paralogs':         project.get('paralogs', 'all'),
        'length_filter_sd': float(project.get('length_filter_sd', 2.0)),
        'trim':             project.get('trim', ''),
        'divergence_filter_sd': project.get('divergence_filter_sd', ''),
        # Alignment
        'aligner':          project.get('aligner', 'mafft'),
        'mafft_mode':       project.get('mafft_mode', 'auto'),
        # Ranking
        'metric':           project.get('metric', 'zscore'),
        # Gene sets
        'gene_sets':        project.get('gene_sets', ''),
        'organism':         project.get('organism', ''),
        'min_size':         int(project.get('min_size', 5)),
        'max_size':         int(project.get('max_size', 500)),
        'permutations':     int(project.get('permutations', 1000)),
        'seed':             int(project.get('seed', 42)),
        'fdr_threshold':    float(project.get('fdr_threshold', 0.25)),
        # Report
        'top_n':            int(project.get('top_n', 20)),
        'obo_cache':        project.get('obo_cache', ''),
        'gene_list_fdr_threshold': project.get('gene_list_fdr_threshold', ''),
        # Two-list mode
        'list1':            project.get('list1', ''),
        'list2':            project.get('list2', ''),
        'list1_name':       project.get('list1_name', 'List 1'),
        'list2_name':       project.get('list2_name', 'List 2'),
        'sources':          project.get('sources', 'GO:BP GO:MF GO:CC KEGG REAC'),
        # Run-mode flags (previously CLI-only; can now also be set here)
        'orthofinder_dir':  project.get('orthofinder_dir', ''),
        'resume':           _parse_bool(project.get('resume', '')),
        'fork_from':        project.get('fork_from', ''),
        'replot':           _parse_bool(project.get('replot', '')),
        'quiet':            _parse_bool(project.get('quiet', '')),
        'all_metrics':      _parse_bool(project.get('all_metrics', '')),
    }
    return parameters


def read_species_list(path: str) -> list[str]:
    """Parse a species list file into a list of species IDs.

    Lines starting with '#' and blank lines are ignored.

    Raises:
        FileNotFoundError: If the file does not exist.
        ValueError: If the file is empty after filtering.
    """
    species = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            species.append(line)
    if not species:
        raise ValueError(
            "Species list file '%s' is empty or contains only comments." % path
        )
    return species


def parse_obo_names(obo_path: str) -> dict[str, str]:
    """Parse a GO OBO file and return a dict mapping GO ID → term name.

    Obsolete terms are excluded.

    Args:
        obo_path: Path to the go-basic.obo file.

    Returns:
        Dict of GO ID → human-readable term name.
    """
    names: dict[str, str] = {}
    current_id   = None
    current_name = None
    is_obsolete  = False
    with open(obo_path) as fh:
        for line in fh:
            line = line.rstrip()
            if line == "[Term]":
                if current_id and not is_obsolete and current_name:
                    names[current_id] = current_name
                current_id   = None
                current_name = None
                is_obsolete  = False
            elif line.startswith("id: GO:"):
                current_id = line[4:]
            elif line.startswith("name: "):
                current_name = line[6:]
            elif line == "is_obsolete: true":
                is_obsolete = True
    if current_id and not is_obsolete and current_name:
        names[current_id] = current_name
    return names


def parse_fasta_records(path: str) -> list[SeqRecord]:
    """Parse a FASTA file into a list of SeqRecord objects, tolerating any
    junk/comment lines that precede the first '>' record header.

    Equivalent to `list(Bio.SeqIO.parse(path, "fasta"))` for every
    well-formed FASTA file, with one difference: newer versions of
    Biopython's strict "fasta" parser (Biopython >= 1.85) raise
    ValueError if a file contains any non-header lines before its first
    '>' record — a pattern that occurs in practice with sequences
    exported by some external databases/tools that prepend a metadata
    banner, and that older Biopython versions silently tolerated. This
    function makes ENHYDRA tolerant of that same input regardless of
    which Biopython version is installed, by stripping everything before
    the first '>' line itself — rather than depending on Biopython's own
    newer 'fasta-pearson'/'fasta-blast' formats (mentioned in that
    ValueError's own message), which are not available on Biopython
    versions older than 1.85 and would make adopting them a version
    compatibility risk for ENHYDRA's own users.

    This deliberately targets only the specific failure mode above
    (leading junk before the first sequence). It does not attempt to
    interpret or strip comment-like lines (e.g. ';'-prefixed) that might
    appear *within* the sequence portion of the file — that is a
    materially different, more invasive change to how a file's actual
    sequence content is interpreted, not merely how tolerantly its
    (possible) leading preamble is skipped, and every ENHYDRA-internal
    consumer of a FASTA file already writes clean, comment-free files of
    its own, so this narrower scope covers every place ENHYDRA actually
    encounters this problem: files originating from outside the pipeline
    (a user's own input directory, or another tool's output).

    A file containing no '>' line at all (including a completely empty
    file, or one containing only comments) returns an empty list rather
    than raising — callers throughout the pipeline already handle "zero
    sequences found" as a normal, if usually degenerate, case (e.g.
    filtering.filter_length()'s existing `len(lengths) < 2` check, which
    already treats a group with too few sequences to filter as one to
    skip rather than an error).

    Args:
        path: Path to a FASTA file.

    Returns:
        A list of Bio.SeqRecord.SeqRecord objects, in file order. Every
        record's .id/.description/.seq have identical semantics to
        Bio.SeqIO.parse(path, "fasta")'s own records, since parsing
        itself (of the part of the file that follows the first '>') is
        delegated to Biopython entirely unchanged.

    Raises:
        FileNotFoundError: If path does not exist (matching
                           Bio.SeqIO.parse()'s own behaviour for a
                           missing file).
    """
    with open(path) as fh:
        lines = fh.readlines()

    first_header_idx = next(
        (i for i, line in enumerate(lines) if line.startswith(">")), None
    )
    if first_header_idx is None:
        return []

    cleaned = "".join(lines[first_header_idx:])
    return list(SeqIO.parse(StringIO(cleaned), "fasta"))
