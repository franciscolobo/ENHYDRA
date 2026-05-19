from __future__ import annotations

import re


def _parse_config(fh) -> dict:
    """Parse a config file handle into a dict, ignoring comments and blank lines."""
    pattern = re.compile(r'^\s*(\w+)\s*=\s*([^#]*)')
    config = {}
    for line in fh:
        m = pattern.match(line)
        if m:
            config[m.group(1)] = m.group(2).strip()
    return config


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
        # Two-list mode
        'list1':            project.get('list1', ''),
        'list2':            project.get('list2', ''),
        'list1_name':       project.get('list1_name', 'List 1'),
        'list2_name':       project.get('list2_name', 'List 2'),
        'sources':          project.get('sources', 'GO:BP GO:MF GO:CC KEGG REAC'),
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
