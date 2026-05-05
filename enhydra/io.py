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
    """Read project and code config files into a parameters dict.

    Required parameters:
        inputdir, outdir, anchor, max_process (project config)
        mafft, trimal (code config)

    Optional parameters (with defaults):
        min_species   (default: 4)
        paralogs      (default: 'all')
        metric        (default: 'zscore')
        gene_sets     (default: '')
        organism      (default: '')
        min_size      (default: 5)
        max_size      (default: 500)
        permutations  (default: 1000)
        seed          (default: 42)
        fdr_threshold (default: 0.25)
        list1         (default: '')
        list2         (default: '')
        sources       (default: 'GO:BP GO:MF GO:CC KEGG REAC')
    """
    project = _parse_config(fh_project)
    code    = _parse_config(fh_code)

    parameters = {
        # Required
        'inputdir':     project.get('inputdir', ''),
        'outdir':       project.get('outdir', ''),
        'anchor':       project.get('anchor', ''),
        'max_process':  int(project.get('max_process', 1)),
        'mafft':        code.get('mafft', ''),
        'trimal':       code.get('trimal', ''),
        # Filtering
        'min_species':  int(project.get('min_species', 4)),
        'paralogs':     project.get('paralogs', 'all'),
        # Ranking
        'metric':       project.get('metric', 'zscore'),
        # Gene sets
        'gene_sets':    project.get('gene_sets', ''),
        'organism':     project.get('organism', ''),
        'min_size':     int(project.get('min_size', 5)),
        'max_size':     int(project.get('max_size', 500)),
        'permutations': int(project.get('permutations', 1000)),
        'seed':         int(project.get('seed', 42)),
        'fdr_threshold': float(project.get('fdr_threshold', 0.25)),
        # Two-list mode
        'list1':        project.get('list1', ''),
        'list2':        project.get('list2', ''),
        'sources':      project.get('sources', 'GO:BP GO:MF GO:CC KEGG REAC'),
    }
    return parameters


def read_species_list(path: str) -> list[str]:
    """Parse a species list file into a list of species IDs.

    The file should contain one species ID per line. Lines starting with
    '#' and blank lines are ignored.

    Args:
        path: Path to the species list file.

    Returns:
        List of species ID strings.

    Raises:
        FileNotFoundError: If the file does not exist.
        ValueError: If the file is empty after filtering comments and blanks.
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
