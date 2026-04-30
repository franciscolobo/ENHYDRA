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
    project = _parse_config(fh_project)
    code = _parse_config(fh_code)
    parameters = {
        'inputdir':    project.get('inputdir', ''),
        'outdir':      project.get('outdir', ''),
        'min_species': int(project.get('min_species', 0)),
        'anchor':      project.get('anchor', ''),
        'max_process': int(project.get('max_process', 1)),
        'mafft':       code.get('mafft', ''),
        'trimal':      code.get('trimal', ''),
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
