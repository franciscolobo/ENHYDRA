import re


def _parse_config(fh):
    """Parse a config file handle into a dict, ignoring comments and blank lines."""
    pattern = re.compile(r'^\s*(\w+)\s*=\s*([^#]*)')
    config = {}
    for line in fh:
        m = pattern.match(line)
        if m:
            config[m.group(1)] = m.group(2).strip()
    return config


def read_config_file(fh_project, fh_code):
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
