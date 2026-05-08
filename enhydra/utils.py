import os
from .exceptions import EnhydraConfigError, EnhydraIOError, EnhydraToolError


def _check_tool(path: str, name: str, code_config: str):
    """Validate that an external tool path is specified, exists, and is executable."""
    if not path:
        raise EnhydraConfigError(
            "No path to %s was specified in code_config: %s" % (name, code_config))
    if not os.path.isfile(path):
        raise EnhydraToolError(
            "%s executable not found at: %s" % (name, path))
    if not os.access(path, os.X_OK):
        raise EnhydraToolError(
            "No execute permission for %s at: %s" % (name, path))


def check_parameters(parameters, code_config):
    if not parameters['inputdir']:
        raise EnhydraConfigError(
            "No path to inputdir was specified in your project's configuration file.")
    if not os.path.isdir(parameters['inputdir']):
        raise EnhydraIOError(
            "inputdir is not a valid directory: %s" % parameters['inputdir'])
    if not os.access(parameters['inputdir'], os.R_OK):
        raise EnhydraIOError(
            "No read permission for inputdir: %s" % parameters['inputdir'])

    if not parameters['outdir']:
        raise EnhydraConfigError(
            "No path to outdir was specified in your project's configuration file.")

    if not parameters['min_species']:
        raise EnhydraConfigError(
            "min_species was not specified in your project's configuration file.")

    if not parameters['anchor']:
        raise EnhydraConfigError(
            "No anchor species was specified in your project's configuration file.")

    if not parameters['max_process']:
        raise EnhydraConfigError(
            "max_process was not specified in your project's configuration file.")

    # trimal is always required
    _check_tool(parameters['trimal'], "trimal", code_config)

    # Only validate the configured aligner
    aligner = parameters.get('aligner', 'mafft')
    if aligner == 'mafft':
        _check_tool(parameters['mafft'], "mafft", code_config)
    elif aligner == 'muscle':
        _check_tool(parameters['muscle'], "muscle", code_config)
    elif aligner == 'prank':
        _check_tool(parameters['prank'], "prank", code_config)
    else:
        raise EnhydraConfigError(
            "Unknown aligner '%s'. Choose from: mafft, muscle, prank." % aligner)
