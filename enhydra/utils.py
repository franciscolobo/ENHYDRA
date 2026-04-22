import os
from .exceptions import EnhydraConfigError, EnhydraIOError, EnhydraToolError


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

    if not parameters['trimal']:
        raise EnhydraConfigError(
            "No path to trimal was specified in code_config: %s" % code_config)
    if not os.path.isfile(parameters['trimal']):
        raise EnhydraToolError(
            "trimal executable not found at: %s" % parameters['trimal'])
    if not os.access(parameters['trimal'], os.X_OK):
        raise EnhydraToolError(
            "No execute permission for trimal at: %s" % parameters['trimal'])

    if not parameters['mafft']:
        raise EnhydraConfigError(
            "No path to mafft was specified in code_config: %s" % code_config)
    if not os.path.isfile(parameters['mafft']):
        raise EnhydraToolError(
            "mafft executable not found at: %s" % parameters['mafft'])
    if not os.access(parameters['mafft'], os.X_OK):
        raise EnhydraToolError(
            "No execute permission for mafft at: %s" % parameters['mafft'])
