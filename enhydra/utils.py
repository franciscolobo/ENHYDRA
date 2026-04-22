import os


def check_parameters(parameters, code_config):
    if parameters['inputdir'] == "":
        exit("No path to the inputdir was specified in your project's configuration file, please fill the parameter 'inputdir'.")
    if not os.path.isdir(parameters['inputdir']):
        exit("The path to your project's species files isn't a valid file, please check if the path in 'infile' is correct: %s" % parameters['inputdir'])
    if not os.access(parameters['inputdir'], os.R_OK):
        exit("You don't have permission to read in your project's species file: %s, please redefine your permissions." % parameters['inputdir'])

    if parameters['outdir'] == "":
        exit("No path to outdir was specified in project_config, please open this file and fill the parameter 'outdir'")
    if parameters['min_species'] == "":
        exit("No minimum number of species per group was specified in your project's configuration file, please fill the the parameter 'min_species'.")

    if parameters['anchor'] == "":
        exit("No anchor species was specified in project_config, please open this file and fill the parameter 'anchor'")
    if parameters['max_process'] == "":
        exit("Maximum number of process was not specified in project_config, please open this file and fill the parameter 'max_process'")

    if parameters['trimal'] == "":
        exit("No path to trimal was specified in code_config at %s, please open this file and fill the parameter 'trimal'" % code_config)
    if not os.path.isfile(parameters['trimal']):
        exit("The executable of trimal wasn't found in the specified path, please check if the path is correct: %s" % parameters['trimal'])
    if not os.access(parameters['trimal'], os.R_OK):
        exit("You don't have permission to execute the trimal file specified at code_config, please check permissions or replace the file")

    if parameters['mafft'] == "":
        exit("No path to mafft was specified in code_config at %s, please open this file and fill the parameter 'mafft'" % code_config)
    if not os.path.isfile(parameters['mafft']):
        exit("The executable of mafft wasn't found in the specified path, please check if the path is correct: %s" % parameters['mafft'])
    if not os.access(parameters['mafft'], os.R_OK):
        exit("You don't have permission to execute the mafft file specified at code_config, please check permissions or replace the file")
