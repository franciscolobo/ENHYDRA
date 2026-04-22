import re


def read_config_file(fh_project, fh_code):
    parameters = {}
    parameters['inputdir'] = read_line("inputdir", fh_project)
    parameters['outdir'] = read_line("outdir", fh_project)
    parameters['min_species'] = int(read_line("min_species", fh_project))
    parameters['anchor'] = read_line("anchor", fh_project)
    parameters['mafft'] = read_line("mafft", fh_code)
    parameters['trimal'] = read_line("trimal", fh_code)
    parameters['max_process'] = int(read_line("max_process", fh_project))
    return parameters


def read_line(parameter, file):
    for line in file:
        line = line.rstrip()
        match = re.search("^%s" % parameter, line)
        if match:
            line = re.sub('^\s*%s\s*=\s*' % parameter, '', line)
            line = re.sub('#.*$', '', line)
            line = re.sub('\s*$', '', line)
            return line
