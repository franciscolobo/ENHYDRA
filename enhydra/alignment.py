import os


def run_mafft(parameters):
    dirpath = parameters['outdir'] + '/group_filter'
    input_files = os.listdir(dirpath)
    outdir = parameters['outdir'] + '/alignment'
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    for file in input_files:
        input_seq = dirpath + "/" + file
        output_seq = outdir + "/" + file + ".aln"
        os.system('%s --auto --quiet --thread %s %s > %s' % (parameters['mafft'], parameters['max_process'], input_seq, output_seq))


def run_trimal(parameters):
    dirpath = parameters['outdir'] + '/alignment'
    input_files = os.listdir(dirpath)
    outdir = parameters['outdir'] + '/ident_alignment'
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    for file in input_files:
        input_seq = dirpath + "/" + file
        outident = outdir + "/" + file + ".ident"
        os.system('%s -sident -in %s > %s' % (parameters['trimal'], input_seq, outident))
