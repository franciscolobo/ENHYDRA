import os
import subprocess
from .exceptions import EnhydraToolError


def run_mafft(parameters):
    dirpath = parameters['outdir'] + '/group_filter'
    input_files = os.listdir(dirpath)
    outdir = parameters['outdir'] + '/alignment'
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    for file in input_files:
        input_seq = dirpath + "/" + file
        output_seq = outdir + "/" + file + ".aln"
        try:
            with open(output_seq, 'w') as outfile:
                subprocess.run(
                    [parameters['mafft'], '--auto', '--quiet',
                     '--thread', str(parameters['max_process']), input_seq],
                    stdout=outfile,
                    stderr=subprocess.PIPE,
                    check=True
                )
        except subprocess.CalledProcessError as e:
            raise EnhydraToolError(
                "mafft failed on %s:\n%s" % (file, e.stderr.decode())
            )


def run_trimal(parameters):
    dirpath = parameters['outdir'] + '/alignment'
    input_files = os.listdir(dirpath)
    outdir = parameters['outdir'] + '/ident_alignment'
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    for file in input_files:
        input_seq = dirpath + "/" + file
        outident = outdir + "/" + file + ".ident"
        try:
            with open(outident, 'w') as outfile:
                subprocess.run(
                    [parameters['trimal'], '-sident', '-in', input_seq],
                    stdout=outfile,
                    stderr=subprocess.PIPE,
                    check=True
                )
        except subprocess.CalledProcessError as e:
            raise EnhydraToolError(
                "trimal failed on %s:\n%s" % (file, e.stderr.decode())
            )
