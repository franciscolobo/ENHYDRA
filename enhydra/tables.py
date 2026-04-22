import os
import re
from Bio import SeqIO


def make_tables(parameters):
    tables = parameters['outdir'] + "/tables"
    if not os.path.isdir(tables):
        os.mkdir(tables)
    ortho_mean = {}
    dirpath1 = parameters['outdir'] + "/alignment"
    input_files1 = os.listdir(dirpath1)
    group2anchor_path = tables + "/group2anchor.tsv"
    group2anchor = open(group2anchor_path, "a")
    dirpath2 = parameters['outdir'] + "/ident_alignment"
    input_files2 = os.listdir(dirpath2)
    group2mean_path = tables + "/group2mean.tsv"
    anchor2mean_path = tables + "/anchor2mean.tsv"
    anchor2mean = open(anchor2mean_path, "a")
    group2mean = open(group2mean_path, "a")
    for file2 in input_files2:
        aux2 = file2.split(".")
        group_name2 = aux2[0]
        ident_files = dirpath2 + '/' + file2
        ident_file = open(ident_files, "r")
        for line in ident_file:
            line = line.rstrip()
            match = re.search("identity:", line)
            if match:
                values = line.split()
                mean_percent = values[-1]
                ortho_mean[group_name2] = mean_percent
                group2mean.write("%s\t%s\n" % (group_name2, mean_percent))
    for file1 in input_files1:
        aux = file1.split(".")
        group_name = aux[0]
        seq_file = dirpath1 + '/' + file1
        for seq_record in SeqIO.parse(seq_file, "fasta"):
            ids_fields = seq_record.id.split("|")
            specie = ids_fields[0]
            sequence_id = ids_fields[1]
            if specie == parameters['anchor']:
                anchor2mean.write("%s\t%s\n" % (sequence_id, ortho_mean[group_name]))
                group2anchor.write("%s\t%s\n" % (group_name, sequence_id))
