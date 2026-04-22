import os
import statistics
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq


def filter_length(parameters, outlog, file):
    inputfile = parameters['inputdir'] + "/" + file
    if os.stat(inputfile).st_size == 0:
        return
    outfile_s_path = parameters['outdir'] + "/length_stats" + "/" + file + "_lengthstats"
    outfile_f_path = parameters['outdir'] + "/length_filter" + "/" + file + "_lengthfilter"
    outstats = open(outfile_s_path, "a")
    outfile = open(outfile_f_path, "a")
    lengths = []
    length_sum = 0
    length_data = {}
    for seq_record in SeqIO.parse(inputfile, "fasta"):
        seq = Seq(seq_record.seq)
        seq_id = seq_record.id
        lenght = len(seq)
        length_data[seq_id] = lenght
        lengths.append(lenght)
        length_sum = length_sum + lenght
    total_seqs = lengths
    mean = statistics.mean(lengths)
    median = statistics.median(lengths)
    stddev = 0
    if len(lengths) > 1:
        stddev = statistics.stdev(lengths)
    outstats.write("##Overal sequence length stats\n")
    outstats.write("Total seqs: %s\nAverage: %s\nMedian: %s\nSD: %s\n" % (total_seqs, mean, median, stddev))
    outstats.write("##Sequence lengths (sorted form smallest to largest)\n")
    outstats.write("#SequenceID\tLength\tPercentageDifFromAvg\n")
    keys = list(length_data.keys())
    values = list(length_data.values())
    sorted_value_index = np.argsort(values)
    length_data_sorted = {keys[i]: values[i] for i in sorted_value_index}
    for key in length_data_sorted:
        dist = length_data_sorted[key]/mean
        outstats.write("%s\t%s\t%s\n" % (key, length_data_sorted[key], dist))
    outstats.close()
    for seq_record in SeqIO.parse(inputfile, "fasta"):
        seq = Seq(seq_record.seq)
        greater = (mean + (2*stddev))
        smaller = (mean - (2*stddev))
        if (len(seq) < smaller) or (len(seq) > greater):
            outlog.write("Sequence %s in group %s removed for length filter step\n" % (seq_record.id, file))
        else:
            outfile.write(">%s\n%s\n" % (seq_record.id, seq))


def filter_groups(parameters, outlog):
    out = parameters['outdir'] + "/group_filter"
    if not os.path.isdir(out):
        os.mkdir(out)
    in_dir_path = parameters['outdir'] + "/length_filter"
    files = os.listdir(in_dir_path)
    for file in files:
        file_fields = file.split(".")
        path_to_file = in_dir_path + "/" + file
        outfile_path = out + "/" + file
        species_ids = []
        for seq_record in SeqIO.parse(path_to_file, "fasta"):
            ids_fields = seq_record.id.split("|")
            sid = ids_fields[0]
            species_ids.append(sid)
        uniq_ids = list(set(species_ids))
        if parameters['anchor'] not in uniq_ids:
            outlog.write("Group %s not contain any sequence of anchor species %s. Group removed from analysis\n" % (file_fields[0], parameters['anchor']))
            continue
        if len(uniq_ids) < int(parameters['min_species']):
            outlog.write("Number of species in group %s is less than minimum required %s , group removed from analysis\n" % (file_fields[0], parameters['min_species']))
        else:
            os.system('cp %s %s' % (path_to_file, outfile_path))
