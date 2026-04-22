import os
import shutil
import logging
import statistics
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq

logger = logging.getLogger(__name__)


def filter_length(parameters, file):
    inputfile = parameters['inputdir'] + "/" + file
    if os.stat(inputfile).st_size == 0:
        return
    outfile_s_path = parameters['outdir'] + "/length_stats/" + file + "_lengthstats"
    outfile_f_path = parameters['outdir'] + "/length_filter/" + file + "_lengthfilter"
    lengths = []
    length_data = {}
    for seq_record in SeqIO.parse(inputfile, "fasta"):
        seq = Seq(seq_record.seq)
        seq_id = seq_record.id
        length = len(seq)
        length_data[seq_id] = length
        lengths.append(length)
    mean = statistics.mean(lengths)
    median = statistics.median(lengths)
    stddev = statistics.stdev(lengths) if len(lengths) > 1 else 0
    sorted_value_index = np.argsort(list(length_data.values()))
    keys = list(length_data.keys())
    values = list(length_data.values())
    length_data_sorted = {keys[i]: values[i] for i in sorted_value_index}
    with open(outfile_s_path, "a") as outstats:
        outstats.write("##Overall sequence length stats\n")
        outstats.write("Total seqs: %s\nAverage: %s\nMedian: %s\nSD: %s\n" % (len(lengths), mean, median, stddev))
        outstats.write("##Sequence lengths (sorted from smallest to largest)\n")
        outstats.write("#SequenceID\tLength\tPercentageDifFromAvg\n")
        for key, value in length_data_sorted.items():
            outstats.write("%s\t%s\t%s\n" % (key, value, value / mean))
    with open(outfile_f_path, "a") as outfile:
        for seq_record in SeqIO.parse(inputfile, "fasta"):
            seq = Seq(seq_record.seq)
            if (len(seq) < mean - 2 * stddev) or (len(seq) > mean + 2 * stddev):
                logger.warning("Sequence %s in group %s removed by length filter", seq_record.id, file)
            else:
                outfile.write(">%s\n%s\n" % (seq_record.id, seq))


def filter_groups(parameters):
    out = parameters['outdir'] + "/group_filter"
    if not os.path.isdir(out):
        os.mkdir(out)
    in_dir_path = parameters['outdir'] + "/length_filter"
    for file in os.listdir(in_dir_path):
        file_fields = file.split(".")
        path_to_file = in_dir_path + "/" + file
        outfile_path = out + "/" + file
        species_ids = []
        for seq_record in SeqIO.parse(path_to_file, "fasta"):
            sid = seq_record.id.split("|")[0]
            species_ids.append(sid)
        uniq_ids = list(set(species_ids))
        if parameters['anchor'] not in uniq_ids:
            logger.warning("Group %s does not contain anchor species %s. Group removed.", file_fields[0], parameters['anchor'])
            continue
        if len(uniq_ids) < int(parameters['min_species']):
            logger.warning("Group %s has fewer species than minimum required (%s). Group removed.", file_fields[0], parameters['min_species'])
        else:
            shutil.copy(path_to_file, outfile_path)
