import os
import random
import shutil
import logging
import statistics
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq

logger = logging.getLogger(__name__)

PARALOG_MODES = ("all", "remove", "longest")


def filter_length(input_path: str, length_stats_dir: str, length_filter_dir: str):
    """Filter sequences in a single FASTA file by length (mean ± 2SD).

    Args:
        input_path:        Path to the input FASTA file.
        length_stats_dir:  Directory where per-group length stats are written.
        length_filter_dir: Directory where filtered FASTA files are written.
    """
    if os.stat(input_path).st_size == 0:
        return
    file = os.path.basename(input_path)
    outfile_s_path = os.path.join(length_stats_dir, file + "_lengthstats")
    outfile_f_path = os.path.join(length_filter_dir, file + "_lengthfilter")
    lengths = []
    length_data = {}
    for seq_record in SeqIO.parse(input_path, "fasta"):
        length = len(seq_record.seq)
        length_data[seq_record.id] = length
        lengths.append(length)
    mean = statistics.mean(lengths)
    median = statistics.median(lengths)
    stddev = statistics.stdev(lengths) if len(lengths) > 1 else 0
    sorted_idx = np.argsort(list(length_data.values()))
    keys, values = list(length_data.keys()), list(length_data.values())
    length_data_sorted = {keys[i]: values[i] for i in sorted_idx}
    with open(outfile_s_path, "a") as outstats:
        outstats.write("##Overall sequence length stats\n")
        outstats.write("Total seqs: %s\nAverage: %s\nMedian: %s\nSD: %s\n" % (len(lengths), mean, median, stddev))
        outstats.write("##Sequence lengths (sorted from smallest to largest)\n")
        outstats.write("#SequenceID\tLength\tPercentageDifFromAvg\n")
        for key, value in length_data_sorted.items():
            outstats.write("%s\t%s\t%s\n" % (key, value, value / mean))
    with open(outfile_f_path, "a") as outfile:
        for seq_record in SeqIO.parse(input_path, "fasta"):
            seq = seq_record.seq
            if (len(seq) < mean - 2 * stddev) or (len(seq) > mean + 2 * stddev):
                logger.warning("Sequence %s in group %s removed by length filter", seq_record.id, file)
            else:
                outfile.write(">%s\n%s\n" % (seq_record.id, seq))


def _resolve_paralogs_longest(records: list) -> list:
    """Keep one sequence per species — the longest, with random tiebreaking.

    Args:
        records: List of SeqRecord objects.

    Returns:
        Filtered list with at most one record per species.
    """
    best: dict[str, object] = {}  # species_id → SeqRecord
    for record in records:
        species_id = record.id.split("|")[0]
        current = best.get(species_id)
        if current is None:
            best[species_id] = record
        else:
            current_len = len(current.seq)
            new_len = len(record.seq)
            if new_len > current_len or (new_len == current_len and random.random() < 0.5):
                best[species_id] = record
    return list(best.values())


def filter_groups(
    length_filter_dir: str,
    group_filter_dir: str,
    anchor: str,
    min_species: int,
    paralog_mode: str = "all",
):
    """Filter groups lacking the anchor species or below the minimum species count.

    Args:
        length_filter_dir: Directory of length-filtered FASTA files.
        group_filter_dir:  Directory where passing groups are written.
        anchor:            Species ID that must be present in each group.
        min_species:       Minimum number of distinct species required.
        paralog_mode:      How to handle paralogs (multiple sequences per species):
                           'all'     — keep all sequences (default).
                           'remove'  — discard any group containing paralogs.
                           'longest' — keep only the longest sequence per species,
                                       breaking ties at random.
    """
    if paralog_mode not in PARALOG_MODES:
        raise ValueError(
            "Invalid paralog_mode '%s'. Choose from: %s" % (paralog_mode, PARALOG_MODES)
        )

    os.makedirs(group_filter_dir, exist_ok=True)
    for file in os.listdir(length_filter_dir):
        group_name = file.split(".")[0]
        path_to_file = os.path.join(length_filter_dir, file)
        outfile_path = os.path.join(group_filter_dir, file)

        records = list(SeqIO.parse(path_to_file, "fasta"))
        if len(records) < 2:
            logger.warning(
                "Group %s has fewer than 2 sequences. Group removed.", group_name
            )
            continue

        species_ids = [r.id.split("|")[0] for r in records]
        uniq_ids = set(species_ids)

        if anchor not in uniq_ids:
            logger.warning(
                "Group %s does not contain anchor species %s. Group removed.",
                group_name, anchor
            )
            continue

        if len(uniq_ids) < min_species:
            logger.warning(
                "Group %s has fewer species than minimum required (%s). Group removed.",
                group_name, min_species
            )
            continue

        # Paralog handling
        has_paralogs = len(species_ids) > len(uniq_ids)
        if has_paralogs:
            if paralog_mode == "remove":
                logger.warning(
                    "Group %s contains paralogs and will be removed (--paralogs remove).",
                    group_name
                )
                continue
            elif paralog_mode == "longest":
                records = _resolve_paralogs_longest(records)
                logger.info(
                    "Group %s: paralogs resolved by keeping longest sequence per species.",
                    group_name
                )

        with open(outfile_path, "w") as out_fh:
            for record in records:
                out_fh.write(">%s\n%s\n" % (record.id, record.seq))
