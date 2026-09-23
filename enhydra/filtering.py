from __future__ import annotations

import os
import random
import logging
import statistics
from collections import Counter
import numpy as np
from Bio import SeqIO
from tqdm import tqdm

logger = logging.getLogger(__name__)

PARALOG_MODES = ("all", "remove", "longest")


def subset_groups(inputdir: str, subset_dir: str, species: list[str],
                  show_progress: bool = False):
    """Subset orthogroup FASTAs to sequences from a given species list.

    Args:
        inputdir:      Directory of input FASTA files (one per orthogroup).
        subset_dir:    Directory where subsetted FASTAs are written.
        species:       List of species IDs to retain.
        show_progress: Show a tqdm progress bar.
    """
    species_set = set(species)
    os.makedirs(subset_dir, exist_ok=True)
    n_written = 0
    n_empty   = 0
    files = os.listdir(inputdir)
    for filename in tqdm(files, desc="  groups", unit="group",
                         leave=False, disable=not show_progress):
        in_path  = os.path.join(inputdir, filename)
        out_path = os.path.join(subset_dir, filename)
        records  = [
            r for r in SeqIO.parse(in_path, "fasta")
            if r.id.split("|")[0] in species_set
        ]
        if not records:
            n_empty += 1
            continue
        with open(out_path, "w") as fh:
            for r in records:
                fh.write(">%s\n%s\n" % (r.id, r.seq))
        n_written += 1
    logger.info(
        "Subset complete: %d groups written, %d groups had no matching species.",
        n_written, n_empty,
    )


def filter_length(
    input_path: str,
    length_stats_dir: str,
    length_filter_dir: str,
    sd_multiplier: float = 2.0,
):
    """Filter sequences in a single FASTA file by length (mean ± sd_multiplier * SD).

    Groups with fewer than 2 sequences before filtering are skipped entirely —
    a single sequence cannot be aligned and would cause downstream failures.

    Both output files are opened in write mode so that re-runs produce clean
    output rather than appending duplicate entries.

    This function is normally called once per group from a multiprocessing
    pool (see cli.py), so it cannot maintain a running cross-group summary
    itself. Instead, every call — including groups skipped entirely — writes
    a per-group '_lengthstats' file recording exactly what happened to that
    group. Call aggregate_length_filter_stats() once after the pool
    finishes to consolidate all of these into summary QC tables.

    Args:
        input_path:        Path to the input FASTA file.
        length_stats_dir:  Directory where per-group length stats are written.
        length_filter_dir: Directory where filtered FASTA files are written.
        sd_multiplier:     Number of SDs from the mean beyond which sequences
                           are removed (default: 2.0).
    """
    file           = os.path.basename(input_path)
    outfile_s_path = os.path.join(length_stats_dir, file + "_lengthstats")
    outfile_f_path = os.path.join(length_filter_dir, file + "_lengthfilter")

    if os.stat(input_path).st_size == 0:
        with open(outfile_s_path, "w") as outstats:
            outstats.write("##GroupSkipped\n")
            outstats.write("reason\tempty_file\n")
        return

    lengths     = []
    length_data = {}
    for seq_record in SeqIO.parse(input_path, "fasta"):
        length = len(seq_record.seq)
        length_data[seq_record.id] = length
        lengths.append(length)

    if len(lengths) < 2:
        logger.warning(
            "Group %s has only 1 sequence before length filtering — skipping.", file
        )
        with open(outfile_s_path, "w") as outstats:
            outstats.write("##GroupSkipped\n")
            outstats.write("reason\tsingle_sequence\n")
            outstats.write("n_sequences\t%d\n" % len(lengths))
        return

    mean   = statistics.mean(lengths)
    median = statistics.median(lengths)
    stddev = statistics.stdev(lengths)   # safe: len >= 2

    lower_bound = mean - sd_multiplier * stddev
    upper_bound = mean + sd_multiplier * stddev

    sorted_idx = np.argsort(list(length_data.values()))
    keys, values = list(length_data.keys()), list(length_data.values())
    length_data_sorted = {keys[i]: values[i] for i in sorted_idx}

    with open(outfile_s_path, "w") as outstats:
        outstats.write("##Overall sequence length stats\n")
        outstats.write("Total seqs: %s\nAverage: %s\nMedian: %s\nSD: %s\n" % (
            len(lengths), mean, median, stddev))
        outstats.write("##Sequence lengths (sorted from smallest to largest)\n")
        outstats.write("#SequenceID\tLength\tPercentageDifFromAvg\tStatus\n")
        for key, value in length_data_sorted.items():
            if value < lower_bound:
                status = "removed_below_min"
            elif value > upper_bound:
                status = "removed_above_max"
            else:
                status = "kept"
            outstats.write("%s\t%s\t%s\t%s\n" % (key, value, value / mean, status))

    with open(outfile_f_path, "w") as outfile:
        for seq_record in SeqIO.parse(input_path, "fasta"):
            seq = seq_record.seq
            if (len(seq) < lower_bound) or (len(seq) > upper_bound):
                logger.warning(
                    "Sequence %s in group %s removed by length filter",
                    seq_record.id, file
                )
            else:
                outfile.write(">%s\n%s\n" % (seq_record.id, seq))


def aggregate_length_filter_stats(
    length_stats_dir: str,
    length_filter_stats_dir: str,
) -> None:
    """Consolidate per-group '_lengthstats' files into summary QC tables.

    filter_length() runs once per group (typically via a multiprocessing
    pool), so it cannot maintain a running summary across groups the way
    filter_groups() does in its own single-threaded loop. Instead, each
    filter_length() call records its own decisions in its per-group
    '_lengthstats' file, and this function is called once afterward — from
    the main process, after the pool has finished — to scan every
    '_lengthstats' file and consolidate them into two summary tables, in
    the same style as filter_groups()'s drop_reasons.tsv / species_counts.tsv:

    - skipped_groups.tsv: one row per group skipped entirely before any
      per-sequence filtering could run. reason is one of 'empty_file' or
      'single_sequence'.
    - drop_reasons.tsv: one row per individual sequence removed by the
      length filter within an otherwise-processed group. reason is one of
      'removed_below_min' or 'removed_above_max', alongside the sequence's
      length and its percentage difference from the group's mean length.

    Both files are opened in write mode, so re-running this function (e.g.
    after --resume re-triggers length filtering) produces clean output
    rather than appending duplicates.

    Args:
        length_stats_dir:        Directory of per-group '_lengthstats' files
                                 (as written by filter_length()).
        length_filter_stats_dir: Directory where the two consolidated
                                 summary files are written.
    """
    os.makedirs(length_filter_stats_dir, exist_ok=True)

    skipped_groups: list[tuple[str, str, str]] = []
    drop_reasons:   list[tuple[str, str, str, str, str]] = []

    suffix = "_lengthstats"
    for filename in os.listdir(length_stats_dir):
        if not filename.endswith(suffix):
            continue
        # Match the same group-id convention used elsewhere (e.g.
        # filter_groups()'s drop_reasons.tsv) so IDs line up across stages.
        group_name = filename[:-len(suffix)].split(".")[0]
        path       = os.path.join(length_stats_dir, filename)

        with open(path) as fh:
            lines = [l.rstrip("\n") for l in fh]

        if lines and lines[0] == "##GroupSkipped":
            reason = ""
            detail = ""
            for line in lines[1:]:
                if line.startswith("reason\t"):
                    reason = line.split("\t", 1)[1]
                elif line.startswith("n_sequences\t"):
                    detail = "n_sequences=%s" % line.split("\t", 1)[1]
            skipped_groups.append((group_name, reason, detail))
            continue

        # Otherwise: a processed group. Scan the per-sequence table for any
        # rows whose Status is not 'kept'.
        in_table = False
        for line in lines:
            if line.startswith("#SequenceID"):
                in_table = True
                continue
            if not in_table or not line.strip():
                continue
            fields = line.split("\t")
            if len(fields) < 4:
                continue
            seq_id, length, pct_diff, status = fields[0], fields[1], fields[2], fields[3]
            if status != "kept":
                drop_reasons.append((group_name, seq_id, length, pct_diff, status))

    skipped_path = os.path.join(length_filter_stats_dir, "skipped_groups.tsv")
    with open(skipped_path, "w") as fh:
        fh.write("group_id\treason\tdetail\n")
        for group_id, reason, detail in sorted(skipped_groups, key=lambda r: r[0]):
            fh.write("%s\t%s\t%s\n" % (group_id, reason, detail))
    logger.info(
        "Length filter: %d group(s) skipped entirely. Path: %s",
        len(skipped_groups), skipped_path,
    )

    drop_reasons_path = os.path.join(length_filter_stats_dir, "drop_reasons.tsv")
    with open(drop_reasons_path, "w") as fh:
        fh.write("group_id\tsequence_id\tlength\tpct_diff_from_avg\treason\n")
        for row in sorted(drop_reasons, key=lambda r: (r[0], r[1])):
            fh.write("%s\t%s\t%s\t%s\t%s\n" % row)
    logger.info(
        "Length filter: %d sequence(s) removed as outliers. Path: %s",
        len(drop_reasons), drop_reasons_path,
    )


def _resolve_paralogs_longest(records: list) -> list:
    """Keep one sequence per species — the longest, with random tiebreaking."""
    best: dict[str, object] = {}
    for record in records:
        species_id = record.id.split("|")[0]
        current    = best.get(species_id)
        if current is None:
            best[species_id] = record
        else:
            current_len = len(current.seq)
            new_len     = len(record.seq)
            if new_len > current_len or (
                new_len == current_len and random.random() < 0.5
            ):
                best[species_id] = record
    return list(best.values())


def filter_groups(
    length_filter_dir: str,
    group_filter_dir: str,
    anchor: str,
    min_species: int,
    min_sequences: int = 2,
    paralog_mode: str = "all",
    require_anchor: bool = True,
    group_stats_dir: str | None = None,
    show_progress: bool = False,
):
    """Filter groups lacking the anchor species or below the minimum species count.

    In addition to writing the passing group FASTAs to group_filter_dir, this
    also writes two small QC/reporting files into a *separate* group_stats_dir
    (mirroring how filter_length keeps its length_stats_dir separate from its
    length_filter_dir output). This separation matters: group_filter_dir is
    later read wholesale by the aligner (run_aligner/run_mafft etc.), which
    assumes every file it lists is a sequence-group FASTA to align. Writing
    stats files into group_filter_dir itself would get them picked up and
    fed into MAFFT/trimAl as if they were real groups.

    Both stats files are opened in write mode so re-runs produce clean
    output rather than appending duplicates:

    - drop_reasons.tsv: one row per group removed at this step, with
      columns group_id, reason, detail. 'reason' is one of:
      'below_min_sequences', 'missing_anchor', 'below_min_species',
      'paralogs_removed'. Groups dropped earlier (e.g. by filter_length)
      are not included here — this file only covers this step's decisions.
    - species_counts.tsv: one row per species, with the number of
      surviving (written) groups that species appears in, sorted ascending
      by count so an underrepresented genome sorts to the top. Useful for
      spotting a genome that is systematically thin (e.g. due to
      consistently failing the anchor or length checks upstream).

    Args:
        length_filter_dir: Directory of length-filtered FASTA files.
        group_filter_dir:  Directory where passing groups are written.
        anchor:            Species ID used for annotation mapping.
        min_species:       Minimum number of distinct species required.
        min_sequences:     Minimum number of sequences required (default: 2).
        paralog_mode:      How to handle paralogs: 'all', 'remove', 'longest'.
        require_anchor:    Discard groups lacking the anchor species.
        group_stats_dir:   Directory for this step's QC/reporting files
                           (drop_reasons.tsv, species_counts.tsv). Defaults
                           to '<group_filter_dir>_stats' if not given.
        show_progress:     Show a tqdm progress bar.
    """
    if paralog_mode not in PARALOG_MODES:
        raise ValueError(
            "Invalid paralog_mode '%s'. Choose from: %s" % (paralog_mode, PARALOG_MODES)
        )

    if group_stats_dir is None:
        group_stats_dir = group_filter_dir.rstrip("/\\") + "_stats"

    os.makedirs(group_filter_dir, exist_ok=True)
    os.makedirs(group_stats_dir, exist_ok=True)
    files = os.listdir(length_filter_dir)

    drop_reasons:  list[tuple[str, str, str]] = []
    species_counts: Counter = Counter()

    for file in tqdm(files, desc="  groups", unit="group",
                     leave=False, disable=not show_progress):
        group_name    = file.split(".")[0]
        path_to_file  = os.path.join(length_filter_dir, file)
        outfile_path  = os.path.join(group_filter_dir, file)

        records = list(SeqIO.parse(path_to_file, "fasta"))
        if len(records) < min_sequences:
            logger.warning(
                "Group %s has fewer than %d sequences. Group removed.",
                group_name, min_sequences,
            )
            drop_reasons.append((
                group_name, "below_min_sequences",
                "n_sequences=%d;min_required=%d" % (len(records), min_sequences),
            ))
            continue

        species_ids = [r.id.split("|")[0] for r in records]
        uniq_ids    = set(species_ids)

        if require_anchor and anchor not in uniq_ids:
            logger.warning(
                "Group %s does not contain anchor species %s. Group removed.",
                group_name, anchor,
            )
            drop_reasons.append((
                group_name, "missing_anchor", "anchor=%s" % anchor,
            ))
            continue

        if len(uniq_ids) < min_species:
            logger.warning(
                "Group %s has fewer species than minimum required (%s). Group removed.",
                group_name, min_species,
            )
            drop_reasons.append((
                group_name, "below_min_species",
                "n_species=%d;min_required=%d" % (len(uniq_ids), min_species),
            ))
            continue

        has_paralogs = len(species_ids) > len(uniq_ids)
        if has_paralogs:
            if paralog_mode == "remove":
                logger.warning(
                    "Group %s contains paralogs and will be removed (--paralogs remove).",
                    group_name,
                )
                drop_reasons.append((
                    group_name, "paralogs_removed",
                    "n_species=%d;n_sequences=%d" % (len(uniq_ids), len(species_ids)),
                ))
                continue
            elif paralog_mode == "longest":
                records = _resolve_paralogs_longest(records)
                logger.info(
                    "Group %s: paralogs resolved by keeping longest per species.",
                    group_name,
                )

        with open(outfile_path, "w") as out_fh:
            for record in records:
                out_fh.write(">%s\n%s\n" % (record.id, record.seq))

        # Species membership is unaffected by paralog resolution (longest
        # mode keeps exactly one sequence per species already present), so
        # uniq_ids correctly reflects the species composition of what was
        # just written regardless of paralog_mode.
        species_counts.update(uniq_ids)

    drop_reasons_path = os.path.join(group_stats_dir, "drop_reasons.tsv")
    with open(drop_reasons_path, "w") as fh:
        fh.write("group_id\treason\tdetail\n")
        for group_id, reason, detail in sorted(drop_reasons):
            fh.write("%s\t%s\t%s\n" % (group_id, reason, detail))
    logger.info(
        "Group filter drop reasons written: %d group(s) removed. Path: %s",
        len(drop_reasons), drop_reasons_path,
    )

    species_counts_path = os.path.join(group_stats_dir, "species_counts.tsv")
    with open(species_counts_path, "w") as fh:
        fh.write("species_id\tn_groups\n")
        for species_id, count in sorted(species_counts.items(),
                                        key=lambda kv: (kv[1], kv[0])):
            fh.write("%s\t%d\n" % (species_id, count))
    logger.info(
        "Group filter species counts written: %d species across surviving groups. Path: %s",
        len(species_counts), species_counts_path,
    )


def strip_species_from_alignments(
    alignment_dir: str,
    stripped_dir: str,
    exclude: set[str],
    show_progress: bool = False,
):
    """Write copies of alignments with sequences from specified species removed.

    Used in two-list mode to remove an injected anchor species from list 1
    alignments before identity estimation, so that group2mean scores reflect
    only genuine list 1 species.

    Args:
        alignment_dir: Directory of alignment files.
        stripped_dir:  Directory where stripped alignments are written.
        exclude:       Set of species IDs to remove.
        show_progress: Show a tqdm progress bar.
    """
    os.makedirs(stripped_dir, exist_ok=True)
    files = os.listdir(alignment_dir)
    for file in tqdm(files, desc="  stripping", unit="group",
                     leave=False, disable=not show_progress):
        in_path  = os.path.join(alignment_dir, file)
        out_path = os.path.join(stripped_dir, file)
        records  = [
            r for r in SeqIO.parse(in_path, "fasta")
            if r.id.split("|")[0] not in exclude
        ]
        if not records:
            logger.warning(
                "All sequences removed from %s after stripping — skipping.", file
            )
            continue
        with open(out_path, "w") as fh:
            for r in records:
                fh.write(">%s\n%s\n" % (r.id, r.seq))
