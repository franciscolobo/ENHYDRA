from __future__ import annotations

import os
import shutil
import random
import logging
import statistics
from collections import Counter
import numpy as np
from Bio import SeqIO
from tqdm import tqdm

logger = logging.getLogger(__name__)

PARALOG_MODES = ("all", "remove", "longest")

# ---------------------------------------------------------------------------
# Pipeline-stage group ID suffixes
# ---------------------------------------------------------------------------
#
# filter_length() appends this suffix when writing a group's surviving
# sequences to length_filter_dir (see its own outfile_f_path in
# filter_length() below). Every later pipeline stage — filter_groups(),
# the aligner, the divergent-sequence filter, make_tables() — reads that
# filename as-is and never strips the suffix back off, so it becomes a
# permanent part of the group_id used for cross-stage file correlation:
# it appears in group2anchor.tsv, group2mean.tsv, alignment filenames, and
# therefore in every alignment page's own file-lookup key. That internal
# ID must not be changed — it is what ties a group's records together
# across every output file — but it is not the group's real, user-facing
# name from the original input directory, and should not be shown
# verbatim to someone reading a report. See display_group_id() below.
_LENGTH_FILTER_SUFFIX = "_lengthfilter"


def display_group_id(internal_group_id: str) -> str:
    """Strip ENHYDRA's internal pipeline-stage suffix from a group ID.

    Every stage after length filtering identifies a group using the exact
    filename filter_length() wrote to length_filter_dir, which carries a
    '_lengthfilter' suffix appended to the group's original name (see the
    module-level note above _LENGTH_FILTER_SUFFIX). That internal ID is
    required for correlating a group across pipeline stages and output
    files — callers must keep using it for lookups, dict keys, and file
    paths — but it is not the group's real name, and reports should
    display this stripped form instead wherever a group_id is shown as
    plain text to a person (e.g. an alignment page's title, or a report
    tab listing groups), as opposed to used as a lookup key.

    This is intentionally a single, non-iterative suffix strip:
    '_lengthfilter' is the only suffix ever appended anywhere in the
    pipeline (no later stage strips or re-suffixes a group_id), so there
    is nothing to loop over.

    Args:
        internal_group_id: The group_id as used internally for file
                           correlation (e.g. as read from group2anchor.tsv,
                           an alignment filename, or any drop_reasons.tsv).

    Returns:
        The group ID with a trailing '_lengthfilter' suffix removed, if
        present; otherwise the input unchanged (e.g. a group skipped
        before length filtering ever wrote a file, or any ID that never
        carried the suffix to begin with).
    """
    return internal_group_id.removesuffix(_LENGTH_FILTER_SUFFIX)


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
    outfile_f_path = os.path.join(length_filter_dir, file + _LENGTH_FILTER_SUFFIX)

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
        base_name = filename[:-len(suffix)]
        path      = os.path.join(length_stats_dir, filename)

        with open(path) as fh:
            lines = [l.rstrip("\n") for l in fh]

        if lines and lines[0] == "##GroupSkipped":
            # Skipped groups (empty_file / single_sequence) never produce a
            # corresponding file in length_filter_dir — filter_length()
            # returns before writing one — so there is no downstream
            # '_lengthfilter'-suffixed group_id to match here, and nothing
            # later in the pipeline will ever look this ID up under a
            # different name (a skipped group has no alignment, or any
            # later stage, to cross-reference against). The bare
            # input-stage name is recorded as-is.
            reason = ""
            detail = ""
            for line in lines[1:]:
                if line.startswith("reason\t"):
                    reason = line.split("\t", 1)[1]
                elif line.startswith("n_sequences\t"):
                    detail = "n_sequences=%s" % line.split("\t", 1)[1]
            skipped_groups.append((base_name, reason, detail))
            continue

        # Otherwise: a processed group. filter_length() always writes this
        # group's survivors to length_filter_dir under
        # '<base_name>_lengthfilter' (see filter_length()'s own
        # outfile_f_path) — the exact filename every later stage
        # (filter_groups(), the divergent-sequence filter, make_tables(),
        # and therefore alignment page IDs) derives its own group_id from.
        # Recording the bare base_name here instead (as a previous version
        # of this function did) silently broke that cross-stage
        # correlation: any caller keying off this function's group_id to
        # look something up in a later stage (e.g. msa_viewer.py's
        # per-group length-filter removal notice on alignment pages) would
        # never find a match, since every later stage's own group_id
        # carries the '_lengthfilter' suffix and this one didn't.
        group_name = base_name + _LENGTH_FILTER_SUFFIX

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
    - species_counts.tsv: one row per species, with three columns:
      n_groups_before (number of groups in length_filter_dir — this step's
      input — that contain that species, counted regardless of whether the
      group ultimately passes this step), n_groups_after (number of
      surviving/written groups that species appears in), and pct_survived
      (n_groups_after / n_groups_before, as a percentage). Rows are sorted
      ascending by n_groups_after so an underrepresented genome sorts to the
      top. Useful for spotting a genome that is systematically thin (e.g.
      due to consistently failing the anchor or length checks upstream) or
      one that is well represented going in but disproportionately dropped
      at this specific step (e.g. because it drives paralog removal).

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
    species_counts_before: Counter = Counter()
    species_counts_after:  Counter = Counter()

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

        # Counted before any of this step's pass/fail decisions, so
        # n_groups_before reflects every species present in this step's
        # input regardless of whether the group ends up surviving.
        species_counts_before.update(uniq_ids)

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
        species_counts_after.update(uniq_ids)

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
        fh.write("species_id\tn_groups_before\tn_groups_after\tpct_survived\n")
        for species_id, n_before in sorted(
            species_counts_before.items(),
            key=lambda kv: (species_counts_after.get(kv[0], 0), kv[0]),
        ):
            n_after = species_counts_after.get(species_id, 0)
            pct     = round(100.0 * n_after / n_before, 1) if n_before else 0.0
            fh.write("%s\t%d\t%d\t%.1f\n" % (species_id, n_before, n_after, pct))
    logger.info(
        "Group filter species counts written: %d species across this step's input. Path: %s",
        len(species_counts_before), species_counts_path,
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


# ---------------------------------------------------------------------------
# Divergent-sequence filter
# ---------------------------------------------------------------------------
#
# A single highly-divergent sequence within an otherwise well-conserved
# group can drag down (or otherwise distort) the group's mean pairwise
# identity — the exact metric ENHYDRA ranks orthogroups by — even though
# the divergence may reflect a misannotation, a paralog mistakenly grouped
# with true orthologs, or a genuinely fast-evolving lineage that shouldn't
# dominate the whole group's score. This filter identifies such sequences
# using trimAl's own per-sequence "identity to most similar sequence"
# values (from `trimal -sident`), rather than a full pairwise average,
# since a sequence's relationship to its single closest match is a more
# direct signal of "does this sequence belong here at all" than its
# average identity against every other member (which is already diluted
# by the rest of the group).
#
# Sequences whose identity-to-closest-match falls more than
# `sd_multiplier` standard deviations below the group's own mean are
# flagged. Removing a flagged sequence invalidates the existing alignment
# for the survivors — the divergent sequence likely forced spurious gap
# placement elsewhere in the alignment — so surviving sequences are
# de-gapped and handed back to the caller for realignment from scratch,
# rather than simply deleting a row from the existing alignment.

_SIDENT_MOST_SIMILAR_HEADER = (
    "## Identity for most similar pair-wise sequences matrix"
)


def parse_sident_most_similar(sident_path: str) -> dict[str, tuple[float, str]]:
    """Parse the "most similar pairwise sequences" section of trimAl -sident output.

    trimAl's `-sident` output includes several sections; this parses only the
    one giving each sequence's identity to its single closest match, e.g.:

        ## Identity for most similar pair-wise sequences matrix
        seqA    0.3716    seqH
        seqB    0.9983    seqC
        ...

    Each data line is whitespace-delimited: sequence_id, identity (float),
    closest_match_id. Parsing stops at the first blank line after the
    section header (or at EOF), so later sections in the same file are not
    accidentally included.

    Args:
        sident_path: Path to a trimAl -sident output file (as written by
                    alignment.run_trimal()).

    Returns:
        Dict mapping sequence_id to (identity_to_closest_match, closest_id).
        Empty dict if the section header is not found (e.g. trimAl produced
        no output for a degenerate alignment, or the file is some other
        format entirely) — callers should treat this the same as "nothing
        to filter" rather than raise, since a missing/malformed section is
        recoverable by simply skipping the filter for that group.

    Raises:
        FileNotFoundError: If sident_path does not exist.
    """
    with open(sident_path) as fh:
        lines = fh.readlines()

    result: dict[str, tuple[float, str]] = {}
    in_section = False
    for raw_line in lines:
        line = raw_line.rstrip("\n")
        if line.startswith(_SIDENT_MOST_SIMILAR_HEADER):
            in_section = True
            continue
        if not in_section:
            continue
        if not line.strip():
            break   # end of section
        fields = line.split()
        if len(fields) < 3:
            continue
        seq_id, identity_str, closest_id = fields[0], fields[1], fields[2]
        try:
            identity = float(identity_str)
        except ValueError:
            continue
        result[seq_id] = (identity, closest_id)

    return result


def detect_divergent_sequences(
    identities: dict[str, tuple[float, str]],
    sd_multiplier: float = 2.0,
) -> tuple[set[str], float, float]:
    """Flag sequences whose identity-to-closest-match is an outlier on the low end.

    A sequence is flagged if its identity value falls more than
    sd_multiplier standard deviations below the mean of all values in
    `identities`. Only the low tail is considered — an unusually *high*
    identity-to-closest-match is not a sign a sequence doesn't belong.

    Args:
        identities:    {seq_id: (identity_to_closest, closest_seq_id)}, as
                      returned by parse_sident_most_similar(). Must contain
                      at least 2 entries (mean/SD are undefined otherwise —
                      callers are expected to have already gated on group
                      size before calling this).
        sd_multiplier: Number of standard deviations below the mean beyond
                      which a sequence is flagged (default: 2.0).

    Returns:
        Tuple of (flagged_seq_ids, mean, stdev).

    Raises:
        statistics.StatisticsError: If identities has fewer than 2 entries.
    """
    values = [identity for identity, _closest in identities.values()]
    mean   = statistics.mean(values)
    stdev  = statistics.stdev(values)
    lower_bound = mean - sd_multiplier * stdev

    flagged = {
        seq_id for seq_id, (identity, _closest) in identities.items()
        if identity < lower_bound
    }
    return flagged, mean, stdev


def filter_divergent_sequences(
    alignment_dir: str,
    sident_dir: str,
    filtered_dir: str,
    realign_input_dir: str,
    stats_dir: str,
    min_species: int,
    min_sequences: int = 2,
    sd_multiplier: float = 2.0,
    show_progress: bool = False,
) -> set[str]:
    """Detect and remove highly divergent sequences from each group's alignment.

    For every alignment in alignment_dir with a matching trimAl -sident
    output in sident_dir (see alignment.run_trimal()):

    - Groups with fewer sequences than min_species are left untouched and
      copied as-is into filtered_dir. With too few sequences, a single
      outlier can drag the mean/SD enough to escape detection entirely
      (the same lesson already documented on test_filtering.py's
      test_outlier_removed for the length filter), so detection is not
      attempted at all below this size rather than attempted unreliably.
    - Groups with a -sident file that has fewer than 2 parseable entries
      (e.g. trimAl produced no usable output) are likewise left untouched,
      since mean/SD cannot be computed.
    - Groups with no sequence flagged by detect_divergent_sequences() are
      left untouched and copied as-is into filtered_dir.
    - Groups with one or more flagged sequences have those sequences
      removed. The alignment itself is now stale for the survivors — a
      divergent sequence typically forces spurious gaps elsewhere in the
      alignment — so the survivors' original (de-gapped) sequences are
      written to realign_input_dir instead of being written to
      filtered_dir directly. Filenames in realign_input_dir use the bare
      group_id (no extension), matching group_filter_dir's convention, so
      that running these through alignment.run_aligner() with
      alignment_dir=filtered_dir produces the correctly-named
      '<group_id>.aln' output directly in filtered_dir with no separate
      renaming step required. Callers are expected to make that
      run_aligner() call themselves after this function returns, only for
      groups that actually needed it (i.e. only if realign_input_dir ends
      up non-empty).
    - Groups that fall below min_species or min_sequences after removing
      their flagged sequence(s) are dropped entirely: nothing is written
      to either filtered_dir or realign_input_dir for that group.

    Every filtering decision is logged to stats_dir/drop_reasons.tsv
    (columns: group_id, sequence_id, identity_to_closest, pct_diff_from_avg,
    reason), opened in write mode so re-runs produce clean output rather
    than appending duplicates. 'reason' is one of:
      - 'removed_divergent_sequence': this sequence was pruned; the group
        survives (to realignment). identity_to_closest and
        pct_diff_from_avg (identity_to_closest / group_mean, matching the
        same ratio convention used in filter_length's own drop_reasons)
        are populated.
      - 'below_min_species_after_divergence_filter': the entire group was
        dropped because removing its flagged sequence(s) left it under
        min_species or min_sequences. sequence_id holds a comma-separated
        list of every sequence that had been flagged (the group, not any
        single sequence, is what was ultimately dropped), and
        identity_to_closest/pct_diff_from_avg are left blank.

    Args:
        alignment_dir:      Directory of alignments to check. Expected to
                            be the point in the pipeline after anchor
                            stripping (if applicable) but before column
                            trimming, so that trimming subsequently
                            operates on the corrected alignment.
        sident_dir:         Directory of trimAl -sident output for the same
                            alignments, e.g. as produced by a prior call to
                            alignment.run_trimal(alignment_dir=alignment_dir,
                            ident_dir=sident_dir, ...). Filenames are
                            expected as '<alignment_filename>.ident'.
        filtered_dir:       Directory where the final, corrected alignments
                            are written: untouched groups go here directly;
                            pruned groups land here once the caller has run
                            them back through the aligner.
        realign_input_dir:  Directory where surviving (de-gapped) sequences
                            of pruned groups are written, awaiting
                            realignment by the caller.
        stats_dir:          Directory where drop_reasons.tsv is written.
        min_species:        Minimum number of distinct species required
                            after pruning, and also the minimum sequence
                            count required before attempting divergence
                            detection at all (see above).
        min_sequences:      Minimum number of sequences required after
                            pruning (default: 2).
        sd_multiplier:      Number of SDs below the mean identity-to-
                            closest-match beyond which a sequence is
                            flagged (default: 2.0).
        show_progress:      Show a tqdm progress bar.

    Returns:
        Set of group_ids dropped entirely by this step.
    """
    os.makedirs(filtered_dir, exist_ok=True)
    os.makedirs(realign_input_dir, exist_ok=True)
    os.makedirs(stats_dir, exist_ok=True)

    drop_reasons:   list[tuple[str, str, str, str, str]] = []
    dropped_groups: set[str] = set()

    files = [f for f in os.listdir(alignment_dir)
             if os.path.isfile(os.path.join(alignment_dir, f))]

    for file in tqdm(files, desc="  groups", unit="group",
                     leave=False, disable=not show_progress):
        group_id    = file.split(".")[0]
        aln_path    = os.path.join(alignment_dir, file)
        sident_path = os.path.join(sident_dir, file + ".ident")

        n_records = sum(1 for _ in SeqIO.parse(aln_path, "fasta"))

        if n_records < min_species:
            shutil.copyfile(aln_path, os.path.join(filtered_dir, file))
            continue

        if not os.path.isfile(sident_path):
            logger.warning(
                "No -sident output found for '%s' (expected: %s) — "
                "skipping divergence filter for this group.",
                file, sident_path,
            )
            shutil.copyfile(aln_path, os.path.join(filtered_dir, file))
            continue

        identities = parse_sident_most_similar(sident_path)
        if len(identities) < 2:
            shutil.copyfile(aln_path, os.path.join(filtered_dir, file))
            continue

        flagged, mean, _stdev = detect_divergent_sequences(
            identities, sd_multiplier=sd_multiplier,
        )

        if not flagged:
            shutil.copyfile(aln_path, os.path.join(filtered_dir, file))
            continue

        records   = list(SeqIO.parse(aln_path, "fasta"))
        survivors = [r for r in records if r.id not in flagged]
        species_ids = {r.id.split("|")[0] for r in survivors}

        # Checked before logging per-sequence removal reasons: if the group
        # ends up dropped entirely, only the single group-level drop row is
        # written (below) — logging both would double-count the same
        # sequences under two different reasons for a group that never
        # actually reaches realignment.
        if len(survivors) < min_sequences or len(species_ids) < min_species:
            logger.warning(
                "Group %s: removing %d divergent sequence(s) leaves only "
                "%d sequence(s) / %d species — below minimum. Group dropped.",
                group_id, len(flagged), len(survivors), len(species_ids),
            )
            drop_reasons.append((
                group_id, ",".join(sorted(flagged)), "", "",
                "below_min_species_after_divergence_filter",
            ))
            dropped_groups.add(group_id)
            continue

        for seq_id in sorted(flagged):
            identity_val, _closest = identities[seq_id]
            pct_diff = (identity_val / mean) if mean else 0.0
            drop_reasons.append((
                group_id, seq_id, "%.4f" % identity_val,
                "%.4f" % pct_diff, "removed_divergent_sequence",
            ))

        realign_out_path = os.path.join(realign_input_dir, group_id)
        with open(realign_out_path, "w") as fh:
            for r in survivors:
                seq = str(r.seq).replace("-", "").replace(".", "")
                fh.write(">%s\n%s\n" % (r.id, seq))

    drop_reasons_path = os.path.join(stats_dir, "drop_reasons.tsv")
    with open(drop_reasons_path, "w") as fh:
        fh.write(
            "group_id\tsequence_id\tidentity_to_closest\tpct_diff_from_avg\treason\n"
        )
        for row in sorted(drop_reasons, key=lambda r: (r[0], r[1])):
            fh.write("%s\t%s\t%s\t%s\t%s\n" % row)
    logger.info(
        "Divergence filter: %d sequence(s) removed, %d group(s) dropped entirely. Path: %s",
        sum(1 for r in drop_reasons if r[4] == "removed_divergent_sequence"),
        len(dropped_groups), drop_reasons_path,
    )

    return dropped_groups
