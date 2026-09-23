import os
import logging
from Bio import SeqIO
from tqdm import tqdm

logger = logging.getLogger(__name__)

_TABLE_FILES = ("group2mean.tsv", "anchor2mean.tsv", "group2anchor.tsv")


def tables_complete(tables_dir: str) -> bool:
    """Return True only if all three output tables exist and are non-empty."""
    return all(
        os.path.isfile(os.path.join(tables_dir, f)) and
        os.path.getsize(os.path.join(tables_dir, f)) > 0
        for f in _TABLE_FILES
    )


def make_tables(
    alignment_dir: str,
    ident_dir: str,
    tables_dir: str,
    anchor: str,
    require_anchor: bool = True,
    show_progress: bool = False,
):
    """Generate ranked output tables from alignments and identity reports.

    Iterates identity files as the primary source and derives the corresponding
    alignment filename for each by stripping the trailing '.ident' suffix.
    A warning is emitted for any mismatch in either direction:

    - Identity file present but alignment file absent → group skipped entirely.
    - Alignment file present but identity file absent → warning only (trimAl
      may have failed silently on that group).

    This single-loop design prevents the silent row omissions that occurred
    with the previous two-loop approach, where a group could appear in
    group2mean without a matching entry in anchor2mean/group2anchor if the
    two directory listings were not perfectly aligned.

    In addition to the three ranked tables, this also writes
    'tables_dir/drop_reasons.tsv', recording every group excluded (fully or
    partially) at this step, with columns group_id, reason, detail. Nothing
    downstream reads tables_dir wholesale (gsea.py, differential.py, etc. all
    reference the three named table files directly by name), so this file can
    sit alongside them safely — unlike group_filter_dir, there's no risk of
    it being picked up by a later step that lists the directory contents.
    'reason' is one of:
      - 'missing_identity'    : alignment exists but no matching .ident file.
      - 'missing_alignment'   : .ident file exists but no matching alignment.
      - 'no_average_identity' : .ident file has no parseable AverageIdentity
                                line (trimAl may have failed silently).
      - 'no_anchor_sequence'  : group has valid alignment + identity and
                                contributes to group2mean.tsv, but has no
                                anchor-species sequence, so it is absent
                                from anchor2mean.tsv / group2anchor.tsv.
                                This is the reason group2mean.tsv and
                                anchor2mean.tsv row counts can diverge.
                                Only recorded here when require_anchor=True
                                (see below) — when the anchor is not expected
                                to be present, this is not an anomaly and
                                would otherwise flood the file with noise.
    This file is opened in write mode, so re-running make_tables() (e.g.
    under --resume) produces clean output rather than appending duplicates.

    Args:
        alignment_dir:  Directory of alignment files (.aln).
        ident_dir:      Directory of trimAl identity report files (.aln.ident).
        tables_dir:     Directory where output tables are written.
        anchor:         Anchor species ID used to map group → gene ID.
        require_anchor: Whether the anchor species is expected to be present
                        in this list's alignments. In two-list differential
                        mode, list2 is not expected to contain the anchor
                        (it belongs to list1) — a missing anchor there is
                        expected behaviour, not a drop, since the group
                        still contributes to group2mean.tsv either way and
                        list2's anchor2mean.tsv/group2anchor.tsv are never
                        read downstream. When False, missing anchors are
                        logged at debug level only and are not recorded in
                        drop_reasons.tsv. Defaults to True (matches previous
                        behaviour for single-list mode and list1).
        show_progress:  False is verbose, True shows progress bar only.
    """
    os.makedirs(tables_dir, exist_ok=True)

    ident_files = os.listdir(ident_dir)
    aln_files   = set(os.listdir(alignment_dir))

    drop_reasons: list[tuple[str, str, str]] = []

    # Pre-check: warn about alignment files that have no matching identity file.
    # These indicate groups where trimAl may have failed silently.
    ident_stems = {f[:-len(".ident")] for f in ident_files if f.endswith(".ident")}
    for aln_file in sorted(aln_files):
        if aln_file not in ident_stems:
            logger.warning(
                "No identity file found for alignment '%s' — "
                "trimAl may have failed silently on this group.", aln_file,
            )
            drop_reasons.append((
                aln_file.split(".")[0], "missing_identity",
                "alignment=%s" % aln_file,
            ))

    with open(os.path.join(tables_dir, "group2mean.tsv"),   "w") as group2mean, \
         open(os.path.join(tables_dir, "anchor2mean.tsv"),  "w") as anchor2mean, \
         open(os.path.join(tables_dir, "group2anchor.tsv"), "w") as group2anchor:

        for ident_file in tqdm(ident_files, desc="  tables",
                               unit="group", leave=False, disable=not show_progress):
            group_name = ident_file.split(".")[0]
            ident_path = os.path.join(ident_dir, ident_file)

            # Derive the alignment filename by stripping the trailing .ident suffix.
            aln_file = (ident_file[:-len(".ident")]
                        if ident_file.endswith(".ident")
                        else ident_file)
            aln_path = os.path.join(alignment_dir, aln_file)

            if not os.path.isfile(aln_path):
                logger.warning(
                    "Alignment file not found for group '%s' (expected: %s) — "
                    "skipping group.", group_name, aln_path,
                )
                drop_reasons.append((
                    group_name, "missing_alignment",
                    "identity=%s" % ident_file,
                ))
                continue

            # Parse identity score from trimAl -sident output.
            mean_percent = None
            with open(ident_path) as fh:
                for line in fh:
                    line = line.rstrip()
                    if line.startswith("## AverageIdentity"):
                        mean_percent = line.split()[-1]
                        break

            if mean_percent is None:
                logger.warning(
                    "No AverageIdentity line found in '%s' — skipping group.",
                    ident_path,
                )
                drop_reasons.append((
                    group_name, "no_average_identity", "",
                ))
                continue

            group2mean.write("%s\t%s\n" % (group_name, mean_percent))

            # Map anchor gene ID(s) from the alignment file.
            # In 'all' paralog mode there may be more than one anchor sequence;
            # all are written to preserve the original behaviour.
            anchor_found = False
            for seq_record in SeqIO.parse(aln_path, "fasta"):
                fields  = seq_record.id.split("|")
                species = fields[0]
                gene_id = fields[1]
                if species == anchor:
                    anchor2mean.write("%s\t%s\n" % (gene_id, mean_percent))
                    group2anchor.write("%s\t%s\n" % (group_name, gene_id))
                    anchor_found = True

            if not anchor_found:
                if require_anchor:
                    logger.warning(
                        "No anchor sequence found in alignment for group '%s' — "
                        "group contributes to group2mean but not anchor2mean or "
                        "group2anchor.", group_name,
                    )
                    drop_reasons.append((
                        group_name, "no_anchor_sequence", "anchor=%s" % anchor,
                    ))
                else:
                    # Expected: this list is not required to contain the
                    # anchor species (e.g. list2 in two-list differential
                    # mode). The group still contributes to group2mean.tsv
                    # above; only anchor2mean/group2anchor are unaffected,
                    # and neither of those is read from this list downstream.
                    logger.debug(
                        "No anchor sequence found in alignment for group '%s' "
                        "(expected — anchor not required for this list).",
                        group_name,
                    )

    drop_reasons_path = os.path.join(tables_dir, "drop_reasons.tsv")
    with open(drop_reasons_path, "w") as fh:
        fh.write("group_id\treason\tdetail\n")
        for group_id, reason, detail in sorted(drop_reasons):
            fh.write("%s\t%s\t%s\n" % (group_id, reason, detail))
    logger.info(
        "Tables step drop reasons written: %d entr%s. Path: %s",
        len(drop_reasons), "y" if len(drop_reasons) == 1 else "ies",
        drop_reasons_path,
    )
