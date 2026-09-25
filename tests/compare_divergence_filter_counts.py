"""
Compare per-group sequence counts before and after the divergent-sequence
filter, to confirm the filter is actually removing sequences (and to see
which groups were affected).

Compares:
    <outdir>/alignment/                       (before)
    <outdir>/alignment_divergence_filtered/    (after)

In two-list mode, point --before/--after at 'outdir/list1/alignment' and
'outdir/list1/alignment_divergence_filtered' (and separately for list2).

Every group filename in 'after' is expected to also exist in 'before'
(the filter never renames a group, only removes sequences and/or drops
the group entirely). Groups present in 'before' but absent from 'after'
were dropped entirely by the filter (see divergence_filter_stats/
drop_reasons.tsv for the reason) and are reported separately.

Usage:
    python tests/compare_divergence_filter_counts.py \\
        --before  outdir/alignment \\
        --after   outdir/alignment_divergence_filtered
"""

from __future__ import annotations

import os
import sys
import argparse
from Bio import SeqIO


def _count_seqs(path: str) -> int:
    return sum(1 for _ in SeqIO.parse(path, "fasta"))


def compare_counts(before_dir: str, after_dir: str) -> None:
    before_files = {f for f in os.listdir(before_dir)
                    if os.path.isfile(os.path.join(before_dir, f))}
    after_files  = {f for f in os.listdir(after_dir)
                    if os.path.isfile(os.path.join(after_dir, f))}

    changed  = []   # (group_id, n_before, n_after)
    unchanged_count = 0
    dropped  = sorted(before_files - after_files)

    for filename in sorted(before_files & after_files):
        group_id = filename.split(".")[0]
        n_before = _count_seqs(os.path.join(before_dir, filename))
        n_after  = _count_seqs(os.path.join(after_dir, filename))
        if n_after < n_before:
            changed.append((group_id, n_before, n_after))
        elif n_after > n_before:
            # Should not happen — flagged loudly since it would indicate
            # a bug (the filter should only ever remove sequences).
            print("WARNING: %s has MORE sequences after filtering "
                  "(%d -> %d) — this should never happen."
                  % (group_id, n_before, n_after))
        else:
            unchanged_count += 1

    print("=" * 60)
    print("Divergent-sequence filter: before/after comparison")
    print("=" * 60)
    print("Groups compared:      %d" % len(before_files & after_files))
    print("Groups unchanged:     %d" % unchanged_count)
    print("Groups with removals: %d" % len(changed))
    print("Groups dropped entirely: %d" % len(dropped))
    print()

    if changed:
        print("--- Groups with sequences removed ---")
        print("%-30s %10s %10s %10s" % ("group_id", "n_before", "n_after", "n_removed"))
        for group_id, n_before, n_after in changed:
            print("%-30s %10d %10d %10d"
                  % (group_id, n_before, n_after, n_before - n_after))
        print()

    if dropped:
        print("--- Groups dropped entirely (see divergence_filter_stats/"
              "drop_reasons.tsv for reason) ---")
        for filename in dropped:
            print("  %s" % filename)
        print()

    if not changed and not dropped:
        print("No groups were affected by the divergence filter. If you "
              "expected removals, check that --divergence-filter-sd was "
              "actually set for this run (see enhydra.log for "
              "'Step 3c: Filtering highly divergent sequences').")


def main():
    parser = argparse.ArgumentParser(
        description="Compare per-group sequence counts before/after the "
                    "divergent-sequence filter."
    )
    parser.add_argument("--before", required=True,
                        help="Path to the pre-filter alignment directory "
                             "(outdir/alignment or outdir/listN/alignment).")
    parser.add_argument("--after", required=True,
                        help="Path to the post-filter alignment directory "
                             "(outdir/alignment_divergence_filtered or "
                             "outdir/listN/alignment_divergence_filtered).")
    args = parser.parse_args()

    if not os.path.isdir(args.before):
        sys.exit("--before directory not found: %s" % args.before)
    if not os.path.isdir(args.after):
        sys.exit(
            "--after directory not found: %s\n"
            "This usually means divergence_filter_sd was not set for this "
            "run, so the step never ran." % args.after
        )

    compare_counts(args.before, args.after)


if __name__ == "__main__":
    main()
