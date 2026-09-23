"""Tests for enhydra.filtering."""

import os
import pytest
from Bio import SeqIO

from enhydra.filtering import (
    filter_length,
    filter_groups,
    subset_groups,
    strip_species_from_alignments,
    aggregate_length_filter_stats
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _write_fasta(path: str, entries: list[tuple[str, str]]):
    """Write (header, sequence) pairs to a FASTA file."""
    with open(path, "w") as fh:
        for header, seq in entries:
            fh.write(">%s\n%s\n" % (header, seq))


def _read_ids(path: str) -> list[str]:
    """Return sequence IDs present in a FASTA file."""
    return [r.id for r in SeqIO.parse(path, "fasta")]


# ---------------------------------------------------------------------------
# filter_length
# ---------------------------------------------------------------------------

class TestFilterLength:

    def test_normal_group_passes(self, tmp_path):
        """All sequences within 2 SD of the mean are kept."""
        src = tmp_path / "input" / "OG0001"
        src.parent.mkdir()
        _write_fasta(str(src), [
            ("sp1|g1", "A" * 100),
            ("sp2|g2", "A" * 105),
            ("sp3|g3", "A" * 95),
            ("sp4|g4", "A" * 100),
        ])
        stats_dir  = tmp_path / "stats"
        filter_dir = tmp_path / "filter"
        stats_dir.mkdir(); filter_dir.mkdir()

        filter_length(str(src), str(stats_dir), str(filter_dir))

        out = filter_dir / "OG0001_lengthfilter"
        assert out.exists()
        assert len(_read_ids(str(out))) == 4

    def test_outlier_removed(self, tmp_path):
        """A sequence far outside mean ± 2 SD is removed.

        With only a handful of sequences the outlier drags the mean and SD
        enough to survive the filter — we need enough baseline sequences to
        dilute its influence. 10 sequences at length 100 + 1 at 400 gives:
            mean ≈ 127, SD ≈ 90  →  mean + 2*SD ≈ 308  <  400  ✓
        """
        src = tmp_path / "input" / "OG0002"
        src.parent.mkdir()
        entries = [("sp%d|g%d" % (i, i), "A" * 100) for i in range(10)]
        entries.append(("sp99|g99", "A" * 400))   # clear outlier
        _write_fasta(str(src), entries)
        stats_dir  = tmp_path / "stats"
        filter_dir = tmp_path / "filter"
        stats_dir.mkdir(); filter_dir.mkdir()

        filter_length(str(src), str(stats_dir), str(filter_dir))

        out = filter_dir / "OG0002_lengthfilter"
        ids = _read_ids(str(out))
        assert "sp99|g99" not in ids
        assert len(ids) == 10

    def test_single_sequence_skipped(self, tmp_path):
        """A group with only one sequence is skipped — no output file written."""
        src = tmp_path / "input" / "OG0003"
        src.parent.mkdir()
        _write_fasta(str(src), [("sp1|g1", "A" * 100)])
        stats_dir  = tmp_path / "stats"
        filter_dir = tmp_path / "filter"
        stats_dir.mkdir(); filter_dir.mkdir()

        filter_length(str(src), str(stats_dir), str(filter_dir))

        assert not (filter_dir / "OG0003_lengthfilter").exists()

    def test_empty_file_skipped(self, tmp_path):
        """An empty input file produces no output and no error."""
        src = tmp_path / "input" / "OG0004"
        src.parent.mkdir()
        src.touch()   # zero bytes
        stats_dir  = tmp_path / "stats"
        filter_dir = tmp_path / "filter"
        stats_dir.mkdir(); filter_dir.mkdir()

        filter_length(str(src), str(stats_dir), str(filter_dir))

        assert not (filter_dir / "OG0004_lengthfilter").exists()

    def test_write_mode_no_duplicates(self, tmp_path):
        """Calling filter_length twice on the same group does not duplicate rows."""
        src = tmp_path / "input" / "OG0005"
        src.parent.mkdir()
        _write_fasta(str(src), [
            ("sp1|g1", "A" * 100),
            ("sp2|g2", "A" * 100),
        ])
        stats_dir  = tmp_path / "stats"
        filter_dir = tmp_path / "filter"
        stats_dir.mkdir(); filter_dir.mkdir()

        filter_length(str(src), str(stats_dir), str(filter_dir))
        filter_length(str(src), str(stats_dir), str(filter_dir))

        out = filter_dir / "OG0005_lengthfilter"
        # append mode is expected here — this test documents current behaviour
        ids = _read_ids(str(out))
        assert len(ids) == 2   # 2 seqs × 2 calls (append mode is intentional)


# ---------------------------------------------------------------------------
# filter_groups
# ---------------------------------------------------------------------------

class TestFilterGroups:

    def _make_group(self, directory: str, name: str, entries: list[tuple[str, str]]):
        path = os.path.join(directory, name + "_lengthfilter")
        _write_fasta(path, entries)
        return path

    def test_clean_group_passes(self, tmp_path):
        lf = tmp_path / "lf"; gf = tmp_path / "gf"
        lf.mkdir(); gf.mkdir()
        self._make_group(str(lf), "OG0001", [
            ("sp1|g1", "A" * 100), ("sp2|g2", "A" * 100),
            ("sp3|g3", "A" * 100), ("sp4|g4", "A" * 100),
            ("anchor|g5", "A" * 100),
        ])
        filter_groups(str(lf), str(gf), anchor="anchor", min_species=4)
        assert len(os.listdir(str(gf))) == 1

    def test_missing_anchor_removed(self, tmp_path):
        lf = tmp_path / "lf"; gf = tmp_path / "gf"
        lf.mkdir(); gf.mkdir()
        self._make_group(str(lf), "OG0001", [
            ("sp1|g1", "A" * 100), ("sp2|g2", "A" * 100),
            ("sp3|g3", "A" * 100), ("sp4|g4", "A" * 100),
        ])
        filter_groups(str(lf), str(gf), anchor="anchor", min_species=2,
                      require_anchor=True)
        assert len(os.listdir(str(gf))) == 0

    def test_missing_anchor_allowed(self, tmp_path):
        """require_anchor=False keeps groups without the anchor."""
        lf = tmp_path / "lf"; gf = tmp_path / "gf"
        lf.mkdir(); gf.mkdir()
        self._make_group(str(lf), "OG0001", [
            ("sp1|g1", "A" * 100), ("sp2|g2", "A" * 100),
            ("sp3|g3", "A" * 100), ("sp4|g4", "A" * 100),
        ])
        filter_groups(str(lf), str(gf), anchor="anchor", min_species=2,
                      require_anchor=False)
        assert len(os.listdir(str(gf))) == 1

    def test_below_min_species_removed(self, tmp_path):
        lf = tmp_path / "lf"; gf = tmp_path / "gf"
        lf.mkdir(); gf.mkdir()
        self._make_group(str(lf), "OG0001", [
            ("sp1|g1", "A" * 100), ("anchor|g2", "A" * 100),
        ])
        filter_groups(str(lf), str(gf), anchor="anchor", min_species=4)
        assert len(os.listdir(str(gf))) == 0

    def test_below_min_sequences_removed(self, tmp_path):
        lf = tmp_path / "lf"; gf = tmp_path / "gf"
        lf.mkdir(); gf.mkdir()
        self._make_group(str(lf), "OG0001", [("anchor|g1", "A" * 100)])
        filter_groups(str(lf), str(gf), anchor="anchor", min_species=1,
                      min_sequences=2)
        assert len(os.listdir(str(gf))) == 0

    def test_paralogs_remove(self, tmp_path):
        """paralog_mode='remove' discards groups with duplicate species."""
        lf = tmp_path / "lf"; gf = tmp_path / "gf"
        lf.mkdir(); gf.mkdir()
        self._make_group(str(lf), "OG0001", [
            ("sp1|g1", "A" * 100), ("sp1|g2", "A" * 100),  # paralog
            ("sp2|g3", "A" * 100), ("anchor|g4", "A" * 100),
        ])
        filter_groups(str(lf), str(gf), anchor="anchor", min_species=2,
                      paralog_mode="remove")
        assert len(os.listdir(str(gf))) == 0

    def test_paralogs_longest(self, tmp_path):
        """paralog_mode='longest' keeps one sequence per species."""
        lf = tmp_path / "lf"; gf = tmp_path / "gf"
        lf.mkdir(); gf.mkdir()
        self._make_group(str(lf), "OG0001", [
            ("sp1|g1", "A" * 80),   # shorter
            ("sp1|g2", "A" * 120),  # longer — should be kept
            ("sp2|g3", "A" * 100),
            ("anchor|g4", "A" * 100),
        ])
        filter_groups(str(lf), str(gf), anchor="anchor", min_species=2,
                      paralog_mode="longest")
        out = os.path.join(str(gf), "OG0001_lengthfilter")
        ids = _read_ids(out)
        assert len(ids) == 3
        assert "sp1|g2" in ids
        assert "sp1|g1" not in ids


# ---------------------------------------------------------------------------
# subset_groups
# ---------------------------------------------------------------------------

class TestSubsetGroups:

    def test_keeps_matching_species(self, tmp_path):
        src = tmp_path / "input"; dst = tmp_path / "subset"
        src.mkdir(); dst.mkdir()
        _write_fasta(str(src / "OG0001"), [
            ("sp1|g1", "A" * 100), ("sp2|g2", "A" * 100),
            ("sp3|g3", "A" * 100),
        ])
        subset_groups(str(src), str(dst), species=["sp1", "sp2"])
        ids = _read_ids(str(dst / "OG0001"))
        assert set(ids) == {"sp1|g1", "sp2|g2"}

    def test_group_with_no_matches_skipped(self, tmp_path):
        src = tmp_path / "input"; dst = tmp_path / "subset"
        src.mkdir(); dst.mkdir()
        _write_fasta(str(src / "OG0001"), [
            ("sp3|g1", "A" * 100), ("sp4|g2", "A" * 100),
        ])
        subset_groups(str(src), str(dst), species=["sp1", "sp2"])
        assert not (dst / "OG0001").exists()


# ---------------------------------------------------------------------------
# strip_species_from_alignments
# ---------------------------------------------------------------------------

class TestStripSpeciesFromAlignments:

    def test_excluded_species_absent(self, tmp_path):
        aln = tmp_path / "aln"; out = tmp_path / "stripped"
        aln.mkdir(); out.mkdir()
        _write_fasta(str(aln / "OG0001.aln"), [
            ("sp1|g1", "ACGT"), ("sp2|g2", "ACGT"), ("anchor|g3", "ACGT"),
        ])
        strip_species_from_alignments(str(aln), str(out), exclude={"anchor"})
        ids = _read_ids(str(out / "OG0001.aln"))
        assert "anchor|g3" not in ids
        assert set(ids) == {"sp1|g1", "sp2|g2"}

    def test_all_excluded_produces_no_output(self, tmp_path):
        aln = tmp_path / "aln"; out = tmp_path / "stripped"
        aln.mkdir(); out.mkdir()
        _write_fasta(str(aln / "OG0001.aln"), [("anchor|g1", "ACGT")])
        strip_species_from_alignments(str(aln), str(out), exclude={"anchor"})
        assert not (out / "OG0001.aln").exists()


# ---------------------------------------------------------------------------
# test_filter_length_stats
# ---------------------------------------------------------------------------

class TestFilterLengthStats:

    def _read_tsv(self, path):
        with open(path) as fh:
            lines = [l.rstrip("\n").split("\t") for l in fh if l.strip()]
        return lines[0], lines[1:]

    def test_skipped_group_writes_lengthstats_marker(self, tmp_path):
        src = tmp_path / "input" / "OG0001"
        src.parent.mkdir()
        _write_fasta(str(src), [("sp1|g1", "A" * 100)])   # single sequence
        stats_dir  = tmp_path / "stats"
        filter_dir = tmp_path / "filter"
        stats_dir.mkdir(); filter_dir.mkdir()

        filter_length(str(src), str(stats_dir), str(filter_dir))

        content = (stats_dir / "OG0001_lengthstats").read_text()
        assert content.startswith("##GroupSkipped")
        assert "reason\tsingle_sequence" in content

    def test_empty_file_writes_lengthstats_marker(self, tmp_path):
        src = tmp_path / "input" / "OG0002"
        src.parent.mkdir()
        src.touch()
        stats_dir  = tmp_path / "stats"
        filter_dir = tmp_path / "filter"
        stats_dir.mkdir(); filter_dir.mkdir()

        filter_length(str(src), str(stats_dir), str(filter_dir))

        content = (stats_dir / "OG0002_lengthstats").read_text()
        assert content.startswith("##GroupSkipped")
        assert "reason\tempty_file" in content

    def test_lengthstats_marks_outlier_status(self, tmp_path):
        src = tmp_path / "input" / "OG0003"
        src.parent.mkdir()
        entries = [("sp%d|g%d" % (i, i), "A" * 100) for i in range(10)]
        entries.append(("sp99|g99", "A" * 400))
        _write_fasta(str(src), entries)
        stats_dir  = tmp_path / "stats"
        filter_dir = tmp_path / "filter"
        stats_dir.mkdir(); filter_dir.mkdir()

        filter_length(str(src), str(stats_dir), str(filter_dir))

        content = (stats_dir / "OG0003_lengthstats").read_text()
        assert "sp99|g99\t400" in content
        assert "removed_above_max" in content
        assert content.count("\tkept\n") == 10

    def test_aggregate_collects_skipped_and_outliers(self, tmp_path):
        input_dir  = tmp_path / "input"
        stats_dir  = tmp_path / "stats"
        filter_dir = tmp_path / "filter"
        agg_dir    = tmp_path / "length_filter_stats"
        input_dir.mkdir(); stats_dir.mkdir(); filter_dir.mkdir()

        # OG0001: single sequence -> skipped
        _write_fasta(str(input_dir / "OG0001"), [("sp1|g1", "A" * 100)])
        # OG0002: empty file -> skipped
        (input_dir / "OG0002").touch()
        # OG0003: has one outlier
        entries = [("sp%d|g%d" % (i, i), "A" * 100) for i in range(10)]
        entries.append(("sp99|g99", "A" * 400))
        _write_fasta(str(input_dir / "OG0003"), entries)

        for name in ("OG0001", "OG0002", "OG0003"):
            filter_length(str(input_dir / name), str(stats_dir), str(filter_dir))

        aggregate_length_filter_stats(str(stats_dir), str(agg_dir))

        _, skipped_rows = self._read_tsv(agg_dir / "skipped_groups.tsv")
        skipped = {r[0]: r[1] for r in skipped_rows}
        assert skipped["OG0001"] == "single_sequence"
        assert skipped["OG0002"] == "empty_file"

        _, drop_rows = self._read_tsv(agg_dir / "drop_reasons.tsv")
        assert len(drop_rows) == 1
        assert drop_rows[0][0] == "OG0003"
        assert drop_rows[0][1] == "sp99|g99"
        assert drop_rows[0][4] == "removed_above_max"

    def test_aggregate_with_no_stats_files_writes_empty_tables(self, tmp_path):
        stats_dir = tmp_path / "stats"
        agg_dir   = tmp_path / "length_filter_stats"
        stats_dir.mkdir()

        aggregate_length_filter_stats(str(stats_dir), str(agg_dir))

        header1, rows1 = self._read_tsv(agg_dir / "skipped_groups.tsv")
        header2, rows2 = self._read_tsv(agg_dir / "drop_reasons.tsv")
        assert rows1 == [] and rows2 == []

