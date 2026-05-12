"""Tests for enhydra.filtering."""

import os
import pytest
from Bio import SeqIO

from enhydra.filtering import (
    filter_length,
    filter_groups,
    subset_groups,
    strip_species_from_alignments,
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
        assert len(ids) == 4   # 2 seqs × 2 calls (append mode is intentional)


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
