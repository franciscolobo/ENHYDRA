"""Tests for enhydra.filtering."""

import os
import pytest
from Bio import SeqIO

from enhydra.filtering import (
    filter_length,
    filter_groups,
    subset_groups,
    strip_species_from_alignments,
    aggregate_length_filter_stats,
    parse_sident_most_similar,
    detect_divergent_sequences,
    filter_divergent_sequences,
    display_group_id,
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


def _write_sident(path: str, entries: list[tuple[str, float, str]],
                  extra_trailer: list[str] | None = None):
    """Write a minimal trimAl -sident-style file containing just the
    'most similar pairwise sequences' section that parse_sident_most_similar()
    reads, optionally followed by a blank line and further trailer lines
    (to exercise the 'section ends at blank line' behaviour).

    entries = [(seq_id, identity, closest_seq_id), ...]
    """
    with open(path, "w") as fh:
        fh.write("## Identity for most similar pair-wise sequences matrix\n")
        for seq_id, identity, closest_id in entries:
            fh.write("%s\t%.4f\t%s\n" % (seq_id, identity, closest_id))
        if extra_trailer:
            fh.write("\n")
            for line in extra_trailer:
                fh.write(line + "\n")


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
        # group_id must carry the '_lengthfilter' suffix — this is the
        # exact filename filter_length() wrote to length_filter_dir for
        # this group, which every downstream stage (group_filter,
        # alignment, tables, and therefore alignment page IDs) derives
        # its own group_id from. See
        # test_drop_reasons_group_id_matches_length_filter_dir_filename
        # below for the stronger, filename-derived version of this check.
        assert drop_rows[0][0] == "OG0003_lengthfilter"
        assert drop_rows[0][1] == "sp99|g99"
        assert drop_rows[0][4] == "removed_above_max"

    def test_drop_reasons_group_id_matches_length_filter_dir_filename(self, tmp_path):
        """The group_id recorded in drop_reasons.tsv must exactly match the
        actual filename filter_length() wrote to length_filter_dir for that
        group — the filename every downstream pipeline stage (group_filter,
        alignment, make_tables, and therefore the alignment-page IDs used
        by msa_viewer.build_alignment_pages()) derives its own group_id
        from. A prior version of aggregate_length_filter_stats() stripped
        this suffix back off, silently breaking that correlation (see
        commit notes): a length-filter removal notice on an alignment page
        would never appear, since the group_id used to look it up never
        matched the group_id used to store it.

        This test derives its expectation directly from the real on-disk
        filename rather than hardcoding the '_lengthfilter' suffix, so it
        stays correct even if that naming convention ever changes, as long
        as it changes consistently across the pipeline.
        """
        input_dir  = tmp_path / "input"
        stats_dir  = tmp_path / "stats"
        filter_dir = tmp_path / "filter"
        agg_dir    = tmp_path / "length_filter_stats"
        input_dir.mkdir(); stats_dir.mkdir(); filter_dir.mkdir()

        entries = [("sp%d|g%d" % (i, i), "A" * 100) for i in range(10)]
        entries.append(("sp99|g99", "A" * 400))
        _write_fasta(str(input_dir / "OG0003"), entries)

        filter_length(str(input_dir / "OG0003"), str(stats_dir), str(filter_dir))
        aggregate_length_filter_stats(str(stats_dir), str(agg_dir))

        filter_dir_files = os.listdir(str(filter_dir))
        assert len(filter_dir_files) == 1
        actual_downstream_group_id = filter_dir_files[0]

        _, drop_rows = self._read_tsv(agg_dir / "drop_reasons.tsv")
        assert len(drop_rows) == 1
        assert drop_rows[0][0] == actual_downstream_group_id

    def test_aggregate_with_no_stats_files_writes_empty_tables(self, tmp_path):
        stats_dir = tmp_path / "stats"
        agg_dir   = tmp_path / "length_filter_stats"
        stats_dir.mkdir()

        aggregate_length_filter_stats(str(stats_dir), str(agg_dir))

        header1, rows1 = self._read_tsv(agg_dir / "skipped_groups.tsv")
        header2, rows2 = self._read_tsv(agg_dir / "drop_reasons.tsv")
        assert rows1 == [] and rows2 == []


# ---------------------------------------------------------------------------
# parse_sident_most_similar
# ---------------------------------------------------------------------------

class TestParseSidentMostSimilar:

    def test_parses_real_trimal_sample(self, tmp_path):
        """End-to-end against an actual `trimal -sident` transcript,
        including the sections before and after the one we care about,
        to make sure parsing isn't accidentally order- or content-sensitive."""
        path = tmp_path / "OG0001.aln.ident"
        path.write_text(
            "## MaxIdentity 1.0000\n"
            "#> MaxIdentity Get the maximum identity value for any pair of "
            "sequences in the alignment\n"
            "\n"
            "## AverageIdentity 0.8305\n"
            "#> AverageIdentity Average identity between all sequences\n"
            "\n"
            "## Identity sequences matrix\n"
            "NC_011601|WP_000190964.1    1.0000    0.3667\n"
            "NC_017626|WP_024261658.1    0.3667    1.0000\n"
            "\n"
            "## AverageMostSimilarIdentity 0.9148\n"
            "#> AverageMostSimilarIdentity Average identity between most "
            "similar pair-wise sequences\n"
            "\n"
            "## Identity for most similar pair-wise sequences matrix\n"
            "NC_011601|WP_000190964.1    0.3716 NC_017626|WP_001529578.1\n"
            "NC_017626|WP_024261658.1    0.9983 NC_011748|WP_001390303.1\n"
            "NC_011748|WP_001390303.1    1.0000 NC_018650|WP_001390303.1\n"
            "NC_018650|WP_001390303.1    1.0000 NC_011748|WP_001390303.1\n"
            "NC_018658|WP_001390303.1    1.0000 NC_011748|WP_001390303.1\n"
            "NC_018661|WP_001390303.1    1.0000 NC_011748|WP_001390303.1\n"
            "NZ_HF572917|WP_001390303.1  1.0000 NC_011748|WP_001390303.1\n"
            "NC_017626|WP_001529578.1    0.9484 NC_017626|WP_024261658.1\n"
        )
        result = parse_sident_most_similar(str(path))
        assert len(result) == 8
        identity, closest = result["NC_011601|WP_000190964.1"]
        assert identity == pytest.approx(0.3716)
        assert closest == "NC_017626|WP_001529578.1"
        identity2, closest2 = result["NC_017626|WP_001529578.1"]
        assert identity2 == pytest.approx(0.9484)
        assert closest2 == "NC_017626|WP_024261658.1"

    def test_stops_at_blank_line_after_section(self, tmp_path):
        """Content after the section's terminating blank line must not leak in,
        even if it happens to look like more data rows."""
        path = tmp_path / "OG0002.aln.ident"
        _write_sident(
            str(path),
            entries=[("seqA", 0.95, "seqB"), ("seqB", 0.95, "seqA")],
            extra_trailer=["seqC\t0.1234\tseqD"],
        )
        result = parse_sident_most_similar(str(path))
        assert set(result.keys()) == {"seqA", "seqB"}
        assert "seqC" not in result

    def test_missing_section_header_returns_empty_dict(self, tmp_path):
        path = tmp_path / "OG0003.aln.ident"
        path.write_text("## AverageIdentity 0.5\nsome unrelated content\n")
        assert parse_sident_most_similar(str(path)) == {}

    def test_malformed_line_skipped(self, tmp_path):
        """A data line with too few whitespace-separated fields is skipped
        rather than raising."""
        path = tmp_path / "OG0004.aln.ident"
        path.write_text(
            "## Identity for most similar pair-wise sequences matrix\n"
            "seqA\n"
            "seqB\t0.9\tseqA\n"
        )
        result = parse_sident_most_similar(str(path))
        assert set(result.keys()) == {"seqB"}

    def test_non_numeric_identity_skipped(self, tmp_path):
        path = tmp_path / "OG0005.aln.ident"
        path.write_text(
            "## Identity for most similar pair-wise sequences matrix\n"
            "seqA\tNaN_or_garbage\tseqB\n"
            "seqB\t0.9\tseqA\n"
        )
        result = parse_sident_most_similar(str(path))
        assert set(result.keys()) == {"seqB"}

    def test_missing_file_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            parse_sident_most_similar(str(tmp_path / "nonexistent.ident"))


# ---------------------------------------------------------------------------
# detect_divergent_sequences
# ---------------------------------------------------------------------------

class TestDetectDivergentSequences:

    def test_flags_low_outlier(self):
        """Mirrors the real-world example: one clear low outlier among
        otherwise tightly-clustered values."""
        identities = {
            "NC_011601|WP_000190964.1":   (0.3716, "x"),
            "NC_017626|WP_024261658.1":   (0.9983, "x"),
            "NC_011748|WP_001390303.1":   (1.0000, "x"),
            "NC_018650|WP_001390303.1":   (1.0000, "x"),
            "NC_018658|WP_001390303.1":   (1.0000, "x"),
            "NC_018661|WP_001390303.1":   (1.0000, "x"),
            "NZ_HF572917|WP_001390303.1": (1.0000, "x"),
            "NC_017626|WP_001529578.1":   (0.9484, "x"),
        }
        flagged, mean, stdev = detect_divergent_sequences(identities, sd_multiplier=2.0)
        assert flagged == {"NC_011601|WP_000190964.1"}
        assert mean == pytest.approx(0.9148, abs=1e-4)

    def test_no_outliers_returns_empty_set(self):
        identities = {
            "a": (0.95, "x"), "b": (0.96, "x"),
            "c": (0.94, "x"), "d": (0.97, "x"),
        }
        flagged, _mean, _stdev = detect_divergent_sequences(identities, sd_multiplier=2.0)
        assert flagged == set()

    def test_high_identity_not_flagged(self):
        """Only the low tail is ever flagged — an unusually high
        identity-to-closest-match is not evidence of not belonging."""
        identities = {
            "a": (0.5, "x"), "b": (0.5, "x"),
            "c": (0.5, "x"), "d": (0.999, "x"),
        }
        flagged, _mean, _stdev = detect_divergent_sequences(identities, sd_multiplier=1.0)
        assert "d" not in flagged

    def test_all_identical_values_no_flags(self):
        """stdev == 0: nothing is strictly less than the mean, so nothing
        is flagged regardless of sd_multiplier."""
        identities = {"a": (0.9, "x"), "b": (0.9, "x"), "c": (0.9, "x")}
        flagged, mean, stdev = detect_divergent_sequences(identities, sd_multiplier=2.0)
        assert flagged == set()
        assert stdev == 0.0

    def test_returns_mean_and_stdev(self):
        identities = {"a": (0.8, "x"), "b": (1.0, "x")}
        _flagged, mean, stdev = detect_divergent_sequences(identities)
        assert mean == pytest.approx(0.9)
        assert stdev == pytest.approx(0.14142, abs=1e-4)

    def test_raises_on_single_entry(self):
        import statistics
        identities = {"a": (0.9, "x")}
        with pytest.raises(statistics.StatisticsError):
            detect_divergent_sequences(identities)

    def test_sd_multiplier_affects_sensitivity(self):
        """A looser multiplier should flag fewer (or equal) sequences than
        a stricter one, for the same data."""
        identities = {
            "a": (0.5, "x"), "b": (0.9, "x"),
            "c": (0.95, "x"), "d": (1.0, "x"),
        }
        strict, _, _ = detect_divergent_sequences(identities, sd_multiplier=0.5)
        loose,  _, _ = detect_divergent_sequences(identities, sd_multiplier=3.0)
        assert loose <= strict


# ---------------------------------------------------------------------------
# filter_divergent_sequences
# ---------------------------------------------------------------------------

class TestFilterDivergentSequences:

    def _dirs(self, tmp_path):
        aln       = tmp_path / "alignment"
        sident    = tmp_path / "sident"
        filtered  = tmp_path / "filtered"
        realign   = tmp_path / "realign_input"
        stats     = tmp_path / "stats"
        aln.mkdir(); sident.mkdir()
        return str(aln), str(sident), str(filtered), str(realign), str(stats)

    def test_group_below_min_species_left_untouched(self, tmp_path):
        """Too few sequences to attempt detection at all — copied as-is
        even if the sident data would otherwise flag something."""
        aln, sident, filtered, realign, stats = self._dirs(tmp_path)
        _write_fasta(os.path.join(aln, "OG0001.aln"), [
            ("sp1|g1", "ACGTACGT"), ("sp2|g2", "TTTTTTTT"),
        ])
        _write_sident(os.path.join(sident, "OG0001.aln.ident"), [
            ("sp1|g1", 0.99, "sp2|g2"), ("sp2|g2", 0.01, "sp1|g1"),
        ])
        dropped = filter_divergent_sequences(
            aln, sident, filtered, realign, stats,
            min_species=4, min_sequences=2, sd_multiplier=2.0,
        )
        assert dropped == set()
        ids = _read_ids(os.path.join(filtered, "OG0001.aln"))
        assert set(ids) == {"sp1|g1", "sp2|g2"}
        assert os.listdir(realign) == []

    def test_missing_sident_file_untouched_with_warning(self, tmp_path, caplog):
        aln, sident, filtered, realign, stats = self._dirs(tmp_path)
        _write_fasta(os.path.join(aln, "OG0001.aln"), [
            ("sp1|g1", "ACGT"), ("sp2|g2", "ACGT"),
            ("sp3|g3", "ACGT"), ("sp4|g4", "ACGT"),
        ])
        # No .ident file written for OG0001.

        with caplog.at_level("WARNING", logger="enhydra.filtering"):
            dropped = filter_divergent_sequences(
                aln, sident, filtered, realign, stats,
                min_species=2, min_sequences=2,
            )
        assert dropped == set()
        assert any("OG0001" in r.message for r in caplog.records)
        ids = _read_ids(os.path.join(filtered, "OG0001.aln"))
        assert len(ids) == 4

    def test_insufficient_sident_entries_untouched(self, tmp_path):
        """A sident file with fewer than 2 parseable entries can't support
        mean/SD — group is copied unchanged."""
        aln, sident, filtered, realign, stats = self._dirs(tmp_path)
        _write_fasta(os.path.join(aln, "OG0001.aln"), [
            ("sp1|g1", "ACGT"), ("sp2|g2", "ACGT"),
            ("sp3|g3", "ACGT"), ("sp4|g4", "ACGT"),
        ])
        _write_sident(os.path.join(sident, "OG0001.aln.ident"), [
            ("sp1|g1", 0.9, "sp2|g2"),
        ])
        dropped = filter_divergent_sequences(
            aln, sident, filtered, realign, stats,
            min_species=2, min_sequences=2,
        )
        assert dropped == set()
        assert len(_read_ids(os.path.join(filtered, "OG0001.aln"))) == 4

    def test_no_outliers_group_copied_unchanged_byte_for_byte(self, tmp_path):
        aln, sident, filtered, realign, stats = self._dirs(tmp_path)
        aln_path = os.path.join(aln, "OG0001.aln")
        _write_fasta(aln_path, [
            ("sp1|g1", "AC-GT"), ("sp2|g2", "AC-GT"),
            ("sp3|g3", "AC-GT"), ("sp4|g4", "AC-GT"),
        ])
        _write_sident(os.path.join(sident, "OG0001.aln.ident"), [
            ("sp1|g1", 0.96, "sp2|g2"), ("sp2|g2", 0.97, "sp1|g1"),
            ("sp3|g3", 0.95, "sp1|g1"), ("sp4|g4", 0.98, "sp1|g1"),
        ])
        filter_divergent_sequences(
            aln, sident, filtered, realign, stats,
            min_species=2, min_sequences=2, sd_multiplier=2.0,
        )
        out_path = os.path.join(filtered, "OG0001.aln")
        assert open(aln_path).read() == open(out_path).read()
        assert os.listdir(realign) == []

    def test_outlier_removed_survivors_written_for_realignment(self, tmp_path):
        """Mirrors the real-world 8-sequence example: 1 outlier removed,
        7 survivors written de-gapped to realign_input_dir under the bare
        group_id (no extension), and NOT written to filtered_dir directly."""
        aln, sident, filtered, realign, stats = self._dirs(tmp_path)
        entries = [
            ("NC_011601|WP_000190964.1",   "AC-GT"),
            ("NC_017626|WP_024261658.1",   "ACGGT"),
            ("NC_011748|WP_001390303.1",   "AC-GT"),
            ("NC_018650|WP_001390303.1",   "AC-GT"),
            ("NC_018658|WP_001390303.1",   "AC-GT"),
            ("NC_018661|WP_001390303.1",   "AC-GT"),
            ("NZ_HF572917|WP_001390303.1", "AC-GT"),
            ("NC_017626|WP_001529578.1",   "ACGGT"),
        ]
        _write_fasta(os.path.join(aln, "OG0001.aln"), entries)
        sident_entries = [
            ("NC_011601|WP_000190964.1",   0.3716, "NC_017626|WP_001529578.1"),
            ("NC_017626|WP_024261658.1",   0.9983, "NC_011748|WP_001390303.1"),
            ("NC_011748|WP_001390303.1",   1.0000, "NC_018650|WP_001390303.1"),
            ("NC_018650|WP_001390303.1",   1.0000, "NC_011748|WP_001390303.1"),
            ("NC_018658|WP_001390303.1",   1.0000, "NC_011748|WP_001390303.1"),
            ("NC_018661|WP_001390303.1",   1.0000, "NC_011748|WP_001390303.1"),
            ("NZ_HF572917|WP_001390303.1", 1.0000, "NC_011748|WP_001390303.1"),
            ("NC_017626|WP_001529578.1",   0.9484, "NC_017626|WP_024261658.1"),
        ]
        _write_sident(os.path.join(sident, "OG0001.aln.ident"), sident_entries)

        dropped = filter_divergent_sequences(
            aln, sident, filtered, realign, stats,
            min_species=3, min_sequences=2, sd_multiplier=2.0,
        )

        assert dropped == set()
        # Not written to filtered_dir yet — that only happens once the
        # caller has realigned it.
        assert not os.path.isfile(os.path.join(filtered, "OG0001.aln"))

        realign_path = os.path.join(realign, "OG0001")
        assert os.path.isfile(realign_path)
        survivor_records = list(SeqIO.parse(realign_path, "fasta"))
        survivor_ids = {r.id for r in survivor_records}
        assert survivor_ids == {e[0] for e in entries
                                if e[0] != "NC_011601|WP_000190964.1"}
        # De-gapped: no '-' characters remain in any written sequence.
        assert all("-" not in str(r.seq) for r in survivor_records)

    def test_group_dropped_when_below_min_species_after_pruning(self, tmp_path):
        aln, sident, filtered, realign, stats = self._dirs(tmp_path)
        _write_fasta(os.path.join(aln, "OG0001.aln"), [
            ("sp1|g1", "ACGT"), ("sp2|g2", "ACGT"), ("sp3|g3", "ACGT"),
        ])
        _write_sident(os.path.join(sident, "OG0001.aln.ident"), [
            ("sp1|g1", 0.1, "sp2|g2"),   # clear outlier
            ("sp2|g2", 0.95, "sp3|g3"),
            ("sp3|g3", 0.95, "sp2|g2"),
        ])
        # min_species=3: after removing sp1's sequence only 2 species remain.
        dropped = filter_divergent_sequences(
            aln, sident, filtered, realign, stats,
            min_species=3, min_sequences=2, sd_multiplier=1.0,
        )
        assert dropped == {"OG0001"}
        assert not os.path.isfile(os.path.join(filtered, "OG0001.aln"))
        assert not os.path.isfile(os.path.join(realign, "OG0001"))

    def test_group_dropped_when_below_min_sequences_but_species_ok(self, tmp_path):
        """Distinguishes the min_sequences check from the min_species check:
        a paralog means the species is still represented after pruning, but
        the raw sequence count drops below min_sequences."""
        aln, sident, filtered, realign, stats = self._dirs(tmp_path)
        _write_fasta(os.path.join(aln, "OG0001.aln"), [
            ("sp1|g1", "ACGT"),   # sp1, normal
            ("sp1|g2", "TTTT"),   # sp1 paralog, divergent outlier
            ("sp2|g3", "ACGT"),
            ("sp3|g4", "ACGT"),
        ])
        _write_sident(os.path.join(sident, "OG0001.aln.ident"), [
            ("sp1|g1", 0.98, "sp2|g3"),
            ("sp1|g2", 0.05, "sp1|g1"),   # clear outlier
            ("sp2|g3", 0.97, "sp1|g1"),
            ("sp3|g4", 0.96, "sp1|g1"),
        ])
        # species count after pruning (sp1 via g1, sp2, sp3) = 3, meets
        # min_species=3, but sequence count after pruning = 3 < min_sequences=4.
        dropped = filter_divergent_sequences(
            aln, sident, filtered, realign, stats,
            min_species=3, min_sequences=4, sd_multiplier=1.0,
        )
        assert dropped == {"OG0001"}

    def test_drop_reasons_tsv_content(self, tmp_path):
        aln, sident, filtered, realign, stats = self._dirs(tmp_path)
        _write_fasta(os.path.join(aln, "OG0001.aln"), [
            ("sp1|g1", "ACGT"), ("sp2|g2", "ACGT"),
            ("sp3|g3", "ACGT"), ("sp4|g4", "ACGT"),
        ])
        _write_sident(os.path.join(sident, "OG0001.aln.ident"), [
            ("sp1|g1", 0.1, "sp2|g2"),
            ("sp2|g2", 0.95, "sp3|g3"),
            ("sp3|g3", 0.95, "sp2|g2"),
            ("sp4|g4", 0.96, "sp2|g2"),
        ])
        filter_divergent_sequences(
            aln, sident, filtered, realign, stats,
            min_species=2, min_sequences=2, sd_multiplier=1.0,
        )
        drop_path = os.path.join(stats, "drop_reasons.tsv")
        assert os.path.isfile(drop_path)
        with open(drop_path) as fh:
            lines = [l.rstrip("\n").split("\t") for l in fh if l.strip()]
        header, rows = lines[0], lines[1:]
        assert header == ["group_id", "sequence_id", "identity_to_closest",
                          "pct_diff_from_avg", "reason"]
        assert len(rows) == 1
        assert rows[0][0] == "OG0001"
        assert rows[0][1] == "sp1|g1"
        assert rows[0][4] == "removed_divergent_sequence"
        assert float(rows[0][2]) == pytest.approx(0.1, abs=1e-4)

    def test_drop_reasons_tsv_for_dropped_group(self, tmp_path):
        aln, sident, filtered, realign, stats = self._dirs(tmp_path)
        _write_fasta(os.path.join(aln, "OG0001.aln"), [
            ("sp1|g1", "ACGT"), ("sp2|g2", "ACGT"), ("sp3|g3", "ACGT"),
        ])
        _write_sident(os.path.join(sident, "OG0001.aln.ident"), [
            ("sp1|g1", 0.1, "sp2|g2"),
            ("sp2|g2", 0.95, "sp3|g3"),
            ("sp3|g3", 0.95, "sp2|g2"),
        ])
        filter_divergent_sequences(
            aln, sident, filtered, realign, stats,
            min_species=3, min_sequences=2, sd_multiplier=1.0,
        )
        with open(os.path.join(stats, "drop_reasons.tsv")) as fh:
            rows = [l.rstrip("\n").split("\t") for l in fh if l.strip()][1:]
        assert len(rows) == 1
        assert rows[0][0] == "OG0001"
        assert rows[0][1] == "sp1|g1"   # comma-joined list of flagged seqs
        assert rows[0][2] == ""
        assert rows[0][3] == ""
        assert rows[0][4] == "below_min_species_after_divergence_filter"

    def test_multiple_groups_mixed_outcomes(self, tmp_path):
        """One untouched, one pruned-and-realignable, one dropped — all in
        a single call, verifying they don't interfere with each other."""
        aln, sident, filtered, realign, stats = self._dirs(tmp_path)

        # OG_clean: no outliers.
        _write_fasta(os.path.join(aln, "OG_clean.aln"), [
            ("sp1|g1", "ACGT"), ("sp2|g2", "ACGT"),
            ("sp3|g3", "ACGT"), ("sp4|g4", "ACGT"),
        ])
        # Identical identity values (zero variance) so this group is
        # guaranteed clean under any sd_multiplier — a near-clean but
        # non-identical spread (e.g. 0.94-0.97) can still trip a strict
        # multiplier like the 1.0 used below, which isn't what this test
        # is trying to exercise.
        _write_sident(os.path.join(sident, "OG_clean.aln.ident"), [
            ("sp1|g1", 0.96, "sp2|g2"), ("sp2|g2", 0.96, "sp1|g1"),
            ("sp3|g3", 0.96, "sp1|g1"), ("sp4|g4", 0.96, "sp1|g1"),
        ])

        # OG_prune: one outlier, group survives pruning.
        _write_fasta(os.path.join(aln, "OG_prune.aln"), [
            ("sp1|g1", "ACGT"), ("sp2|g2", "ACGT"),
            ("sp3|g3", "ACGT"), ("sp4|g4", "ACGT"), ("sp5|g5", "ACGT"),
        ])
        _write_sident(os.path.join(sident, "OG_prune.aln.ident"), [
            ("sp1|g1", 0.1, "sp2|g2"),
            ("sp2|g2", 0.95, "sp3|g3"), ("sp3|g3", 0.96, "sp2|g2"),
            ("sp4|g4", 0.94, "sp2|g2"), ("sp5|g5", 0.97, "sp2|g2"),
        ])

        # OG_drop: one outlier, group falls below min_species after pruning.
        _write_fasta(os.path.join(aln, "OG_drop.aln"), [
            ("sp1|g1", "ACGT"), ("sp2|g2", "ACGT"), ("sp3|g3", "ACGT"),
        ])
        _write_sident(os.path.join(sident, "OG_drop.aln.ident"), [
            ("sp1|g1", 0.1, "sp2|g2"),
            ("sp2|g2", 0.95, "sp3|g3"), ("sp3|g3", 0.96, "sp2|g2"),
        ])

        dropped = filter_divergent_sequences(
            aln, sident, filtered, realign, stats,
            min_species=3, min_sequences=2, sd_multiplier=1.0,
        )

        assert dropped == {"OG_drop"}
        assert os.path.isfile(os.path.join(filtered, "OG_clean.aln"))
        assert not os.path.isfile(os.path.join(filtered, "OG_prune.aln"))
        assert os.path.isfile(os.path.join(realign, "OG_prune"))
        assert not os.path.isfile(os.path.join(realign, "OG_drop"))
        assert not os.path.isfile(os.path.join(filtered, "OG_drop.aln"))

    def test_returns_set_of_dropped_group_ids(self, tmp_path):
        aln, sident, filtered, realign, stats = self._dirs(tmp_path)
        _write_fasta(os.path.join(aln, "OG0001.aln"), [
            ("sp1|g1", "ACGT"), ("sp2|g2", "ACGT"), ("sp3|g3", "ACGT"),
        ])
        _write_sident(os.path.join(sident, "OG0001.aln.ident"), [
            ("sp1|g1", 0.1, "sp2|g2"),
            ("sp2|g2", 0.95, "sp3|g3"), ("sp3|g3", 0.96, "sp2|g2"),
        ])
        dropped = filter_divergent_sequences(
            aln, sident, filtered, realign, stats,
            min_species=3, min_sequences=2, sd_multiplier=1.0,
        )
        assert isinstance(dropped, set)
        assert dropped == {"OG0001"}


# ---------------------------------------------------------------------------
# display_group_id
# ---------------------------------------------------------------------------

class TestDisplayGroupId:

    def test_strips_lengthfilter_suffix(self):
        assert display_group_id("OG0001_lengthfilter") == "OG0001"

    def test_no_suffix_returned_unchanged(self):
        """A group_id that never carried the suffix (e.g. one recorded in
        skipped_groups.tsv, which never reaches length_filter_dir) must be
        returned exactly as given."""
        assert display_group_id("OG0001") == "OG0001"

    def test_only_strips_trailing_suffix_not_mid_string_occurrence(self):
        """The literal substring '_lengthfilter' appearing somewhere other
        than as the trailing suffix must not be stripped — only an exact
        trailing match counts."""
        assert display_group_id("OG_lengthfilter_extra") == "OG_lengthfilter_extra"

    def test_idempotent(self):
        """Applying twice must not strip anything further (there is only
        ever one suffix to strip)."""
        once  = display_group_id("OG0001_lengthfilter")
        twice = display_group_id(once)
        assert once == twice == "OG0001"

    def test_matches_actual_length_filter_dir_filename(self, tmp_path):
        """End-to-end: applying display_group_id() to the real filename
        filter_length() writes to length_filter_dir must recover the
        original input group name."""
        src = tmp_path / "input" / "OG0007"
        src.parent.mkdir()
        _write_fasta(str(src), [("sp1|g1", "A" * 100), ("sp2|g2", "A" * 100)])
        stats_dir  = tmp_path / "stats"
        filter_dir = tmp_path / "filter"
        stats_dir.mkdir(); filter_dir.mkdir()

        filter_length(str(src), str(stats_dir), str(filter_dir))

        written_files = os.listdir(str(filter_dir))
        assert len(written_files) == 1
        assert display_group_id(written_files[0]) == "OG0007"
