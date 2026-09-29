"""Tests for enhydra.stage_markers."""

import os
import json

import pytest

from enhydra.stage_markers import (
    MARKER_FILENAME,
    write_stage_marker,
    check_stage,
    aggregate_stage_markers,
    list_data_files,
)


# ---------------------------------------------------------------------------
# write_stage_marker
# ---------------------------------------------------------------------------

class TestWriteStageMarker:

    def test_writes_marker_file(self, tmp_path):
        d = tmp_path / "group_filter"
        write_stage_marker(str(d), "group_filter", {"min_species": 4})
        marker_path = d / MARKER_FILENAME
        assert marker_path.is_file()

    def test_creates_directory_if_missing(self, tmp_path):
        d = tmp_path / "not_yet_created"
        write_stage_marker(str(d), "alignment", {"aligner": "mafft"})
        assert d.is_dir()

    def test_marker_content(self, tmp_path):
        d = tmp_path / "tables"
        write_stage_marker(str(d), "tables", {"anchor": "hsapiens"})
        with open(d / MARKER_FILENAME) as fh:
            marker = json.load(fh)
        assert marker["stage"] == "tables"
        assert marker["parameters"] == {"anchor": "hsapiens"}
        assert "completed_at" in marker
        assert "enhydra_version" in marker

    def test_overwrites_existing_marker(self, tmp_path):
        d = tmp_path / "length_filter"
        write_stage_marker(str(d), "length_filter", {"length_filter_sd": 2.0})
        write_stage_marker(str(d), "length_filter", {"length_filter_sd": 1.5})
        with open(d / MARKER_FILENAME) as fh:
            marker = json.load(fh)
        assert marker["parameters"] == {"length_filter_sd": 1.5}


# ---------------------------------------------------------------------------
# check_stage
# ---------------------------------------------------------------------------

def _touch(path):
    with open(path, "w") as fh:
        fh.write("x")


class TestCheckStageNotRun:

    def test_nonexistent_directory(self, tmp_path):
        d = tmp_path / "nope"
        result = check_stage([(str(d), None)], str(d), {"a": 1})
        assert result.status == "not_run"

    def test_empty_directory(self, tmp_path):
        d = tmp_path / "empty"
        d.mkdir()
        result = check_stage([(str(d), None)], str(d), {"a": 1})
        assert result.status == "not_run"

    def test_directory_with_only_marker_file_is_not_run(self, tmp_path):
        """A directory containing nothing but its own (perhaps stale)
        marker, with no real output alongside it, must not be mistaken
        for 'has output' — the marker file itself is excluded from the
        generic has-output check."""
        d = tmp_path / "weird"
        d.mkdir()
        write_stage_marker(str(d), "weird", {"a": 1})
        result = check_stage([(str(d), None)], str(d), {"a": 1})
        assert result.status == "not_run"

    def test_sentinel_files_missing(self, tmp_path):
        d = tmp_path / "tables"
        d.mkdir()
        _touch(d / "group2mean.tsv")
        # anchor2mean.tsv and group2anchor.tsv missing.
        result = check_stage(
            [(str(d), ["group2mean.tsv", "anchor2mean.tsv", "group2anchor.tsv"])],
            str(d), {"anchor": "hsapiens"},
        )
        assert result.status == "not_run"

    def test_sentinel_file_empty_counts_as_missing(self, tmp_path):
        d = tmp_path / "tables"
        d.mkdir()
        for f in ("group2mean.tsv", "anchor2mean.tsv", "group2anchor.tsv"):
            (d / f).touch()   # zero-byte
        result = check_stage(
            [(str(d), ["group2mean.tsv", "anchor2mean.tsv", "group2anchor.tsv"])],
            str(d), {"anchor": "hsapiens"},
        )
        assert result.status == "not_run"


class TestCheckStageOk:

    def test_matching_parameters(self, tmp_path):
        d = tmp_path / "group_filter"
        d.mkdir()
        _touch(d / "OG0001")
        params = {"min_species": 4, "anchor": "hsapiens"}
        write_stage_marker(str(d), "group_filter", params)
        result = check_stage([(str(d), None)], str(d), params)
        assert result.status == "ok"
        assert result.changed_parameters is None

    def test_sentinel_based_ok(self, tmp_path):
        d = tmp_path / "tables"
        d.mkdir()
        for f in ("group2mean.tsv", "anchor2mean.tsv", "group2anchor.tsv"):
            (d / f).write_text("data")
        params = {"anchor": "hsapiens", "require_anchor_in_tables": True}
        write_stage_marker(str(d), "tables", params)
        result = check_stage(
            [(str(d), ["group2mean.tsv", "anchor2mean.tsv", "group2anchor.tsv"])],
            str(d), params,
        )
        assert result.status == "ok"

    def test_multi_dir_spec_all_present_and_matching(self, tmp_path):
        filtered_dir = tmp_path / "alignment_divergence_filtered"
        stats_dir    = tmp_path / "divergence_filter_stats"
        filtered_dir.mkdir()
        stats_dir.mkdir()
        _touch(filtered_dir / "OG0001.aln")
        (stats_dir / "drop_reasons.tsv").write_text("group_id\treason\n")
        params = {"divergence_filter_sd": 2.0, "min_species": 4, "min_sequences": 2}
        write_stage_marker(str(filtered_dir), "divergence_filter", params)
        result = check_stage(
            [(str(filtered_dir), None), (str(stats_dir), ["drop_reasons.tsv"])],
            str(filtered_dir), params,
        )
        assert result.status == "ok"


class TestCheckStageMismatch:

    def test_single_changed_parameter(self, tmp_path):
        d = tmp_path / "group_filter"
        d.mkdir()
        _touch(d / "OG0001")
        write_stage_marker(str(d), "group_filter", {"min_species": 4, "anchor": "hsapiens"})
        result = check_stage(
            [(str(d), None)], str(d), {"min_species": 8, "anchor": "hsapiens"},
        )
        assert result.status == "mismatch"
        assert result.changed_parameters == {"min_species": (4, 8)}

    def test_multiple_changed_parameters(self, tmp_path):
        d = tmp_path / "group_filter"
        d.mkdir()
        _touch(d / "OG0001")
        write_stage_marker(
            str(d), "group_filter",
            {"min_species": 4, "paralog_mode": "all"},
        )
        result = check_stage(
            [(str(d), None)], str(d),
            {"min_species": 8, "paralog_mode": "remove"},
        )
        assert result.status == "mismatch"
        assert result.changed_parameters == {
            "min_species":  (4, 8),
            "paralog_mode": ("all", "remove"),
        }

    def test_added_parameter_key_counts_as_mismatch(self, tmp_path):
        """A parameter key present now but absent from the old marker
        (e.g. a new tunable added to a stage in a later ENHYDRA version)
        must be detected, not silently ignored."""
        d = tmp_path / "group_filter"
        d.mkdir()
        _touch(d / "OG0001")
        write_stage_marker(str(d), "group_filter", {"min_species": 4})
        result = check_stage(
            [(str(d), None)], str(d), {"min_species": 4, "require_anchor": True},
        )
        assert result.status == "mismatch"
        assert result.changed_parameters == {"require_anchor": ("<unset>", True)}

    def test_removed_parameter_key_counts_as_mismatch(self, tmp_path):
        d = tmp_path / "group_filter"
        d.mkdir()
        _touch(d / "OG0001")
        write_stage_marker(
            str(d), "group_filter", {"min_species": 4, "require_anchor": True},
        )
        result = check_stage([(str(d), None)], str(d), {"min_species": 4})
        assert result.status == "mismatch"
        assert result.changed_parameters == {"require_anchor": (True, "<unset>")}


class TestCheckStageMissingMarker:

    def test_output_present_no_marker_file(self, tmp_path):
        """Simulates a run from before this feature existed: real output
        is present but no marker was ever written."""
        d = tmp_path / "alignment"
        d.mkdir()
        _touch(d / "OG0001.aln")
        result = check_stage([(str(d), None)], str(d), {"aligner": "mafft"})
        assert result.status == "missing_marker"

    def test_corrupted_marker_file_treated_as_missing(self, tmp_path):
        d = tmp_path / "alignment"
        d.mkdir()
        _touch(d / "OG0001.aln")
        with open(d / MARKER_FILENAME, "w") as fh:
            fh.write("{not valid json")
        result = check_stage([(str(d), None)], str(d), {"aligner": "mafft"})
        assert result.status == "missing_marker"

    def test_marker_missing_parameters_key_treated_as_missing(self, tmp_path):
        d = tmp_path / "alignment"
        d.mkdir()
        _touch(d / "OG0001.aln")
        with open(d / MARKER_FILENAME, "w") as fh:
            json.dump({"stage": "alignment"}, fh)   # no 'parameters' key
        result = check_stage([(str(d), None)], str(d), {"aligner": "mafft"})
        assert result.status == "missing_marker"


# ---------------------------------------------------------------------------
# aggregate_stage_markers
# ---------------------------------------------------------------------------

class TestAggregateStageMarkers:

    def test_collects_all_markers(self, tmp_path):
        root = tmp_path / "outdir"
        root.mkdir()
        write_stage_marker(str(root / "length_filter"), "length_filter",
                          {"length_filter_sd": 2.0})
        write_stage_marker(str(root / "group_filter"), "group_filter",
                          {"min_species": 4})
        summary = aggregate_stage_markers(str(root))
        assert set(summary.keys()) == {"length_filter", "group_filter"}
        assert summary["length_filter"]["parameters"] == {"length_filter_sd": 2.0}
        assert "completed_at" in summary["length_filter"]

    def test_ignores_directories_without_markers(self, tmp_path):
        root = tmp_path / "outdir"
        root.mkdir()
        (root / "some_scratch_dir").mkdir()
        write_stage_marker(str(root / "tables"), "tables", {"anchor": "a"})
        summary = aggregate_stage_markers(str(root))
        assert set(summary.keys()) == {"tables"}

    def test_nonexistent_root_returns_empty_dict(self, tmp_path):
        assert aggregate_stage_markers(str(tmp_path / "nope")) == {}

    def test_no_markers_at_all_returns_empty_dict(self, tmp_path):
        root = tmp_path / "outdir"
        root.mkdir()
        (root / "some_dir").mkdir()
        _touch(root / "some_dir" / "file.txt")
        assert aggregate_stage_markers(str(root)) == {}

    def test_corrupted_marker_skipped_with_warning(self, tmp_path, caplog):
        root = tmp_path / "outdir"
        root.mkdir()
        bad_dir = root / "bad_stage"
        bad_dir.mkdir()
        with open(bad_dir / MARKER_FILENAME, "w") as fh:
            fh.write("{not valid json")
        write_stage_marker(str(root / "good_stage"), "good_stage", {"x": 1})

        with caplog.at_level("WARNING", logger="enhydra.stage_markers"):
            summary = aggregate_stage_markers(str(root))

        assert "good_stage" in summary
        assert "bad_stage" not in summary
        assert any("Could not parse stage marker" in r.message
                  for r in caplog.records)

    def test_uses_marker_stage_field_not_directory_name(self, tmp_path):
        """The summary is keyed by the marker's own recorded 'stage' name,
        not by the directory's filename, in case they ever diverge."""
        root = tmp_path / "outdir"
        root.mkdir()
        write_stage_marker(
            str(root / "alignment_divergence_filtered"),
            "divergence_filter", {"divergence_filter_sd": 2.0},
        )
        summary = aggregate_stage_markers(str(root))
        assert "divergence_filter" in summary
        assert "alignment_divergence_filtered" not in summary


# ---------------------------------------------------------------------------
# list_data_files
# ---------------------------------------------------------------------------
#
# Regression coverage for a real crash: write_stage_marker() writes
# MARKER_FILENAME directly into a stage's own output directory, and any
# *later* stage that reads that same directory wholesale via plain
# os.listdir() (to process "every file in it" as pipeline data, as
# opposed to merely checking whether the directory has any output at all)
# will pick up the marker file too and typically choke on it — e.g.
# passing '.stage_marker.json' to trimAl as if it were an alignment file.
# list_data_files() is the drop-in os.listdir() replacement every such
# call site in the package now uses instead.

class TestListDataFiles:

    def test_excludes_marker_file(self, tmp_path):
        d = tmp_path / "alignment"
        d.mkdir()
        _touch(d / "OG0001.aln")
        _touch(d / "OG0002.aln")
        write_stage_marker(str(d), "alignment", {"aligner": "mafft"})

        files = list_data_files(str(d))

        assert MARKER_FILENAME not in files
        assert set(files) == {"OG0001.aln", "OG0002.aln"}

    def test_no_marker_present_returns_everything(self, tmp_path):
        d = tmp_path / "alignment"
        d.mkdir()
        _touch(d / "OG0001.aln")
        assert list_data_files(str(d)) == ["OG0001.aln"]

    def test_empty_directory_returns_empty_list(self, tmp_path):
        d = tmp_path / "empty"
        d.mkdir()
        assert list_data_files(str(d)) == []

    def test_directory_with_only_marker_returns_empty_list(self, tmp_path):
        d = tmp_path / "alignment"
        d.mkdir()
        write_stage_marker(str(d), "alignment", {"aligner": "mafft"})
        assert list_data_files(str(d)) == []

    def test_does_not_filter_files_merely_containing_dots(self, tmp_path):
        """Only the exact MARKER_FILENAME is excluded — a real data file
        whose own name happens to contain dots (e.g. the '.fa' extension
        bug fixed previously) must not be filtered out."""
        d = tmp_path / "alignment"
        d.mkdir()
        _touch(d / "OG0001.fa_lengthfilter.aln")
        write_stage_marker(str(d), "alignment", {"aligner": "mafft"})
        files = list_data_files(str(d))
        assert files == ["OG0001.fa_lengthfilter.aln"]
