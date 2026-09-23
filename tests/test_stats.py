import os
import json
import pytest

from enhydra.stats import aggregate_pipeline_stats, compute_differential_stats


def _write_tsv(path, header, rows):
    with open(path, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for row in rows:
            fh.write("\t".join(str(c) for c in row) + "\n")


class TestAggregatePipelineStats:

    def _make_listdir(self, tmp_path):
        d = tmp_path / "list1"
        for sub in ("length_filter", "length_filter_stats",
                    "group_filter", "group_filter_stats", "tables"):
            (d / sub).mkdir(parents=True)
        return d

    def test_basic_funnel_counts(self, tmp_path):
        d = self._make_listdir(tmp_path)
        for f in ("OG1", "OG2", "OG3"):
            (d / "length_filter" / f).touch()
        for f in ("OG1", "OG2"):
            (d / "group_filter" / f).touch()
        _write_tsv(d / "tables" / "group2mean.tsv", [], [("OG1", 0.9), ("OG2", 0.8)])
        _write_tsv(d / "tables" / "anchor2mean.tsv", [], [("gA", 0.9)])

        summary = aggregate_pipeline_stats(str(d), n_input=5)

        assert summary["n_input"] == 5
        assert summary["n_after_length_filter"] == 3
        assert summary["n_after_group_filter"] == 2
        assert summary["n_final_groups"] == 2
        assert summary["n_anchor_mapped"] == 1

    def test_drop_reason_breakdown(self, tmp_path):
        d = self._make_listdir(tmp_path)
        _write_tsv(d / "length_filter_stats" / "skipped_groups.tsv",
                  ["group_id", "reason", "detail"],
                  [("OG1", "empty_file", ""), ("OG2", "single_sequence", "")])
        _write_tsv(d / "length_filter_stats" / "drop_reasons.tsv",
                  ["group_id", "sequence_id", "length", "pct_diff_from_avg", "reason"],
                  [("OG3", "sp1|g1", 400, 2.5, "removed_above_max")])
        _write_tsv(d / "group_filter_stats" / "drop_reasons.tsv",
                  ["group_id", "reason", "detail"],
                  [("OG4", "missing_anchor", "")])
        _write_tsv(d / "tables" / "drop_reasons.tsv",
                  ["group_id", "reason", "detail"],
                  [("OG5", "no_average_identity", "")])
        _write_tsv(d / "tables" / "group2mean.tsv", [], [])

        summary = aggregate_pipeline_stats(str(d), n_input=10)
        ds = summary["drop_summary"]
        assert ds["length_filter"]["groups_skipped"] == {
            "empty_file": 1, "single_sequence": 1,
        }
        assert ds["length_filter"]["sequences_removed"] == {"removed_above_max": 1}
        assert ds["group_filter"] == {"missing_anchor": 1}
        assert ds["tables"] == {"no_average_identity": 1}

    def test_missing_stats_files_produce_empty_breakdowns(self, tmp_path):
        d = self._make_listdir(tmp_path)
        summary = aggregate_pipeline_stats(str(d), n_input=0)
        assert summary["drop_summary"]["length_filter"]["groups_skipped"] == {}
        assert summary["drop_summary"]["group_filter"] == {}

    def test_writes_json_file(self, tmp_path):
        d = self._make_listdir(tmp_path)
        aggregate_pipeline_stats(str(d), n_input=1)
        out = d / "pipeline_stats.json"
        assert out.is_file()
        with open(out) as fh:
            data = json.load(fh)
        assert data["n_input"] == 1

    def test_custom_output_path(self, tmp_path):
        d = self._make_listdir(tmp_path)
        custom = tmp_path / "custom.json"
        aggregate_pipeline_stats(str(d), n_input=1, output_path=str(custom))
        assert custom.is_file()
        assert not (d / "pipeline_stats.json").is_file()


class TestComputeDifferentialStats:

    def test_overlap_counts(self, tmp_path):
        t1 = tmp_path / "t1"; t2 = tmp_path / "t2"
        t1.mkdir(); t2.mkdir()
        _write_tsv(t1 / "group2mean.tsv", [], [("OG1", 0.9), ("OG2", 0.8), ("OG3", 0.7)])
        _write_tsv(t2 / "group2mean.tsv", [], [("OG2", 0.5), ("OG3", 0.6), ("OG4", 0.4)])

        summary = compute_differential_stats(str(t1), str(t2),
                                             list1_name="Pathogenic",
                                             list2_name="Non-pathogenic")
        assert summary["n_groups_list1"] == 3
        assert summary["n_groups_list2"] == 3
        assert summary["n_common_groups"] == 2
        assert summary["n_list1_only"] == 1
        assert summary["n_list2_only"] == 1
        assert summary["list1_name"] == "Pathogenic"
        assert summary["list2_name"] == "Non-pathogenic"

    def test_writes_json_when_output_path_given(self, tmp_path):
        t1 = tmp_path / "t1"; t2 = tmp_path / "t2"
        t1.mkdir(); t2.mkdir()
        _write_tsv(t1 / "group2mean.tsv", [], [("OG1", 0.9)])
        _write_tsv(t2 / "group2mean.tsv", [], [("OG1", 0.5)])
        out = tmp_path / "differential_stats.json"

        compute_differential_stats(str(t1), str(t2), output_path=str(out))
        assert out.is_file()

    def test_no_file_written_without_output_path(self, tmp_path):
        t1 = tmp_path / "t1"; t2 = tmp_path / "t2"
        t1.mkdir(); t2.mkdir()
        _write_tsv(t1 / "group2mean.tsv", [], [("OG1", 0.9)])
        _write_tsv(t2 / "group2mean.tsv", [], [("OG1", 0.5)])

        summary = compute_differential_stats(str(t1), str(t2))
        assert summary["n_common_groups"] == 1
        assert not (tmp_path / "differential_stats.json").is_file()
