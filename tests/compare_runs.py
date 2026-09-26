"""Tests for enhydra.compare_runs."""

import os
import json

import pandas as pd
import pytest

from enhydra.compare_runs import (
    _discover_run,
    _load_run_metric_df,
    _merge_pair,
    _sankey_counts,
    _spearman_matrix,
    _diff_rows,
    _build_run_params_diff_html,
    _render_correlation_scatter_svg,
    _render_sankey_svg,
    _pairwise_correlation_block,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _write_gsea_csv(path: str, rows: list[tuple]):
    df = pd.DataFrame(rows, columns=["Term", "NES", "NOM p-val", "FDR q-val"])
    df.to_csv(path, index=False)


def _write_run_parameters(path: str, record: dict):
    with open(path, "w") as fh:
        json.dump(record, fh)


def _minimal_run_parameters(two_list=False, all_metrics=False, metrics_run=None):
    metrics_run = metrics_run if metrics_run is not None else ["zscore"]
    record = {
        "enhydra_version": "0.1.0",
        "input": {"two_list_mode": two_list, "anchor": "hsapiens"},
        "ranking": {"all_metrics": all_metrics, "metrics_run": metrics_run},
        "filtering": {"min_species": 4},
    }
    if two_list:
        record["input"]["list1"] = {"name": "List 1", "n_species": 3}
        record["input"]["list2"] = {"name": "List 2", "n_species": 2}
    return record


# ---------------------------------------------------------------------------
# _discover_run
# ---------------------------------------------------------------------------

class TestDiscoverRunSingleList:

    def test_simple_single_metric(self, tmp_path):
        outdir = tmp_path / "run1"
        outdir.mkdir()
        _write_run_parameters(
            str(outdir / "run_parameters.json"),
            _minimal_run_parameters(two_list=False, all_metrics=False,
                                    metrics_run=["zscore"]),
        )
        res_dir = outdir / "enrichment"
        res_dir.mkdir()
        _write_gsea_csv(str(res_dir / "gseapy.gene_set.prerank.report.csv"),
                        [("GO:0001", 1.5, 0.01, 0.1)])
        (outdir / "tables").mkdir()

        info = _discover_run(str(outdir))
        assert info["two_list_mode"] is False
        assert info["all_metrics"] is False
        assert set(info["metric_paths"].keys()) == {"zscore"}
        assert info["metric_paths"]["zscore"]["tables_dir2"] is None

    def test_all_metrics(self, tmp_path):
        outdir = tmp_path / "run2"
        outdir.mkdir()
        _write_run_parameters(
            str(outdir / "run_parameters.json"),
            _minimal_run_parameters(all_metrics=True,
                                    metrics_run=["identity", "zscore", "rank"]),
        )
        for m in ("identity", "zscore", "rank"):
            d = outdir / ("enrichment_%s" % m)
            d.mkdir()
            _write_gsea_csv(str(d / "gseapy.gene_set.prerank.report.csv"),
                            [("GO:0001", 1.0, 0.05, 0.2)])

        info = _discover_run(str(outdir))
        assert set(info["metric_paths"].keys()) == {"identity", "zscore", "rank"}

    def test_missing_csv_for_declared_metric_skipped(self, tmp_path):
        outdir = tmp_path / "run3"
        outdir.mkdir()
        _write_run_parameters(
            str(outdir / "run_parameters.json"),
            _minimal_run_parameters(metrics_run=["zscore"]),
        )
        # No enrichment/ dir created at all.
        info = _discover_run(str(outdir))
        assert info["metric_paths"] == {}


class TestDiscoverRunTwoList:

    def test_two_list_all_metrics(self, tmp_path):
        outdir = tmp_path / "run_diff"
        outdir.mkdir()
        _write_run_parameters(
            str(outdir / "run_parameters.json"),
            _minimal_run_parameters(two_list=True, all_metrics=True,
                                    metrics_run=["identity", "zscore"]),
        )
        for m in ("identity", "zscore"):
            d = outdir / ("differential_%s" % m) / "enrichment"
            d.mkdir(parents=True)
            _write_gsea_csv(str(d / "gseapy.gene_set.prerank.report.csv"),
                            [("GO:0001", 0.5, 0.2, 0.3)])
        (outdir / "list1" / "tables").mkdir(parents=True)
        (outdir / "list2" / "tables").mkdir(parents=True)

        info = _discover_run(str(outdir))
        assert info["two_list_mode"] is True
        assert set(info["metric_paths"].keys()) == {"identity", "zscore"}
        assert info["metric_paths"]["identity"]["tables_dir1"].endswith(
            os.path.join("list1", "tables")
        )
        assert info["metric_paths"]["identity"]["tables_dir2"].endswith(
            os.path.join("list2", "tables")
        )


class TestDiscoverRunFallback:

    def test_missing_run_parameters_sniffs_single_list(self, tmp_path):
        outdir = tmp_path / "old_run"
        outdir.mkdir()
        res_dir = outdir / "enrichment"
        res_dir.mkdir()
        _write_gsea_csv(str(res_dir / "gseapy.gene_set.prerank.report.csv"),
                        [("GO:0001", 1.0, 0.05, 0.2)])

        info = _discover_run(str(outdir))
        assert info["run_parameters"] is None
        assert info["two_list_mode"] is False
        assert set(info["metric_paths"].keys()) == {"unknown"}

    def test_missing_run_parameters_sniffs_all_metrics(self, tmp_path):
        outdir = tmp_path / "old_run_am"
        outdir.mkdir()
        for m in ("identity", "rank"):
            d = outdir / ("enrichment_%s" % m)
            d.mkdir()
            _write_gsea_csv(str(d / "gseapy.gene_set.prerank.report.csv"),
                            [("GO:0001", 1.0, 0.05, 0.2)])

        info = _discover_run(str(outdir))
        assert info["all_metrics"] is True
        assert set(info["metric_paths"].keys()) == {"identity", "rank"}

    def test_missing_run_parameters_sniffs_two_list(self, tmp_path):
        outdir = tmp_path / "old_run_diff"
        outdir.mkdir()
        (outdir / "list1").mkdir()
        (outdir / "list2").mkdir()
        d = outdir / "differential" / "enrichment"
        d.mkdir(parents=True)
        _write_gsea_csv(str(d / "gseapy.gene_set.prerank.report.csv"),
                        [("GO:0001", 1.0, 0.05, 0.2)])

        info = _discover_run(str(outdir))
        assert info["two_list_mode"] is True
        assert set(info["metric_paths"].keys()) == {"unknown"}


# ---------------------------------------------------------------------------
# _load_run_metric_df
# ---------------------------------------------------------------------------

class TestLoadRunMetricDf:

    def test_loads_expected_columns(self, tmp_path):
        d = tmp_path / "enrichment"
        d.mkdir()
        _write_gsea_csv(str(d / "gseapy.gene_set.prerank.report.csv"),
                        [("GO:0001", 1.2, 0.01, 0.05), ("GO:0002", -0.8, 0.2, 0.4)])
        df = _load_run_metric_df({"results_dir": str(d)})
        assert list(df.columns) == ["Term", "NES", "FDR q-val"]
        assert len(df) == 2

    def test_missing_file_returns_none(self, tmp_path):
        df = _load_run_metric_df({"results_dir": str(tmp_path / "nope")})
        assert df is None

    def test_missing_required_column_returns_none(self, tmp_path):
        d = tmp_path / "enrichment"
        d.mkdir()
        pd.DataFrame({"Term": ["GO:0001"], "NES": [1.0]}).to_csv(
            str(d / "gseapy.gene_set.prerank.report.csv"), index=False
        )
        df = _load_run_metric_df({"results_dir": str(d)})
        assert df is None


# ---------------------------------------------------------------------------
# _merge_pair / _sankey_counts / _spearman_matrix
# ---------------------------------------------------------------------------

def _df(rows):
    return pd.DataFrame(rows, columns=["Term", "NES", "FDR q-val"])


class TestMergePair:

    def test_inner_join_on_term(self):
        a = _df([("GO:1", 1.0, 0.1), ("GO:2", 0.5, 0.3)])
        b = _df([("GO:1", 0.9, 0.2), ("GO:3", -1.0, 0.05)])
        merged = _merge_pair(a, b)
        assert set(merged["Term"]) == {"GO:1"}
        assert "FDR_a" in merged.columns and "FDR_b" in merged.columns

    def test_no_overlap_returns_empty(self):
        a = _df([("GO:1", 1.0, 0.1)])
        b = _df([("GO:2", 1.0, 0.1)])
        merged = _merge_pair(a, b)
        assert merged.empty


class TestSankeyCounts:

    def test_all_four_categories(self):
        a = _df([("GO:1", 1.0, 0.01), ("GO:2", 1.0, 0.01),
                 ("GO:3", 1.0, 0.5),  ("GO:4", 1.0, 0.5)])
        b = _df([("GO:1", 1.0, 0.01), ("GO:2", 1.0, 0.5),
                 ("GO:3", 1.0, 0.01), ("GO:4", 1.0, 0.5)])
        merged = _merge_pair(a, b)
        counts = _sankey_counts(merged, fdr_threshold=0.25)
        assert counts == {
            (True, True): 1,    # GO:1 sig in both
            (True, False): 1,   # GO:2 lost
            (False, True): 1,   # GO:3 gained
            (False, False): 1,  # GO:4 stable not-sig
        }

    def test_empty_merge_all_zero(self):
        merged = pd.DataFrame(columns=["Term", "NES_a", "FDR_a", "NES_b", "FDR_b"])
        counts = _sankey_counts(merged, fdr_threshold=0.25)
        assert counts == {
            (True, True): 0, (True, False): 0,
            (False, True): 0, (False, False): 0,
        }

    def test_counts_sum_to_common_terms(self):
        a = _df([("GO:%d" % i, 1.0, 0.01 if i % 2 == 0 else 0.5) for i in range(6)])
        b = _df([("GO:%d" % i, 1.0, 0.01 if i % 3 == 0 else 0.5) for i in range(6)])
        merged = _merge_pair(a, b)
        counts = _sankey_counts(merged, fdr_threshold=0.25)
        assert sum(counts.values()) == len(merged)


class TestSpearmanMatrix:

    def test_perfectly_correlated_runs(self):
        dfs = {
            "A": _df([("GO:1", 1.0, 0.1), ("GO:2", 2.0, 0.1), ("GO:3", 3.0, 0.1)]),
            "B": _df([("GO:1", 1.5, 0.1), ("GO:2", 2.5, 0.1), ("GO:3", 3.5, 0.1)]),
        }
        mat = _spearman_matrix(dfs)
        assert mat.loc["A", "B"] == pytest.approx(1.0)
        assert mat.loc["A", "A"] == pytest.approx(1.0)

    def test_anticorrelated_runs(self):
        dfs = {
            "A": _df([("GO:1", 1.0, 0.1), ("GO:2", 2.0, 0.1), ("GO:3", 3.0, 0.1)]),
            "B": _df([("GO:1", 3.0, 0.1), ("GO:2", 2.0, 0.1), ("GO:3", 1.0, 0.1)]),
        }
        mat = _spearman_matrix(dfs)
        assert mat.loc["A", "B"] == pytest.approx(-1.0)

    def test_too_few_common_terms_is_nan(self):
        dfs = {
            "A": _df([("GO:1", 1.0, 0.1)]),
            "B": _df([("GO:2", 1.0, 0.1)]),
        }
        mat = _spearman_matrix(dfs)
        assert pd.isna(mat.loc["A", "B"])

    def test_matrix_is_symmetric(self):
        dfs = {
            "A": _df([("GO:1", 1.0, 0.1), ("GO:2", 0.2, 0.1), ("GO:3", -1.0, 0.1)]),
            "B": _df([("GO:1", 0.5, 0.1), ("GO:2", 1.5, 0.1), ("GO:3", -0.5, 0.1)]),
            "C": _df([("GO:1", -1.0, 0.1), ("GO:2", -0.2, 0.1), ("GO:3", 1.0, 0.1)]),
        }
        mat = _spearman_matrix(dfs)
        assert mat.loc["A", "C"] == pytest.approx(mat.loc["C", "A"])


# ---------------------------------------------------------------------------
# Run parameters diff table
# ---------------------------------------------------------------------------

class TestDiffRows:

    def test_identical_values_not_flagged(self):
        specs = [("min_species", "Minimum species")]
        records = [{"filtering": {"min_species": 4}},
                   {"filtering": {"min_species": 4}}]
        html_out = _diff_rows(
            specs, lambda r: (r or {}).get("filtering") or {}, records, ["A", "B"],
        )
        assert "diff-row" not in html_out

    def test_differing_values_flagged(self):
        specs = [("min_species", "Minimum species")]
        records = [{"filtering": {"min_species": 4}},
                   {"filtering": {"min_species": 8}}]
        html_out = _diff_rows(
            specs, lambda r: (r or {}).get("filtering") or {}, records, ["A", "B"],
        )
        assert 'class="diff-row"' in html_out

    def test_none_record_shows_not_set(self):
        specs = [("min_species", "Minimum species")]
        records = [{"filtering": {"min_species": 4}}, None]
        html_out = _diff_rows(
            specs, lambda r: (r or {}).get("filtering") or {}, records, ["A", "B"],
        )
        assert "not set" in html_out
        # None vs a real value must count as differing.
        assert 'class="diff-row"' in html_out


class TestBuildRunParamsDiffHtml:

    def test_two_list_sections_included_when_any_run_is_two_list(self):
        records = [
            _minimal_run_parameters(two_list=True),
            _minimal_run_parameters(two_list=False),
        ]
        html_out = _build_run_params_diff_html(records, ["A", "B"])
        assert "List 1" in html_out
        assert "List 2" in html_out

    def test_two_list_sections_omitted_when_all_single_list(self):
        records = [
            _minimal_run_parameters(two_list=False),
            _minimal_run_parameters(two_list=False),
        ]
        html_out = _build_run_params_diff_html(records, ["A", "B"])
        assert "Two-list mode" not in html_out

    def test_handles_missing_record(self):
        records = [_minimal_run_parameters(), None]
        html_out = _build_run_params_diff_html(records, ["A", "B"])
        assert "not set" in html_out


# ---------------------------------------------------------------------------
# SVG rendering smoke tests
# ---------------------------------------------------------------------------

class TestRenderCorrelationScatterSvg:

    def test_produces_valid_svg_markup(self):
        merged = _merge_pair(
            _df([("GO:1", 1.0, 0.1), ("GO:2", 0.5, 0.3), ("GO:3", -0.5, 0.4)]),
            _df([("GO:1", 0.9, 0.2), ("GO:2", 0.4, 0.05), ("GO:3", -0.4, 0.3)]),
        )
        svg = _render_correlation_scatter_svg("RunA", "RunB", merged, 0.25, 0.87)
        assert svg.startswith("<svg")
        assert svg.count("<circle") == 3
        assert "RunA" in svg and "RunB" in svg

    def test_via_pairwise_block_with_insufficient_overlap(self):
        a = _df([("GO:1", 1.0, 0.1)])
        b = _df([("GO:2", 1.0, 0.1)])
        result = _pairwise_correlation_block("Z-score", "RunA", "RunB", a, b, 0.25)
        assert "Fewer than 2 gene sets" in result
        assert "<svg" not in result


class TestRenderSankeySvg:

    def test_produces_valid_svg_with_all_four_flows(self):
        counts = {(True, True): 3, (True, False): 2,
                  (False, True): 1, (False, False): 5}
        svg = _render_sankey_svg("RunA", "RunB", "Z-score", counts)
        assert svg.startswith("<svg")
        assert svg.count("<path") == 4
        assert svg.count("<rect") == 4

    def test_zero_flow_omitted(self):
        counts = {(True, True): 0, (True, False): 2,
                  (False, True): 0, (False, False): 5}
        svg = _render_sankey_svg("RunA", "RunB", "Z-score", counts)
        # Only the two non-zero flows should produce a ribbon path.
        assert svg.count("<path") == 2

    def test_all_zero_returns_message_not_svg(self):
        counts = {(True, True): 0, (True, False): 0,
                  (False, True): 0, (False, False): 0}
        result = _render_sankey_svg("RunA", "RunB", "Z-score", counts)
        assert "<svg" not in result
        assert "No common gene sets" in result
