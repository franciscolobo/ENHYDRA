"""Tests for enhydra.provenance."""

import os
import json

import pytest

from enhydra.provenance import collect_run_parameters, write_run_parameters


def _minimal_resolved_params(**overrides) -> dict:
    """A representative resolved_params dict, as cli.py would build it
    after _resolve() has been applied to every CLI/config value."""
    base = {
        "min_species":          4,
        "min_sequences":        2,
        "paralogs":             "all",
        "length_filter_sd":     2.0,
        "divergence_filter_sd": None,
        "trim":                 "",
        "aligner":              "mafft",
        "mafft_mode":           "auto",
        "metric":               "zscore",
        "permutations":         1000,
        "min_size":             5,
        "max_size":             500,
        "seed":                 42,
        "fdr_threshold":        0.25,
        "top_n":                20,
        "gene_sets":            "gmt/NC_004431_GO_BP.gmt",
        "organism":             None,
        "sources":              ["GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC"],
    }
    base.update(overrides)
    return base


# ---------------------------------------------------------------------------
# collect_run_parameters
# ---------------------------------------------------------------------------

class TestCollectRunParametersSingleList:

    def test_top_level_fields(self, tmp_path):
        record = collect_run_parameters(
            outdir=str(tmp_path / "out"),
            code_config_path="code_config",
            project_config_path="project_config",
            command_line=["enhydra", "code_config", "project_config"],
            inputdir=str(tmp_path / "input"),
            anchor="hsapiens",
            resolved_params=_minimal_resolved_params(),
            metrics_run=("zscore",),
        )
        assert record["enhydra_version"]
        assert "timestamp" in record
        assert record["command_line"] == ["enhydra", "code_config", "project_config"]
        assert record["outdir"] == os.path.abspath(str(tmp_path / "out"))
        assert record["code_config_path"] == os.path.abspath("code_config")
        assert record["project_config_path"] == os.path.abspath("project_config")

    def test_input_section(self, tmp_path):
        record = collect_run_parameters(
            outdir=str(tmp_path / "out"),
            code_config_path="code_config",
            project_config_path="project_config",
            command_line=["enhydra"],
            inputdir=str(tmp_path / "input"),
            anchor="hsapiens",
            resolved_params=_minimal_resolved_params(),
            metrics_run=("zscore",),
        )
        assert record["input"]["inputdir"] == os.path.abspath(str(tmp_path / "input"))
        assert record["input"]["anchor"] == "hsapiens"
        assert record["input"]["two_list_mode"] is False
        assert record["input"]["orthofinder_dir"] is None
        assert "list1" not in record["input"]
        assert "list2" not in record["input"]

    def test_orthofinder_dir_recorded_when_given(self, tmp_path):
        of_dir = tmp_path / "orthofinder_results"
        record = collect_run_parameters(
            outdir=str(tmp_path / "out"),
            code_config_path="code_config",
            project_config_path="project_config",
            command_line=["enhydra"],
            inputdir=str(tmp_path / "input"),
            anchor="hsapiens",
            resolved_params=_minimal_resolved_params(),
            orthofinder_dir=str(of_dir),
            metrics_run=("zscore",),
        )
        assert record["input"]["orthofinder_dir"] == os.path.abspath(str(of_dir))

    def test_filtering_section(self, tmp_path):
        record = collect_run_parameters(
            outdir=str(tmp_path), code_config_path="c", project_config_path="p",
            command_line=["enhydra"], inputdir=str(tmp_path), anchor="a",
            resolved_params=_minimal_resolved_params(
                min_species=6, length_filter_sd=1.5, divergence_filter_sd=2.0,
                trim="strict",
            ),
            metrics_run=("zscore",),
        )
        f = record["filtering"]
        assert f["min_species"] == 6
        assert f["min_sequences"] == 2
        assert f["paralogs"] == "all"
        assert f["length_filter_sd"] == 1.5
        assert f["divergence_filter_sd"] == 2.0
        assert f["trim"] == "strict"

    def test_empty_trim_recorded_as_none(self, tmp_path):
        """trim='' (disabled) should be recorded as null, not an empty string,
        so downstream consumers can treat 'not set' uniformly."""
        record = collect_run_parameters(
            outdir=str(tmp_path), code_config_path="c", project_config_path="p",
            command_line=["enhydra"], inputdir=str(tmp_path), anchor="a",
            resolved_params=_minimal_resolved_params(trim=""),
            metrics_run=("zscore",),
        )
        assert record["filtering"]["trim"] is None

    def test_alignment_section(self, tmp_path):
        record = collect_run_parameters(
            outdir=str(tmp_path), code_config_path="c", project_config_path="p",
            command_line=["enhydra"], inputdir=str(tmp_path), anchor="a",
            resolved_params=_minimal_resolved_params(aligner="muscle", mafft_mode="linsi"),
            metrics_run=("zscore",),
        )
        assert record["alignment"] == {"aligner": "muscle", "mafft_mode": "linsi"}

    def test_gsea_section(self, tmp_path):
        record = collect_run_parameters(
            outdir=str(tmp_path), code_config_path="c", project_config_path="p",
            command_line=["enhydra"], inputdir=str(tmp_path), anchor="a",
            resolved_params=_minimal_resolved_params(
                permutations=500, min_size=10, max_size=200,
                seed=7, fdr_threshold=0.1, top_n=30,
            ),
            metrics_run=("zscore",),
        )
        g = record["gsea"]
        assert g == {
            "permutations": 500, "min_size": 10, "max_size": 200,
            "seed": 7, "fdr_threshold": 0.1, "top_n": 30,
        }

    def test_gene_sets_section_with_local_gmt(self, tmp_path):
        record = collect_run_parameters(
            outdir=str(tmp_path), code_config_path="c", project_config_path="p",
            command_line=["enhydra"], inputdir=str(tmp_path), anchor="a",
            resolved_params=_minimal_resolved_params(
                gene_sets="gmt/foo.gmt", organism=None,
            ),
            metrics_run=("zscore",),
        )
        assert record["gene_sets"]["gene_sets_path"] == "gmt/foo.gmt"
        assert record["gene_sets"]["organism"] is None

    def test_gene_sets_section_with_organism(self, tmp_path):
        record = collect_run_parameters(
            outdir=str(tmp_path), code_config_path="c", project_config_path="p",
            command_line=["enhydra"], inputdir=str(tmp_path), anchor="a",
            resolved_params=_minimal_resolved_params(
                gene_sets="", organism="hsapiens",
            ),
            metrics_run=("zscore",),
        )
        assert record["gene_sets"]["gene_sets_path"] is None
        assert record["gene_sets"]["organism"] == "hsapiens"

    def test_ranking_section_single_metric(self, tmp_path):
        record = collect_run_parameters(
            outdir=str(tmp_path), code_config_path="c", project_config_path="p",
            command_line=["enhydra"], inputdir=str(tmp_path), anchor="a",
            resolved_params=_minimal_resolved_params(),
            all_metrics=False, metrics_run=("zscore",),
        )
        assert record["ranking"] == {
            "all_metrics": False, "metrics_run": ["zscore"],
        }

    def test_ranking_section_all_metrics(self, tmp_path):
        record = collect_run_parameters(
            outdir=str(tmp_path), code_config_path="c", project_config_path="p",
            command_line=["enhydra"], inputdir=str(tmp_path), anchor="a",
            resolved_params=_minimal_resolved_params(),
            all_metrics=True, metrics_run=("identity", "zscore", "rank"),
        )
        assert record["ranking"] == {
            "all_metrics": True,
            "metrics_run": ["identity", "zscore", "rank"],
        }

    def test_missing_resolved_param_key_is_none(self, tmp_path):
        """A resolved_params dict missing an expected key must not raise —
        the corresponding field is simply recorded as null."""
        incomplete = _minimal_resolved_params()
        del incomplete["seed"]
        record = collect_run_parameters(
            outdir=str(tmp_path), code_config_path="c", project_config_path="p",
            command_line=["enhydra"], inputdir=str(tmp_path), anchor="a",
            resolved_params=incomplete,
            metrics_run=("zscore",),
        )
        assert record["gsea"]["seed"] is None


class TestCollectRunParametersTwoList:

    def test_two_list_fields_present(self, tmp_path):
        record = collect_run_parameters(
            outdir=str(tmp_path), code_config_path="c", project_config_path="p",
            command_line=["enhydra", "--list1", "path1.txt", "--list2", "path2.txt"],
            inputdir=str(tmp_path), anchor="anchor_sp",
            resolved_params=_minimal_resolved_params(),
            two_list_mode=True,
            list1_path="path1.txt", list2_path="path2.txt",
            list1_name="Pathogenic", list2_name="Non-pathogenic",
            species1=["sp1", "sp2", "sp3", "anchor_sp"],
            species2=["sp4", "sp5"],
            anchor_injected=True,
            metrics_run=("identity", "zscore", "rank"),
        )
        assert record["input"]["two_list_mode"] is True
        assert record["input"]["list1"] == {
            "path": "path1.txt", "name": "Pathogenic",
            "n_species": 4, "anchor_injected": True,
        }
        assert record["input"]["list2"] == {
            "path": "path2.txt", "name": "Non-pathogenic", "n_species": 2,
        }

    def test_anchor_not_injected(self, tmp_path):
        record = collect_run_parameters(
            outdir=str(tmp_path), code_config_path="c", project_config_path="p",
            command_line=["enhydra"], inputdir=str(tmp_path), anchor="anchor_sp",
            resolved_params=_minimal_resolved_params(),
            two_list_mode=True,
            list1_path="path1.txt", list2_path="path2.txt",
            list1_name="List 1", list2_name="List 2",
            species1=["anchor_sp", "sp2"], species2=["sp3"],
            anchor_injected=False,
            metrics_run=("zscore",),
        )
        assert record["input"]["list1"]["anchor_injected"] is False

    def test_species_lists_recorded_as_counts_only(self, tmp_path):
        """Full species IDs are not duplicated into this record — only
        counts, since the original list files are already recorded by
        path and are the source of truth."""
        record = collect_run_parameters(
            outdir=str(tmp_path), code_config_path="c", project_config_path="p",
            command_line=["enhydra"], inputdir=str(tmp_path), anchor="a",
            resolved_params=_minimal_resolved_params(),
            two_list_mode=True,
            list1_path="p1.txt", list2_path="p2.txt",
            list1_name="L1", list2_name="L2",
            species1=["sp1", "sp2"], species2=["sp3"],
            anchor_injected=False,
            metrics_run=("zscore",),
        )
        serialised = json.dumps(record)
        assert "sp1" not in serialised
        assert "sp2" not in serialised
        assert "sp3" not in serialised

    def test_single_list_mode_omits_two_list_keys(self, tmp_path):
        record = collect_run_parameters(
            outdir=str(tmp_path), code_config_path="c", project_config_path="p",
            command_line=["enhydra"], inputdir=str(tmp_path), anchor="a",
            resolved_params=_minimal_resolved_params(),
            two_list_mode=False,
            metrics_run=("zscore",),
        )
        assert "list1" not in record["input"]
        assert "list2" not in record["input"]


class TestCollectRunParametersJsonSerialisable:

    def test_full_record_is_json_serialisable(self, tmp_path):
        record = collect_run_parameters(
            outdir=str(tmp_path), code_config_path="c", project_config_path="p",
            command_line=["enhydra", "c", "p", "--all-metrics"],
            inputdir=str(tmp_path), anchor="hsapiens",
            resolved_params=_minimal_resolved_params(divergence_filter_sd=2.0),
            two_list_mode=True,
            list1_path="p1.txt", list2_path="p2.txt",
            list1_name="L1", list2_name="L2",
            species1=["a", "b"], species2=["c"],
            anchor_injected=True,
            all_metrics=True, metrics_run=("identity", "zscore", "rank"),
        )
        # Must not raise.
        json.dumps(record)


# ---------------------------------------------------------------------------
# write_run_parameters
# ---------------------------------------------------------------------------

class TestWriteRunParameters:

    def test_default_output_path(self, tmp_path):
        record = {"enhydra_version": "0.1.0"}
        outdir = str(tmp_path / "out")
        os.makedirs(outdir)
        path = write_run_parameters(record, outdir)
        assert path == os.path.join(outdir, "run_parameters.json")
        assert os.path.isfile(path)

    def test_custom_output_path(self, tmp_path):
        record = {"enhydra_version": "0.1.0"}
        outdir = str(tmp_path / "out")
        os.makedirs(outdir)
        custom = str(tmp_path / "custom_params.json")
        path = write_run_parameters(record, outdir, output_path=custom)
        assert path == custom
        assert os.path.isfile(custom)
        assert not os.path.isfile(os.path.join(outdir, "run_parameters.json"))

    def test_written_content_matches_input(self, tmp_path):
        record = {
            "enhydra_version": "0.1.0",
            "filtering": {"min_species": 4},
        }
        outdir = str(tmp_path)
        path = write_run_parameters(record, outdir)
        with open(path) as fh:
            loaded = json.load(fh)
        assert loaded == record

    def test_overwrites_existing_file(self, tmp_path):
        outdir = str(tmp_path)
        write_run_parameters({"enhydra_version": "0.1.0"}, outdir)
        write_run_parameters({"enhydra_version": "0.2.0"}, outdir)
        with open(os.path.join(outdir, "run_parameters.json")) as fh:
            loaded = json.load(fh)
        assert loaded == {"enhydra_version": "0.2.0"}

    def test_returns_path_written(self, tmp_path):
        outdir = str(tmp_path)
        path = write_run_parameters({"enhydra_version": "0.1.0"}, outdir)
        assert os.path.isfile(path)


# ---------------------------------------------------------------------------
# enhydra_version lookup
# ---------------------------------------------------------------------------

class TestEnhydraVersion:

    def test_version_is_a_non_empty_string(self, tmp_path):
        record = collect_run_parameters(
            outdir=str(tmp_path), code_config_path="c", project_config_path="p",
            command_line=["enhydra"], inputdir=str(tmp_path), anchor="a",
            resolved_params=_minimal_resolved_params(),
            metrics_run=("zscore",),
        )
        assert isinstance(record["enhydra_version"], str)
        assert len(record["enhydra_version"]) > 0
