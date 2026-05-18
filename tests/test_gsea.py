"""Tests for enhydra.gsea.

build_gmt_from_gprofiler: GProfiler is mocked at the import site so no
network calls are made.

_check_overlap: tested via a real GMT file written to tmp_path.
"""
from __future__ import annotations

from unittest.mock import MagicMock, patch

import pandas as pd
import pytest

from enhydra.exceptions import EnhydraIOError
from enhydra.gsea import build_gmt_from_gprofiler, _check_overlap


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_profile_df(rows: list[dict]) -> pd.DataFrame:
    """Build a minimal g:Profiler-style DataFrame."""
    return pd.DataFrame(rows, columns=["native", "name", "intersections"])


def _mock_gprofiler(df: pd.DataFrame):
    """Return a mock GProfiler class whose instance.profile() returns df."""
    mock_inst = MagicMock()
    mock_inst.profile.return_value = df
    mock_cls = MagicMock(return_value=mock_inst)
    return mock_cls, mock_inst


def _write_gmt(path: str, lines: list[str]) -> None:
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


def _ranked_df(gene_ids: list[str]) -> pd.DataFrame:
    """Minimal ranked DataFrame matching _check_overlap's expectations."""
    return pd.DataFrame({
        "gene_id":  gene_ids,
        "identity": [float(i) for i in range(len(gene_ids), 0, -1)],
    })


# ---------------------------------------------------------------------------
# build_gmt_from_gprofiler
# ---------------------------------------------------------------------------

class TestBuildGmtFromGprofiler:

    def test_gmt_written_in_correct_format(self, tmp_path):
        """Each output line is: term_id <TAB> term_name <TAB> gene1 <TAB> gene2 ..."""
        df = _make_profile_df([
            {"native": "GO:0001234", "name": "some process",
             "intersections": ["geneA", "geneB"]},
            {"native": "GO:0005678", "name": "other process",
             "intersections": ["geneC"]},
        ])
        gmt_path = str(tmp_path / "out.gmt")
        mock_cls, _ = _mock_gprofiler(df)

        with patch("enhydra.gsea.GProfiler", mock_cls):
            build_gmt_from_gprofiler(["geneA", "geneB", "geneC"], "hsapiens", gmt_path)

        lines = open(gmt_path).read().splitlines()
        assert len(lines) == 2
        assert lines[0] == "GO:0001234\tsome process\tgeneA\tgeneB"
        assert lines[1] == "GO:0005678\tother process\tgeneC"

    def test_returns_gmt_path(self, tmp_path):
        df  = _make_profile_df([
            {"native": "GO:0001", "name": "foo", "intersections": ["g1"]},
        ])
        gmt_path = str(tmp_path / "out.gmt")
        mock_cls, _ = _mock_gprofiler(df)

        with patch("enhydra.gsea.GProfiler", mock_cls):
            result = build_gmt_from_gprofiler(["g1"], "hsapiens", gmt_path)

        assert result == gmt_path

    def test_empty_results_raise_enhydra_io_error(self, tmp_path):
        gmt_path = str(tmp_path / "out.gmt")
        mock_cls, _ = _mock_gprofiler(pd.DataFrame())

        with patch("enhydra.gsea.GProfiler", mock_cls):
            with pytest.raises(EnhydraIOError, match="g:Profiler returned no annotations"):
                build_gmt_from_gprofiler(["geneA"], "invalid_organism", gmt_path)

    def test_rows_with_empty_intersections_skipped(self, tmp_path):
        """Rows where intersections is empty or None must not appear in the GMT."""
        df = _make_profile_df([
            {"native": "GO:0001", "name": "present", "intersections": ["g1"]},
            {"native": "GO:0002", "name": "empty list", "intersections": []},
            {"native": "GO:0003", "name": "none value", "intersections": None},
        ])
        gmt_path = str(tmp_path / "out.gmt")
        mock_cls, _ = _mock_gprofiler(df)

        with patch("enhydra.gsea.GProfiler", mock_cls):
            build_gmt_from_gprofiler(["g1"], "hsapiens", gmt_path)

        lines = [l for l in open(gmt_path).read().splitlines() if l]
        assert len(lines) == 1
        assert lines[0].startswith("GO:0001")

    def test_default_sources_when_none_passed(self, tmp_path):
        """sources=None must result in the default list being sent to g:Profiler."""
        df = _make_profile_df([
            {"native": "GO:0001", "name": "foo", "intersections": ["g1"]},
        ])
        gmt_path = str(tmp_path / "out.gmt")
        mock_cls, mock_inst = _mock_gprofiler(df)

        with patch("enhydra.gsea.GProfiler", mock_cls):
            build_gmt_from_gprofiler(["g1"], "hsapiens", gmt_path, sources=None)

        call_kwargs = mock_inst.profile.call_args[1]
        assert call_kwargs["sources"] == ["GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC"]

    def test_custom_sources_passed_through(self, tmp_path):
        df = _make_profile_df([
            {"native": "GO:0001", "name": "foo", "intersections": ["g1"]},
        ])
        gmt_path = str(tmp_path / "out.gmt")
        mock_cls, mock_inst = _mock_gprofiler(df)

        with patch("enhydra.gsea.GProfiler", mock_cls):
            build_gmt_from_gprofiler(["g1"], "hsapiens", gmt_path, sources=["GO:BP"])

        call_kwargs = mock_inst.profile.call_args[1]
        assert call_kwargs["sources"] == ["GO:BP"]

    def test_gene_ids_forwarded_as_query(self, tmp_path):
        """The gene_ids argument must be forwarded as the 'query' kwarg."""
        gene_ids = ["ENSG00000141510", "ENSG00000012048"]
        df = _make_profile_df([
            {"native": "GO:0001", "name": "foo", "intersections": gene_ids},
        ])
        gmt_path = str(tmp_path / "out.gmt")
        mock_cls, mock_inst = _mock_gprofiler(df)

        with patch("enhydra.gsea.GProfiler", mock_cls):
            build_gmt_from_gprofiler(gene_ids, "hsapiens", gmt_path)

        call_kwargs = mock_inst.profile.call_args[1]
        assert call_kwargs["query"] == gene_ids

    def test_organism_forwarded(self, tmp_path):
        df = _make_profile_df([
            {"native": "GO:0001", "name": "foo", "intersections": ["g1"]},
        ])
        gmt_path = str(tmp_path / "out.gmt")
        mock_cls, mock_inst = _mock_gprofiler(df)

        with patch("enhydra.gsea.GProfiler", mock_cls):
            build_gmt_from_gprofiler(["g1"], "mmusculus", gmt_path)

        call_kwargs = mock_inst.profile.call_args[1]
        assert call_kwargs["organism"] == "mmusculus"

    def test_gprofiler_instantiated_with_correct_kwargs(self, tmp_path):
        """GProfiler must be constructed with return_dataframe=True."""
        df = _make_profile_df([
            {"native": "GO:0001", "name": "foo", "intersections": ["g1"]},
        ])
        gmt_path = str(tmp_path / "out.gmt")
        mock_cls, _ = _mock_gprofiler(df)

        with patch("enhydra.gsea.GProfiler", mock_cls):
            build_gmt_from_gprofiler(["g1"], "hsapiens", gmt_path)

        _, ctor_kwargs = mock_cls.call_args
        assert ctor_kwargs.get("return_dataframe") is True


# ---------------------------------------------------------------------------
# _check_overlap
# ---------------------------------------------------------------------------

class TestCheckOverlap:

    def test_sufficient_overlap_does_not_raise(self, tmp_path):
        gmt = str(tmp_path / "test.gmt")
        _write_gmt(gmt, [
            "GO:0001\tdescription\tgeneA\tgeneB\tgeneC",
        ])
        ranked = _ranked_df(["geneA", "geneB", "geneZ"])
        # Two genes in common (geneA, geneB) — should not raise.
        _check_overlap(ranked, gmt)

    def test_zero_overlap_raises(self, tmp_path):
        gmt = str(tmp_path / "test.gmt")
        _write_gmt(gmt, [
            "GO:0001\tdescription\tgeneX\tgeneY",
        ])
        ranked = _ranked_df(["geneA", "geneB"])
        with pytest.raises(EnhydraIOError, match="Too few genes overlap"):
            _check_overlap(ranked, gmt)

    def test_single_overlap_raises(self, tmp_path):
        """GSEApy requires at least 2 overlapping genes; 1 must raise."""
        gmt = str(tmp_path / "test.gmt")
        _write_gmt(gmt, [
            "GO:0001\tdescription\tgeneA\tgeneX",
        ])
        ranked = _ranked_df(["geneA", "geneB"])
        with pytest.raises(EnhydraIOError, match="Too few genes overlap"):
            _check_overlap(ranked, gmt)

    def test_error_message_includes_example_ids(self, tmp_path):
        """The error must show example IDs from both sides to aid debugging."""
        gmt = str(tmp_path / "test.gmt")
        _write_gmt(gmt, ["GO:0001\tdescription\tgmtGeneA\tgmtGeneB"])
        ranked = _ranked_df(["rankedGeneX", "rankedGeneY"])
        with pytest.raises(EnhydraIOError) as exc_info:
            _check_overlap(ranked, gmt)
        msg = str(exc_info.value)
        # Both sides must be represented in the diagnostic output.
        assert any(g in msg for g in ["rankedGeneX", "rankedGeneY"])
        assert any(g in msg for g in ["gmtGeneA", "gmtGeneB"])

    def test_exactly_two_overlap_does_not_raise(self, tmp_path):
        """The boundary case: exactly 2 overlapping genes is sufficient."""
        gmt = str(tmp_path / "test.gmt")
        _write_gmt(gmt, ["GO:0001\tdescription\tgeneA\tgeneB\tgeneX"])
        ranked = _ranked_df(["geneA", "geneB", "geneZ"])
        _check_overlap(ranked, gmt)   # must not raise

    def test_multiple_gene_sets_in_gmt(self, tmp_path):
        """Overlap is computed across all gene sets, not just the first line."""
        gmt = str(tmp_path / "test.gmt")
        _write_gmt(gmt, [
            "GO:0001\tfoo\tgeneX\tgeneY",
            "GO:0002\tbar\tgeneA\tgeneB",
        ])
        ranked = _ranked_df(["geneA", "geneB"])
        # Overlap via GO:0002 only — still ≥ 2, must not raise.
        _check_overlap(ranked, gmt)
