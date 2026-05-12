"""Tests for enhydra.differential."""

import os
import pytest
import pandas as pd

from enhydra.differential import normalise_scores, compute_differential


# ---------------------------------------------------------------------------
# normalise_scores
# ---------------------------------------------------------------------------

class TestNormaliseScores:

    def _series(self, values):
        return pd.Series(values, index=["OG%04d" % i for i in range(len(values))])

    def test_identity_passthrough(self):
        s = self._series([0.8, 0.6, 0.9, 0.7])
        result = normalise_scores(s, "identity")
        pd.testing.assert_series_equal(result, s)

    def test_zscore_mean_zero_std_one(self):
        s = self._series([0.8, 0.6, 0.9, 0.7])
        result = normalise_scores(s, "zscore")
        assert abs(result.mean()) < 1e-10
        assert abs(result.std() - 1.0) < 1e-10

    def test_rank_range(self):
        s = self._series([0.8, 0.6, 0.9, 0.7])
        result = normalise_scores(s, "rank")
        assert result.min() > 0
        assert result.max() <= 1.0

    def test_rank_highest_identity_lowest_rank_value(self):
        """Rank=1/N should go to the most conserved (highest identity) group."""
        s = self._series([0.9, 0.5, 0.7])
        result = normalise_scores(s, "rank")
        assert result.idxmin() == "OG0000"   # 0.9 gets rank 1 → 1/3

    def test_unknown_metric_raises(self):
        s = self._series([0.8, 0.6])
        with pytest.raises(ValueError, match="Unknown metric"):
            normalise_scores(s, "banana")


# ---------------------------------------------------------------------------
# compute_differential
# ---------------------------------------------------------------------------

def _write_tsv(path: str, rows: list[tuple]):
    with open(path, "w") as fh:
        for row in rows:
            fh.write("\t".join(str(v) for v in row) + "\n")


def _make_tables_dir(tmp_path, name: str,
                     group2mean: list[tuple],
                     group2anchor: list[tuple]) -> str:
    d = tmp_path / name
    d.mkdir()
    _write_tsv(str(d / "group2mean.tsv"),   group2mean)
    _write_tsv(str(d / "group2anchor.tsv"), group2anchor)
    return str(d)


class TestComputeDifferential:

    def test_identity_difference(self, tmp_path):
        t1 = _make_tables_dir(tmp_path, "t1",
            group2mean   = [("OG0001", 0.8), ("OG0002", 0.6)],
            group2anchor = [("OG0001", "geneA"), ("OG0002", "geneB")],
        )
        t2 = _make_tables_dir(tmp_path, "t2",
            group2mean   = [("OG0001", 0.5), ("OG0002", 0.7)],
            group2anchor = [("OG0001", "geneA"), ("OG0002", "geneB")],
        )
        diff_dir = str(tmp_path / "diff")
        compute_differential(t1, t2, diff_dir, metric="identity")

        a2m = pd.read_csv(
            os.path.join(diff_dir, "anchor2mean.tsv"),
            sep="\t", header=None, names=["gene_id", "score"],
        ).set_index("gene_id")["score"]

        assert abs(a2m["geneA"] - 0.3) < 1e-9   # 0.8 - 0.5
        assert abs(a2m["geneB"] - (-0.1)) < 1e-9  # 0.6 - 0.7

    def test_only_common_groups_used(self, tmp_path):
        """Groups present in only one list are excluded from the differential."""
        t1 = _make_tables_dir(tmp_path, "t1",
            group2mean   = [("OG0001", 0.8), ("OG0003", 0.9)],
            group2anchor = [("OG0001", "geneA"), ("OG0003", "geneC")],
        )
        t2 = _make_tables_dir(tmp_path, "t2",
            group2mean   = [("OG0001", 0.5), ("OG0002", 0.7)],
            group2anchor = [("OG0001", "geneA"), ("OG0002", "geneB")],
        )
        diff_dir = str(tmp_path / "diff")
        compute_differential(t1, t2, diff_dir, metric="identity")

        a2m = pd.read_csv(
            os.path.join(diff_dir, "anchor2mean.tsv"),
            sep="\t", header=None, names=["gene_id", "score"],
        )
        assert set(a2m["gene_id"]) == {"geneA"}

    def test_no_common_groups_raises(self, tmp_path):
        t1 = _make_tables_dir(tmp_path, "t1",
            group2mean   = [("OG0001", 0.8)],
            group2anchor = [("OG0001", "geneA")],
        )
        t2 = _make_tables_dir(tmp_path, "t2",
            group2mean   = [("OG0002", 0.5)],
            group2anchor = [("OG0002", "geneB")],
        )
        diff_dir = str(tmp_path / "diff")
        with pytest.raises(ValueError, match="No common orthogroups"):
            compute_differential(t1, t2, diff_dir, metric="identity")

    def test_zscore_metric(self, tmp_path):
        """zscore mode: each list is normalised independently before differencing."""
        t1 = _make_tables_dir(tmp_path, "t1",
            group2mean   = [("OG0001", 0.9), ("OG0002", 0.5),
                            ("OG0003", 0.7), ("OG0004", 0.6)],
            group2anchor = [("OG0001", "gA"), ("OG0002", "gB"),
                            ("OG0003", "gC"), ("OG0004", "gD")],
        )
        t2 = _make_tables_dir(tmp_path, "t2",
            group2mean   = [("OG0001", 0.4), ("OG0002", 0.8),
                            ("OG0003", 0.6), ("OG0004", 0.5)],
            group2anchor = [("OG0001", "gA"), ("OG0002", "gB"),
                            ("OG0003", "gC"), ("OG0004", "gD")],
        )
        diff_dir = str(tmp_path / "diff")
        compute_differential(t1, t2, diff_dir, metric="zscore")

        scores = pd.read_csv(
            os.path.join(diff_dir, "differential_scores.tsv"), sep="\t"
        )
        # Scores should be signed and centred around 0
        assert scores["score"].abs().max() > 0
        assert scores["score"].dtype == float

    def test_differential_scores_tsv_written(self, tmp_path):
        t1 = _make_tables_dir(tmp_path, "t1",
            group2mean   = [("OG0001", 0.8)],
            group2anchor = [("OG0001", "geneA")],
        )
        t2 = _make_tables_dir(tmp_path, "t2",
            group2mean   = [("OG0001", 0.5)],
            group2anchor = [("OG0001", "geneA")],
        )
        diff_dir = str(tmp_path / "diff")
        compute_differential(t1, t2, diff_dir, metric="identity")
        assert os.path.isfile(os.path.join(diff_dir, "differential_scores.tsv"))
        assert os.path.isfile(os.path.join(diff_dir, "anchor2mean.tsv"))
