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

    def test_rank_highest_identity_gets_highest_rank_value(self):
        """Rank=1.0 (N/N) should go to the most conserved (highest identity)
        group, matching normalise_scores()'s own documented convention
        ("1 = most conserved").

        A prior version of this test asserted the opposite direction
        (idxmin() instead of idxmax()), matching a since-fixed bug in the
        'rank' branch of normalise_scores() itself (it used
        ascending=False instead of ascending=True). The test was wrong in
        the same direction as the code, so it passed despite the bug —
        see test_all_three_metrics_agree_in_sign below for the test that
        actually would have caught it via compute_differential().
        """
        s = self._series([0.9, 0.5, 0.7])
        result = normalise_scores(s, "rank")
        assert result.idxmax() == "OG0000"          # 0.9 is highest -> rank 3/3
        assert result["OG0000"] == pytest.approx(1.0)
        assert result["OG0001"] == pytest.approx(1 / 3)   # 0.5 is lowest -> rank 1/3

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

    def test_all_three_metrics_agree_in_sign(self, tmp_path):
        """Regression test for a metric-specific sign-inversion bug.

        For a group that is clearly more conserved *relative to its own
        list's distribution* in list 1 than in list 2, all three ranking
        metrics (identity, zscore, rank) must agree the differential score
        is positive. A silent ascending/descending mistake in any one
        metric's own normalisation (as previously happened for 'rank' in
        normalise_scores(), which used ascending=False instead of
        ascending=True) flips that metric's sign relative to the other two
        without raising any error — GSEA's own statistics are symmetric
        under negation, so the bug never surfaced as a crash, only as
        rank-metric NES signs disagreeing with identity/zscore for
        genuinely-enriched gene sets, and the Cross-metric consensus tab
        misreporting real 2/3 agreement as "Mixed direction".

        OG1 is deliberately an extreme high performer in list1's own
        distribution (its max) while sitting near the *average* of list2's
        distribution — not just uniformly shifted or rescaled relative to
        the rest of its own list, since z-score and (up to ties) rank are
        both invariant to a shared additive/multiplicative transform
        applied to an entire list, and would trivially agree on sign for a
        less carefully constructed case regardless of any sign bug.
        """
        group2mean1 = [
            ("OG1", 0.95), ("OG2", 0.50), ("OG3", 0.55),
            ("OG4", 0.60), ("OG5", 0.45),
        ]
        group2mean2 = [
            ("OG1", 0.50), ("OG2", 0.50), ("OG3", 0.55),
            ("OG4", 0.60), ("OG5", 0.45),
        ]
        group2anchor = [(g, "gene_%s" % g) for g, _ in group2mean1]

        t1 = _make_tables_dir(tmp_path, "t1_sign",
            group2mean=group2mean1, group2anchor=group2anchor)
        t2 = _make_tables_dir(tmp_path, "t2_sign",
            group2mean=group2mean2, group2anchor=group2anchor)

        for metric in ("identity", "zscore", "rank"):
            diff_dir = str(tmp_path / ("diff_sign_%s" % metric))
            compute_differential(t1, t2, diff_dir, metric=metric)
            scores = pd.read_csv(
                os.path.join(diff_dir, "differential_scores.tsv"), sep="\t"
            ).set_index("group_id")["score"]
            assert scores["OG1"] > 0, (
                "metric=%s gave OG1 a non-positive differential score "
                "(%.4f), even though OG1 is the top performer within "
                "list1's own distribution but only average within "
                "list2's — all three metrics must agree this is positive."
                % (metric, scores["OG1"])
            )
