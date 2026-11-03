"""Tests for collapsing identical-result metrics in the consensus tab."""

import pandas as pd

from enhydra.report import (
    group_identical_results,
    build_significance_agreement_table_html,
)


def _df(rows):
    return pd.DataFrame(rows, columns=["Term", "NES", "FDR q-val"])


IDENT = _df([("GO:1", 2.0, 0.01), ("GO:2", -1.5, 0.30), ("GO:3", 1.1, 0.60)])
ZSCORE = _df([("GO:1", 1.8, 0.40), ("GO:2", -1.4, 0.35), ("GO:3", 1.0, 0.70)])


class TestGroupIdenticalResults:

    def test_identical_columns_grouped(self):
        g = group_identical_results(
            {"identity": IDENT, "zscore": ZSCORE, "rank": IDENT.copy()}
        )
        assert g == [["identity", "rank"], ["zscore"]]

    def test_all_distinct(self):
        other = _df([("GO:1", 0.5, 0.9), ("GO:2", 0.1, 0.9), ("GO:3", 0.2, 0.9)])
        g = group_identical_results(
            {"a": IDENT, "b": ZSCORE, "c": other}
        )
        assert g == [["a"], ["b"], ["c"]]

    def test_tiny_float_noise_still_grouped(self):
        noisy = IDENT.copy()
        noisy["NES"] = noisy["NES"] * (1 + 1e-9)
        assert group_identical_results({"a": IDENT, "b": noisy}) == [["a", "b"]]

    def test_different_term_sets_not_grouped(self):
        fewer = IDENT.iloc[:2].copy()
        assert group_identical_results({"a": IDENT, "b": fewer}) == [["a"], ["b"]]

    def test_every_column_in_exactly_one_group(self):
        cols = {"identity": IDENT, "zscore": ZSCORE, "rank": IDENT.copy()}
        flat = [c for g in group_identical_results(cols) for c in g]
        assert sorted(flat) == sorted(cols)


class TestConsensusCountsDistinctAnalyses:

    def test_without_dedup_duplicate_pair_alone_reaches_majority(self):
        """Documents the original problem: identity+rank alone give 2/3."""
        html_, _ = build_significance_agreement_table_html(
            {"identity": IDENT, "zscore": ZSCORE, "rank": IDENT.copy()},
            {"identity": "Identity", "zscore": "Z-score", "rank": "Rank"},
            {}, fdr_threshold=0.25,
        )
        assert "GO:1" in html_  # listed although z-score disagrees

    def test_with_dedup_requires_both_distinct_analyses(self):
        cols = {"identity": IDENT, "zscore": ZSCORE, "rank": IDENT.copy()}
        groups = group_identical_results(cols)
        dedup = {g[0]: cols[g[0]] for g in groups}
        html_, num = build_significance_agreement_table_html(
            dedup, {"identity": "Identity = Rank", "zscore": "Z-score"},
            {}, fdr_threshold=0.25,
        )
        # GO:1 is significant only in the unweighted analysis -> 1 of 2, below
        # the majority of 2, so nothing qualifies.
        assert "No gene sets reached majority significance" in html_
        assert num is None

    def test_with_dedup_gene_set_significant_in_both_is_listed(self):
        z = ZSCORE.copy()
        z.loc[0, "FDR q-val"] = 0.05
        cols = {"identity": IDENT, "zscore": z, "rank": IDENT.copy()}
        groups = group_identical_results(cols)
        dedup = {g[0]: cols[g[0]] for g in groups}
        html_, _ = build_significance_agreement_table_html(
            dedup, {"identity": "Identity = Rank", "zscore": "Z-score"},
            {}, fdr_threshold=0.25,
        )
        assert "GO:1" in html_
        assert "2 / 2" in html_
        assert 'class="sig-row"' in html_
