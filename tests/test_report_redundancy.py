"""Tests for the GO-term redundancy clustering logic in enhydra.report."""

import pandas as pd
import pytest

from enhydra.report import _overlap_coefficient, _cluster_redundant_terms


# ---------------------------------------------------------------------------
# _overlap_coefficient
# ---------------------------------------------------------------------------

class TestOverlapCoefficient:

    def test_identical_sets(self):
        a = {"g1", "g2", "g3"}
        assert _overlap_coefficient(a, set(a)) == pytest.approx(1.0)

    def test_disjoint_sets(self):
        a = {"g1", "g2"}
        b = {"g3", "g4"}
        assert _overlap_coefficient(a, b) == 0.0

    def test_small_set_fully_contained_in_large_set_is_maximal(self):
        """The key property motivating overlap coefficient over Jaccard:
        a small child term fully nested inside a much larger parent term
        must score 1.0 (maximally redundant), even though the two sets
        are very different sizes."""
        child  = {"g%d" % i for i in range(20)}
        parent = {"g%d" % i for i in range(800)}
        assert _overlap_coefficient(child, parent) == pytest.approx(1.0)
        # Jaccard on the same pair would be 20/800 = 0.025 — overlap
        # coefficient must not degrade to anything close to that.
        jaccard = len(child & parent) / len(child | parent)
        assert jaccard < 0.05
        assert _overlap_coefficient(child, parent) > 10 * jaccard

    def test_partial_overlap(self):
        a = {"g1", "g2", "g3", "g4"}
        b = {"g3", "g4", "g5"}
        # intersection = {g3, g4} = 2; min(|a|, |b|) = min(4, 3) = 3
        assert _overlap_coefficient(a, b) == pytest.approx(2 / 3)

    def test_empty_set_returns_zero(self):
        assert _overlap_coefficient(set(), {"g1"}) == 0.0
        assert _overlap_coefficient({"g1"}, set()) == 0.0
        assert _overlap_coefficient(set(), set()) == 0.0

    def test_symmetric(self):
        a = {"g1", "g2", "g3"}
        b = {"g2", "g3", "g4", "g5"}
        assert _overlap_coefficient(a, b) == pytest.approx(_overlap_coefficient(b, a))


# ---------------------------------------------------------------------------
# _cluster_redundant_terms
# ---------------------------------------------------------------------------

def _df(rows):
    """rows = [(term, fdr), ...]"""
    return pd.DataFrame(rows, columns=["Term", "FDR q-val"])


class TestClusterRedundantTerms:

    def test_two_highly_overlapping_terms_cluster(self):
        df = _df([("GO:A", 0.01), ("GO:B", 0.02)])
        gene_sets = {
            "GO:A": ["g1", "g2", "g3", "g4"],
            "GO:B": ["g1", "g2", "g3"],   # overlap coeff with A = 3/3 = 1.0
        }
        result = _cluster_redundant_terms(df, gene_sets, fdr_threshold=0.25)
        assert result["GO:A"]["role"] == "representative"
        assert result["GO:A"]["satellites"] == ["GO:B"]
        assert result["GO:B"]["role"] == "satellite"
        assert result["GO:B"]["representative"] == "GO:A"

    def test_most_significant_term_becomes_representative(self):
        """Processing order is ascending FDR — the more significant term
        must be the representative regardless of input row order."""
        df = _df([("GO:LESS_SIG", 0.10), ("GO:MORE_SIG", 0.01)])
        gene_sets = {
            "GO:LESS_SIG": ["g1", "g2", "g3"],
            "GO:MORE_SIG": ["g1", "g2", "g3"],
        }
        result = _cluster_redundant_terms(df, gene_sets, fdr_threshold=0.25)
        assert result["GO:MORE_SIG"]["role"] == "representative"
        assert result["GO:LESS_SIG"]["role"] == "satellite"
        assert result["GO:LESS_SIG"]["representative"] == "GO:MORE_SIG"

    def test_below_threshold_terms_stay_independent(self):
        df = _df([("GO:A", 0.01), ("GO:B", 0.02)])
        gene_sets = {
            "GO:A": ["g1", "g2", "g3", "g4", "g5", "g6", "g7", "g8"],
            "GO:B": ["g1", "g9"],   # overlap coeff = 1/2 = 0.5, below default 0.75
        }
        result = _cluster_redundant_terms(df, gene_sets, fdr_threshold=0.25)
        assert result["GO:A"]["role"] == "representative"
        assert result["GO:A"]["satellites"] == []
        assert result["GO:B"]["role"] == "representative"
        assert result["GO:B"]["satellites"] == []

    def test_custom_overlap_threshold(self):
        df = _df([("GO:A", 0.01), ("GO:B", 0.02)])
        gene_sets = {
            "GO:A": ["g1", "g2", "g3", "g4"],
            "GO:B": ["g1", "g9"],   # overlap coeff = 1/2 = 0.5
        }
        # Default threshold (0.75): should NOT cluster.
        default_result = _cluster_redundant_terms(df, gene_sets, fdr_threshold=0.25)
        assert default_result["GO:B"]["role"] == "representative"
        # Lowered threshold: should cluster.
        loose_result = _cluster_redundant_terms(
            df, gene_sets, fdr_threshold=0.25, overlap_threshold=0.4,
        )
        assert loose_result["GO:B"]["role"] == "satellite"

    def test_exactly_at_threshold_clusters(self):
        """overlap_threshold is inclusive (>=), not a strict '>' cutoff."""
        df = _df([("GO:A", 0.01), ("GO:B", 0.02)])
        gene_sets = {
            "GO:A": ["g1", "g2", "g3", "g4"],
            "GO:B": ["g1", "g2", "g3"],   # overlap coeff exactly 1.0
        }
        result = _cluster_redundant_terms(
            df, gene_sets, fdr_threshold=0.25, overlap_threshold=1.0,
        )
        assert result["GO:B"]["role"] == "satellite"

    def test_non_significant_terms_excluded_entirely(self):
        df = _df([("GO:SIG", 0.01), ("GO:NOTSIG", 0.5)])
        gene_sets = {
            "GO:SIG":    ["g1", "g2", "g3"],
            "GO:NOTSIG": ["g1", "g2", "g3"],
        }
        result = _cluster_redundant_terms(df, gene_sets, fdr_threshold=0.25)
        assert "GO:SIG" in result
        assert "GO:NOTSIG" not in result

    def test_term_missing_from_gene_sets_becomes_own_representative(self):
        """Graceful degradation: a term absent from the GMT gene_sets dict
        (e.g. a GMT/results mismatch) is treated as having an empty gene
        set, so it can never overlap with anything and always ends up as
        its own representative rather than raising."""
        df = _df([("GO:KNOWN", 0.01), ("GO:MISSING", 0.02)])
        gene_sets = {
            "GO:KNOWN": ["g1", "g2", "g3"],
            # "GO:MISSING" intentionally absent from gene_sets.
        }
        result = _cluster_redundant_terms(df, gene_sets, fdr_threshold=0.25)
        assert result["GO:MISSING"]["role"] == "representative"
        assert result["GO:MISSING"]["satellites"] == []

    def test_three_terms_one_independent(self):
        df = _df([("GO:A", 0.01), ("GO:B", 0.02), ("GO:C", 0.03)])
        gene_sets = {
            "GO:A": ["g1", "g2", "g3", "g4"],
            "GO:B": ["g1", "g2", "g3"],       # redundant with A
            "GO:C": ["h1", "h2", "h3", "h4"],  # disjoint from A/B
        }
        result = _cluster_redundant_terms(df, gene_sets, fdr_threshold=0.25)
        assert result["GO:A"]["role"] == "representative"
        assert set(result["GO:A"]["satellites"]) == {"GO:B"}
        assert result["GO:B"]["representative"] == "GO:A"
        assert result["GO:C"]["role"] == "representative"
        assert result["GO:C"]["satellites"] == []

    def test_no_transitive_clustering_through_a_satellite(self):
        """A greedy algorithm compares each term only against already-
        chosen *representatives*, not against every previously-seen term.
        So if B joins A as a satellite, and C overlaps highly with B but
        only weakly with A directly, C must NOT be folded into A's
        cluster merely because of its relationship with B — there is no
        transitive closure through a satellite. This is a deliberate,
        documented property of _cluster_redundant_terms(), not an
        oversight, and is worth pinning down explicitly since it is easy
        to misread the greedy loop as doing full transitive clustering."""
        df = _df([("GO:A", 0.01), ("GO:B", 0.02), ("GO:C", 0.03)])
        gene_sets = {
            "GO:A": ["g1", "g2", "g3", "g4"],
            "GO:B": ["g1", "g2", "g3"],            # overlap(A,B) = 3/3 = 1.0 -> satellite of A
            "GO:C": ["g3", "x1", "x2"],            # overlap(A,C) = 1/3 = 0.33 (below threshold)
                                                    # overlap(B,C) = 1/3 = 0.33 (also below threshold here,
                                                    # but even if it were high, B is not a
                                                    # representative and so is never checked)
        }
        result = _cluster_redundant_terms(df, gene_sets, fdr_threshold=0.25)
        assert result["GO:B"]["representative"] == "GO:A"
        assert result["GO:C"]["role"] == "representative"

    def test_returns_empty_dict_when_nothing_significant(self):
        df = _df([("GO:A", 0.5), ("GO:B", 0.8)])
        gene_sets = {"GO:A": ["g1"], "GO:B": ["g1"]}
        result = _cluster_redundant_terms(df, gene_sets, fdr_threshold=0.25)
        assert result == {}

    def test_does_not_mutate_input_fdr_column(self):
        """The function must work from a numeric FDR column without
        requiring (or causing) any stringification of the caller's own
        DataFrame — _results_table_html() relies on 'FDR q-val' still
        being numeric immediately after this call, for its own
        subsequent formatting pass."""
        df = _df([("GO:A", 0.01), ("GO:B", 0.5)])
        gene_sets = {"GO:A": ["g1"], "GO:B": ["g1"]}
        _cluster_redundant_terms(df, gene_sets, fdr_threshold=0.25)
        assert pd.api.types.is_numeric_dtype(df["FDR q-val"])
