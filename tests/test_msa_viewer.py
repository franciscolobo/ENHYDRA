"""Tests for enhydra.msa_viewer."""

from enhydra.msa_viewer import parse_trimal_colnumbering


def _write(path, content):
    with open(path, "w") as fh:
        fh.write(content)


class TestParseTrimalColnumbering:

    def test_plain_comma_separated(self, tmp_path):
        f = tmp_path / "OG0001.colnumbering"
        _write(str(f), "0, 1, 2, 4, 5, 6, 9, 10\n")
        assert parse_trimal_colnumbering(str(f)) == {0, 1, 2, 4, 5, 6, 9, 10}

    def test_with_header_line(self, tmp_path):
        """Some trimAl versions prepend a '## ColumnsMap' style header."""
        f = tmp_path / "OG0002.colnumbering"
        _write(str(f), "## ColumnsMap\n0, 1, 2, 3, 5, 6\n")
        assert parse_trimal_colnumbering(str(f)) == {0, 1, 2, 3, 5, 6}

    def test_no_delimiter_whitespace_only(self, tmp_path):
        f = tmp_path / "OG0003.colnumbering"
        _write(str(f), "0 1 2 3 7 8\n")
        assert parse_trimal_colnumbering(str(f)) == {0, 1, 2, 3, 7, 8}

    def test_empty_file_returns_empty_set(self, tmp_path):
        f = tmp_path / "OG0004.colnumbering"
        _write(str(f), "")
        assert parse_trimal_colnumbering(str(f)) == set()

    def test_no_parseable_integers_returns_empty_set(self, tmp_path):
        f = tmp_path / "OG0005.colnumbering"
        _write(str(f), "no columns retained\n")
        assert parse_trimal_colnumbering(str(f)) == set()

    def test_missing_file_raises(self, tmp_path):
        import pytest
        with pytest.raises(FileNotFoundError):
            parse_trimal_colnumbering(str(tmp_path / "nonexistent.colnumbering"))

    def test_single_column_retained(self, tmp_path):
        f = tmp_path / "OG0006.colnumbering"
        _write(str(f), "42\n")
        assert parse_trimal_colnumbering(str(f)) == {42}

    def test_zero_is_a_valid_retained_column(self, tmp_path):
        """Column index 0 must not be treated as falsy/absent."""
        f = tmp_path / "OG0007.colnumbering"
        _write(str(f), "0\n")
        assert parse_trimal_colnumbering(str(f)) == {0}

import math
import pytest
from enhydra.msa_viewer import compute_column_stats


class TestComputeColumnStats:

    def test_all_identical_column(self):
        identities, entropies = compute_column_stats(["A", "A", "A", "A"])
        assert identities == [1.0]
        assert entropies == [0.0]

    def test_all_gap_column(self):
        identities, entropies = compute_column_stats(["-", "-", "-"])
        assert identities == [0.0]
        assert entropies == [0.0]

    def test_maximally_diverse_column(self):
        """4 distinct residues, each once: identity=1/4, entropy=2 bits
        normalised by log2(20)."""
        identities, entropies = compute_column_stats(["A", "C", "D", "E"])
        assert identities[0] == pytest.approx(0.25)
        expected_entropy = 2.0 / math.log2(20)
        assert entropies[0] == pytest.approx(expected_entropy, abs=1e-6)

    def test_three_of_four_same(self):
        identities, entropies = compute_column_stats(["A", "A", "A", "C"])
        assert identities[0] == pytest.approx(0.75)
        p_a, p_c = 0.75, 0.25
        expected_bits = -(p_a * math.log2(p_a) + p_c * math.log2(p_c))
        assert entropies[0] == pytest.approx(expected_bits / math.log2(20), abs=1e-6)

    def test_gaps_excluded_from_calculation(self):
        """Column 'A','C','-','A': gap excluded, identity computed over the
        3 non-gap residues only."""
        identities, entropies = compute_column_stats(["A", "C", "-", "A"])
        assert identities[0] == pytest.approx(2 / 3)
        p_a, p_c = 2 / 3, 1 / 3
        expected_bits = -(p_a * math.log2(p_a) + p_c * math.log2(p_c))
        assert entropies[0] == pytest.approx(expected_bits / math.log2(20), abs=1e-6)

    def test_dot_character_also_treated_as_gap(self):
        identities, _ = compute_column_stats(["A", "A", ".", "A"])
        assert identities[0] == pytest.approx(1.0)

    def test_tie_for_most_common_residue(self):
        """2 A's and 2 C's: identity uses the shared max count (2/4),
        regardless of which residue is 'the' most common."""
        identities, _ = compute_column_stats(["A", "A", "C", "C"])
        assert identities[0] == pytest.approx(0.5)

    def test_lowercase_residues_treated_as_uppercase(self):
        identities, entropies = compute_column_stats(["a", "A", "a", "A"])
        assert identities[0] == pytest.approx(1.0)
        assert entropies[0] == pytest.approx(0.0)

    def test_single_sequence_alignment(self):
        identities, entropies = compute_column_stats(["ACDE"])
        assert identities == [1.0, 1.0, 1.0, 1.0]
        assert entropies == [0.0, 0.0, 0.0, 0.0]

    def test_multi_column_alignment_returns_one_value_per_column(self):
        seqs = ["AAC-", "AAG-", "AAT-"]
        identities, entropies = compute_column_stats(seqs)
        assert len(identities) == 4
        assert len(entropies) == 4
        # Columns 0 and 1 are fully conserved.
        assert identities[0] == pytest.approx(1.0)
        assert identities[1] == pytest.approx(1.0)
        # Column 2 has three distinct residues (C, G, T), no majority.
        assert identities[2] == pytest.approx(1 / 3)
        # Column 3 is all-gap.
        assert identities[3] == 0.0
        assert entropies[3] == 0.0

    def test_empty_sequence_list_returns_empty_lists(self):
        assert compute_column_stats([]) == ([], [])

    def test_zero_width_sequences_return_empty_lists(self):
        assert compute_column_stats(["", "", ""]) == ([], [])

    def test_unequal_sequence_lengths_raises(self):
        with pytest.raises(ValueError, match="same length"):
            compute_column_stats(["ACD", "AC"])

    def test_entropy_values_always_bounded(self):
        """Sanity check across a variety of columns: entropy must stay in
        [0, 1] regardless of residue diversity."""
        test_columns = [
            ["A", "A", "A"],
            ["A", "C", "D", "E", "F", "G", "H", "I"],
            ["A", "-", "-"],
            ["A", "C"],
        ]
        for col in test_columns:
            _, entropies = compute_column_stats(col)
            assert 0.0 <= entropies[0] <= 1.0

    def test_identity_values_always_bounded(self):
        test_columns = [
            ["A", "A", "A"],
            ["A", "C", "D", "E"],
            ["-", "-", "-"],
        ]
        for col in test_columns:
            identities, _ = compute_column_stats(col)
            assert 0.0 <= identities[0] <= 1.0
