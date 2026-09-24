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
