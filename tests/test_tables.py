"""Tests for enhydra.tables."""

import os
import pytest

from enhydra.tables import make_tables


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _write_fasta(path: str, entries: list[tuple[str, str]]):
    with open(path, "w") as fh:
        for header, seq in entries:
            fh.write(">%s\n%s\n" % (header, seq))


def _write_ident(path: str, identity: float):
    """Write a minimal trimAl -sident output with a single AverageIdentity line."""
    with open(path, "w") as fh:
        fh.write("## AverageIdentity %.6f\n" % identity)


def _read_tsv(path: str) -> list[tuple[str, str]]:
    rows = []
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            rows.append(tuple(parts))
    return rows


# ---------------------------------------------------------------------------
# make_tables
# ---------------------------------------------------------------------------

class TestMakeTables:

    def _setup(self, tmp_path, groups: dict[str, dict]) -> tuple[str, str, str]:
        """
        groups = {
            "OG0001": {
                "identity": 0.85,
                "sequences": [("sp1|geneA", "ACGT"), ("anchor|geneX", "ACGT")],
            }, ...
        }
        Returns (alignment_dir, ident_dir, tables_dir).
        """
        aln_dir   = str(tmp_path / "alignment")
        ident_dir = str(tmp_path / "ident")
        tbl_dir   = str(tmp_path / "tables")
        os.makedirs(aln_dir);  os.makedirs(ident_dir)

        for name, data in groups.items():
            aln_file   = os.path.join(aln_dir,   name + ".aln")
            ident_file = os.path.join(ident_dir, name + ".aln.ident")
            _write_fasta(aln_file, data["sequences"])
            _write_ident(ident_file, data["identity"])

        return aln_dir, ident_dir, tbl_dir

    def test_group2mean_written(self, tmp_path):
        aln, ident, tbl = self._setup(tmp_path, {
            "OG0001": {"identity": 0.85,
                       "sequences": [("sp1|g1", "ACGT"), ("anchor|gX", "ACGT")]},
        })
        make_tables(aln, ident, tbl, anchor="anchor")

        rows = _read_tsv(os.path.join(tbl, "group2mean.tsv"))
        assert len(rows) == 1
        assert rows[0][0] == "OG0001"
        assert abs(float(rows[0][1]) - 0.85) < 1e-5

    def test_anchor2mean_written(self, tmp_path):
        aln, ident, tbl = self._setup(tmp_path, {
            "OG0001": {"identity": 0.85,
                       "sequences": [("sp1|g1", "ACGT"), ("anchor|geneX", "ACGT")]},
        })
        make_tables(aln, ident, tbl, anchor="anchor")

        rows = _read_tsv(os.path.join(tbl, "anchor2mean.tsv"))
        assert len(rows) == 1
        assert rows[0][0] == "geneX"
        assert abs(float(rows[0][1]) - 0.85) < 1e-5

    def test_group2anchor_written(self, tmp_path):
        aln, ident, tbl = self._setup(tmp_path, {
            "OG0001": {"identity": 0.85,
                       "sequences": [("sp1|g1", "ACGT"), ("anchor|geneX", "ACGT")]},
        })
        make_tables(aln, ident, tbl, anchor="anchor")

        rows = _read_tsv(os.path.join(tbl, "group2anchor.tsv"))
        assert len(rows) == 1
        assert rows[0] == ("OG0001", "geneX")

    def test_multiple_groups(self, tmp_path):
        aln, ident, tbl = self._setup(tmp_path, {
            "OG0001": {"identity": 0.9,
                       "sequences": [("sp1|g1", "ACGT"), ("anchor|gA", "ACGT")]},
            "OG0002": {"identity": 0.7,
                       "sequences": [("sp1|g2", "ACGT"), ("anchor|gB", "ACGT")]},
        })
        make_tables(aln, ident, tbl, anchor="anchor")

        g2m = _read_tsv(os.path.join(tbl, "group2mean.tsv"))
        a2m = _read_tsv(os.path.join(tbl, "anchor2mean.tsv"))
        g2a = _read_tsv(os.path.join(tbl, "group2anchor.tsv"))
        assert len(g2m) == 2
        assert len(a2m) == 2
        assert len(g2a) == 2

    def test_write_mode_no_duplicates(self, tmp_path):
        """Calling make_tables twice on the same data must not duplicate rows."""
        aln, ident, tbl = self._setup(tmp_path, {
            "OG0001": {"identity": 0.85,
                       "sequences": [("sp1|g1", "ACGT"), ("anchor|gX", "ACGT")]},
        })
        make_tables(aln, ident, tbl, anchor="anchor")
        make_tables(aln, ident, tbl, anchor="anchor")

        rows = _read_tsv(os.path.join(tbl, "group2mean.tsv"))
        assert len(rows) == 1

    def test_group_without_anchor_excluded_from_anchor2mean(self, tmp_path):
        """A group with no anchor sequence contributes to group2mean but not anchor2mean."""
        aln, ident, tbl = self._setup(tmp_path, {
            "OG0001": {"identity": 0.9,
                       "sequences": [("sp1|g1", "ACGT"), ("anchor|gA", "ACGT")]},
            "OG0002": {"identity": 0.7,
                       "sequences": [("sp1|g2", "ACGT"), ("sp2|g3", "ACGT")]},
        })
        make_tables(aln, ident, tbl, anchor="anchor")

        g2m = _read_tsv(os.path.join(tbl, "group2mean.tsv"))
        a2m = _read_tsv(os.path.join(tbl, "anchor2mean.tsv"))
        assert len(g2m) == 2
        assert len(a2m) == 1
