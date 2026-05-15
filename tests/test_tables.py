"""Tests for enhydra.tables."""

import logging
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


def _read_tsv(path: str) -> list[tuple]:
    rows = []
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            rows.append(tuple(parts))
    return rows


# ---------------------------------------------------------------------------
# Fixtures / setup helpers
# ---------------------------------------------------------------------------

def _setup(tmp_path, groups: dict) -> tuple[str, str, str]:
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
    os.makedirs(aln_dir)
    os.makedirs(ident_dir)

    for name, data in groups.items():
        _write_fasta(os.path.join(aln_dir,   name + ".aln"), data["sequences"])
        _write_ident(os.path.join(ident_dir, name + ".aln.ident"), data["identity"])

    return aln_dir, ident_dir, tbl_dir


# ---------------------------------------------------------------------------
# Core output tests
# ---------------------------------------------------------------------------

class TestMakeTables:

    def test_group2mean_written(self, tmp_path):
        aln, ident, tbl = _setup(tmp_path, {
            "OG0001": {"identity": 0.85,
                       "sequences": [("sp1|g1", "ACGT"), ("anchor|gX", "ACGT")]},
        })
        make_tables(aln, ident, tbl, anchor="anchor")

        rows = _read_tsv(os.path.join(tbl, "group2mean.tsv"))
        assert len(rows) == 1
        assert rows[0][0] == "OG0001"
        assert abs(float(rows[0][1]) - 0.85) < 1e-5

    def test_anchor2mean_written(self, tmp_path):
        aln, ident, tbl = _setup(tmp_path, {
            "OG0001": {"identity": 0.85,
                       "sequences": [("sp1|g1", "ACGT"), ("anchor|geneX", "ACGT")]},
        })
        make_tables(aln, ident, tbl, anchor="anchor")

        rows = _read_tsv(os.path.join(tbl, "anchor2mean.tsv"))
        assert len(rows) == 1
        assert rows[0][0] == "geneX"
        assert abs(float(rows[0][1]) - 0.85) < 1e-5

    def test_group2anchor_written(self, tmp_path):
        aln, ident, tbl = _setup(tmp_path, {
            "OG0001": {"identity": 0.85,
                       "sequences": [("sp1|g1", "ACGT"), ("anchor|geneX", "ACGT")]},
        })
        make_tables(aln, ident, tbl, anchor="anchor")

        rows = _read_tsv(os.path.join(tbl, "group2anchor.tsv"))
        assert len(rows) == 1
        assert rows[0] == ("OG0001", "geneX")

    def test_multiple_groups(self, tmp_path):
        aln, ident, tbl = _setup(tmp_path, {
            "OG0001": {"identity": 0.9,
                       "sequences": [("sp1|g1", "ACGT"), ("anchor|gA", "ACGT")]},
            "OG0002": {"identity": 0.7,
                       "sequences": [("sp1|g2", "ACGT"), ("anchor|gB", "ACGT")]},
        })
        make_tables(aln, ident, tbl, anchor="anchor")

        assert len(_read_tsv(os.path.join(tbl, "group2mean.tsv")))   == 2
        assert len(_read_tsv(os.path.join(tbl, "anchor2mean.tsv")))  == 2
        assert len(_read_tsv(os.path.join(tbl, "group2anchor.tsv"))) == 2

    def test_write_mode_no_duplicates(self, tmp_path):
        """Calling make_tables twice on the same data must not duplicate rows."""
        aln, ident, tbl = _setup(tmp_path, {
            "OG0001": {"identity": 0.85,
                       "sequences": [("sp1|g1", "ACGT"), ("anchor|gX", "ACGT")]},
        })
        make_tables(aln, ident, tbl, anchor="anchor")
        make_tables(aln, ident, tbl, anchor="anchor")

        rows = _read_tsv(os.path.join(tbl, "group2mean.tsv"))
        assert len(rows) == 1

    def test_group_without_anchor_excluded_from_anchor2mean(self, tmp_path):
        """A group with no anchor sequence contributes to group2mean only."""
        aln, ident, tbl = _setup(tmp_path, {
            "OG0001": {"identity": 0.9,
                       "sequences": [("sp1|g1", "ACGT"), ("anchor|gA", "ACGT")]},
            "OG0002": {"identity": 0.7,
                       "sequences": [("sp1|g2", "ACGT"), ("sp2|g3",    "ACGT")]},
        })
        make_tables(aln, ident, tbl, anchor="anchor")

        assert len(_read_tsv(os.path.join(tbl, "group2mean.tsv")))  == 2
        assert len(_read_tsv(os.path.join(tbl, "anchor2mean.tsv"))) == 1


# ---------------------------------------------------------------------------
# Cross-validation / mismatch tests
# ---------------------------------------------------------------------------

class TestMakeTablesMismatch:

    def test_ident_without_alignment_skipped(self, tmp_path):
        """A group whose identity file has no matching alignment is skipped
        entirely — it must not appear in any output table."""
        aln_dir   = str(tmp_path / "alignment")
        ident_dir = str(tmp_path / "ident")
        tbl_dir   = str(tmp_path / "tables")
        os.makedirs(aln_dir)
        os.makedirs(ident_dir)

        # Identity file present; alignment file deliberately absent.
        _write_ident(os.path.join(ident_dir, "OG0001.aln.ident"), 0.85)

        make_tables(aln_dir, ident_dir, tbl_dir, anchor="anchor")

        assert _read_tsv(os.path.join(tbl_dir, "group2mean.tsv"))   == []
        assert _read_tsv(os.path.join(tbl_dir, "anchor2mean.tsv"))  == []
        assert _read_tsv(os.path.join(tbl_dir, "group2anchor.tsv")) == []

    def test_ident_without_alignment_warns(self, tmp_path, caplog):
        """The missing-alignment case must emit a WARNING."""
        aln_dir   = str(tmp_path / "alignment")
        ident_dir = str(tmp_path / "ident")
        tbl_dir   = str(tmp_path / "tables")
        os.makedirs(aln_dir)
        os.makedirs(ident_dir)

        _write_ident(os.path.join(ident_dir, "OG0001.aln.ident"), 0.85)

        with caplog.at_level(logging.WARNING, logger="enhydra.tables"):
            make_tables(aln_dir, ident_dir, tbl_dir, anchor="anchor")

        messages = [r.message for r in caplog.records]
        assert any("OG0001" in m for m in messages)

    def test_alignment_without_ident_warns(self, tmp_path, caplog):
        """An alignment file that has no matching identity file must emit a WARNING.
        The pre-check catches this before the main loop so the group is never
        silently added to any table."""
        aln_dir   = str(tmp_path / "alignment")
        ident_dir = str(tmp_path / "ident")
        tbl_dir   = str(tmp_path / "tables")
        os.makedirs(aln_dir)
        os.makedirs(ident_dir)

        # Alignment present; identity file deliberately absent.
        _write_fasta(
            os.path.join(aln_dir, "OG0001.aln"),
            [("sp1|g1", "ACGT"), ("anchor|gX", "ACGT")],
        )

        with caplog.at_level(logging.WARNING, logger="enhydra.tables"):
            make_tables(aln_dir, ident_dir, tbl_dir, anchor="anchor")

        messages = [r.message for r in caplog.records]
        assert any("OG0001" in m for m in messages)
        # Group must not appear in any table.
        assert _read_tsv(os.path.join(tbl_dir, "group2mean.tsv")) == []

    def test_valid_and_mismatched_groups_together(self, tmp_path):
        """A run mixing a valid group and a group with a missing alignment
        produces correct output for the valid group and skips the broken one."""
        aln_dir   = str(tmp_path / "alignment")
        ident_dir = str(tmp_path / "ident")
        tbl_dir   = str(tmp_path / "tables")
        os.makedirs(aln_dir)
        os.makedirs(ident_dir)

        # OG0001: valid — both files present.
        _write_fasta(
            os.path.join(aln_dir,   "OG0001.aln"),
            [("sp1|g1", "ACGT"), ("anchor|gA", "ACGT")],
        )
        _write_ident(os.path.join(ident_dir, "OG0001.aln.ident"), 0.9)

        # OG0002: broken — identity file present, alignment absent.
        _write_ident(os.path.join(ident_dir, "OG0002.aln.ident"), 0.7)

        make_tables(aln_dir, ident_dir, tbl_dir, anchor="anchor")

        g2m = _read_tsv(os.path.join(tbl_dir, "group2mean.tsv"))
        a2m = _read_tsv(os.path.join(tbl_dir, "anchor2mean.tsv"))

        assert len(g2m) == 1
        assert g2m[0][0] == "OG0001"
        assert len(a2m) == 1
        assert a2m[0][0] == "gA"
