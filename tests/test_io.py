"""Tests for enhydra.io."""

import io
import pytest

from enhydra.io import (
    read_config_file, read_species_list, parse_obo_names, parse_fasta_records,
)


# ---------------------------------------------------------------------------
# read_config_file
# ---------------------------------------------------------------------------

class TestReadConfigFile:

    def _make_fh(self, text: str):
        return io.StringIO(text)

    def test_required_fields_parsed(self):
        project = self._make_fh(
            "inputdir    = /data/input\n"
            "outdir      = /data/output\n"
            "anchor      = NC_004431\n"
            "max_process = 8\n"
        )
        code = self._make_fh(
            "mafft  = /usr/bin/mafft\n"
            "trimal = /usr/bin/trimal\n"
        )
        p = read_config_file(project, code)
        assert p["inputdir"]    == "/data/input"
        assert p["outdir"]      == "/data/output"
        assert p["anchor"]      == "NC_004431"
        assert p["max_process"] == 8
        assert p["mafft"]       == "/usr/bin/mafft"
        assert p["trimal"]      == "/usr/bin/trimal"

    def test_defaults_applied(self):
        project = self._make_fh(
            "inputdir    = /data/input\n"
            "outdir      = /data/output\n"
            "anchor      = NC_004431\n"
            "max_process = 4\n"
        )
        code = self._make_fh("")
        p = read_config_file(project, code)
        assert p["min_species"]      == 4
        assert p["min_sequences"]    == 2
        assert p["paralogs"]         == "all"
        assert p["metric"]           == "zscore"
        assert p["permutations"]     == 1000
        assert p["fdr_threshold"]    == 0.25
        assert p["aligner"]          == "mafft"
        assert p["length_filter_sd"] == 2.0

    def test_comments_ignored(self):
        project = self._make_fh(
            "inputdir = /data/input  # trailing comment\n"
            "# full line comment\n"
            "outdir = /data/output\n"
            "anchor = NC_004431\n"
            "max_process = 2\n"
        )
        code = self._make_fh("")
        p = read_config_file(project, code)
        assert p["inputdir"] == "/data/input"

    def test_optional_tool_paths_empty_by_default(self):
        project = self._make_fh(
            "inputdir = /in\noutdir = /out\nanchor = sp1\nmax_process = 1\n"
        )
        code = self._make_fh("")
        p = read_config_file(project, code)
        assert p["muscle"] == ""
        assert p["prank"]  == ""

    def test_cli_override_values_parsed(self):
        """list1/list2 and sources are parsed from project config."""
        project = self._make_fh(
            "inputdir = /in\noutdir = /out\nanchor = sp1\nmax_process = 1\n"
            "list1 = pathogenic.txt\nlist2 = non_pathogenic.txt\n"
            "sources = GO:BP KEGG\n"
        )
        code = self._make_fh("")
        p = read_config_file(project, code)
        assert p["list1"]   == "pathogenic.txt"
        assert p["sources"] == "GO:BP KEGG"


# ---------------------------------------------------------------------------
# read_species_list
# ---------------------------------------------------------------------------

class TestReadSpeciesList:

    def test_reads_species(self, tmp_path):
        f = tmp_path / "species.txt"
        f.write_text("NC_001\nNC_002\nNC_003\n")
        assert read_species_list(str(f)) == ["NC_001", "NC_002", "NC_003"]

    def test_ignores_comments_and_blanks(self, tmp_path):
        f = tmp_path / "species.txt"
        f.write_text("# header\nNC_001\n\nNC_002\n")
        assert read_species_list(str(f)) == ["NC_001", "NC_002"]

    def test_empty_file_raises(self, tmp_path):
        f = tmp_path / "species.txt"
        f.write_text("# only a comment\n\n")
        with pytest.raises(ValueError, match="empty"):
            read_species_list(str(f))

    def test_missing_file_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            read_species_list(str(tmp_path / "nonexistent.txt"))


# ---------------------------------------------------------------------------
# parse_obo_names
# ---------------------------------------------------------------------------

class TestParseOboNames:

    def _write_obo(self, path, terms):
        """
        terms = [{"id": "GO:0001", "name": "foo", "obsolete": False}, ...]
        """
        lines = ["format-version: 1.2\n\n"]
        for t in terms:
            lines.append("[Term]\n")
            lines.append("id: %s\n" % t["id"])
            lines.append("name: %s\n" % t["name"])
            lines.append("namespace: biological_process\n")
            if t.get("obsolete"):
                lines.append("is_obsolete: true\n")
            lines.append("\n")
        with open(path, "w") as fh:
            fh.writelines(lines)

    def test_parses_term_names(self, tmp_path):
        obo = str(tmp_path / "go.obo")
        self._write_obo(obo, [
            {"id": "GO:0001234", "name": "some process", "obsolete": False},
            {"id": "GO:0005678", "name": "other process", "obsolete": False},
        ])
        names = parse_obo_names(obo)
        assert names["GO:0001234"] == "some process"
        assert names["GO:0005678"] == "other process"

    def test_obsolete_terms_excluded(self, tmp_path):
        obo = str(tmp_path / "go.obo")
        self._write_obo(obo, [
            {"id": "GO:0001234", "name": "active term",   "obsolete": False},
            {"id": "GO:0009999", "name": "obsolete term", "obsolete": True},
        ])
        names = parse_obo_names(obo)
        assert "GO:0009999" not in names
        assert "GO:0001234" in names

    def test_empty_obo_returns_empty_dict(self, tmp_path):
        obo = str(tmp_path / "go.obo")
        with open(obo, "w") as fh:
            fh.write("format-version: 1.2\n")
        assert parse_obo_names(obo) == {}


# ---------------------------------------------------------------------------
# parse_fasta_records
# ---------------------------------------------------------------------------
#
# Regression coverage for a real crash: newer Biopython versions
# (>= 1.85) raise ValueError from Bio.SeqIO.parse(path, "fasta") if a
# file contains any lines before its first '>' record header — a pattern
# some external tools/databases produce (e.g. a metadata banner), and
# which older Biopython versions silently tolerated. parse_fasta_records()
# strips any such leading lines itself before delegating to Bio.SeqIO,
# so it never depends on Biopython's own newer 'fasta-pearson'/
# 'fasta-blast' formats (which don't exist on older Biopython installs).

class TestParseFastaRecords:

    def test_well_formed_fasta_unaffected(self, tmp_path):
        f = tmp_path / "OG0001"
        f.write_text(">sp1|g1\nACDE\n>sp2|g2\nACDF\n")
        records = parse_fasta_records(str(f))
        assert [r.id for r in records] == ["sp1|g1", "sp2|g2"]
        assert [str(r.seq) for r in records] == ["ACDE", "ACDF"]

    def test_leading_hash_comment_stripped(self, tmp_path):
        f = tmp_path / "OG0002"
        f.write_text("# generated by some external tool\n>sp1|g1\nACDE\n")
        records = parse_fasta_records(str(f))
        assert len(records) == 1
        assert records[0].id == "sp1|g1"
        assert str(records[0].seq) == "ACDE"

    def test_leading_semicolon_comment_stripped(self, tmp_path):
        f = tmp_path / "OG0003"
        f.write_text("; Pearson-style comment\n>sp1|g1\nACDE\n")
        records = parse_fasta_records(str(f))
        assert len(records) == 1
        assert records[0].id == "sp1|g1"

    def test_multiple_leading_comment_and_blank_lines_stripped(self, tmp_path):
        f = tmp_path / "OG0004"
        f.write_text(
            "# banner line 1\n"
            "# banner line 2\n"
            "\n"
            ";another style of comment\n"
            "\n"
            ">sp1|g1\nACDE\n>sp2|g2\nACDF\n"
        )
        records = parse_fasta_records(str(f))
        assert [r.id for r in records] == ["sp1|g1", "sp2|g2"]

    def test_leading_comments_do_not_truncate_first_record(self, tmp_path):
        """The whole first record (header + full multi-line sequence) must
        survive intact — only the junk strictly before the first '>' is
        removed."""
        f = tmp_path / "OG0005"
        f.write_text(
            "# comment\n"
            ">sp1|g1\nACDEFGHIKL\nMNPQRSTVWY\n>sp2|g2\nACDE\n"
        )
        records = parse_fasta_records(str(f))
        assert str(records[0].seq) == "ACDEFGHIKLMNPQRSTVWY"
        assert str(records[1].seq) == "ACDE"

    def test_completely_empty_file_returns_empty_list(self, tmp_path):
        f = tmp_path / "OG0006"
        f.write_text("")
        assert parse_fasta_records(str(f)) == []

    def test_comments_only_no_sequences_returns_empty_list(self, tmp_path):
        f = tmp_path / "OG0007"
        f.write_text("# just a comment\n; another comment\n")
        assert parse_fasta_records(str(f)) == []

    def test_missing_file_raises_file_not_found(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            parse_fasta_records(str(tmp_path / "nonexistent"))

    def test_no_leading_junk_identical_record_count_and_content(self, tmp_path):
        """Sanity check against Bio.SeqIO.parse directly for a
        well-formed file, confirming no behavioural drift for the common
        case."""
        from Bio import SeqIO
        f = tmp_path / "OG0008"
        f.write_text(">sp1|g1\nACDE\n>sp2|g2\nACDF\n>sp3|g3\nACDG\n")
        via_helper  = parse_fasta_records(str(f))
        via_biopython = list(SeqIO.parse(str(f), "fasta"))
        assert [r.id for r in via_helper] == [r.id for r in via_biopython]
        assert [str(r.seq) for r in via_helper] == [str(r.seq) for r in via_biopython]

    def test_enhydra_header_convention_species_pipe_gene_preserved(self, tmp_path):
        """Confirms downstream code's r.id.split('|')[0] species extraction
        still works correctly on records recovered from a file with
        leading comments."""
        f = tmp_path / "OG0009"
        f.write_text("# comment\n>Hsapiens|ENSG00000141510\nACDE\n")
        records = parse_fasta_records(str(f))
        assert records[0].id.split("|")[0] == "Hsapiens"
        assert records[0].id.split("|")[1] == "ENSG00000141510"
