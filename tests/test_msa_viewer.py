"""Tests for enhydra.msa_viewer."""

import os
from enhydra.msa_viewer import (
    parse_trimal_colnumbering,
    compute_column_stats,
    render_alignment_page,
    build_alignment_pages,
)

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

class TestRenderAlignmentPage:

    def _records(self):
        return [
            ("Species_a|GENE_A", "ACD-EFG"),
            ("Species_b|GENE_B", "ACDEEFG"),
            ("Species_c|GENE_C", "ACD-EFG"),
        ]

    def test_creates_file(self, tmp_path):
        out = tmp_path / "OG0001.html"
        render_alignment_page("OG0001", self._records(), str(out))
        assert out.is_file()

    def test_contains_group_and_metadata(self, tmp_path):
        out = tmp_path / "OG0001.html"
        render_alignment_page(
            "OG0001", self._records(), str(out),
            anchor_gene_id="ENSG00000141510",
            anchor_species="Species_a",
            mean_identity=0.85,
        )
        content = out.read_text()
        assert "OG0001" in content
        assert "ENSG00000141510" in content
        assert "Species_a" in content
        assert "0.8500" in content

    def test_group_id_suffix_stripped_from_title_and_metadata(self, tmp_path):
        """The pipeline-internal '_lengthfilter' suffix (see
        filtering.display_group_id()) is a file-correlation artifact, not
        part of the group's real name — it must never appear in the
        rendered page's title or 'Group ID' metadata row.

        Callers (build_alignment_pages()) are still expected to pass the
        raw, suffixed group_id here — that's the correct value for this
        function's own log messages and for correlating with whatever
        file this alignment came from — but the page's own human-facing
        text must show the stripped form. This is handled internally by
        render_alignment_page() itself, not by callers stripping the
        suffix beforehand.
        """
        out = tmp_path / "OG0001_lengthfilter.html"
        render_alignment_page("OG0001_lengthfilter", self._records(), str(out))
        content = out.read_text()
        assert "OG0001_lengthfilter" not in content
        assert "<title>Alignment: OG0001</title>" in content
        assert "<dt>Group ID</dt><dd>OG0001</dd>" in content

    def test_group_id_without_suffix_unaffected(self, tmp_path):
        """A group_id that never carried the pipeline suffix (e.g. any ID
        not derived from a length-filtered group) must render unchanged —
        display_group_id() is a documented no-op in that case."""
        out = tmp_path / "OG0002.html"
        render_alignment_page("OG0002", self._records(), str(out))
        content = out.read_text()
        assert "<title>Alignment: OG0002</title>" in content
        assert "<dt>Group ID</dt><dd>OG0002</dd>" in content

    def test_no_anchor_species_no_warning(self, tmp_path):
        out = tmp_path / "OG0001.html"
        render_alignment_page("OG0001", self._records(), str(out))
        content = out.read_text()
        assert '<div class="warning">' not in content

    def test_anchor_pinned_to_top(self, tmp_path):
        out = tmp_path / "OG0001.html"
        render_alignment_page(
            "OG0001", self._records(), str(out), anchor_species="Species_c",
        )
        content = out.read_text()
        # The pinned anchor's row-label should appear before the others.
        pos_c = content.index("Species_c|GENE_C")
        pos_a = content.index("Species_a|GENE_A")
        pos_b = content.index("Species_b|GENE_B")
        assert pos_c < pos_a
        assert pos_c < pos_b
    def test_anchor_not_found_adds_warning_and_preserves_order(self, tmp_path):
        out = tmp_path / "OG0001.html"
        render_alignment_page(
            "OG0001", self._records(), str(out), anchor_species="Nonexistent_species",
        )
        content = out.read_text()
        assert '<div class="warning">' in content
        assert "Nonexistent_species" in content
        pos_a = content.index("Species_a|GENE_A")
        pos_b = content.index("Species_b|GENE_B")
        pos_c = content.index("Species_c|GENE_C")
        assert pos_a < pos_b < pos_c   # original order unchanged

    def test_masked_columns_get_trimmed_class(self, tmp_path):
        out = tmp_path / "OG0001.html"
        # Records are 7 columns wide; mask out column index 3 (the gap column).
        retained = {0, 1, 2, 4, 5, 6}
        render_alignment_page(
            "OG0001", self._records(), str(out), retained_columns=retained,
        )
        content = out.read_text()
        assert 'class="aa trimmed"' in content
        assert "trimmed by trimAl" in content
        assert "Trimmed columns" in content
        assert "1 of 7 removed" in content

    def test_no_retained_columns_means_no_masking(self, tmp_path):
        out = tmp_path / "OG0001.html"
        render_alignment_page("OG0001", self._records(), str(out))
        content = out.read_text()
        assert 'class="aa trimmed"' not in content
        assert "Trimmed columns" not in content

    def test_all_aa_spans_present_including_gaps(self, tmp_path):
        """7 columns x 3 sequences = 21 residue spans, gaps included."""
        out = tmp_path / "OG0001.html"
        render_alignment_page("OG0001", self._records(), str(out))
        content = out.read_text()
        assert content.count('class="aa"') + content.count('class="aa trimmed"') == 21

    def test_block_width_wraps_into_multiple_blocks(self, tmp_path):
        records = [("Sp_a|G1", "A" * 130), ("Sp_b|G2", "A" * 130)]
        out = tmp_path / "OG0002.html"
        render_alignment_page("OG0002", records, str(out), block_width=60)
        content = out.read_text()
        # 130 columns / 60 per block -> 3 blocks (60, 60, 10)
        assert content.count('class="aln-block"') == 3
        assert "Columns 1\u201360" in content
        assert "Columns 121\u2013130" in content

    def test_raises_on_empty_records(self, tmp_path):
        import pytest
        out = tmp_path / "OG0003.html"
        with pytest.raises(ValueError, match="must not be empty"):
            render_alignment_page("OG0003", [], str(out))

    def test_raises_on_unequal_length_sequences(self, tmp_path):
        """Propagated from compute_column_stats()."""
        import pytest
        out = tmp_path / "OG0004.html"
        with pytest.raises(ValueError, match="same length"):
            render_alignment_page(
                "OG0004",
                [("Sp_a|G1", "ACDE"), ("Sp_b|G2", "AC")],
                str(out),
            )

    def test_mean_identity_none_shows_na(self, tmp_path):
        out = tmp_path / "OG0005.html"
        render_alignment_page("OG0005", self._records(), str(out), mean_identity=None)
        content = out.read_text()
        assert "N/A" in content

    def test_legend_present(self, tmp_path):
        out = tmp_path / "OG0006.html"
        render_alignment_page("OG0006", self._records(), str(out))
        content = out.read_text()
        assert "Aliphatic" in content
        assert "Aromatic" in content

    def test_out_path_accepts_str(self, tmp_path):
        """Confirms the function works with plain str paths (not just Path)."""
        out_str = str(tmp_path / "OG0007.html")
        render_alignment_page("OG0007", self._records(), out_str)
        import os
        assert os.path.isfile(out_str)

class TestBuildAlignmentPages:

    def _write_fasta(self, path, entries):
        with open(path, "w") as fh:
            for header, seq in entries:
                fh.write(">%s\n%s\n" % (header, seq))

    def _write_tsv(self, path, rows):
        with open(path, "w") as fh:
            for row in rows:
                fh.write("\t".join(row) + "\n")

    def _setup(self, tmp_path, groups: dict):
        """groups = {group_id: {"entries": [(header, seq)], "mean": float,
                                "anchor_gene": str | None}}"""
        aln_dir = tmp_path / "alignment"
        tables_dir = tmp_path / "tables"
        aln_dir.mkdir()
        tables_dir.mkdir()

        mean_rows, anchor_rows = [], []
        for gid, data in groups.items():
            self._write_fasta(str(aln_dir / (gid + ".aln")), data["entries"])
            mean_rows.append((gid, str(data["mean"])))
            if data.get("anchor_gene"):
                anchor_rows.append((gid, data["anchor_gene"]))

        self._write_tsv(str(tables_dir / "group2mean.tsv"), mean_rows)
        self._write_tsv(str(tables_dir / "group2anchor.tsv"), anchor_rows)
        return str(aln_dir), str(tables_dir)

    def test_default_group_ids_from_group2anchor(self, tmp_path):
        aln_dir, tables_dir = self._setup(tmp_path, {
            "OG0001": {"entries": [("Species_a|G1", "ACD"), ("Species_b|G2", "ACD")],
                      "mean": 0.9, "anchor_gene": "G1"},
            "OG0002": {"entries": [("Species_a|G3", "ACD"), ("Species_b|G4", "ACD")],
                      "mean": 0.8, "anchor_gene": "G3"},
        })
        outdir = str(tmp_path / "pages")
        pages = build_alignment_pages(aln_dir, tables_dir, outdir)
        assert set(pages.keys()) == {"OG0001", "OG0002"}
        for path in pages.values():
            assert os.path.isfile(path)
            assert os.path.isabs(path)

    def test_explicit_group_ids_restricts_rendering(self, tmp_path):
        aln_dir, tables_dir = self._setup(tmp_path, {
            "OG0001": {"entries": [("Species_a|G1", "ACD")], "mean": 0.9, "anchor_gene": "G1"},
            "OG0002": {"entries": [("Species_a|G2", "ACD")], "mean": 0.8, "anchor_gene": "G2"},
        })
        outdir = str(tmp_path / "pages")
        pages = build_alignment_pages(aln_dir, tables_dir, outdir, group_ids=["OG0001"])
        assert set(pages.keys()) == {"OG0001"}

    def test_group_missing_from_group2anchor_renders_without_anchor_gene(self, tmp_path):
        """Simulates a list2-style call: group_id explicitly requested but
        absent from this list's own (mostly-empty) group2anchor.tsv."""
        aln_dir, tables_dir = self._setup(tmp_path, {
            "OG0001": {"entries": [("Species_a|G1", "ACD")], "mean": 0.9, "anchor_gene": None},
        })
        outdir = str(tmp_path / "pages")
        pages = build_alignment_pages(aln_dir, tables_dir, outdir, group_ids=["OG0001"])
        content = open(pages["OG0001"]).read()
        assert "Anchor gene" not in content

    def test_missing_alignment_file_skipped(self, tmp_path):
        aln_dir, tables_dir = self._setup(tmp_path, {
            "OG0001": {"entries": [("Species_a|G1", "ACD")], "mean": 0.9, "anchor_gene": "G1"},
        })
        outdir = str(tmp_path / "pages")
        pages = build_alignment_pages(
            aln_dir, tables_dir, outdir, group_ids=["OG0001", "OG9999"],
        )
        assert set(pages.keys()) == {"OG0001"}

    def test_mean_identity_loaded_correctly(self, tmp_path):
        aln_dir, tables_dir = self._setup(tmp_path, {
            "OG0001": {"entries": [("Species_a|G1", "ACD")], "mean": 0.7321, "anchor_gene": "G1"},
        })
        outdir = str(tmp_path / "pages")
        pages = build_alignment_pages(aln_dir, tables_dir, outdir)
        content = open(pages["OG0001"]).read()
        assert "0.7321" in content

    def test_colnumbering_applies_masking(self, tmp_path):
        aln_dir, tables_dir = self._setup(tmp_path, {
            "OG0001": {"entries": [("Species_a|G1", "ACDEF")], "mean": 0.9, "anchor_gene": "G1"},
        })
        colnum_dir = tmp_path / "colnumbering"
        colnum_dir.mkdir()
        with open(str(colnum_dir / "OG0001.aln.colnumbering"), "w") as fh:
            fh.write("#ColumnsMap\t0, 1, 3, 4\n")   # column 2 trimmed

        outdir = str(tmp_path / "pages")
        pages = build_alignment_pages(
            aln_dir, tables_dir, outdir, colnumbering_dir=str(colnum_dir),
        )
        content = open(pages["OG0001"]).read()
        assert 'class="aa trimmed"' in content

    def test_no_colnumbering_dir_means_no_masking(self, tmp_path):
        aln_dir, tables_dir = self._setup(tmp_path, {
            "OG0001": {"entries": [("Species_a|G1", "ACDEF")], "mean": 0.9, "anchor_gene": "G1"},
        })
        outdir = str(tmp_path / "pages")
        pages = build_alignment_pages(aln_dir, tables_dir, outdir)
        content = open(pages["OG0001"]).read()
        assert 'class="aa trimmed"' not in content

    def test_anchor_gene_lookup_overrides_tables_dir(self, tmp_path):
        aln_dir, tables_dir = self._setup(tmp_path, {
            "OG0001": {"entries": [("Species_a|G1", "ACD")], "mean": 0.9, "anchor_gene": None},
        })
        outdir = str(tmp_path / "pages")
        pages = build_alignment_pages(
            aln_dir, tables_dir, outdir,
            group_ids=["OG0001"],
            anchor_gene_lookup={"OG0001": "ENSG00000141510"},
        )
        content = open(pages["OG0001"]).read()
        assert "ENSG00000141510" in content

    def test_anchor_species_pinning_applied(self, tmp_path):
        aln_dir, tables_dir = self._setup(tmp_path, {
            "OG0001": {"entries": [("Species_b|G2", "ACD"), ("Species_a|G1", "ACD")],
                      "mean": 0.9, "anchor_gene": "G1"},
        })
        outdir = str(tmp_path / "pages")
        pages = build_alignment_pages(
            aln_dir, tables_dir, outdir, anchor_species="Species_a",
        )
        content = open(pages["OG0001"]).read()
        assert content.index("Species_a|G1") < content.index("Species_b|G2")

    def test_outdir_created_if_missing(self, tmp_path):
        aln_dir, tables_dir = self._setup(tmp_path, {
            "OG0001": {"entries": [("Species_a|G1", "ACD")], "mean": 0.9, "anchor_gene": "G1"},
        })
        outdir = str(tmp_path / "does_not_exist_yet" / "pages")
        pages = build_alignment_pages(aln_dir, tables_dir, outdir)
        assert os.path.isdir(outdir)
        assert os.path.isfile(pages["OG0001"])

    def test_empty_group_ids_returns_empty_dict(self, tmp_path):
        aln_dir, tables_dir = self._setup(tmp_path, {})
        outdir = str(tmp_path / "pages")
        pages = build_alignment_pages(aln_dir, tables_dir, outdir, group_ids=[])
        assert pages == {}

    def test_group_id_with_suffix_displays_stripped_on_page(self, tmp_path):
        """End-to-end through the batch API: a group_id carrying the
        '_lengthfilter' suffix (as it would in a real pipeline run) must
        still be used correctly as the file-lookup/dict key (the returned
        pages dict, and the alignment/table file lookups), while the
        rendered page content itself shows only the stripped display name.
        """
        aln_dir, tables_dir = self._setup(tmp_path, {
            "OG0001_lengthfilter": {
                "entries": [("Species_a|G1", "ACD"), ("Species_b|G2", "ACD")],
                "mean": 0.9, "anchor_gene": "G1",
            },
        })
        outdir = str(tmp_path / "pages")
        pages = build_alignment_pages(aln_dir, tables_dir, outdir)

        # Lookup key / returned dict must remain the raw internal ID.
        assert set(pages.keys()) == {"OG0001_lengthfilter"}

        content = open(pages["OG0001_lengthfilter"]).read()
        assert "OG0001_lengthfilter" not in content
        assert "<dt>Group ID</dt><dd>OG0001</dd>" in content
