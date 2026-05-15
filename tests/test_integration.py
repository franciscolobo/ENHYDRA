"""Integration tests for the ENHYDRA pipeline.

Exercises the full filter → tables → GSEA → report path using pre-canned
fixture files in place of MAFFT and trimAl output. No external tools are
required.

Run selectively:
    pytest -m integration             # integration tests only
    pytest -m "not integration"       # unit tests only (fast CI)
"""
from __future__ import annotations

import os
import shutil

import pytest
from Bio import SeqIO

from enhydra.filtering import filter_length, filter_groups
from enhydra.tables import make_tables
from enhydra.gsea import run_gsea
from enhydra.differential import compute_differential
from enhydra.report import build_report, build_multi_metric_report
from enhydra.plotting import make_single_list_plots

pytestmark = pytest.mark.integration


# ---------------------------------------------------------------------------
# Constants (mirror generate_test_data.py)
# ---------------------------------------------------------------------------

ANCHOR      = "hsapiens"
SPECIES     = [ANCHOR, "mmusculus", "drerio", "ggallus", "xtropicalis"]
BASE_LEN    = 400
MUT_RATE    = 0.15
AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"

# Ensembl IDs used as anchor gene IDs — must match the GMT gene sets below.
ANCHOR_GENE_IDS = [
    "ENSG00000141510",  # OG000001  TP53
    "ENSG00000012048",  # OG000002  BRCA1
    "ENSG00000139618",  # OG000003  BRCA2
    "ENSG00000157764",  # OG000004  BRAF  — filtered (no anchor in group)
    "ENSG00000133703",  # OG000005  KRAS  — filtered (only 2 species)
    "ENSG00000171862",  # OG000006  PTEN
]

# Identity values assigned to the four groups that pass all filters.
IDENTITIES = {
    "OG000001": 0.85,
    "OG000002": 0.75,
    "OG000003": 0.65,
    "OG000006": 0.90,
}

# Two gene sets that overlap with the four passing anchor gene IDs.
_GMT_LINES = [
    "GO:0000001\ttest process alpha\t"
    "ENSG00000141510\tENSG00000012048\tENSG00000139618\tENSG00000171862",
    "GO:0000002\ttest process beta\t"
    "ENSG00000139618\tENSG00000171862",
]

# list2 identities (inverted rankings relative to list1) used in
# the differential mode fixture.
_LIST2_IDENTITIES = {
    "OG000001": 0.60,
    "OG000002": 0.70,
    "OG000003": 0.80,
    "OG000006": 0.65,
}


# ---------------------------------------------------------------------------
# Sequence helpers
# ---------------------------------------------------------------------------

def _rng_seq(length: int, seed: int) -> str:
    import random
    rng = random.Random(seed)
    return "".join(rng.choices(AMINO_ACIDS, k=length))


def _mutate(seq: str, rate: float, seed: int) -> str:
    import random
    rng = random.Random(seed)
    aa = list(seq)
    for i in range(len(aa)):
        if rng.random() < rate:
            aa[i] = rng.choice(AMINO_ACIDS)
    return "".join(aa)


def _write_fasta(path: str, entries: list[tuple[str, str]]) -> None:
    with open(path, "w") as fh:
        for header, seq in entries:
            fh.write(">%s\n%s\n" % (header, seq))


def _write_ident(path: str, identity: float) -> None:
    with open(path, "w") as fh:
        fh.write("## AverageIdentity %.6f\n" % identity)


# ---------------------------------------------------------------------------
# Session-scoped fixtures
# (each fixture runs once; downstream fixtures reuse its output path)
# ---------------------------------------------------------------------------

@pytest.fixture(scope="session")
def workspace(tmp_path_factory):
    """Root temp directory shared across all integration fixtures."""
    return tmp_path_factory.mktemp("integration")


@pytest.fixture(scope="session")
def input_dir(workspace):
    """Synthetic FASTA input groups mirroring generate_test_data.py.

    OG000001–3 : all 5 species, anchor present  → pass all filters
    OG000004   : no anchor species              → removed by filter_groups
    OG000005   : only 2 species                 → removed by filter_groups
    OG000006   : drerio sequence is 3× longer   → drerio removed by filter_length;
                 remaining 4 species pass
    """
    d = workspace / "input"
    d.mkdir()
    non_anchor = [sp for sp in SPECIES if sp != ANCHOR]

    # OG000001–3: clean groups
    for i in range(1, 4):
        name = "OG%06d" % i
        base = _rng_seq(BASE_LEN, seed=hash(name) & 0xFFFF)
        entries = [("%s|%s" % (ANCHOR, ANCHOR_GENE_IDS[i - 1]), base)]
        for sp in non_anchor:
            entries.append((
                "%s|%s_%s" % (sp, name, sp[:3].upper()),
                _mutate(base, MUT_RATE, seed=hash(name + sp) & 0xFFFF),
            ))
        _write_fasta(str(d / name), entries)

    # OG000004: no anchor
    name = "OG000004"
    base = _rng_seq(BASE_LEN, seed=40)
    _write_fasta(str(d / name), [
        ("%s|%s_%s" % (sp, name, sp[:3].upper()),
         _mutate(base, MUT_RATE, seed=hash(name + sp) & 0xFFFF))
        for sp in non_anchor
    ])

    # OG000005: only 2 species
    _write_fasta(str(d / "OG000005"), [
        ("%s|%s" % (ANCHOR, ANCHOR_GENE_IDS[4]), _rng_seq(BASE_LEN, seed=50)),
        ("mmusculus|OG000005_MMU",                _rng_seq(BASE_LEN, seed=51)),
    ])

    # OG000006: drerio is 3× longer than the rest
    name = "OG000006"
    base = _rng_seq(BASE_LEN, seed=60)
    entries = []
    for sp in SPECIES:
        gid = ANCHOR_GENE_IDS[5] if sp == ANCHOR else "%s_%s" % (name, sp[:3].upper())
        seq = (_rng_seq(BASE_LEN * 3, seed=61)
               if sp == "drerio"
               else _mutate(base, MUT_RATE, seed=hash(name + sp) & 0xFFFF))
        entries.append(("%s|%s" % (sp, gid), seq))
    _write_fasta(str(d / name), entries)

    return str(d)


@pytest.fixture(scope="session")
def filtered(workspace, input_dir):
    """Run filter_length (serially) then filter_groups.

    Returns:
        Tuple of (length_filter_dir, group_filter_dir).
    """
    stats_dir = workspace / "length_stats"
    lf_dir    = workspace / "length_filter"
    gf_dir    = workspace / "group_filter"
    for p in (stats_dir, lf_dir, gf_dir):
        p.mkdir()

    for fname in os.listdir(input_dir):
        filter_length(
            os.path.join(input_dir, fname),
            str(stats_dir), str(lf_dir),
        )

    filter_groups(
        length_filter_dir=str(lf_dir),
        group_filter_dir=str(gf_dir),
        anchor=ANCHOR,
        min_species=4,
    )
    return str(lf_dir), str(gf_dir)


@pytest.fixture(scope="session")
def alignment_dir(workspace):
    """Pre-canned alignment files — bypasses MAFFT entirely.

    Files are minimal valid FASTA sequences of uniform length (no real
    alignment required since make_tables only parses headers and the
    AverageIdentity is supplied separately via the ident fixture).
    """
    d = workspace / "alignment"
    d.mkdir()
    non_anchor = [sp for sp in SPECIES if sp != ANCHOR]

    # OG000001–3: all 5 species
    for i in range(1, 4):
        name = "OG%06d" % i
        base = _rng_seq(BASE_LEN, seed=hash(name) & 0xFFFF)
        entries = [("%s|%s" % (ANCHOR, ANCHOR_GENE_IDS[i - 1]), base)]
        for sp in non_anchor:
            entries.append((
                "%s|%s_%s" % (sp, name, sp[:3].upper()),
                _mutate(base, MUT_RATE, seed=hash(name + sp) & 0xFFFF),
            ))
        _write_fasta(str(d / (name + ".aln")), entries)

    # OG000006: 4 species (drerio already removed by filter_length)
    name = "OG000006"
    base = _rng_seq(BASE_LEN, seed=60)
    entries = [("%s|%s" % (ANCHOR, ANCHOR_GENE_IDS[5]), base)]
    for sp in [s for s in SPECIES if s not in (ANCHOR, "drerio")]:
        entries.append((
            "%s|%s_%s" % (sp, name, sp[:3].upper()),
            _mutate(base, MUT_RATE, seed=hash(name + sp) & 0xFFFF),
        ))
    _write_fasta(str(d / (name + ".aln")), entries)

    return str(d)


@pytest.fixture(scope="session")
def ident_dir(workspace):
    """Pre-canned trimAl identity files — bypasses trimAl entirely."""
    d = workspace / "ident"
    d.mkdir()
    for og, identity in IDENTITIES.items():
        _write_ident(str(d / (og + ".aln.ident")), identity)
    return str(d)


@pytest.fixture(scope="session")
def tables_dir(workspace, alignment_dir, ident_dir):
    """Run make_tables on the pre-canned alignment + ident fixtures."""
    d = str(workspace / "tables")
    make_tables(alignment_dir, ident_dir, d, anchor=ANCHOR)
    return d


@pytest.fixture(scope="session")
def gmt_file(workspace):
    """Tiny GMT whose gene IDs match the four passing anchor sequences."""
    path = workspace / "test.gmt"
    path.write_text("\n".join(_GMT_LINES) + "\n")
    return str(path)


@pytest.fixture(scope="session")
def gsea_dir(workspace, tables_dir, gmt_file):
    """Run GSEApy prerank on the fixture data."""
    d = str(workspace / "enrichment")
    run_gsea(
        anchor2mean_path=os.path.join(tables_dir, "anchor2mean.tsv"),
        results_dir=d,
        gene_sets=gmt_file,
        permutations=100,
        min_size=2,
        max_size=500,
        seed=42,
    )
    return d


@pytest.fixture(scope="session")
def diff_dir(workspace, tables_dir):
    """Compute differential scores between list1 (IDENTITIES) and list2
    (_LIST2_IDENTITIES), which have deliberately inverted rankings.

    OG000001: list1=0.85 > list2=0.60  → positive differential score
    OG000003: list1=0.65 < list2=0.80  → negative differential score
    """
    tbl2 = workspace / "tables2"
    tbl2.mkdir()

    with open(str(tbl2 / "group2mean.tsv"), "w") as fh:
        for og, val in _LIST2_IDENTITIES.items():
            fh.write("%s\t%.6f\n" % (og, val))

    # group2anchor is the same for both lists (same anchor gene IDs)
    shutil.copy(
        os.path.join(tables_dir, "group2anchor.tsv"),
        str(tbl2 / "group2anchor.tsv"),
    )

    d = str(workspace / "differential")
    compute_differential(
        tables_dir1=tables_dir,
        tables_dir2=str(tbl2),
        diff_dir=d,
        metric="zscore",
    )
    return d


# ---------------------------------------------------------------------------
# Tests: filtering
# ---------------------------------------------------------------------------

class TestFilteringIntegration:

    def test_og4_removed_no_anchor(self, filtered):
        _, gf_dir = filtered
        assert not any("OG000004" in n for n in os.listdir(gf_dir))

    def test_og5_removed_too_few_species(self, filtered):
        _, gf_dir = filtered
        assert not any("OG000005" in n for n in os.listdir(gf_dir))

    def test_og1_to_3_pass(self, filtered):
        _, gf_dir = filtered
        names = os.listdir(gf_dir)
        for i in range(1, 4):
            assert any("OG%06d" % i in n for n in names)

    def test_og6_drerio_removed_by_length_filter(self, filtered):
        lf_dir, _ = filtered
        lf_path = os.path.join(lf_dir, "OG000006_lengthfilter")
        ids = [r.id for r in SeqIO.parse(lf_path, "fasta")]
        assert not any(i.startswith("drerio") for i in ids)

    def test_og6_passes_group_filter_with_four_species(self, filtered):
        _, gf_dir = filtered
        assert any("OG000006" in n for n in os.listdir(gf_dir))


# ---------------------------------------------------------------------------
# Tests: tables
# ---------------------------------------------------------------------------

class TestTablesIntegration:

    def _read_tsv(self, path):
        with open(path) as fh:
            return [line.strip().split("\t") for line in fh if line.strip()]

    def test_group2mean_has_four_groups(self, tables_dir):
        rows = self._read_tsv(os.path.join(tables_dir, "group2mean.tsv"))
        assert len(rows) == 4

    def test_group2mean_values_match_fixtures(self, tables_dir):
        rows   = self._read_tsv(os.path.join(tables_dir, "group2mean.tsv"))
        result = {r[0]: float(r[1]) for r in rows}
        for og, expected in IDENTITIES.items():
            assert abs(result[og] - expected) < 1e-5

    def test_anchor2mean_has_correct_gene_ids(self, tables_dir):
        rows     = self._read_tsv(os.path.join(tables_dir, "anchor2mean.tsv"))
        gene_ids = {r[0] for r in rows}
        expected = {
            ANCHOR_GENE_IDS[0],   # OG000001
            ANCHOR_GENE_IDS[1],   # OG000002
            ANCHOR_GENE_IDS[2],   # OG000003
            ANCHOR_GENE_IDS[5],   # OG000006
        }
        assert gene_ids == expected

    def test_group2anchor_maps_correctly(self, tables_dir):
        rows    = self._read_tsv(os.path.join(tables_dir, "group2anchor.tsv"))
        mapping = {r[0]: r[1] for r in rows}
        assert mapping["OG000001"] == ANCHOR_GENE_IDS[0]
        assert mapping["OG000006"] == ANCHOR_GENE_IDS[5]


# ---------------------------------------------------------------------------
# Tests: GSEA
# ---------------------------------------------------------------------------

class TestGSEAIntegration:

    def test_results_csv_written(self, gsea_dir):
        assert os.path.isfile(
            os.path.join(gsea_dir, "gseapy.gene_set.prerank.report.csv")
        )

    def test_results_have_expected_columns(self, gsea_dir):
        import pandas as pd
        df = pd.read_csv(
            os.path.join(gsea_dir, "gseapy.gene_set.prerank.report.csv")
        )
        for col in ("Term", "NES", "NOM p-val", "FDR q-val"):
            assert col in df.columns

    def test_both_gene_sets_tested(self, gsea_dir):
        import pandas as pd
        df = pd.read_csv(
            os.path.join(gsea_dir, "gseapy.gene_set.prerank.report.csv")
        )
        assert set(df["Term"]) == {"GO:0000001", "GO:0000002"}


# ---------------------------------------------------------------------------
# Tests: differential
# ---------------------------------------------------------------------------

class TestDifferentialIntegration:

    def test_anchor2mean_written(self, diff_dir):
        assert os.path.isfile(os.path.join(diff_dir, "anchor2mean.tsv"))

    def test_differential_scores_written(self, diff_dir):
        assert os.path.isfile(os.path.join(diff_dir, "differential_scores.tsv"))

    def test_all_common_groups_present(self, diff_dir):
        import pandas as pd
        df = pd.read_csv(os.path.join(diff_dir, "differential_scores.tsv"), sep="\t")
        assert set(df["group_id"]) == set(IDENTITIES.keys())

    def test_score_signs_reflect_rankings(self, diff_dir):
        """OG000001 is more conserved in list1 than list2 → positive score.
        OG000003 is more conserved in list2 than list1 → negative score.
        This holds regardless of the normalisation details of zscore mode.
        """
        import pandas as pd
        df = (pd.read_csv(os.path.join(diff_dir, "differential_scores.tsv"), sep="\t")
                .set_index("group_id"))
        assert df.loc["OG000001", "score"] > 0
        assert df.loc["OG000003", "score"] < 0


# ---------------------------------------------------------------------------
# Tests: report
# ---------------------------------------------------------------------------

class TestReportIntegration:

    def test_single_metric_report_written(self, workspace, gsea_dir, tables_dir):
        plots_dir   = str(workspace / "plots_single")
        report_path = str(workspace / "report_single.html")
        make_single_list_plots(
            anchor2mean_path=os.path.join(tables_dir, "anchor2mean.tsv"),
            results_dir=gsea_dir,
            plots_dir=plots_dir,
        )
        build_report(
            results_dir=gsea_dir,
            plots_dir=plots_dir,
            report_path=report_path,
            mode="single",
        )
        assert os.path.isfile(report_path)
        content = open(report_path).read()
        assert "ENHYDRA" in content
        assert "GO:0000001" in content
        assert "GO:0000002" in content

    def test_multi_metric_report_has_all_tab_labels(self, workspace, gsea_dir, tables_dir):
        from enhydra.report import METRIC_LABELS
        metric_data = {}
        for m in ("identity", "zscore", "rank"):
            plots_dir = str(workspace / ("plots_multi_" + m))
            make_single_list_plots(
                anchor2mean_path=os.path.join(tables_dir, "anchor2mean.tsv"),
                results_dir=gsea_dir,
                plots_dir=plots_dir,
            )
            metric_data[m] = {"results_dir": gsea_dir, "plots_dir": plots_dir}

        report_path = str(workspace / "report_multi.html")
        build_multi_metric_report(
            metric_data=metric_data,
            report_path=report_path,
            mode="single",
        )
        assert os.path.isfile(report_path)
        content = open(report_path).read()
        for label in METRIC_LABELS.values():
            assert label in content
