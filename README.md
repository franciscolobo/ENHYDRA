# ENHYDRA

**ENHYDRA** is a bioinformatics pipeline for Gene Set Enrichment Analysis (GSEA) in an evolutionary genomics context. Starting from groups of homologs, it filters, aligns, and ranks them by sequence identity, then applies the GSEA algorithm to detect functional enrichment patterns across a phylogenetic comparative framework.

---

## Overview

Traditional GSEA is designed for expression data. ENHYDRA extends this concept to evolutionary genomics by ranking homolog groups according to their average pairwise sequence identity across species, and testing whether gene sets (e.g. Gene Ontology terms) are enriched at the top or bottom of this ranking — revealing which biological functions are evolving faster or slower across a set of species.

Functional annotations are mapped through an **anchor genome** — a reference species for which GO annotation data is available.

### Pipeline steps

1. **[Optional] OrthoFinder preprocessing** — converts OrthoFinder 3 output into ENHYDRA-compatible input.
2. **Length filtering** — sequences deviating more than a configurable number of standard deviations from the group mean length are removed.
3. **Group filtering** — groups lacking the anchor species, falling below the minimum species count, or failing paralog criteria are discarded.
4. **Alignment** — surviving groups are aligned (MAFFT/MUSCLE/PRANK).
5. **[Optional] Divergent sequence filtering** — sequences that are highly divergent from the rest of their group are identified and removed, and the group is realigned without them. See [Divergent sequence filtering](#divergent-sequence-filtering).
6. **[Optional] Column trimming** — alignment columns are trimmed with trimAl before identity estimation.
7. **Identity estimation** — trimAl computes per-alignment average sequence identity.
8. **Table generation** — outputs ranked group-to-identity and anchor gene-to-identity tables.
9. **GSEA** — GSEApy prerank is run on the ranked gene list using a local GMT file.

---

## Dependencies

| Tool / Library | Purpose |
|---|---|
| Python ≥ 3.9 | Runtime |
| [Biopython](https://biopython.org/) | FASTA parsing |
| [NumPy](https://numpy.org/) | Sequence length sorting |
| [pandas](https://pandas.pydata.org/) | Ranked list handling |
| [GSEApy](https://gseapy.readthedocs.io/) | GSEA algorithm |
| [gprofiler-official](https://pypi.org/project/gprofiler-official/) | g:Profiler API access |
| [MAFFT](https://mafft.cbrc.jp/alignment/software/) | Multiple sequence alignment |
| [MUSCLE5](https://drive5.com/muscle/) | Multiple sequence alignment |
| [PRANK](http://wasabiapp.org/software/prank/) | Multiple sequence alignment |
| [trimAl](http://trimal.cgenomics.org/) | Alignment identity computation, column trimming, and divergent sequence detection |

Install Python dependencies:

```bash
pip install enhydra
```

MAFFT, MUSCLE, PRANK and trimAl must be installed separately and their paths provided in the code configuration file.

---

## Input

ENHYDRA accepts two types of input:

### Mode A — Pre-computed homolog groups
A directory of FASTA files, one per homolog group. Sequence identifiers must follow the format:

```
>speciesID|geneID
```

### Mode B — OrthoFinder 3 output
Pass the OrthoFinder output directory via `--orthofinder-dir`. ENHYDRA will automatically convert `Orthogroup_Sequences/` into the required format. Sequence headers in OrthoFinder output must follow the `speciesname_geneID` convention (automatically produced when input proteome files are named `<speciesID>.fa`).

---

## Preparing input proteomes

For bacterial datasets, ENHYDRA provides a helper script to download proteomes from NCBI and prepare them for OrthoFinder:

```bash
# Download proteomes (chromosome + plasmids merged per strain)
python tests/prepare_ecoli_proteomes.py strain_table.tsv proteomes/

# Remove pseudogenes from existing proteome FASTAs
python tests/remove_pseudogenes.py proteomes/
```

The strain table should be a TSV file with columns `genomeID`, `Pathogenicity`, and `plasmidIDs` (pipe-separated, or `-` if none):

```
genomeID    Pathogenicity   plasmidIDs
NC_004431   pathogenic      -
NC_002655   pathogenic      NC_007414
CP019778    non-pathogenic  -
```

Pseudogene removal is enabled by default (sequences flagged `[pseudo=true]` in efetch output). To keep pseudogenes:

```bash
python tests/prepare_ecoli_proteomes.py strain_table.tsv proteomes/ --keep-pseudogenes
```

---

## Building GMT files

ENHYDRA requires a GMT file for GSEA. The recommended approach is to run InterProScan on the anchor proteome and build the GMT from its output:

```bash
# Run InterProScan (requires separate installation)
interproscan.sh -i proteomes/NC_004431.fa -f TSV \
    -o NC_004431.interpro.tsv -goterms -pa

# Build GMT (one per namespace)
python tests/build_gmt_interproscan.py \
    --interproscan NC_004431.interpro.tsv \
    --anchor       NC_004431 \
    --outdir       gmt/ \
    --cache        obo_cache/ \
    --namespaces   GO_BP
```

To build all three GO namespaces at once:

```bash
python tests/build_gmt_interproscan.py \
    --interproscan NC_004431.interpro.tsv \
    --anchor       NC_004431 \
    --outdir       gmt/ \
    --cache        obo_cache/ \
    --namespaces   GO_BP GO_MF GO_CC
```

The GO OBO file is downloaded automatically and cached for reuse across runs and anchors.

Alternatively, annotations can be fetched from g:Profiler at runtime using `--organism` instead of `--gene-sets` (see [Usage](#usage)). This requires a stable internet connection and a g:Profiler-recognised organism code.

---

## Configuration

ENHYDRA uses two configuration files.

### Project configuration file

```
inputdir    = /path/to/input/fasta/directory
outdir      = /path/to/output/directory
min_species = 40
anchor      = NC_004431
max_process = 8
```

| Parameter | Description |
|---|---|
| `inputdir` | Directory containing input FASTA files (one per homolog group) |
| `outdir` | Output directory (must not already exist unless `--resume` is used) |
| `min_species` | Minimum number of distinct species required to retain a group |
| `length_filter_sd` | Number of standard deviations from a group's mean sequence length beyond which a sequence is removed (default: 2.0) |
| `divergence_filter_sd` | Number of standard deviations below a group's mean identity-to-closest-match (trimAl `-sident`) beyond which a sequence is removed as divergent. Empty by default (filter disabled). See [Divergent sequence filtering](#divergent-sequence-filtering) |
| `aligner` | Alignment tool to use: `mafft` (default), `muscle`, or `prank` |
| `anchor` | Species ID of the anchor genome used for annotation mapping |
| `max_process` | Number of parallel processes for alignment and length filtering |

### Code configuration file

```
mafft  = /usr/local/bin/mafft
trimal = /usr/local/bin/trimal
muscle = /usr/local/bin/muscle   # optional, required only if aligner = muscle
prank  = /usr/local/bin/prank    # optional, required only if aligner = prank
```

---

## Usage

```bash
enhydra <code_config> <project_config> --gene-sets <gmt_file> [options]
```

### Basic usage

```bash
enhydra code_config project_config --gene-sets gmt/NC_004431_GO_BP.gmt
```

### With OrthoFinder input

```bash
enhydra code_config project_config \
    --orthofinder-dir OrthoFinder/Results_Jan01 \
    --gene-sets gmt/NC_004431_GO_BP.gmt
```

### Resume an interrupted run

```bash
enhydra code_config project_config \
    --gene-sets gmt/NC_004431_GO_BP.gmt \
    --resume
```

Alignment, identity estimation, and table generation are skipped if their
output already exists. Plots and the HTML report are always regenerated so
that changes to `--fdr-threshold` or `--top-n` take effect without
rerunning the full pipeline.

### Re-run GSEA and regenerate plots on a completed run

```bash
enhydra code_config project_config \
    --gene-sets gmt/NC_004431_GO_BP.gmt \
    --replot \
    --fdr-threshold 0.10 \
    --top-n 30
```

Skips alignment, identity estimation, and table generation (treating them
as already complete), then re-runs GSEA, plots, and the HTML report with
the new parameters.  Combine with `--all-metrics` to switch between
single-metric and multi-metric output on an existing run.

### Run all three ranking metrics in one pass

```bash
enhydra code_config project_config \
    --gene-sets gmt/NC_004431_GO_BP.gmt \
    --all-metrics
```

Runs GSEA for identity, z-score, and rank in a single pipeline execution
and produces a single tabbed HTML report (`report.html`) with one tab per
metric. See [Ranking metrics](#ranking-metrics) for details.

### Paralog handling

```bash
# Keep all sequences including paralogs (default)
enhydra code_config project_config --gene-sets gmt/NC_004431_GO_BP.gmt

# Remove any group containing paralogs
enhydra code_config project_config --gene-sets gmt/NC_004431_GO_BP.gmt \
    --paralogs remove

# Keep only the longest sequence per species
enhydra code_config project_config --gene-sets gmt/NC_004431_GO_BP.gmt \
    --paralogs longest
```

### Removing highly divergent sequences

```bash
enhydra code_config project_config --gene-sets gmt/NC_004431_GO_BP.gmt \
    --divergence-filter-sd 2.0
```

Removes sequences whose identity to their closest match falls more than
2 standard deviations below their group's mean, then realigns the group
without them. See [Divergent sequence filtering](#divergent-sequence-filtering)
for details. Disabled by default.

### Full options

```
positional arguments:
  code_config           Path to the code configuration file.
  project_config        Path to the project configuration file.

options:
  --orthofinder-dir     Path to an OrthoFinder 3 output directory.
  --resume              Resume a previously interrupted run, skipping steps
                        whose output already exists. Plots and the HTML
                        report are always regenerated.
  --replot              Re-run GSEA, plots, and the HTML report on an
                        existing output directory without repeating
                        alignment or identity estimation. Use this to apply
                        a new --fdr-threshold, --top-n, --permutations, or
                        --all-metrics to a completed run. Implies --resume
                        for all pipeline steps preceding GSEA.
  --quiet               Suppress INFO/WARNING on the console; show progress
                        bars instead. All messages are still written to
                        enhydra.log.
  --paralogs {all,remove,longest}
                        How to handle paralogs (default: all).
  --min-species         Minimum number of distinct species required to
                        retain a group (default: 4, or the value of
                        min_species in the project config).
  --length-filter-sd    Number of standard deviations from a group's mean
                        sequence length beyond which a sequence is removed
                        (default: 2.0, or the value of length_filter_sd in
                        the project config).
  --divergence-filter-sd
                        Remove sequences whose identity to their closest
                        match (trimAl -sident) falls more than this many
                        standard deviations below their group's mean, then
                        realign the group without them. Disabled by default.
                        See Divergent sequence filtering below.
  --trim                Trim alignment columns with trimAl before identity
                        estimation. Accepts a number between 0 and 1 (used
                        as trimAl's -gt gap threshold), or one of: strict,
                        strictplus, automated. Disabled by default.
  --gene-sets           Path to a local .gmt file.
  --organism            g:Profiler organism name (alternative to --gene-sets).
  --sources             g:Profiler data sources (default: GO:BP GO:MF GO:CC KEGG REAC).
  --all-metrics         Run GSEA for all three ranking metrics (identity,
                        zscore, rank) in a single pass and produce a tabbed
                        HTML report. When set, --metric is ignored.
  --metric {identity,zscore,rank}
                        Ranking metric to use when --all-metrics is not set
                        (default: zscore). See Ranking metrics below.
  --permutations        Number of GSEA permutations (default: 1000).
  --min-size            Minimum gene set size (default: 5).
  --max-size            Maximum gene set size (default: 500).
  --seed                Random seed (default: 42).
  --fdr-threshold       FDR threshold for significance and plot filtering
                        (default: 0.25).
  --top-n               Maximum number of gene sets shown in bar plots
                        (default: 20).

two-list differential mode:
  --list1               Path to a file listing species IDs for group 1.
  --list2               Path to a file listing species IDs for group 2.
```

---

## Output

### Single-metric run

```
outdir/
├── enhydra.log
├── length_stats/
├── length_filter/
├── group_filter/
├── alignment/
├── divergence_sident/                    (only if divergence_filter_sd is set)
├── alignment_divergence_filtered/        (only if divergence_filter_sd is set)
├── divergence_realign_input/             (only if divergence_filter_sd is set)
├── divergence_filter_stats/              (only if divergence_filter_sd is set)
├── ident_alignment/
├── tables/
│   ├── group2mean.tsv
│   ├── anchor2mean.tsv
│   └── group2anchor.tsv
├── enrichment/
│   └── gseapy.gene_set.prerank.report.csv
├── plots/
│   ├── identity_distribution.png
│   └── gsea_barplot.png/.svg
├── alignments/
│   └── <group_id>.html                   (per-group alignment viewer pages)
└── report.html
```

### Multi-metric run (`--all-metrics`)

Alignment, filtering, and table generation run once. GSEA and plots are
produced separately for each metric in suffixed directories. A single tabbed
`report.html` collects all three.

```
outdir/
├── enhydra.log
├── length_stats/
├── length_filter/
├── group_filter/
├── alignment/
├── alignment_divergence_filtered/   (only if divergence_filter_sd is set)
├── ident_alignment/
├── tables/
├── enrichment_identity/   enrichment_zscore/   enrichment_rank/
├── plots_identity/        plots_zscore/        plots_rank/
└── report.html            ← single tabbed report
```

### Two-list differential run

```
outdir/
├── list1/   list2/        ← per-list filter/align/tables trees
├── differential/
│   ├── anchor2mean.tsv
│   ├── differential_scores.tsv
│   ├── enrichment/
│   ├── plots/
│   └── report.html
└── enhydra.log
```

Each of `list1/` and `list2/` follows the same per-list layout shown above
(`length_filter/`, `group_filter/`, `alignment/`, and — if
`divergence_filter_sd` is set — `alignment_divergence_filtered/` and
`divergence_filter_stats/`, etc.), independently for each list.

With `--all-metrics` the `differential/` directory is replaced by
`differential_identity/`, `differential_zscore/`, and `differential_rank/`,
and a single tabbed `report.html` is written at the root of `outdir/`.

---

## Interpreting GSEA results

The results table contains one row per tested gene set with the following key columns:

| Column | Description |
|---|---|
| `Name` | GO term ID |
| `Term` | GO term name |
| `NES` | Normalised enrichment score (positive = more conserved, negative = faster evolving) |
| `NOM p-val` | Nominal p-value |
| `FDR q-val` | FDR-corrected p-value |
| `Lead_genes` | Genes in the leading edge |

Significant gene sets are conventionally those with FDR < 0.25. To extract them:

```bash
awk -F',' 'NR==1 || $6 < 0.25' \
    outdir/enrichment/gseapy.gene_set.prerank.report.csv
```

Or in R:

```r
library(tidyverse)
read_csv("outdir/enrichment/gseapy.gene_set.prerank.report.csv") |>
    filter(`FDR q-val` < 0.25) |>
    arrange(`FDR q-val`) |>
    select(Name, Term, NES, `NOM p-val`, `FDR q-val`, Lead_genes)
```

The HTML report (`report.html`) provides an interactive table with per-column
filtering and clickable GO IDs that open individual enrichment plots. When
generated with `--all-metrics`, the report presents results for all three
metrics in separate tabs so they can be compared side by side.

The report also includes:

- An **Alignments** tab, showing a rendered, colour-coded alignment page for
  every orthogroup that fed GSEA, nested by GO term. If a group had a
  sequence removed by the divergence filter, its page shows the corrected
  (realigned) alignment along with a notice identifying which sequence was
  removed and its identity-to-closest-match value.
- A **Filtering summary** tab, showing a funnel of how many groups survived
  each pipeline stage and a breakdown of why groups or sequences were
  dropped at each step, with links to the affected group/sequence IDs. If
  `--divergence-filter-sd` was used, this includes the number of groups
  affected and sequences removed by that step.

---

## Ranking metrics

Three ranking metrics are available, selectable via `--metric` or all run
together with `--all-metrics`:

| Metric | Description |
|---|---|
| `identity` | Raw mean pairwise sequence identity. No normalisation applied. |
| `zscore` | Z-score normalised identity: (score − mean) / SD across all groups. Positive = more conserved than average. **Default.** |
| `rank` | Normalised rank: groups ranked by identity; most conserved group receives score 1.0, least conserved receives 1/*N*. |

All three metrics rank genes such that the most conserved appear at the top
of the GSEA input list, consistent with a positive NES indicating functional
conservation.

---

## Divergent sequence filtering

A single highly divergent sequence in an otherwise well-conserved homolog
group can distort that group's mean identity score — the exact metric
ENHYDRA ranks orthogroups by — even though the divergence may reflect a
misannotation, a paralog mistakenly grouped with true orthologs, or one
fast-evolving lineage dominating an otherwise conserved group's score.

When `divergence_filter_sd` (or `--divergence-filter-sd`) is set, ENHYDRA
runs trimAl's `-sident` analysis on each alignment and examines every
sequence's identity to its single closest match in the group. A sequence is
flagged and removed if this value falls more than `divergence_filter_sd`
standard deviations below the group's own mean. Because removing a sequence
invalidates the existing alignment for the survivors — the divergent
sequence often forces spurious gap placement elsewhere — the remaining
sequences are **realigned from scratch** rather than simply deleted from the
existing alignment, and identity is then computed on the corrected
alignment.

Notes:

- Groups with fewer sequences than `min_species` are exempt from this
  check: with too few sequences, a single outlier can drag the group's own
  mean and SD enough to escape detection (the same limitation that applies
  to the length filter — see below).
- If removing the flagged sequence(s) would leave a group below
  `min_species` or `min_sequences`, the entire group is dropped rather than
  kept with too few members.
- This is a single-pass filter: each group is checked once, pruned at most
  once, and realigned at most once. It does not iterate to re-check the
  realigned result for further outliers.
- The anchor species itself is eligible for removal if it is the divergent
  outlier; in that case the group survives (with its identity recomputed
  from the remaining species) but is excluded from `anchor2mean.tsv` for
  that run.
- Disabled by default. Recommended when input orthogroups may contain
  misannotated, mis-clustered, or unusually fast-evolving sequences (e.g.
  automatically generated OrthoFinder groups on divergent taxa).

---

## Notes on sequence identity and paralogs

ENHYDRA uses **protein-level** average pairwise sequence identity as the ranking metric. Synonymous mutations at the nucleotide level are invisible at the protein level, making protein identity a better proxy for functional conservation and a stronger signal for positive Darwinian selection.

Paralog handling is controlled by `--paralogs`:
- For **bacterial datasets**, `--paralogs all` (default) is recommended — duplicate proteins from multi-locus coding or horizontal gene transfer are biologically meaningful and correctly increase a group's average identity score.
- For **eukaryotic datasets**, `--paralogs longest` or `--paralogs remove` is recommended to avoid bias from within-species gene family expansions.

---

## Citation

> *(Citation information will be added upon publication.)*

---

## License

> *(License information to be added.)*
