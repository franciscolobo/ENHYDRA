# ENHYDRA

**ENHYDRA** is a Python pipeline that applies Gene Set Enrichment Analysis (GSEA) to evolutionary genomics. It starts from groups of homologous proteins (orthogroups), aligns them, ranks each group by its average pairwise sequence identity across species, and then asks which gene sets (for example Gene Ontology terms) are concentrated at the conserved or fast-evolving end of that ranking.

The result is a list of biological functions that are evolving unusually slowly or unusually fast across the species you supplied, delivered as an interactive HTML report.

---

## Table of contents

1. [How it works](#how-it-works)
2. [Pipeline steps](#pipeline-steps)
3. [Installation](#installation)
4. [Input data](#input-data)
5. [Gene sets (GMT files)](#gene-sets-gmt-files)
6. [Configuration](#configuration)
7. [Usage](#usage)
8. [Two-list differential mode](#two-list-differential-mode)
9. [Command-line reference](#command-line-reference)
10. [Output](#output)
11. [Interpreting results](#interpreting-results)
12. [Ranking metrics](#ranking-metrics)
13. [Divergent sequence filtering](#divergent-sequence-filtering)
14. [Paralogs and sequence identity](#paralogs-and-sequence-identity)
15. [Resuming, replotting and forking runs](#resuming-replotting-and-forking-runs)
16. [Comparing runs](#comparing-runs)
17. [Reproducibility](#reproducibility)
18. [Helper scripts and example datasets](#helper-scripts-and-example-datasets)
19. [Testing](#testing)
20. [Project layout](#project-layout)
21. [Troubleshooting](#troubleshooting)
22. [Citation and license](#citation-and-license)

---

## How it works

Classic GSEA ranks genes by a measurement such as differential expression and tests whether a gene set sits near the top or bottom of the list. ENHYDRA keeps the same statistics but changes what is ranked.

1. Each orthogroup receives a score: the mean pairwise protein sequence identity of its alignment.
2. Orthogroups are mapped to genes of an **anchor genome**, a reference species that has functional annotation (GO terms, for instance).
3. The anchor genes are ranked from most to least conserved and passed to GSEApy's `prerank` algorithm together with a GMT file of gene sets.
4. A gene set with a positive normalised enrichment score (NES) is concentrated among conserved groups. A gene set with a negative NES is concentrated among fast-evolving groups.

Protein-level identity is used rather than nucleotide identity because synonymous changes are invisible at the protein level, which makes it a better proxy for functional constraint.

---

## Pipeline steps

| # | Step | Description |
|---|---|---|
| 0 | OrthoFinder preprocessing *(optional)* | Converts OrthoFinder 3 `Orthogroup_Sequences/` into ENHYDRA's input format. |
| 1 | Length filtering | Removes sequences whose length is more than `length_filter_sd` standard deviations from their group's mean. |
| 2 | Group filtering | Drops groups that lack the anchor species, have fewer than `min_species` species, or fail the paralog rule. |
| 3 | Alignment | Aligns surviving groups with MAFFT, MUSCLE or PRANK, in parallel. |
| 3b | Anchor stripping *(two-list mode only)* | Removes an injected anchor from alignments before identity is computed. |
| 3c | Divergent sequence filtering *(optional)* | Removes outlier sequences using trimAl `-sident` and realigns the group without them. |
| 4 | Column trimming *(optional)* | Trims alignment columns with trimAl before identity estimation. |
| 5 | Identity estimation | trimAl computes average pairwise identity for each alignment. |
| 6 | Table generation | Writes ranked group-to-identity, anchor-gene-to-identity and group-to-anchor tables. |
| 7 | GSEA | Runs GSEApy prerank against a local GMT file or g:Profiler annotations, for one metric or all three. |
| 8 | Plots and report | Produces figures and a self-contained interactive `report.html`, including per-group alignment viewers. |

Every stage writes a small marker file recording the parameters it used, which is what makes safe resuming possible (see [Resuming, replotting and forking runs](#resuming-replotting-and-forking-runs)).

---

## Installation

### Python requirements

- Python 3.9 or newer
- [Biopython](https://biopython.org/), [NumPy](https://numpy.org/), [pandas](https://pandas.pydata.org/), [matplotlib](https://matplotlib.org/), [tqdm](https://tqdm.github.io/)
- [GSEApy](https://gseapy.readthedocs.io/) for the GSEA algorithm
- [gprofiler-official](https://pypi.org/project/gprofiler-official/) for optional g:Profiler annotation access

```bash
pip install enhydra
```

To install from a source checkout instead:

```bash
git clone <repository-url>
cd ENHYDRA
pip install -e .
```

### External tools

These must be installed separately. Their paths are given in the [code configuration file](#code-configuration-file).

| Tool | Required | Purpose |
|---|---|---|
| [MAFFT](https://mafft.cbrc.jp/alignment/software/) | Yes (default aligner) | Multiple sequence alignment |
| [trimAl](http://trimal.cgenomics.org/) | Yes | Identity computation, column trimming, divergent sequence detection |
| [MUSCLE 5](https://drive5.com/muscle/) | Only if `aligner = muscle` | Alternative aligner |
| [PRANK](http://wasabiapp.org/software/prank/) | Only if `aligner = prank` | Alternative aligner |
| [OrthoFinder 3](https://github.com/OrthoFinder/OrthoFinder) | Optional | Building orthogroups from proteomes |
| [InterProScan](https://www.ebi.ac.uk/interpro/about/interproscan/) | Optional | Annotating the anchor proteome to build a GMT file |

Find the paths of tools already on your system with `which mafft` and `which trimal`.

The HTML report loads jQuery and DataTables from a CDN on first generation and caches them, so the first report build needs internet access.

---

## Input data

ENHYDRA accepts homolog groups in one of two ways.

### Mode A: pre-computed homolog groups

A directory of FASTA files, one file per homolog group. Every sequence header must follow:

```
>speciesID|geneID
```

The `speciesID` is what ENHYDRA counts as a species. The `anchor` setting must match one of these IDs, and the `geneID` of the anchor's sequences must match the gene IDs in your GMT file.

### Mode B: OrthoFinder 3 output

Pass the OrthoFinder output directory with `--orthofinder-dir` (or `orthofinder_dir` in the project config). ENHYDRA converts `Orthogroup_Sequences/` automatically. You can point it at either the top-level OrthoFinder directory or a specific `Results_<date>` directory; if several runs exist, the most recent is used. Headers in OrthoFinder output must follow `speciesname_geneID`, which OrthoFinder produces when the input proteomes are named `<speciesID>.fa`.

---

## Gene sets (GMT files)

ENHYDRA needs gene sets keyed by the same gene IDs as the anchor's sequences. There are two sources.

### Local GMT file (recommended)

Supply a GMT file with `--gene-sets` (or `gene_sets` in the project config). A GMT file is tab-separated: `term_id`, `description`, then one gene ID per remaining column. ENHYDRA checks the overlap between the ranked genes and the GMT file up front and stops with a diagnostic message showing example IDs from both sides if the overlap is negligible.

Three helper scripts in `tests/` build GMT files:

**From InterProScan output** (bacteria and other genomes without curated annotation):

```bash
interproscan.sh -i proteomes/NC_004431.fa -f TSV \
    -o NC_004431.interpro.tsv -goterms -pa

python tests/build_gmt_interproscan.py \
    --interproscan NC_004431.interpro.tsv \
    --anchor       NC_004431 \
    --outdir       gmt/ \
    --cache        obo_cache/ \
    --namespaces   GO_BP GO_MF GO_CC
```

**From a GAF annotation file** (organisms with curated GO annotation; optionally propagates annotations up the GO hierarchy):

```bash
python tests/hymenoptera/build_gmt_gaf.py --gaf Amel.gaf --anchor Apis_mellifera ...
```

**From NCBI `gene2go` and `gene2refseq`** (NCBI-annotated genomes):

```bash
python tests/build_gmt.py --anchor NC_004431 --email you@example.org \
    --outdir gmt/ --cache ncbi_cache/ --namespaces GO_BP
```

`gene2refseq.gz` is about 1.2 GB; it is downloaded once and cached. The GO OBO file used by the InterProScan and GAF builders is also downloaded once and cached for reuse across runs and anchors.

### g:Profiler at runtime

Use `--organism` instead of `--gene-sets` to fetch annotations from g:Profiler, with `--sources` choosing the databases (default `GO:BP GO:MF GO:CC KEGG REAC`). This needs a stable internet connection and a g:Profiler-recognised organism name. `--gene-sets` and `--organism` are mutually exclusive, and one of them is required.

---

## Configuration

ENHYDRA reads two plain-text configuration files made of `key = value` lines. Blank lines and anything after `#` are ignored. Any project-config value can be overridden on the command line, and an explicit command-line flag always wins.

### Code configuration file

Paths to external executables. This file describes the machine, so it can be shared across projects.

```
mafft  = /usr/local/bin/mafft
trimal = /usr/local/bin/trimal
muscle = /usr/local/bin/muscle   # only needed if aligner = muscle
prank  = /usr/local/bin/prank    # only needed if aligner = prank
```

### Project configuration file

Minimal example:

```
inputdir    = /path/to/input/fasta/directory
outdir      = /path/to/output/directory
min_species = 40
anchor      = NC_004431
max_process = 8
```

All recognised keys:

| Key | Default | Description |
|---|---|---|
| `inputdir` | *required* | Directory of input FASTA files, one per homolog group. |
| `outdir` | *required* | Output directory. Must not already exist unless resuming (`--resume`, `--replot`) or forking (`--fork-from`). |
| `anchor` | *required* | Species ID of the anchor genome used for annotation mapping. |
| `max_process` | `1` | Parallel processes for length filtering and alignment. |
| `min_species` | `4` | Minimum number of distinct species a group needs to be kept. |
| `min_sequences` | `2` | Minimum number of sequences a group needs to be kept. |
| `paralogs` | `all` | Paralog policy: `all`, `remove` or `longest`. |
| `length_filter_sd` | `2.0` | Standard-deviation cutoff for the length filter. |
| `divergence_filter_sd` | *(off)* | Standard-deviation cutoff for the divergent sequence filter. Empty disables it. |
| `trim` | *(off)* | Column trimming: a number between 0 and 1 (trimAl `-gt`), or `strict`, `strictplus`, `automated`. |
| `aligner` | `mafft` | `mafft`, `muscle` or `prank`. |
| `mafft_mode` | `auto` | MAFFT strategy. `auto` uses MAFFT's `--auto`; any other value is passed as `--<value>` (for example `localpair`). |
| `metric` | `zscore` | Ranking metric: `identity`, `zscore` or `rank`. |
| `gene_sets` | | Path to a local GMT file. |
| `organism` | | g:Profiler organism name (alternative to `gene_sets`). |
| `sources` | `GO:BP GO:MF GO:CC KEGG REAC` | g:Profiler data sources. |
| `min_size` / `max_size` | `5` / `500` | Gene set size limits. |
| `permutations` | `1000` | GSEA permutations. |
| `seed` | `42` | Random seed. |
| `fdr_threshold` | `0.25` | FDR cutoff for significance and plot filtering. |
| `top_n` | `20` | Maximum gene sets shown in bar plots. |
| `gene_list_fdr_threshold` | same as `fdr_threshold` | FDR below which a gene set's full and leading-edge gene lists are embedded in the report. Raise it (1.0 or higher) to embed lists for nearly all gene sets, at the cost of a much larger report. |
| `obo_cache` | | Directory for the cached GO OBO file used to name terms in reports. |
| `list1`, `list2` | | Species-list files for [two-list mode](#two-list-differential-mode). |
| `list1_name`, `list2_name` | `List 1`, `List 2` | Display names for the two lists. |
| `orthofinder_dir` | | OrthoFinder output directory (same as `--orthofinder-dir`). |
| `resume`, `replot`, `quiet`, `all_metrics` | `false` | Boolean run-mode flags (`true`/`yes`/`on`/`1` and their opposites). |
| `fork_from` | | Same as `--fork-from`. |

---

## Usage

```bash
enhydra <code_config> <project_config> (--gene-sets <file.gmt> | --organism <name>) [options]
```

### Basic run

```bash
enhydra code_config project_config --gene-sets gmt/NC_004431_GO_BP.gmt
```

### OrthoFinder input

```bash
enhydra code_config project_config \
    --orthofinder-dir OrthoFinder/Results_Jan01 \
    --gene-sets gmt/NC_004431_GO_BP.gmt
```

### All three ranking metrics in one pass

```bash
enhydra code_config project_config \
    --gene-sets gmt/NC_004431_GO_BP.gmt \
    --all-metrics
```

Alignment, filtering and table generation run once; GSEA and plots run for each metric; a single tabbed `report.html` collects the results.

### Paralog handling

```bash
enhydra code_config project_config --gene-sets gmt/X.gmt --paralogs remove    # drop any group with paralogs
enhydra code_config project_config --gene-sets gmt/X.gmt --paralogs longest   # keep the longest sequence per species
```

### Column trimming

```bash
enhydra code_config project_config --gene-sets gmt/X.gmt --trim strict
enhydra code_config project_config --gene-sets gmt/X.gmt --trim 0.8
```

### Divergent sequence removal

```bash
enhydra code_config project_config --gene-sets gmt/X.gmt --divergence-filter-sd 2.0
```

See [Divergent sequence filtering](#divergent-sequence-filtering).

### Quiet mode

`--quiet` replaces console INFO/WARNING messages with progress bars. Everything is still written to `enhydra.log` in the output directory.

---

## Two-list differential mode

Two-list mode asks a comparative question: which gene sets are more conserved in one set of species than in another? Examples include pathogenic versus non-pathogenic strains, or ants versus bees.

Provide two text files, each listing one species ID per line (blank lines and `#` comments are ignored):

```bash
enhydra code_config project_config \
    --gene-sets gmt/X.gmt \
    --list1 pathogenic.txt --list2 non_pathogenic.txt \
    --metric zscore
```

How it works:

- Each list is filtered, aligned and scored independently (the `list1/` and `list2/` directories).
- Both lists are restricted to the orthogroups they share, and a signed differential score is computed per group: `identity_1 − identity_2` for `identity`, or the difference of the normalised scores for `zscore` and `rank`.
- Groups are mapped to anchor genes through list 1, so the GMT file must be built for the anchor.
- If the anchor is not already in list 1 it is added so groups can be mapped to genes, then stripped from the alignments before identity is computed so it does not influence the scores.
- **Sign convention:** a positive score, and therefore a positive NES, means *more conserved in list 1*. A negative one means more conserved in list 2.

Both lists must be given together, and both must be derived from the same set of homolog groups. Two-list mode works with `--all-metrics` too.

---

## Command-line reference

```
positional arguments:
  code_config             Path to the code configuration file.
  project_config          Path to the project configuration file.

run control:
  --orthofinder-dir DIR   OrthoFinder 3 output directory.
  --resume                Reuse a stage's existing output only if its marker shows it was
                          produced with identical parameters; otherwise recompute that stage
                          and everything downstream.
  --fork-from OUTDIR      Seed a new outdir from a completed run, reusing every stage whose
                          parameters are unchanged. Implies --resume. outdir must not exist.
  --replot                Re-run GSEA, plots and report without repeating alignment or
                          identity estimation. Implies --resume.
  --quiet                 Progress bars on the console; full messages go to enhydra.log.

filtering and alignment:
  --min-species N         Minimum species per group (default 4).
  --paralogs {all,remove,longest}
  --length-filter-sd X    Length-filter cutoff in SDs (default 2.0).
  --divergence-filter-sd X
                          Divergent sequence cutoff in SDs (default: disabled).
  --trim VALUE            0-1 gap threshold, or strict | strictplus | automated
                          (default: no trimming).

gene sets (mutually exclusive, one required):
  --gene-sets FILE        Local GMT file.
  --organism NAME         g:Profiler organism name.
  --sources S [S ...]     g:Profiler sources (default: GO:BP GO:MF GO:CC KEGG REAC).
  --obo-cache DIR         Cache directory for the GO OBO file.

GSEA and ranking:
  --metric {identity,zscore,rank}
                          Ranking metric (default zscore). Ignored with --all-metrics.
  --all-metrics           Run all three metrics and build a tabbed report.
  --permutations N        GSEA permutations (default 1000).
  --min-size N            Minimum gene set size (default 5).
  --max-size N            Maximum gene set size (default 500).
  --seed N                Random seed (default 42).

report:
  --fdr-threshold X       FDR cutoff for significance and plots (default 0.25).
  --gene-list-fdr-threshold X
                          FDR below which gene lists are embedded in the report
                          (default: same as --fdr-threshold).
  --top-n N               Maximum gene sets in bar plots (default 20).

two-list differential mode:
  --list1 FILE            Species IDs for group 1.
  --list2 FILE            Species IDs for group 2.
```

---

## Output

### Single-metric run

```
outdir/
├── enhydra.log
├── run_parameters.json                   effective parameters for this run
├── length_stats/  length_filter/  group_filter/  group_filter_stats/
├── alignment/
├── divergence_sident/                    only if divergence_filter_sd is set
├── alignment_divergence_filtered/        only if divergence_filter_sd is set
├── divergence_realign_input/             only if divergence_filter_sd is set
├── divergence_filter_stats/              only if divergence_filter_sd is set
├── alignment_trimmed/                    only if trim is set
├── ident_alignment/                      trimAl identity reports
├── tables/
│   ├── group2mean.tsv                    group → mean identity
│   ├── anchor2mean.tsv                   anchor gene → mean identity (the GSEA input)
│   ├── group2anchor.tsv                  group → anchor gene
│   └── drop_reasons.tsv                  groups excluded at this step, with reasons
├── enrichment/
│   └── gseapy.gene_set.prerank.report.csv
├── plots/
│   ├── identity_distribution.png/.pdf
│   └── gsea_barplot.png/.pdf/.svg
├── alignments/
│   └── <group_id>.html                   per-group alignment viewer pages
└── report.html
```

Exact stage directory names can vary slightly with the options used; the `Filtering summary` tab in the report always reflects what actually ran.

### Multi-metric run (`--all-metrics`)

Shared stages run once. GSEA and plots are written to suffixed directories (`enrichment_identity/`, `enrichment_zscore/`, `enrichment_rank/`, and the matching `plots_*/`), and one tabbed `report.html` sits at the root.

### Two-list run

```
outdir/
├── list1/   list2/          per-list filter, alignment and tables trees
├── differential/
│   ├── anchor2mean.tsv      differential score per anchor gene (GSEA input)
│   ├── differential_scores.tsv
│   ├── enrichment/
│   ├── plots/
│   └── report.html
├── run_parameters.json
└── enhydra.log
```

With `--all-metrics`, `differential/` becomes `differential_identity/`, `differential_zscore/` and `differential_rank/`, and one tabbed `report.html` is written at the root.

---

## Interpreting results

The GSEA table (`gseapy.gene_set.prerank.report.csv`) has one row per tested gene set. The key columns are:

| Column | Meaning |
|---|---|
| `Name` | Gene set ID (for example a GO term ID). |
| `Term` | Human-readable term name. |
| `NES` | Normalised enrichment score. Positive means more conserved; negative means faster evolving. In two-list mode, positive means more conserved in list 1. |
| `NOM p-val` | Nominal p-value. |
| `FDR q-val` | FDR-corrected p-value. |
| `Lead_genes` | Genes in the leading edge. |

Gene sets with FDR below 0.25 are conventionally treated as significant in GSEA; tighten this with `--fdr-threshold`. To pull the significant rows in R:

```r
library(tidyverse)
read_csv("outdir/enrichment/gseapy.gene_set.prerank.report.csv") |>
    filter(`FDR q-val` < 0.25) |>
    arrange(`FDR q-val`) |>
    select(Name, Term, NES, `NOM p-val`, `FDR q-val`, Lead_genes)
```

### The HTML report

`report.html` is the main deliverable. It contains:

- An interactive results table with per-column filtering and clickable gene set IDs that open the enrichment plot.
- Redundancy clustering of significant gene sets, so overlapping GO terms can be read together.
- An **Alignments** tab with a colour-coded alignment page for every orthogroup that fed GSEA, nested by gene set. If the divergence filter removed a sequence from a group, its page shows the corrected alignment and names the removed sequence and its identity to its closest match.
- A **Filtering summary** tab showing how many groups survived each stage and why groups or sequences were dropped, with links to the affected IDs.
- A **Run parameters** tab listing the effective parameters of the run.
- With `--all-metrics`, one tab per metric plus a cross-metric significance-agreement table.

---

## Ranking metrics

Select a metric with `--metric`, or run all three with `--all-metrics`.

| Metric | Description |
|---|---|
| `identity` | Raw mean pairwise sequence identity, with no normalisation. |
| `zscore` | `(score − mean) / SD` across all groups. Positive means more conserved than average. **Default.** |
| `rank` | Normalised rank. The most conserved group scores 1.0 and the least conserved scores 1/*N*. |

All three order genes so that the most conserved sit at the top of the GSEA input, which is why a positive NES always means conservation. The metrics usually agree on strong signals; running all three is a quick check that a result does not depend on how the scores were scaled.

---

## Divergent sequence filtering

A single highly divergent sequence can drag down a whole group's mean identity, even when the divergence reflects a misannotation, a mis-clustered paralog, or one fast-evolving lineage rather than the group's typical behaviour. This filter is off by default and is worth enabling for automatically generated orthogroups on divergent taxa.

When `divergence_filter_sd` is set, ENHYDRA runs trimAl `-sident` on each alignment and looks at each sequence's identity to its single closest match in the group. A sequence is removed when that value falls more than `divergence_filter_sd` standard deviations below the group's own mean. Because a divergent sequence often forces spurious gap placement, the survivors are **realigned from scratch** rather than simply deleted from the old alignment.

Details to be aware of:

- Groups with fewer sequences than `min_species` are exempt, since a tiny group gives an unreliable mean and SD.
- If removal would leave a group below `min_species` or `min_sequences`, the whole group is dropped.
- The filter runs once per group: one check, at most one removal step, at most one realignment. It does not iterate.
- The anchor can itself be removed if it is the outlier. The group then survives with identity recomputed from the remaining species but is left out of `anchor2mean.tsv`.
- It runs after alignment (and after anchor stripping in two-list mode) and before column trimming.

To confirm the filter is doing something, `tests/compare_divergence_filter_counts.py --before outdir/alignment --after outdir/alignment_divergence_filtered` reports per-group sequence counts before and after.

---

## Paralogs and sequence identity

Paralog handling is set with `--paralogs`:

- `all` (default) keeps every sequence. This is recommended for bacterial datasets, where multiple loci or horizontally transferred copies are biologically meaningful and legitimately raise a group's average identity.
- `longest` keeps only the longest sequence per species, and `remove` discards any group containing paralogs. Either is recommended for eukaryotic datasets, to avoid bias from within-species gene family expansions.

---

## Resuming, replotting and forking runs

**Resume** (`--resume`). Alignment is the slow step, so interrupted runs can continue. ENHYDRA writes a hidden marker file (`.stage_marker.json`) into each stage directory recording the parameters used. On resume, a stage is skipped only if its marker matches the current parameters. If it does not match, or the marker is missing or unreadable, that stage and everything downstream are recomputed with a warning. This prevents the silent reuse of stale output, for example resuming with a different `--min-species`. Plots and the report are always regenerated.

**Replot** (`--replot`). Treats alignment, identity estimation and tables as complete and re-runs only GSEA, plots and the report. Use it to change `--fdr-threshold`, `--top-n`, `--permutations`, `--metric` or `--all-metrics` on a finished run.

```bash
enhydra code_config project_config --gene-sets gmt/X.gmt \
    --replot --fdr-threshold 0.10 --top-n 30
```

**Fork** (`--fork-from`). Changes a parameter from some stage onward without recomputing what comes before it. ENHYDRA seeds a *new* output directory from a completed run, reuses every stage whose parameters are unchanged, and recomputes the first changed stage and everything after it. For example, switching `--trim` from `strict` to `strictplus` reuses length filtering, group filtering and alignment.

```bash
# project_config now points outdir at a fresh directory
enhydra code_config project_config --gene-sets gmt/X.gmt \
    --fork-from old_outdir --trim strictplus
```

The new `outdir` must not exist, and the source run must have used the same `inputdir`.

---

## Comparing runs

`enhydra.compare_runs` is a standalone tool that compares two or more completed runs and writes a single self-contained HTML report.

```bash
python -m enhydra.compare_runs \
    --runs outdir_A outdir_B outdir_C \
    --names "min_species=4" "min_species=8" "min_species=12" \
    --output comparison_report.html
```

The report has four tabs:

- **Run parameters**: an N-column diff of each run's `run_parameters.json`, with differing rows flagged.
- **Significance agreement**: majority agreement on significant gene sets across runs, for each ranking metric present in all of them.
- **Rank correlation**: Spearman correlation of GO term NES between runs, as a scatter for two runs or a correlation matrix with pairwise scatters for more.
- **Significance transitions**: Sankey-style diagrams showing gene sets that gain, lose or keep significance between each pair of runs.

Comparison is based on GO term identity, which stays meaningful across runs with different anchors, `min_species` or filtering, unlike orthogroup IDs. Runs must share at least one ranking metric to be compared on the data tabs, since comparing different metrics would conflate "the runs differ" with "the scoring differs".

---

## Reproducibility

Every run writes `run_parameters.json` to its output directory. It records the command line, the configuration file paths, the input, the anchor, and the *effective* values of every parameter after command-line overrides were applied, rather than the raw config file contents. It also records the mode (single-list or two-list), the metrics run, and, in two-list mode, the species in each list and whether the anchor was injected. Per-stage markers are summarised in the same file, and the report's Run parameters tab and `compare_runs` both read it. GSEA runs are seeded (`--seed`, default 42).

Software versions of the aligner, trimAl and Python packages are not currently recorded, so keep your own record of those.

---

## Helper scripts and example datasets

Data-preparation and example scripts live in `tests/`. They are not part of the installed package.

| Script | Purpose |
|---|---|
| `prepare_ecoli_proteomes.py` | Downloads *E. coli* proteomes from NCBI (chromosome and plasmids merged per strain) and prepares them for OrthoFinder. |
| `remove_pseudogenes.py` | Strips `[pseudo=true]` sequences from existing proteome FASTAs. |
| `reconstruct_orthogroups_E_coli.py` | Rebuilds per-orthogroup FASTA files (all, pathogenic, non-pathogenic) from OrthoFinder output. |
| `build_gmt_interproscan.py`, `build_gmt.py` | Build GMT files from InterProScan output or NCBI annotation. |
| `hymenoptera/prepare_hgd.py` | Converts Hymenoptera Genome Database proteomes and homolog tables into ENHYDRA input. |
| `hymenoptera/build_gmt_gaf.py` | Builds GMT files from a GAF annotation file, with optional GO hierarchy propagation. |
| `reconstruct_orthogroups_Hymenoptera.py` | Subsets groups by clade (ants, bees, sawflies, wasps) for two-list comparisons. |
| `generate_test_data.py` | Generates a small synthetic dataset that exercises each filter. |
| `compare_divergence_filter_counts.py` | Shows what the divergent sequence filter removed. |

### Bacterial example (*E. coli*)

```bash
python tests/prepare_ecoli_proteomes.py tests/E_coli_genomes.tsv proteomes/
```

The strain table is a TSV with columns `genomeID`, `Pathogenicity` and `plasmidIDs` (pipe-separated, or `-` for none):

```
genomeID    Pathogenicity   plasmidIDs
NC_004431   pathogenic      -
NC_002655   pathogenic      NC_007414
CP019778    non-pathogenic  -
```

Pseudogenes are removed by default; pass `--keep-pseudogenes` to retain them. Then run OrthoFinder on `proteomes/`, build a GMT file for the anchor strain, and run ENHYDRA either on all strains or in two-list mode with `tests/E_coli_pathogenic.txt` and `tests/E_coli_non_pathogenic.txt`.

### Synthetic example

`tests/generate_test_data.py` builds a six-group dataset with anchor `hsapiens` and four other species. It contains clean groups plus groups designed to be removed for a missing anchor, too few species, or an outlier sequence length. This is a quick way to confirm the install works. `tests/test_config_file` shows matching code and project configs.

---

## Testing

```bash
pip install pytest
pytest tests/
```

The suite covers configuration and I/O, filtering, statistics, GSEA, differential scoring, tables, the alignment viewer, provenance, stage markers, resume and fork behaviour, report redundancy clustering and the integration path. The integration tests use pre-canned fixture files in place of MAFFT and trimAl output, so no external tools are needed to run the tests. To run subsets:

```bash
pytest -m integration          # integration tests only
pytest -m "not integration"    # unit tests only (fast)
```

---

## Project layout

```
enhydra/
├── cli.py            command-line entry point and pipeline orchestration
├── io.py             config parsing, FASTA/OBO/species-list readers
├── utils.py          parameter validation, trim-argument resolution
├── orthofinder.py    OrthoFinder 3 preprocessing
├── filtering.py      length, group, paralog and divergent-sequence filters
├── alignment.py      MAFFT / MUSCLE / PRANK and trimAl wrappers
├── tables.py         ranked identity tables
├── gsea.py           GSEApy prerank and g:Profiler integration
├── differential.py   score normalisation and two-list differential scores
├── plotting.py       figures
├── report.py         interactive HTML report
├── msa_viewer.py     per-group alignment viewer pages
├── stats.py          pipeline funnel and drop-reason statistics
├── stage_markers.py  per-stage resume markers
├── provenance.py     run_parameters.json
├── compare_runs.py   multi-run comparison report
└── exceptions.py     EnhydraConfigError, EnhydraIOError, EnhydraToolError
tests/                unit and integration tests, data-prep scripts, example inputs
```

---

## Troubleshooting

| Symptom | Likely cause and fix |
|---|---|
| `Configuration error: ... executable not found` | A path in the code config is wrong or not executable. Check it with `which mafft`. |
| `outdir` already exists | ENHYDRA refuses to overwrite. Choose a new `outdir`, or use `--resume`, `--replot` or `--fork-from`. |
| GSEA fails with almost no gene overlap | Anchor gene IDs in your FASTA headers do not match the GMT gene IDs. Compare the example IDs in the error message. |
| `--fork-from` rejected | The destination exists, or the source used a different `inputdir`. |
| Stage recomputed on `--resume` | Its marker is missing or its parameters changed. This is intentional; the log names the changed parameter. |
| Two-list mode finds no common groups | The two lists were not derived from the same homolog groups. |
| Report fails to build offline | The first build downloads jQuery and DataTables. Run once with internet access. |

Everything is logged to `enhydra.log` in the output directory.

---

## Citation and license

Citation information will be added upon publication.

License information is to be added.

