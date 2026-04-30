# ENHYDRA

**ENHYDRA** is a bioinformatics pipeline for Gene Set Enrichment Analysis (GSEA) in an evolutionary genomics context. Starting from groups of homologs, it filters, aligns, and ranks them by sequence identity, then applies the GSEA algorithm to detect functional enrichment patterns across a phylogenetic comparative framework.

---

## Overview

Traditional GSEA is designed for expression data. ENHYDRA extends this concept to evolutionary genomics by ranking homolog groups according to their average pairwise sequence identity across species, and testing whether gene sets (e.g. Gene Ontology terms) are enriched at the top or bottom of this ranking — revealing which biological functions are evolving faster or slower across a set of species.

Functional annotations are mapped through an **anchor genome** — a reference species for which GO annotation data is available.

### Pipeline steps

1. **[Optional] OrthoFinder preprocessing** — converts OrthoFinder 3 output into ENHYDRA-compatible input.
2. **Length filtering** — sequences deviating more than 2 standard deviations from the group mean length are removed.
3. **Group filtering** — groups lacking the anchor species, falling below the minimum species count, or failing paralog criteria are discarded.
4. **Alignment** — surviving groups are aligned with MAFFT.
5. **Identity estimation** — trimAl computes per-alignment average sequence identity.
6. **Table generation** — outputs ranked group-to-identity and anchor gene-to-identity tables.
7. **GSEA** — GSEApy prerank is run on the ranked gene list using a local GMT file.

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
| [trimAl](http://trimal.cgenomics.org/) | Alignment identity computation |

Install Python dependencies:

```bash
pip install enhydra
```

MAFFT and trimAl must be installed separately and their paths provided in the code configuration file.

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
| `anchor` | Species ID of the anchor genome used for annotation mapping |
| `max_process` | Number of parallel processes for the length filtering step |

### Code configuration file

```
mafft  = /usr/local/bin/mafft
trimal = /usr/local/bin/trimal
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

### Full options

```
positional arguments:
  code_config           Path to the code configuration file.
  project_config        Path to the project configuration file.

options:
  --orthofinder-dir     Path to an OrthoFinder 3 output directory.
  --resume              Resume a previously interrupted run.
  --paralogs {all,remove,longest}
                        How to handle paralogs (default: all).
  --gene-sets           Path to a local .gmt file.
  --organism            g:Profiler organism name (alternative to --gene-sets).
  --sources             g:Profiler data sources (default: GO:BP GO:MF GO:CC KEGG REAC).
  --permutations        Number of GSEA permutations (default: 1000).
  --min-size            Minimum gene set size (default: 5).
  --max-size            Maximum gene set size (default: 500).
  --seed                Random seed (default: 42).
```

---

## Output

```
outdir/
├── enhydra.log                   # Run log
├── length_stats/                 # Per-group length distribution statistics
├── length_filter/                # FASTA files after length filtering
├── group_filter/                 # FASTA files after group filtering
├── alignment/                    # MAFFT alignments
├── ident_alignment/              # trimAl identity reports
├── tables/
│   ├── group2mean.tsv            # Homolog group → mean alignment identity
│   ├── anchor2mean.tsv           # Anchor gene ID → mean alignment identity
│   └── group2anchor.tsv          # Homolog group → anchor gene ID
└── enrichment/
    └── gseapy.gene_set.prerank.report.csv   # GSEA results
```

### Interpreting GSEA results

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
read_csv("outdir/enrichment/gseapy.gene_set.prerank.report.csv") %>%
    filter(`FDR q-val` < 0.25) %>%
    arrange(`FDR q-val`) %>%
    select(Name, Term, NES, `NOM p-val`, `FDR q-val`, Lead_genes)
```

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
