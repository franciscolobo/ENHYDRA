# ENHYDRA

**ENHYDRA** is a bioinformatics pipeline for Gene Set Enrichment Analysis (GSEA) in an evolutionary genomics context. Starting from groups of homologs, it filters, aligns, and ranks them by sequence identity, then applies the GSEA algorithm to detect functional enrichment patterns across a phylogenetic comparative framework.

---

## Overview

Traditional GSEA is designed for expression data. ENHYDRA extends this concept to evolutionary genomics by ranking homolog groups according to their average pairwise sequence identity across groups of homologs, and testing whether gene sets (e.g. Gene Ontology terms) are enriched at the top or bottom of this ranking. Functional annotations are mapped through an **anchor genome** — a reference species for which annotation data is available.

### Pipeline steps

1. **Length filtering** — sequences deviating more than 2 standard deviations from the group mean length are removed.
2. **Group filtering** — groups lacking the anchor species or falling below the minimum species count are discarded.
3. **Alignment** — surviving groups are aligned with MAFFT or MUSCLE.
4. **Identity estimation** — trimAl is used to compute per-alignment average sequence identity.
5. **Table generation** — outputs ranked group-to-identity and anchor gene-to-identity tables for downstream GSEA.

---

## Dependencies

| Tool / Library | Purpose |
|---|---|
| Python ≥ 3.x | Runtime |
| [Biopython](https://biopython.org/) | FASTA parsing (`Bio.SeqIO`, `Bio.Seq`) |
| [NumPy](https://numpy.org/) | Sorting sequence length data |
| [MAFFT](https://mafft.cbrc.jp/alignment/software/) | Multiple sequence alignment |
| [trimAl](http://trimal.cgenomics.org/) | Alignment trimming & identity computation |

Install Python dependencies with:

```bash
pip install biopython numpy
```

MAFFT, MUSCLE and trimAl must be installed separately and their paths provided in the code configuration file (see below).

---

## Input

ENHYDRA expects a directory of **FASTA files**, one per homolog group. Each file should contain sequences from multiple species. Sequence identifiers must follow the format:

```
>speciesID|geneID
```

where `speciesID` is used to identify species (including the anchor) and `geneID` is used to map functional annotations.

---

## Configuration

ENHYDRA uses two configuration files.

### Project configuration file

```
inputdir    = /path/to/input/fasta/directory
outdir      = /path/to/output/directory
min_species = 4
anchor      = Hsapiens
max_process = 8
```

| Parameter | Description |
|---|---|
| `inputdir` | Path to directory containing input FASTA files (one per homolog group) |
| `outdir` | Path to output directory (must not already exist) |
| `min_species` | Minimum number of distinct species required to retain a group |
| `anchor` | Species ID of the anchor genome used for annotation mapping |
| `max_process` | Number of parallel processes for the length filtering step |

### Code configuration file

```
mafft  = /usr/local/bin/mafft
trimal = /usr/local/bin/trimal
```

| Parameter | Description |
|---|---|
| `mafft` | Path to the MAFFT executable |
| `trimal` | Path to the trimAl executable |

---

## Usage

```bash
python enhydra.py <code_config> <project_config>
```

**Example:**

```bash
python enhydra.py code_config.txt project_config.txt
```

---

## Output

All output is written to the directory specified by `outdir`. The structure is as follows:

```
outdir/
├── outlog                        # Run log with filtered sequences and groups
├── length_stats/                 # Per-group length distribution statistics
│   └── <group>_lengthstats
├── length_filter/                # FASTA files after length filtering
│   └── <group>_lengthfilter
├── group_filter/                 # FASTA files after species/anchor filtering
├── alignment/                    # MAFFT alignments (.aln)
├── ident_alignment/              # trimAl identity reports (.ident)
└── tables/
    ├── group2mean.tsv            # Homolog group → mean alignment identity
    ├── anchor2mean.tsv           # Anchor gene ID → mean alignment identity
    └── group2anchor.tsv          # Homolog group → anchor gene ID
```

The files in `tables/` are the primary inputs for the GSEA step:
- **`anchor2mean.tsv`** provides the ranked gene list (anchor gene → identity score).
- **`group2anchor.tsv`** maps homolog groups to anchor genes for annotation lookup.
- **`group2mean.tsv`** provides group-level identity scores.

---

## How filtering works

### Length filter

For each homolog group, the mean (μ) and standard deviation (σ) of sequence lengths are computed. Sequences outside the range **[μ − 2σ, μ + 2σ]** are removed. Filtering is parallelised across groups using Python's `multiprocessing.Pool`.

### Group filter

After length filtering, groups are retained only if:
- They contain **at least one sequence from the anchor species**.
- They contain sequences from **at least `min_species` distinct species**.

---

## Citation

> *(Citation information will be added upon publication.)*

---

## License

> *(License information to be added.)*
