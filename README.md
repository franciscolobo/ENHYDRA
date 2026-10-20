## Replace the existing "### Code configuration file" subsection with:

### Code configuration file

The code configuration file records paths to external tools (MAFFT, trimAl,
and optionally MUSCLE/PRANK). Unlike the project configuration file, it
rarely needs to change between analyses, so it is set up **once per
machine** rather than passed on every invocation.

Set it up automatically:

```bash
enhydra --init-code-config
```

This scans `$PATH` for `mafft`, `trimal`, `muscle`, and `prank`, and writes
the result to `~/.enhydra/code_config`:

```
mafft  = /usr/local/bin/mafft
trimal = /usr/local/bin/trimal
```

If a required tool (`mafft` or `trimal`) isn't found on `$PATH`, the file is
still written with a blank value for that key, and a warning tells you to
edit it in by hand. `muscle`/`prank` are optional — only required if you set
`aligner = muscle` or `aligner = prank` in your project config — and are
simply omitted if not found.

To write the file somewhere else, or to overwrite an existing one:

```bash
enhydra --init-code-config /path/to/code_config
enhydra --init-code-config --force   # overwrite ~/.enhydra/code_config
```

Once this file exists (at the default location, or pointed to by the
`ENHYDRA_CODE_CONFIG` environment variable), every `enhydra` invocation
finds it automatically — no need to specify it again. To use a different
code configuration for a single run, pass `--code-config` explicitly:

```bash
enhydra project_config --code-config /path/to/other_code_config --gene-sets ...
```

---

## Replace the entire "## Usage" section (through the end of "Full options") with:

## Usage

```bash
enhydra <project_config> --gene-sets <gmt_file> [options]
```

A code configuration file (tool paths) is resolved automatically — see
[Code configuration file](#code-configuration-file) above for the one-time
setup. Pass `--code-config <file>` to override it for a single run.

### Basic usage

```bash
enhydra project_config --gene-sets gmt/NC_004431_GO_BP.gmt
```

### With OrthoFinder input

```bash
enhydra project_config \
    --orthofinder-dir OrthoFinder/Results_Jan01 \
    --gene-sets gmt/NC_004431_GO_BP.gmt
```

### Resume an interrupted run

```bash
enhydra project_config \
    --gene-sets gmt/NC_004431_GO_BP.gmt \
    --resume
```

Alignment, identity estimation, and table generation are skipped if their
output already exists. Plots and the HTML report are always regenerated so
that changes to `--fdr-threshold` or `--top-n` take effect without
rerunning the full pipeline.

### Re-run GSEA and regenerate plots on a completed run

```bash
enhydra project_config \
    --gene-sets gmt/NC_004431_GO_BP.gmt \
    --replot \
    --fdr-threshold 0.10 \
    --top-n 30
```

Skips alignment, identity estimation, and table generation (treating them
as already complete), then re-runs GSEA, plots, and the HTML report with
the new parameters. Combine with `--all-metrics` to switch between
single-metric and multi-metric output on an existing run.

### Run all three ranking metrics in one pass

```bash
enhydra project_config \
    --gene-sets gmt/NC_004431_GO_BP.gmt \
    --all-metrics
```

Runs GSEA for identity, z-score, and rank in a single pipeline execution
and produces a single tabbed HTML report (`report.html`) with one tab per
metric. See [Ranking metrics](#ranking-metrics) for details.

### Paralog handling

```bash
# Keep all sequences including paralogs (default)
enhydra project_config --gene-sets gmt/NC_004431_GO_BP.gmt

# Remove any group containing paralogs
enhydra project_config --gene-sets gmt/NC_004431_GO_BP.gmt \
    --paralogs remove

# Keep only the longest sequence per species
enhydra project_config --gene-sets gmt/NC_004431_GO_BP.gmt \
    --paralogs longest
```

### Removing highly divergent sequences

```bash
enhydra project_config --gene-sets gmt/NC_004431_GO_BP.gmt \
    --divergence-filter-sd 2.0
```

Removes sequences whose identity to their closest match falls more than
2 standard deviations below their group's mean, then realigns the group
without them. See [Divergent sequence filtering](#divergent-sequence-filtering)
for details. Disabled by default.

### Full options

```
positional arguments:
  project_config        Path to the project configuration file. Not
                         required when using --init-code-config.

options:
  --code-config         Path to the code configuration file. Optional —
                        resolved automatically (see Code configuration
                        file above) if not given.
  --init-code-config [FILE]
                        Auto-detect mafft/trimal/muscle/prank on $PATH and
                        write a code configuration file, then exit.
                        Defaults to ~/.enhydra/code_config if FILE is
                        omitted.
  --force               With --init-code-config, overwrite an existing
                        file rather than raising an error.
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
  --fork-from OUTDIR    Seed a new output directory from a prior completed
                        run, reusing any pipeline stage whose parameters
                        are unchanged and recomputing the rest. Implies
                        --resume.
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
  --gene-list-fdr-threshold
                        (Advanced) FDR threshold below which a gene set's
                        full gene list / leading-edge gene list is
                        embedded in the report as a clickable link, rather
                        than shown as a count only. Defaults to the same
                        value as --fdr-threshold.
  --top-n               Maximum number of gene sets shown in bar plots
                        (default: 20).
  --obo-cache           Directory containing a cached go-basic.obo file,
                        for resolving GO term names.
  --package-dir         After a successful run, copy only the files needed
                        to view the report into this directory, as a
                        single self-contained, shareable folder.
  --package-overwrite   With --package-dir, overwrite an existing
                        directory rather than raising an error.

two-list differential mode:
  --list1               Path to a file listing species IDs for group 1.
  --list2               Path to a file listing species IDs for group 2.
```
