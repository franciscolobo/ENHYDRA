import os
from enhydra.report import _build_alignment_tree_html, _gmt_gene_sets, _resolve_term_names
from enhydra.msa_viewer import build_alignment_pages

# --- adjust these three paths to match your real run ---
GMT_PATH    = "test_data/hymenoptera/GO_annotation/Apis_mellifera_GO_BP.gmt"   # from the `find` command above
TABLES_DIR1 = "validation/hymenoptera_output_bee_vs_ant_paralogs_small/list1/tables"
ALN_DIR1    = "validation/hymenoptera_output_bee_vs_ant_paralogs_small/list1/alignment"
ALIGNMENTS1 = "validation/hymenoptera_output_bee_vs_ant_paralogs_small/list1/alignments"

TABLES_DIR2 = "validation/hymenoptera_output_bee_vs_ant_paralogs_small/list2/tables"
ALIGNMENTS2 = "validation/hymenoptera_output_bee_vs_ant_paralogs_small/list2/alignments"
# ---------------------------------------------------------

# Reuse pages already rendered by the real pipeline run (Commit 5), rather
# than re-rendering — build_alignment_pages() just needs to know where they
# already are, so pass group_ids=None to have it re-derive from group2anchor.tsv,
# or simplest: just reconstruct the {group_id: path} map directly from disk.
alignment_pages1 = {
    os.path.splitext(f)[0]: os.path.abspath(os.path.join(ALIGNMENTS1, f))
    for f in os.listdir(ALIGNMENTS1) if f.endswith(".html")
}
alignment_pages2 = {
    os.path.splitext(f)[0]: os.path.abspath(os.path.join(ALIGNMENTS2, f))
    for f in os.listdir(ALIGNMENTS2) if f.endswith(".html")
} if os.path.isdir(ALIGNMENTS2) else None

gmt_gene_sets = _gmt_gene_sets(GMT_PATH)
obo_names = {}  # empty is fine for this smoke test; GO IDs will show without names

print("Loaded %d gene sets from GMT" % len(gmt_gene_sets))
print("list1 alignment pages: %d" % len(alignment_pages1))
if alignment_pages2 is not None:
    print("list2 alignment pages: %d" % len(alignment_pages2))

tree_html = _build_alignment_tree_html(
    gmt_gene_sets=gmt_gene_sets,
    tables_dir1=TABLES_DIR1,
    obo_names=obo_names,
    alignment_pages1=alignment_pages1,
    alignment_pages2=alignment_pages2,
    report_dir=os.path.abspath("validation/hymenoptera_output_bee_vs_ant_paralogs_small"),
    label1="Pathogenic",   # adjust to your actual list1_name/list2_name
    label2="Non-pathogenic",
)

# Count how many <details> term blocks got built, as a quick sanity number
n_terms = tree_html.count('class="aln-tree-term"')
n_rows  = tree_html.count('class="aln-tree-row"')
print("Tree contains %d GO term blocks, %d group rows total" % (n_terms, n_rows))

out_path = "/tmp/alignment_tree_check.html"
with open(out_path, "w") as fh:
    fh.write("<html><body>%s</body></html>" % tree_html)
print("Wrote standalone snippet to: %s" % out_path)
