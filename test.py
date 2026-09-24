from enhydra.msa_viewer import build_alignment_pages

pages = build_alignment_pages(
    alignment_dir="validation/hymenoptera/hymenoptera_output_bee_vs_ant_orthologs/list1/alignment",
    tables_dir="validation/hymenoptera/hymenoptera_output_bee_vs_ant_orthologs/list1/tables",
    outdir="/tmp/alignment_pages_smoke_test",
    show_progress=True,
)

print("Rendered %d pages" % len(pages))

