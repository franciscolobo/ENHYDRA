"""
Generate a synthetic ENHYDRA test dataset.

Groups produced and what they test:
  OG000001  5 species, anchor present  → passes all filters (clean group)
  OG000002  5 species, anchor present  → passes all filters (clean group)
  OG000003  5 species, anchor present  → passes all filters (clean group)
  OG000004  5 species, NO anchor       → removed by filter_groups (anchor check)
  OG000005  2 species, anchor present  → removed by filter_groups (min_species=4)
  OG000006  5 species, anchor present, one seq too long/short
                                       → seq removed by filter_length

Anchor species : hsapiens
Other species  : mmusculus, drerio, ggallus, xtropicalis
Sequence type  : amino acid (single-letter codes)
"""

import os
import random

random.seed(42)

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
ANCHOR      = "hsapiens"
SPECIES     = [ANCHOR, "mmusculus", "drerio", "ggallus", "xtropicalis"]
BASE_LENGTH = 400   # aa residues
MUTATION_RATE = 0.15

# Real human Ensembl gene IDs — used as anchor gene IDs so g:Profiler
# can map them to GO terms. These are well-annotated housekeeping genes.
ANCHOR_GENE_IDS = [
    "ENSG00000141510",  # TP53
    "ENSG00000012048",  # BRCA1
    "ENSG00000139618",  # BRCA2
    "ENSG00000157764",  # BRAF
    "ENSG00000133703",  # KRAS
    "ENSG00000171862",  # PTEN
]


def random_sequence(length: int) -> str:
    return "".join(random.choices(AMINO_ACIDS, k=length))


def mutate(seq: str, rate: float) -> str:
    """Introduce point mutations at a given rate."""
    aa = list(seq)
    for i in range(len(aa)):
        if random.random() < rate:
            aa[i] = random.choice(AMINO_ACIDS)
    return "".join(aa)


def write_group(outdir: str, name: str, entries: list[tuple[str, str, str]]):
    """Write a FASTA file. entries = [(species, gene_id, sequence)]"""
    path = os.path.join(outdir, name)
    with open(path, "w") as fh:
        for species, gene_id, seq in entries:
            fh.write(">%s|%s\n%s\n" % (species, gene_id, seq))
    print("Written: %s (%d sequences)" % (path, len(entries)))


def make_group(name: str, species_list: list[str], base_seq: str,
               anchor_gene_id: str,
               mutation_rate: float = MUTATION_RATE) -> list[tuple]:
    """Generate one homolog group with per-species mutations."""
    entries = []
    for sp in species_list:
        # Anchor uses the real Ensembl ID; others get a placeholder
        if sp == ANCHOR:
            gene_id = anchor_gene_id
        else:
            gene_id = "%s_%s" % (name, sp[:3].upper())
        seq = mutate(base_seq, mutation_rate)
        entries.append((sp, gene_id, seq))
    return entries


def main(outdir: str = "test_data/input"):
    os.makedirs(outdir, exist_ok=True)

    # OG000001–OG000003: clean groups, all species, anchor present
    for i in range(1, 4):
        name = "OG%06d" % i
        base = random_sequence(BASE_LENGTH)
        entries = make_group(name, SPECIES, base, anchor_gene_id=ANCHOR_GENE_IDS[i - 1])
        write_group(outdir, name, entries)

    # OG000004: no anchor species → removed by filter_groups
    name = "OG000004"
    non_anchor = [sp for sp in SPECIES if sp != ANCHOR]
    base = random_sequence(BASE_LENGTH)
    entries = make_group(name, non_anchor, base, anchor_gene_id=ANCHOR_GENE_IDS[3])
    write_group(outdir, name, entries)

    # OG000005: only 2 species (anchor + 1) → removed by filter_groups
    name = "OG000005"
    base = random_sequence(BASE_LENGTH)
    entries = make_group(name, [ANCHOR, "mmusculus"], base, anchor_gene_id=ANCHOR_GENE_IDS[4])
    write_group(outdir, name, entries)

    # OG000006: anchor present, 5 species, but one seq is an outlier length
    name = "OG000006"
    base = random_sequence(BASE_LENGTH)
    entries = make_group(name, SPECIES, base, anchor_gene_id=ANCHOR_GENE_IDS[5])
    # Replace drerio sequence with one that is 3x longer (clear outlier)
    outlier_seq = random_sequence(BASE_LENGTH * 3)
    entries = [(sp, gid, outlier_seq if sp == "drerio" else seq)
               for sp, gid, seq in entries]
    write_group(outdir, name, entries)

    print("\nTest data written to '%s'" % outdir)
    print("\nExpected behaviour:")
    print("  OG000001–3 : pass all filters → aligned → identity computed")
    print("  OG000004   : removed (no anchor species)")
    print("  OG000005   : removed (only 2 species, min_species=4)")
    print("  OG000006   : drerio seq removed by length filter; "
          "group passes with 4 remaining species")


if __name__ == "__main__":
    import sys
    outdir = sys.argv[1] if len(sys.argv) > 1 else "test_data/input"
    main(outdir)
