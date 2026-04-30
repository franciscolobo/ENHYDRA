"""Quick script to remove pseudogene sequences from existing proteome FASTAs.
Filters out any sequence whose header contains [pseudo=true].
Overwrites files in place.

Usage:
    python tests/remove_pseudogenes.py <proteomes_dir>
"""

import os
import sys

proteomes_dir = sys.argv[1]

for filename in os.listdir(proteomes_dir):
    if not filename.endswith(".fa"):
        continue
    path = os.path.join(proteomes_dir, filename)
    kept, removed = [], 0
    current_header, current_seq = None, []

    def flush():
        global removed
        if current_header is None:
            return
        if "[pseudo=true]" in current_header:
            removed += 1
        else:
            kept.append(current_header)
            kept.extend(current_seq)

    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                flush()
                current_header = line
                current_seq = []
            else:
                current_seq.append(line)
    flush()

    with open(path, "w") as fh:
        fh.write("\n".join(kept) + "\n")

    if removed:
        print("%s: removed %d pseudogene(s), kept %d sequences" %
              (filename, removed, len([l for l in kept if l.startswith(">")])))
