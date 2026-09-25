from __future__ import annotations

import os
import re
import html
import math
import logging
from collections import Counter

from Bio import SeqIO
from tqdm import tqdm

from .filtering import display_group_id

logger = logging.getLogger(__name__)

__all__ = [
    "parse_trimal_colnumbering",
    "compute_column_stats",
    "render_alignment_page",
    "build_alignment_pages",
    "load_group_anchor",
]

_GAP_CHARS = frozenset("-.")
# Reference max entropy: a uniform distribution over the 20 standard amino
# acids. Used as a fixed denominator so normalised entropy values stay
# comparable across different columns and alignments, rather than being
# rescaled relative to however many distinct residues happen to appear in
# one particular column.
_MAX_ENTROPY_BITS = math.log2(20)


def parse_trimal_colnumbering(path: str) -> set[int]:
    """Parse a trimAl -colnumbering sidecar file into retained column indices.

    The sidecar file is expected to contain a single line as isolated by
    alignment._extract_colnumbering_line(): trimAl's '#ColumnsMap' label
    (if present) followed by a comma-separated list of 0-based column
    indices from the *original* (pre-trimming) alignment that survive
    trimming, e.g.:

        #ColumnsMap\t0, 1, 2, 3, 4, 5, 6, 7, ...

    Older or different trimAl builds may omit the '#ColumnsMap' label
    entirely and just list the numbers, or use a different label —
    parsing is therefore defensive: any '#ColumnsMap' prefix is stripped
    if present, then every integer remaining in the line is extracted.
    Because the sidecar is expected to contain only this one isolated
    line (not the full alignment), there is no risk of this regex
    picking up unrelated digits such as sequence header gene IDs.

    Args:
        path: Path to a sidecar file written by
              alignment.run_trimal_columns() with colnumbering_dir set.

    Returns:
        Set of 0-based column indices (relative to the original,
        untrimmed alignment) that were retained after trimming. A file
        with no parseable integers (e.g. trimAl produced no output for a
        degenerate alignment) returns an empty set.

    Raises:
        FileNotFoundError: If path does not exist.
    """
    with open(path) as fh:
        content = fh.read()

    match = re.search(r"#ColumnsMap\s*(.*)", content, re.DOTALL)
    if match:
        content = match.group(1)

    return {int(tok) for tok in re.findall(r"\d+", content)}


def compute_column_stats(seqs: list[str]) -> tuple[list[float], list[float]]:
    """Compute per-column identity fraction and normalised Shannon entropy.

    For each alignment column:
      - Identity fraction: count of the most common non-gap residue,
        divided by the total number of non-gap residues in that column.
        A tie for most common still just uses that shared count (which
        residue "wins" the tie doesn't affect the fraction).
      - Normalised entropy: Shannon entropy (base 2) of the non-gap
        residue frequency distribution in that column, divided by
        log2(20) so values fall in [0, 1] and are comparable across
        columns/alignments regardless of how many distinct residues
        happen to appear in any single column. Clipped to 1.0 as a
        defensive measure in the unusual case where more than 20
        distinct symbols appear in one column (e.g. heavy use of
        ambiguity codes), which would otherwise push the raw ratio
        above 1.0.

    Both '-' and '.' are treated as gap characters and excluded from
    both calculations. All other characters (including ambiguity codes
    such as 'X') are treated as ordinary residues — this is a
    deliberate simplification, not an oversight. Residues are
    case-normalised (uppercased) before counting.

    A column with no non-gap residues at all (all sequences have a gap
    there) has no conservation signal to measure. Rather than return
    NaN or None — which would force every caller to special-case
    missing values (e.g. when rendering a bounded color scale) — both
    stats are defined as 0.0 for an all-gap column.

    Args:
        seqs: Aligned sequences (equal length, gap-padded). An empty
              list returns ([], []).

    Returns:
        Tuple of (identity_fractions, entropies), each a list with one
        entry per alignment column, in column order.

    Raises:
        ValueError: If the sequences are not all the same length (i.e.
                    not a valid alignment).
    """
    if not seqs:
        return [], []

    width = len(seqs[0])
    lengths = {len(s) for s in seqs}
    if len(lengths) > 1:
        raise ValueError(
            "All sequences must be the same length to compute column "
            "statistics (this must be an alignment, not raw unaligned "
            "sequences). Got lengths: %s" % sorted(lengths)
        )

    identities: list[float] = []
    entropies: list[float] = []

    for col_idx in range(width):
        residues = [
            s[col_idx].upper() for s in seqs if s[col_idx] not in _GAP_CHARS
        ]
        if not residues:
            identities.append(0.0)
            entropies.append(0.0)
            continue

        counts = Counter(residues)
        n = len(residues)

        identities.append(max(counts.values()) / n)

        entropy_bits = -sum(
            (count / n) * math.log2(count / n) for count in counts.values()
        )
        entropies.append(min(entropy_bits / _MAX_ENTROPY_BITS, 1.0))

    return identities, entropies


# ---------------------------------------------------------------------------
# Alignment page rendering
# ---------------------------------------------------------------------------

# Zappo colour scheme (as used in Jalview): amino acids grouped by
# physicochemical property. A recognised standard rather than an invented
# palette, so rendering stays cross-referenceable against other MSA viewers.
_ZAPPO_GROUPS = [
    ("Aliphatic / hydrophobic", "ILVAM", "#ffafaf"),
    ("Aromatic",                "FWY",   "#ffc800"),
    ("Positive",                "KRH",   "#6464ff"),
    ("Negative",                "DE",    "#ff0000"),
    ("Hydrophilic",             "STNQ",  "#00ff00"),
    ("Conformationally special", "PG",   "#ff00ff"),
    ("Cysteine",                 "C",    "#ffff00"),
]
_OTHER_RESIDUE_COLOR = "#e0e0e0"   # ambiguity codes (X, B, Z, ...), etc.
_GAP_COLOR           = "#ffffff"

# Dark end of the conservation-track gradient. Reuses the report's own
# brand colour so alignment pages read as part of the same product.
_CONSERVATION_DARK = (0x1a, 0x3a, 0x5c)
_CONSERVATION_LIGHT = (0xff, 0xff, 0xff)

_LABEL_TRUNCATE_LEN = 26


def _build_aa_color_map() -> dict[str, tuple[str, str]]:
    mapping: dict[str, tuple[str, str]] = {}
    for group_name, residues, color in _ZAPPO_GROUPS:
        for residue in residues:
            mapping[residue] = (color, group_name)
    return mapping


_AA_COLOR_MAP = _build_aa_color_map()


def _conservation_color(value: float) -> str:
    """Map a conservation value in [0, 1] to a white -> dark-teal gradient.

    1.0 (fully conserved) is darkest; 0.0 is white. Callers wanting the
    entropy track to follow the same "darker = more conserved" direction
    as the identity track should pass (1 - normalised_entropy), since raw
    entropy increases with diversity, not conservation.
    """
    value = max(0.0, min(1.0, value))
    r = round(_CONSERVATION_LIGHT[0] + value * (_CONSERVATION_DARK[0] - _CONSERVATION_LIGHT[0]))
    g = round(_CONSERVATION_LIGHT[1] + value * (_CONSERVATION_DARK[1] - _CONSERVATION_LIGHT[1]))
    b = round(_CONSERVATION_LIGHT[2] + value * (_CONSERVATION_DARK[2] - _CONSERVATION_LIGHT[2]))
    return "rgb(%d,%d,%d)" % (r, g, b)


def _truncate_label(label: str, max_len: int = _LABEL_TRUNCATE_LEN) -> str:
    if len(label) <= max_len:
        return label
    return label[: max_len - 1] + "\u2026"


def _residue_span(residue: str, masked: bool) -> str:
    if residue in _GAP_CHARS:
        color, text = _GAP_COLOR, ""
    else:
        upper = residue.upper()
        color, _group = _AA_COLOR_MAP.get(upper, (_OTHER_RESIDUE_COLOR, "other"))
        text = html.escape(upper)
    cls   = "aa trimmed" if masked else "aa"
    title = ' title="trimmed by trimAl"' if masked else ""
    return '<span class="%s" style="background-color:%s"%s>%s</span>' % (
        cls, color, title, text,
    )


def _track_cell(value: float, tooltip_label: str, masked: bool) -> str:
    color = _conservation_color(value)
    cls   = "track-cell trimmed" if masked else "track-cell"
    return '<span class="%s" style="background-color:%s" title="%s"></span>' % (
        cls, color, html.escape(tooltip_label),
    )


def _render_row(label: str, full_label: str, cells_html: str) -> str:
    return (
        '<div class="row"><span class="row-label" title="%s">%s</span>'
        '<span class="row-body">%s</span></div>'
    ) % (html.escape(full_label), html.escape(label), cells_html)


def _render_block(
    ordered_records: list[tuple[str, str]],
    identities: list[float],
    entropies: list[float],
    retained_columns: set[int] | None,
    start: int,
    end: int,
) -> str:
    rows = []
    for seq_id, seq in ordered_records:
        cells = "".join(
            _residue_span(
                seq[i],
                masked=(retained_columns is not None and i not in retained_columns),
            )
            for i in range(start, end)
        )
        rows.append(_render_row(_truncate_label(seq_id), seq_id, cells))

    identity_cells = "".join(
        _track_cell(
            identities[i], "Identity: %.2f" % identities[i],
            masked=(retained_columns is not None and i not in retained_columns),
        )
        for i in range(start, end)
    )
    rows.append(_render_row("Identity", "Identity (per-column)", identity_cells))

    entropy_cells = "".join(
        _track_cell(
            1.0 - entropies[i], "Entropy: %.2f" % entropies[i],
            masked=(retained_columns is not None and i not in retained_columns),
        )
        for i in range(start, end)
    )
    rows.append(_render_row("Entropy", "Entropy (per-column)", entropy_cells))

    ruler = '<div class="block-ruler">Columns %d\u2013%d</div>' % (start + 1, end)
    return '<div class="aln-block">%s%s</div>' % (ruler, "".join(rows))


def _build_legend_html() -> str:
    items = []
    for group_name, residues, color in _ZAPPO_GROUPS:
        items.append(
            '<span class="legend-item">'
            '<span class="legend-swatch" style="background-color:%s"></span>'
            '%s (%s)</span>' % (color, html.escape(group_name), residues)
        )
    items.append(
        '<span class="legend-item">'
        '<span class="legend-swatch" style="background-color:#ffffff;'
        'border:1px solid #ccc;"></span>Gap</span>'
    )
    items.append(
        '<span class="legend-item">'
        '<span class="legend-swatch trimmed" style="background-color:#999;">'
        '</span>Trimmed by trimAl</span>'
    )
    cons_steps = [0.0, 0.25, 0.5, 0.75, 1.0]
    cons_swatches = "".join(
        '<span class="cons-swatch" style="background-color:%s"></span>'
        % _conservation_color(v)
        for v in cons_steps
    )
    items.append(
        '<span class="legend-item">Conservation: %s Low \u2192 High '
        '(applies to both tracks; entropy track is inverted so darker '
        'always means more conserved)</span>' % cons_swatches
    )
    return '<div class="legend">%s</div>' % "".join(items)


_LEGEND_HTML = _build_legend_html()


_PAGE_TEMPLATE = """\
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8"/>
<title>{title}</title>
<style>
body {{ font-family: Arial, sans-serif; margin: 0; padding: 0;
        background: #f5f5f5; color: #222; }}
header {{ background: #1a3a5c; color: white; padding: 16px 28px; }}
header h1 {{ margin: 0; font-size: 1.3em; }}
main {{ max-width: 100%; margin: 0 auto; padding: 20px 28px 40px; }}
.meta {{ background: white; border-radius: 6px; padding: 14px 20px;
         margin-bottom: 16px; box-shadow: 0 1px 3px rgba(0,0,0,0.08);
         font-size: 13px; }}
.meta dl {{ display: grid; grid-template-columns: max-content 1fr;
            gap: 4px 14px; margin: 0; }}
.meta dt {{ font-weight: bold; color: #1a3a5c; }}
.meta dd {{ margin: 0; }}
.warning {{ background: #fff4e5; border-left: 3px solid #d9822b;
            padding: 8px 14px; margin-bottom: 16px; font-size: 13px; }}
.warning ul {{ margin: 6px 0 0; padding-left: 20px; }}
.legend {{ background: white; border-radius: 6px; padding: 10px 18px;
           margin-bottom: 20px; box-shadow: 0 1px 3px rgba(0,0,0,0.08);
           font-size: 12px; }}
.legend-item {{ margin-right: 18px; white-space: nowrap; }}
.legend-swatch, .cons-swatch {{ display: inline-block; width: 12px; height: 12px;
                                border-radius: 2px; margin-right: 4px;
                                vertical-align: middle; }}
.cons-swatch {{ width: 14px; margin-right: 0; border-radius: 0; }}
.aln-block {{ background: white; border-radius: 6px; padding: 10px 16px 14px;
              margin-bottom: 14px; box-shadow: 0 1px 3px rgba(0,0,0,0.06);
              overflow-x: auto; }}
.block-ruler {{ font-family: 'Courier New', monospace; font-size: 11px;
                color: #888; margin: 2px 0 4px 200px; }}
.row {{ white-space: nowrap; font-size: 0; line-height: 0; border-bottom: 2px solid #ffffff; }}
.row-label {{ display: inline-block; width: 195px; font-family: 'Courier New', monospace;
              font-size: 12px; line-height: 18px; vertical-align: top;
              white-space: nowrap; overflow: hidden; text-overflow: ellipsis;
              padding-right: 6px; box-sizing: border-box; }}
.row-body {{ display: inline-block; vertical-align: top; }}
.aa {{ display: inline-block; vertical-align: top; width: 14px; height: 18px;
       line-height: 18px; text-align: center; font-family: 'Courier New', monospace;
       font-size: 12px; box-sizing: border-box; color: #222; }}
.track-cell {{ display: inline-block; vertical-align: top; width: 14px; height: 9px;
               box-sizing: border-box; }}
.trimmed {{ opacity: 0.35;
            background-image: repeating-linear-gradient(
                45deg, rgba(0,0,0,0.3) 0, rgba(0,0,0,0.3) 2px,
                transparent 2px, transparent 4px); }}
footer {{ text-align: center; padding: 16px; font-size: 0.85em; color: #888; }}
</style>
</head>
<body>
<header><h1>{title}</h1></header>
<main>
<div class="meta">
  <dl>{meta_html}</dl>
</div>
{warning_html}
{legend_html}
{blocks_html}
</main>
<footer>Generated by ENHYDRA</footer>
</body>
</html>"""


def render_alignment_page(
    group_id: str,
    records: list[tuple[str, str]],
    out_path: str,
    anchor_gene_id: str | None = None,
    anchor_species: str | None = None,
    mean_identity: float | None = None,
    retained_columns: set[int] | None = None,
    removed_divergent_sequences: list[tuple[str, float]] | None = None,
    removed_length_filtered_sequences: list[tuple[str, int, float, str]] | None = None,
    block_width: int = 60,
) -> None:
    """Render a static, self-contained HTML page for one alignment.

    Sequences are colour-coded by amino acid property using the Zappo
    scheme, wrapped into fixed-width blocks (block_width columns per
    block — long alignments would otherwise be unreadable as a single
    unwrapped line). Each block carries two per-column conservation
    tracks beneath the sequences: identity fraction and normalised
    Shannon entropy (see compute_column_stats()). Both tracks use the
    same white -> dark-teal colour direction, with darker always meaning
    more conserved (the entropy track is coloured by 1 - entropy so the
    two tracks read the same visual direction despite being different
    statistics); the actual entropy value is still shown in that cell's
    hover tooltip.

    If retained_columns is given, any column *not* in that set is
    rendered with reduced opacity and a hatched overlay (both the
    residue row and both conservation tracks for that column), to show
    at a glance which columns trimAl removed before identity estimation,
    without needing to render a second, fully separate alignment. Column
    indices in retained_columns are 0-based and relative to *this*
    alignment's own column positions (as produced by
    alignment.run_trimal_columns() on this same file) — see
    msa_viewer.parse_trimal_colnumbering().

    Note this masking-overlay approach does NOT extend to the divergent-
    sequence filter (see removed_divergent_sequences below): unlike
    column trimming, removing a divergent sequence is followed by
    realigning the survivors, which changes column positions for every
    remaining sequence. There is no stable column correspondence between
    the pre- and post-removal alignments to overlay a mask against, so
    `records` here is expected to already be the *final*, post-removal
    (and, if applicable, post-realignment) alignment — the same one
    identity was actually computed from — with removed_divergent_sequences
    only used to surface a textual notice of what was removed, not to
    render the removed sequence's row.

    removed_length_filtered_sequences is handled the same textual-notice
    way, for a different reason: length filtering runs before alignment
    even happens, so a sequence removed there was never part of any
    alignment in the first place — there is no row to mask or overlay,
    only a fact to disclose.

    If anchor_species is given, the first record whose ID's species
    field (text before '|') matches it is moved to the top of the
    alignment as a pinned reference row; this only applies to
    single-list runs / list1 in two-list mode, since list2 alignments
    have no anchor sequence to pin. If anchor_species is given but no
    matching record is found, a visible warning banner is added to the
    page (rather than failing silently) and the original record order
    is preserved unchanged.

    This function assumes the parent directory of out_path already
    exists — it is not created here. Callers rendering many pages in a
    batch (see build_alignment_pages()) are expected to create the
    output directory once up front, rather than have every call
    redundantly check for it.

    Args:
        group_id:         Orthogroup identifier as used internally by the
                          pipeline, and as the key build_alignment_pages()
                          uses to look up this group's alignment, table,
                          and colnumbering files. This may carry a
                          pipeline-stage suffix such as '_lengthfilter'
                          (see filtering.display_group_id()) that is not
                          part of the group's real, user-facing name. The
                          page's title and 'Group ID' metadata row show
                          the stripped display form instead — the raw
                          value passed here is never shown directly.
        records:          List of (sequence_id, aligned_sequence) tuples.
                          sequence_id is expected in 'species|gene_id'
                          form (ENHYDRA's convention) so the species
                          field can be used for anchor pinning; the full
                          ID is always shown in the row label's hover
                          tooltip regardless of length. Expected to be
                          the final, post-divergence-filter alignment —
                          see the note above.
        out_path:         Destination HTML file path.
        anchor_gene_id:   Anchor gene ID to display in the metadata
                          panel (display only — does not affect pinning
                          or row lookup).
        anchor_species:   Species ID to pin to the top of the alignment,
                          if present. None disables pinning entirely
                          (e.g. for list2 pages).
        mean_identity:    This group's mean alignment identity score
                          (from group2mean.tsv), shown in the metadata
                          panel. None displays as "N/A".
        retained_columns: Set of 0-based column indices retained after
                          trimAl column trimming, or None if trimming
                          was not applied to this run (no masking).
        removed_divergent_sequences: List of (sequence_id,
                          identity_to_closest_match) tuples for
                          sequences that the divergent-sequence filter
                          (filtering.filter_divergent_sequences())
                          removed from this group before `records` was
                          computed. When non-empty, a warning banner
                          lists each removed sequence and its identity
                          value, and a summary line is added to the
                          metadata panel. None or an empty list omits
                          the banner entirely (e.g. divergence filtering
                          was not enabled for this run, or this
                          particular group was unaffected).
        removed_length_filtered_sequences: List of (sequence_id, length,
                          pct_diff_from_avg, reason) tuples for sequences
                          removed from this group by the length filter
                          (filtering.filter_length()) before alignment
                          — reason is 'removed_below_min' or
                          'removed_above_max'. When non-empty, a separate
                          warning banner lists each removed sequence,
                          its length, and how it compared to the group's
                          mean length at that stage. Unlike the
                          divergent-sequence notice, this never implies
                          realignment — these sequences were excluded
                          before any alignment was attempted. None or an
                          empty list omits the banner entirely.
        block_width:      Number of alignment columns rendered per
                          visual block (default: 60).

    Raises:
        ValueError: If records is empty, or if the sequences it
                    contains are not all the same length (propagated
                    from compute_column_stats()).
    """
    if not records:
        raise ValueError(
            "records must not be empty (group '%s' has no sequences to "
            "render)." % group_id
        )

    # group_id is the internal, file-correlation form (see the Args note
    # above) — disp_id is what actually gets shown to a person reading
    # this page. Computed once here rather than at each of the two call
    # sites below so both stay in sync automatically if either changes.
    disp_id = display_group_id(group_id)

    seqs = [seq for _, seq in records]
    identities, entropies = compute_column_stats(seqs)
    width = len(seqs[0])

    ordered = list(records)
    warnings: list[str] = []
    if anchor_species is not None:
        anchor_idx = next(
            (i for i, (sid, _) in enumerate(ordered)
             if sid.split("|", 1)[0] == anchor_species),
            None,
        )
        if anchor_idx is not None:
            ordered.insert(0, ordered.pop(anchor_idx))
        else:
            logger.warning(
                "Anchor species '%s' not found in alignment for group "
                "'%s' — rendering without pinning.", anchor_species, group_id,
            )
            warnings.append(
                '<div class="warning">Anchor species \u2018%s\u2019 was not '
                'found in this alignment \u2014 rows are shown in their '
                'original order.</div>' % html.escape(anchor_species)
            )

    if removed_length_filtered_sequences:
        _reason_labels = {
            "removed_below_min": "too short",
            "removed_above_max": "too long",
        }
        items = "".join(
            "<li>%s \u2014 length %d (%s; %.2f\u00d7 the group's mean "
            "length at that step)</li>"
            % (html.escape(seq_id), length, _reason_labels.get(reason, reason),
               pct_diff)
            for seq_id, length, pct_diff, reason in removed_length_filtered_sequences
        )
        warnings.append(
            '<div class="warning">%d sequence(s) were removed from this '
            'group by the length filter before alignment (their length '
            'deviated too far from the group mean at that stage), and '
            'never appear in the alignment shown below:'
            '<ul>%s</ul></div>'
            % (len(removed_length_filtered_sequences), items)
        )

    if removed_divergent_sequences:
        items = "".join(
            "<li>%s (identity to closest match: %.4f)</li>"
            % (html.escape(seq_id), identity_val)
            for seq_id, identity_val in removed_divergent_sequences
        )
        warnings.append(
            '<div class="warning">%d sequence(s) were removed from this '
            'group by the divergence filter (identity to closest match '
            'fell too far below the group mean), and the remaining '
            'sequences shown below were realigned without them:'
            '<ul>%s</ul></div>'
            % (len(removed_divergent_sequences), items)
        )
    warning_html = "".join(warnings)

    n_masked = None
    if retained_columns is not None:
        n_masked = width - len([i for i in range(width) if i in retained_columns])

    meta_rows = [
        ("Group ID", html.escape(disp_id)),
        ("Sequences", str(len(records))),
        ("Alignment length", "%d columns" % width),
    ]
    if anchor_gene_id is not None:
        meta_rows.append(("Anchor gene", html.escape(anchor_gene_id)))
    if anchor_species is not None:
        meta_rows.append(("Anchor species", html.escape(anchor_species)))
    meta_rows.append((
        "Mean identity",
        "%.4f" % mean_identity if mean_identity is not None else "N/A",
    ))
    if n_masked is not None:
        meta_rows.append((
            "Trimmed columns",
            "%d of %d removed by trimAl (shown masked below)" % (n_masked, width),
        ))
    if removed_length_filtered_sequences:
        meta_rows.append((
            "Length-filtered sequences",
            "%d (see notice above; excluded before alignment)"
            % len(removed_length_filtered_sequences),
        ))
    if removed_divergent_sequences:
        meta_rows.append((
            "Divergent sequences removed",
            "%d (see notice above; alignment shown is post-removal, "
            "realigned)" % len(removed_divergent_sequences),
        ))
    meta_html = "".join(
        "<dt>%s</dt><dd>%s</dd>" % (k, v) for k, v in meta_rows
    )

    blocks_html = "".join(
        _render_block(ordered, identities, entropies, retained_columns,
                      block_start, min(block_start + block_width, width))
        for block_start in range(0, width, block_width)
    )

    html_content = _PAGE_TEMPLATE.format(
        title="Alignment: %s" % html.escape(disp_id),
        meta_html=meta_html,
        warning_html=warning_html,
        legend_html=_LEGEND_HTML,
        blocks_html=blocks_html,
    )
    with open(out_path, "w", encoding="utf-8") as fh:
        fh.write(html_content)


# ---------------------------------------------------------------------------
# Batch generation
# ---------------------------------------------------------------------------

def _load_group_mean(tables_dir: str) -> dict[str, float]:
    path = os.path.join(tables_dir, "group2mean.tsv")
    result: dict[str, float] = {}
    if not os.path.isfile(path):
        return result
    with open(path) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                continue
            try:
                result[fields[0]] = float(fields[1])
            except ValueError:
                continue
    return result


def load_group_anchor(tables_dir: str) -> dict[str, str]:
    """Load a list's group2anchor.tsv as a {group_id: anchor_gene_id} dict.

    Public (unlike _load_group_mean) because cli.py needs this mapping
    directly when building an anchor_gene_lookup override for list2's
    alignment pages in two-list mode — list2's own group2anchor.tsv is
    expected to be empty by design (list2 never contains the anchor
    species), so list1's mapping is passed through instead.

    Returns an empty dict if the file does not exist.
    """
    return _load_group_anchor(tables_dir)


def _load_group_anchor(tables_dir: str) -> dict[str, str]:
    path = os.path.join(tables_dir, "group2anchor.tsv")
    result: dict[str, str] = {}
    if not os.path.isfile(path):
        return result
    with open(path) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                continue
            result[fields[0]] = fields[1]
    return result


def _load_divergence_removed(
    path: str | None,
) -> dict[str, list[tuple[str, float]]]:
    """Load a divergence filter's drop_reasons.tsv into per-group removal lists.

    Only rows with reason == 'removed_divergent_sequence' are included.
    The file's other possible reason, 'below_min_species_after_
    divergence_filter', describes a group dropped entirely — such a group
    never has a rendered alignment page in the first place (it has no
    surviving alignment file at all), so there is nothing for that reason
    to annotate here.

    Args:
        path: Path to filtering.filter_divergent_sequences()'s
             drop_reasons.tsv, or None if divergence filtering was not
             enabled for this run.

    Returns:
        Dict mapping group_id to a list of (sequence_id,
        identity_to_closest_match) tuples. Empty dict if path is None,
        the file does not exist, or the file's header is missing any of
        the expected columns (defensive against a mismatched or
        corrupted stats file — silently rendering pages without the
        notice is preferable to crashing report generation over it).
    """
    if not path or not os.path.isfile(path):
        return {}

    result: dict[str, list[tuple[str, float]]] = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            reason_idx   = header.index("reason")
            group_idx    = header.index("group_id")
            seq_idx      = header.index("sequence_id")
            identity_idx = header.index("identity_to_closest")
        except ValueError:
            logger.warning(
                "Divergence drop-reasons file has an unexpected header — "
                "skipping removal notices in alignment pages: %s", path,
            )
            return {}

        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) <= max(reason_idx, group_idx, seq_idx, identity_idx):
                continue
            if fields[reason_idx] != "removed_divergent_sequence":
                continue
            try:
                identity_val = float(fields[identity_idx])
            except ValueError:
                continue
            result.setdefault(fields[group_idx], []).append(
                (fields[seq_idx], identity_val)
            )

    return result


_LENGTH_FILTER_REMOVAL_REASONS = frozenset({
    "removed_below_min", "removed_above_max",
})


def _load_length_filter_removed(
    path: str | None,
) -> dict[str, list[tuple[str, int, float, str]]]:
    """Load the length filter's drop_reasons.tsv into per-group removal lists.

    Unlike the divergence filter, length filtering is not optional — it
    always runs, and aggregate_length_filter_stats() always writes this
    file unconditionally (even if empty of rows) — so, unlike
    _load_divergence_removed(), callers do not need to treat a missing
    path as "step didn't run"; it is simply "nothing to report" either
    way. path is still accepted as None/missing defensively (e.g. a
    directory layout from an older run), in which case this returns an
    empty dict exactly as it would for a file with zero drop rows.

    Only rows whose reason is one of 'removed_below_min' or
    'removed_above_max' are included (filtering.filter_length()'s other
    possible reasons — 'empty_file', 'single_sequence' — describe a group
    skipped entirely before any alignment was attempted, and are recorded
    in a separate skipped_groups.tsv, not this file; such a group has no
    alignment page to annotate in the first place).

    Args:
        path: Path to length_filter_stats/drop_reasons.tsv (as written by
             filtering.aggregate_length_filter_stats()), or None.

    Returns:
        Dict mapping group_id to a list of (sequence_id, length,
        pct_diff_from_avg, reason) tuples, where pct_diff_from_avg
        follows the same value/mean ratio convention used throughout
        this codebase's other drop-reasons files (e.g. the divergence
        filter's own pct_diff_from_avg column) — a ratio around 1.0, not
        a percentage. Empty dict if path is None, the file does not
        exist, or its header is missing any expected column.
    """
    if not path or not os.path.isfile(path):
        return {}

    result: dict[str, list[tuple[str, int, float, str]]] = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            group_idx    = header.index("group_id")
            seq_idx      = header.index("sequence_id")
            length_idx   = header.index("length")
            pct_diff_idx = header.index("pct_diff_from_avg")
            reason_idx   = header.index("reason")
        except ValueError:
            logger.warning(
                "Length filter drop-reasons file has an unexpected header "
                "— skipping removal notices in alignment pages: %s", path,
            )
            return {}

        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) <= max(group_idx, seq_idx, length_idx,
                                  pct_diff_idx, reason_idx):
                continue
            reason = fields[reason_idx]
            if reason not in _LENGTH_FILTER_REMOVAL_REASONS:
                continue
            try:
                length   = int(float(fields[length_idx]))
                pct_diff = float(fields[pct_diff_idx])
            except ValueError:
                continue
            result.setdefault(fields[group_idx], []).append(
                (fields[seq_idx], length, pct_diff, reason)
            )

    return result


def build_alignment_pages(
    alignment_dir: str,
    tables_dir: str,
    outdir: str,
    group_ids: list[str] | None = None,
    anchor_species: str | None = None,
    anchor_gene_lookup: dict[str, str] | None = None,
    colnumbering_dir: str | None = None,
    divergence_drop_reasons_path: str | None = None,
    length_filter_drop_reasons_path: str | None = None,
    block_width: int = 60,
    show_progress: bool = False,
) -> dict[str, str]:
    """Render one alignment page per group, returning {group_id: page path}.

    Alignment filenames are reconstructed as '<group_id>.aln', matching
    how make_tables() itself derives group_id (everything before the
    first '.' in a group's identity-report filename) — this holds for
    every alignment produced by run_aligner(), since group filenames
    never contain literal periods. This also holds for groups realigned
    by the divergent-sequence filter (filtering.filter_divergent_
    sequences()): its realignment input files are named with the bare
    group_id (no extension) precisely so that re-running them through
    run_aligner() reproduces this same '<group_id>.aln' naming with no
    separate rename step.

    Args:
        alignment_dir:      Directory of alignment files for this list.
                            Callers should point this at whichever
                            directory holds the *final* alignment used
                            for identity estimation — i.e. if the
                            divergent-sequence filter was enabled for
                            this run, this should be that filter's own
                            output directory (filtering.
                            filter_divergent_sequences()'s filtered_dir),
                            not the raw pre-filter alignment_dir, since
                            rendering the raw alignment would show a
                            since-removed sequence with no indication it
                            was excluded from the actual score, and (for
                            groups that were realigned) would show
                            entirely different column positions than the
                            alignment identity was actually computed
                            from. If divergence filtering was not
                            enabled, this is simply the aligner's own
                            output directory, as before.
        tables_dir:          This list's own tables/ directory, used to
                            load group2mean.tsv (identity scores) and,
                            unless anchor_gene_lookup overrides it,
                            group2anchor.tsv (anchor gene IDs).
        outdir:              Directory where rendered pages are written,
                            one per group as '<group_id>.html'. Created
                            if it does not already exist (unlike
                            render_alignment_page(), which assumes its
                            own parent directory is already present).
        group_ids:           Explicit set of group IDs to render. If
                            None, defaults to every group_id found in
                            tables_dir/group2anchor.tsv — i.e. every
                            group that reached anchor2mean.tsv in
                            single-list mode. In two-list mode, callers
                            should instead pass the group_id column of
                            differential_scores.tsv, since that is the
                            actual set of groups that fed GSEA for that
                            run (see msa_viewer module docs / project
                            notes on 'bucket 3').
        anchor_species:      Species ID to pin at the top of every
                            rendered alignment. None disables pinning
                            (e.g. for list2 pages, which never contain
                            the anchor species).
        anchor_gene_lookup:  Optional {group_id: anchor_gene_id} mapping
                            that takes priority over tables_dir's own
                            group2anchor.tsv. Intended for rendering a
                            list2 page while still showing the anchor
                            gene ID that was resolved from list1 for
                            that same group, since list2's own
                            group2anchor.tsv is expected to be empty by
                            design (list2 never contains the anchor
                            species). Single-list callers should omit
                            this and let it fall back to tables_dir.
        colnumbering_dir:    Directory of trimAl -colnumbering sidecar
                            files for this list (as produced by
                            alignment.run_trimal_columns()), or None if
                            column trimming was not used for this run.
                            When given, a group without a matching
                            sidecar file (e.g. trimAl produced no
                            output for a degenerate alignment) is
                            rendered without masking rather than being
                            skipped. Column indices are relative to
                            whatever alignment_dir points at above —
                            i.e. the post-divergence-filter alignment,
                            if that filter ran, since column trimming
                            in the pipeline always runs after it.
        divergence_drop_reasons_path: Path to the divergent-sequence
                            filter's drop_reasons.tsv for this list (see
                            filtering.filter_divergent_sequences()), or
                            None if that filter was not enabled for this
                            run. When given, groups with one or more
                            removed sequences get a visible notice on
                            their rendered page (see
                            render_alignment_page()'s
                            removed_divergent_sequences parameter) rather
                            than silently showing only the corrected
                            alignment with no explanation of the
                            discrepancy from what the person might expect
                            to see (e.g. the original group_filter
                            sequence count).
        length_filter_drop_reasons_path: Path to the length filter's own
                            aggregated drop_reasons.tsv (see filtering.
                            aggregate_length_filter_stats()) for this
                            list. Unlike divergence_drop_reasons_path,
                            this step always runs, so this path is
                            expected to be given on every call in
                            practice — passed as None only degrades
                            gracefully (no notice shown), it does not
                            indicate anything about whether length
                            filtering itself ran. When given, groups with
                            one or more sequences removed at that stage
                            get a separate visible notice (see
                            render_alignment_page()'s
                            removed_length_filtered_sequences parameter)
                            distinct from the divergence-filter notice,
                            since these sequences were excluded before
                            alignment even began and never implied any
                            realignment.
        block_width:         Passed through to render_alignment_page().
        show_progress:       Show a tqdm progress bar.

    Returns:
        Dict mapping group_id to the absolute path of its rendered HTML
        page. Only successfully rendered groups are included; groups
        requested (explicitly or via group2anchor.tsv) but missing an
        alignment file are logged as a warning and skipped, not raised
        as an error, since a partially-resumed run or an upstream
        filtering edge case can legitimately produce this.
    """
    os.makedirs(outdir, exist_ok=True)

    if group_ids is None:
        group_ids = sorted(_load_group_anchor(tables_dir).keys())
    else:
        group_ids = sorted(set(group_ids))

    mean_lookup = _load_group_mean(tables_dir)
    anchor_lookup = anchor_gene_lookup
    if anchor_lookup is None:
        anchor_lookup = _load_group_anchor(tables_dir)
    divergence_removed = _load_divergence_removed(divergence_drop_reasons_path)
    length_filter_removed = _load_length_filter_removed(
        length_filter_drop_reasons_path
    )

    pages: dict[str, str] = {}
    n_missing_alignment = 0

    for group_id in tqdm(group_ids, desc="  alignments", unit="group",
                         leave=False, disable=not show_progress):
        aln_path = os.path.join(alignment_dir, group_id + ".aln")
        if not os.path.isfile(aln_path):
            logger.warning(
                "No alignment file found for group '%s' (expected: %s) — "
                "skipping alignment page.", group_id, aln_path,
            )
            n_missing_alignment += 1
            continue

        records = [(r.id, str(r.seq)) for r in SeqIO.parse(aln_path, "fasta")]

        retained_columns = None
        if colnumbering_dir is not None:
            colnum_path = os.path.join(
                colnumbering_dir, group_id + ".aln.colnumbering"
            )
            if os.path.isfile(colnum_path):
                retained_columns = parse_trimal_colnumbering(colnum_path)

        out_path = os.path.abspath(os.path.join(outdir, group_id + ".html"))
        render_alignment_page(
            group_id=group_id,
            records=records,
            out_path=out_path,
            anchor_gene_id=anchor_lookup.get(group_id),
            anchor_species=anchor_species,
            mean_identity=mean_lookup.get(group_id),
            retained_columns=retained_columns,
            removed_divergent_sequences=divergence_removed.get(group_id),
            removed_length_filtered_sequences=length_filter_removed.get(group_id),
            block_width=block_width,
        )
        pages[group_id] = out_path

    logger.info(
        "Alignment pages written: %d. Groups skipped (no alignment found): %d. Path: %s",
        len(pages), n_missing_alignment, outdir,
    )
    return pages
