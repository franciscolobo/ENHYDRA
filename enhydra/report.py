from __future__ import annotations
from .io import parse_obo_names as _parse_obo_names

import json
import os
import base64
import logging
import urllib.request
import ssl
import pandas as pd

logger = logging.getLogger(__name__)

_DATATABLES_JS_URL  = "https://cdn.datatables.net/1.13.6/js/jquery.dataTables.min.js"
_JQUERY_URL         = "https://code.jquery.com/jquery-3.7.0.min.js"
_DATATABLES_CSS_URL = "https://cdn.datatables.net/1.13.6/css/jquery.dataTables.min.css"

METRIC_LABELS = {
    "identity": "Identity",
    "zscore":   "Z-score",
    "rank":     "Rank",
}

_METRIC_DESCS = {
    "identity": (
        "Raw mean pairwise sequence identity averaged across all species pairs "
        "in each orthogroup. No normalisation is applied."
    ),
    "zscore": (
        "Z-score normalised identity: each group\u2019s score is expressed in "
        "standard deviations from the mean across all groups. "
        "Positive\u00a0=\u00a0more conserved than average; "
        "negative\u00a0=\u00a0faster evolving."
    ),
    "rank": (
        "Normalised rank: groups are ranked by identity and scores are divided "
        "by N so that the most-conserved group receives a score of 1.0 and the "
        "fastest-evolving group receives 1/N."
    ),
}

_TEMPLATE = """\
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8"/>
<meta name="viewport" content="width=device-width, initial-scale=1.0"/>
<title>{title}</title>
<style>
{dt_css}
body {{ font-family: Arial, sans-serif; margin: 0; padding: 0;
        background: #f5f5f5; color: #222; }}
header {{ background: #1a3a5c; color: white; padding: 24px 40px; }}
header h1 {{ margin: 0; font-size: 1.8em; }}
header p  {{ margin: 4px 0 0; font-size: 0.95em; opacity: 0.85; }}
main {{ max-width: 1300px; margin: 32px auto; padding: 0 24px; }}
section {{ background: white; border-radius: 8px;
           box-shadow: 0 1px 4px rgba(0,0,0,0.1);
           padding: 28px 32px; margin-bottom: 28px; }}
h2 {{ margin-top: 0; color: #1a3a5c; border-bottom: 2px solid #e0e0e0;
      padding-bottom: 8px; }}
.plot-block {{ margin: 20px 0; text-align: center; }}
.plot-block img {{ max-width: 100%; border: 1px solid #e0e0e0; border-radius: 4px; }}
.plot-caption {{ font-size: 0.9em; color: #555; margin-bottom: 6px; }}
.plot-grid {{ display: grid; grid-template-columns: 1fr 1fr; gap: 20px; }}
tr.sig-row {{ background-color: #eaf3fb !important; font-weight: bold; }}
a.go-link {{ color: #1a3a5c; text-decoration: underline dotted; cursor: pointer; }}
.col-tip {{ display: inline-block; width: 14px; height: 14px; line-height: 14px;
            font-size: 10px; text-align: center; border-radius: 50%;
            background: #aaa; color: white; cursor: help; margin-left: 3px;
            position: relative; }}
.col-tip .tip-text {{ display: none; position: absolute; bottom: 120%; left: 50%;
                      transform: translateX(-50%); background: #333; color: #fff;
                      padding: 6px 10px; border-radius: 4px; font-size: 11px;
                      white-space: normal; width: 220px; z-index: 999;
                      font-weight: normal; line-height: 1.4; }}
.col-tip:hover .tip-text {{ display: block; }}
thead tr.filter-row th input {{
    width: 100%; box-sizing: border-box; font-size: 11px;
    padding: 3px; border: 1px solid #ccc; border-radius: 3px; }}
thead tr.filter-row th {{ padding: 4px 8px; }}
.svg-tooltip {{ position: fixed; background: #333; color: #fff;
                padding: 8px 12px; border-radius: 4px; font-size: 12px;
                max-width: 320px; line-height: 1.5; pointer-events: none;
                z-index: 2000; display: none; white-space: normal; }}
#modal-overlay {{ display: none; position: fixed; top: 0; left: 0;
                  width: 100%; height: 100%; background: rgba(0,0,0,0.6);
                  z-index: 1000; justify-content: center; align-items: center; }}
#modal-overlay.active {{ display: flex; }}
#modal-box {{ background: white; border-radius: 8px; padding: 24px;
              max-width: 700px; width: 90%; position: relative; }}
#modal-title {{ font-size: 1.1em; font-weight: bold; color: #1a3a5c;
                margin-bottom: 12px; }}
#modal-img {{ width: 100%; border: 1px solid #e0e0e0; border-radius: 4px; }}
#modal-close {{ position: absolute; top: 12px; right: 16px; font-size: 1.4em;
                cursor: pointer; color: #555; background: none; border: none; }}
footer {{ text-align: center; padding: 20px; font-size: 0.85em; color: #888; }}
</style>
</head>
<body>
<header>
  <h1>{title}</h1>
  <p>ENHYDRA &mdash; Gene Set Enrichment Analysis for evolutionary genomics</p>
</header>

<div id="svg-tooltip" class="svg-tooltip"></div>

<div id="modal-overlay">
  <div id="modal-box">
    <button id="modal-close" title="Close">&times;</button>
    <div id="modal-title"></div>
    <img id="modal-img" src="" alt="Enrichment plot"/>
  </div>
</div>

<main>
<section>
  <h2>Plots</h2>
  <div class="plot-grid">
{plots_html}
  </div>
</section>

<section>
  <h2>Enrichment results</h2>
  <p>Significant gene sets are highlighted in blue.
     Click a GO ID to view its enrichment plot.
     Use the filter boxes below each column header to filter by that column.</p>
{table_html}
</section>
</main>
<footer>Generated by ENHYDRA</footer>

<script>{jquery_js}</script>
<script>{dt_js}</script>
<script>
{plot_data_js}
$(document).ready(function() {{

    // SVG bar tooltips
    var tip = document.getElementById('svg-tooltip');
    document.querySelectorAll('[data-tip]').forEach(function(el) {{
        el.addEventListener('mousemove', function(e) {{
            tip.innerHTML = this.getAttribute('data-tip');
            tip.style.display = 'block';
            tip.style.left = (e.clientX + 15) + 'px';
            tip.style.top  = (e.clientY + 15) + 'px';
        }});
        el.addEventListener('mouseleave', function() {{
            tip.style.display = 'none';
        }});
    }});

    // Per-column filters
    var numericCols = {numeric_col_indices};
    var colFilters = {{}};

    $.fn.dataTable.ext.search.push(function(settings, data) {{
        for (var i in colFilters) {{
            var f = colFilters[i];
            if (f.text !== undefined) {{
                if (data[i].toLowerCase().indexOf(f.text) === -1) return false;
            }} else {{
                var cell = parseFloat(data[i]);
                if (isNaN(cell)) return false;
                if (f.op === '<'  && !(cell <  f.num)) return false;
                if (f.op === '<=' && !(cell <= f.num)) return false;
                if (f.op === '>'  && !(cell >  f.num)) return false;
                if (f.op === '>=' && !(cell >= f.num)) return false;
                if ((f.op === '=' || f.op === '==') && !(cell === f.num)) return false;
                if (f.op === '!=' && !(cell !== f.num)) return false;
            }}
        }}
        return true;
    }});

    var table = $('#results-table').DataTable({{
        pageLength: 25,
        orderCellsTop: true,
        order: [[4, 'asc']],
        columnDefs: [{{ targets: numericCols, type: 'num' }}]
    }});

    $('#results-table thead tr.filter-row th').each(function(i) {{
        var isNumeric = numericCols.indexOf(i) !== -1;
        var placeholder = isNumeric ? "e.g. < 0.05" : "Filter...";
        var input = $('<input type="text" placeholder="' + placeholder + '"/>');
        $(this).html(input);
        input.on('keyup change', function() {{
            var val = this.value.trim();
            if (val === '') {{
                delete colFilters[i];
            }} else if (isNumeric) {{
                var m = val.match(/^([<>=!]=?)\s*([\d.eE+\-]+)$/);
                if (m) colFilters[i] = {{ op: m[1], num: parseFloat(m[2]) }};
                else delete colFilters[i];
            }} else {{
                colFilters[i] = {{ text: val.toLowerCase() }};
            }}
            table.draw();
        }});
    }});

    // GO enrichment plot modal
    $(document).on('click', '.go-link', function(e) {{
        e.preventDefault();
        var goId = $(this).data('goid');
        var uri  = enrichmentPlots[goId];
        if (uri) {{
            $('#modal-title').text(goId);
            $('#modal-img').attr('src', uri);
            $('#modal-overlay').addClass('active');
        }}
    }});
    $('#modal-close, #modal-overlay').on('click', function(e) {{
        if (e.target === this) $('#modal-overlay').removeClass('active');
    }});

}});
</script>
</body>
</html>"""

_MULTI_TEMPLATE = """\
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8"/>
<meta name="viewport" content="width=device-width, initial-scale=1.0"/>
<title>{title}</title>
<style>
{dt_css}
body {{ font-family: Arial, sans-serif; margin: 0; padding: 0;
        background: #f5f5f5; color: #222; }}
header {{ background: #1a3a5c; color: white; padding: 24px 40px; }}
header h1 {{ margin: 0; font-size: 1.8em; }}
header p  {{ margin: 4px 0 0; font-size: 0.95em; opacity: 0.85; }}
main {{ max-width: 1300px; margin: 32px auto; padding: 0 24px; }}
section {{ background: white; border-radius: 8px;
           box-shadow: 0 1px 4px rgba(0,0,0,0.1);
           padding: 28px 32px; margin-bottom: 28px; }}
h2 {{ margin-top: 0; color: #1a3a5c; border-bottom: 2px solid #e0e0e0;
      padding-bottom: 8px; }}
h3 {{ color: #2c5282; margin: 24px 0 12px; }}
.tab-nav {{ display: flex; gap: 0; border-bottom: 3px solid #1a3a5c;
            margin-bottom: 28px; flex-wrap: wrap; }}
.tab-btn {{ padding: 11px 32px; border: none; border-radius: 6px 6px 0 0;
            background: #e2eaf3; cursor: pointer; font-size: 14px;
            font-weight: 600; color: #555; margin-right: 3px;
            transition: background .15s, color .15s; }}
.tab-btn:hover:not(.active) {{ background: #c8d8ec; color: #1a3a5c; }}
.tab-btn.active {{ background: #1a3a5c; color: white; }}
.tab-panel {{ display: none; }}
.tab-panel.active {{ display: block; }}
.metric-desc {{ font-size: 0.9em; color: #444; margin: 0 0 20px;
                padding: 10px 14px; background: #f0f5fa;
                border-left: 3px solid #1a3a5c; border-radius: 0 4px 4px 0; }}
.plot-block {{ margin: 20px 0; text-align: center; }}
.plot-block img {{ max-width: 100%; border: 1px solid #e0e0e0; border-radius: 4px; }}
.plot-caption {{ font-size: 0.9em; color: #555; margin-bottom: 6px; }}
.plot-grid {{ display: grid; grid-template-columns: 1fr 1fr; gap: 20px; }}
tr.sig-row {{ background-color: #eaf3fb !important; font-weight: bold; }}
a.go-link {{ color: #1a3a5c; text-decoration: underline dotted; cursor: pointer; }}
.col-tip {{ display: inline-block; width: 14px; height: 14px; line-height: 14px;
            font-size: 10px; text-align: center; border-radius: 50%;
            background: #aaa; color: white; cursor: help; margin-left: 3px;
            position: relative; }}
.col-tip .tip-text {{ display: none; position: absolute; bottom: 120%; left: 50%;
                      transform: translateX(-50%); background: #333; color: #fff;
                      padding: 6px 10px; border-radius: 4px; font-size: 11px;
                      white-space: normal; width: 220px; z-index: 999;
                      font-weight: normal; line-height: 1.4; }}
.col-tip:hover .tip-text {{ display: block; }}
thead tr.filter-row th input {{
    width: 100%; box-sizing: border-box; font-size: 11px;
    padding: 3px; border: 1px solid #ccc; border-radius: 3px; }}
thead tr.filter-row th {{ padding: 4px 8px; }}
.svg-tooltip {{ position: fixed; background: #333; color: #fff;
                padding: 8px 12px; border-radius: 4px; font-size: 12px;
                max-width: 320px; line-height: 1.5; pointer-events: none;
                z-index: 2000; display: none; white-space: normal; }}
#modal-overlay {{ display: none; position: fixed; top: 0; left: 0;
                  width: 100%; height: 100%; background: rgba(0,0,0,0.6);
                  z-index: 1000; justify-content: center; align-items: center; }}
#modal-overlay.active {{ display: flex; }}
#modal-box {{ background: white; border-radius: 8px; padding: 24px;
              max-width: 700px; width: 90%; position: relative; }}
#modal-title {{ font-size: 1.1em; font-weight: bold; color: #1a3a5c;
                margin-bottom: 12px; }}
#modal-img {{ width: 100%; border: 1px solid #e0e0e0; border-radius: 4px; }}
#modal-close {{ position: absolute; top: 12px; right: 16px; font-size: 1.4em;
                cursor: pointer; color: #555; background: none; border: none; }}
footer {{ text-align: center; padding: 20px; font-size: 0.85em; color: #888; }}
</style>
</head>
<body>
<header>
  <h1>{title}</h1>
  <p>ENHYDRA &mdash; Gene Set Enrichment Analysis for evolutionary genomics</p>
</header>
<div id="svg-tooltip" class="svg-tooltip"></div>
<div id="modal-overlay">
  <div id="modal-box">
    <button id="modal-close" title="Close">&times;</button>
    <div id="modal-title"></div>
    <img id="modal-img" src="" alt="Enrichment plot"/>
  </div>
</div>
<main>
<section>
  <h2>Results by ranking metric</h2>
  <nav class="tab-nav" role="tablist">
{tab_buttons}
  </nav>
{tab_panels}
</section>
</main>
<footer>Generated by ENHYDRA</footer>
<script>{jquery_js}</script>
<script>{dt_js}</script>
<script>
// Per-metric enrichment plots: {{ metric: {{ goId: base64uri }} }}
var enrichmentPlotsMap = {enrichment_plots_map};
// Per-metric numeric column indices for DataTables type detection
var numericColsMap = {numeric_cols_map};
var dtInstances   = {{}};
var colFiltersMap  = {{}};

$.fn.dataTable.ext.search.push(function(settings, data) {{
    var metric = settings.nTable.id.replace('results-table-', '');
    var cf = colFiltersMap[metric] || {{}};
    for (var i in cf) {{
        var f = cf[i];
        if (f.text !== undefined) {{
            if (data[i].toLowerCase().indexOf(f.text) === -1) return false;
        }} else {{
            var v = parseFloat(data[i]);
            if (isNaN(v)) return false;
            if (f.op === '<'  && !(v <  f.num)) return false;
            if (f.op === '<=' && !(v <= f.num)) return false;
            if (f.op === '>'  && !(v >  f.num)) return false;
            if (f.op === '>=' && !(v >= f.num)) return false;
            if ((f.op === '=' || f.op === '==') && v !== f.num) return false;
            if (f.op === '!=' && v === f.num)   return false;
        }}
    }}
    return true;
}});

function initTable(metric) {{
    if (dtInstances[metric]) return;
    colFiltersMap[metric] = {{}};
    var numericCols = numericColsMap[metric] || [];
    var dt = $('#results-table-' + metric).DataTable({{
        pageLength: 25,
        orderCellsTop: true,
        order: [[4, 'asc']],
        columnDefs: [{{ targets: numericCols, type: 'num' }}],
    }});
    $('#results-table-' + metric + ' thead tr.filter-row th').each(function(i) {{
        var isNum = numericCols.indexOf(i) !== -1;
        var inp   = $('<input type="text" placeholder="' +
                      (isNum ? 'e.g. < 0.05' : 'Filter...') + '"/>');
        $(this).html(inp);
        inp.on('keyup change', (function(col) {{
            return function() {{
                var val = $.trim(this.value);
                if (!val) {{
                    delete colFiltersMap[metric][col];
                }} else if (isNum) {{
                    var m = val.match(/^([<>=!]=?)\s*([\d.eE+\-]+)$/);
                    if (m) colFiltersMap[metric][col] = {{ op: m[1], num: parseFloat(m[2]) }};
                    else   delete colFiltersMap[metric][col];
                }} else {{
                    colFiltersMap[metric][col] = {{ text: val.toLowerCase() }};
                }}
                dt.draw();
            }};
        }})(i));
    }});
    dtInstances[metric] = dt;
}}

$(document).ready(function() {{
    var svgTip = document.getElementById('svg-tooltip');
    document.querySelectorAll('[data-tip]').forEach(function(el) {{
        el.addEventListener('mousemove', function(e) {{
            svgTip.innerHTML     = this.getAttribute('data-tip');
            svgTip.style.display = 'block';
            svgTip.style.left    = (e.clientX + 15) + 'px';
            svgTip.style.top     = (e.clientY + 15) + 'px';
        }});
        el.addEventListener('mouseleave', function() {{
            svgTip.style.display = 'none';
        }});
    }});

    document.querySelectorAll('.tab-btn').forEach(function(btn) {{
        btn.addEventListener('click', function() {{
            var metric = this.dataset.metric;
            document.querySelectorAll('.tab-btn').forEach(function(b) {{
                b.classList.remove('active');
            }});
            document.querySelectorAll('.tab-panel').forEach(function(p) {{
                p.classList.remove('active');
            }});
            this.classList.add('active');
            document.getElementById('tab-' + metric).classList.add('active');
            initTable(metric);
        }});
    }});

    $(document).on('click', '.go-link', function(e) {{
        e.preventDefault();
        var goId   = $(this).data('goid');
        var metric = $(this).data('metric');
        var plots  = metric ? enrichmentPlotsMap[metric] : enrichmentPlots;
        var uri    = plots ? plots[goId] : undefined;
        if (uri) {{
            $('#modal-title').text(goId);
            $('#modal-img').attr('src', uri);
            $('#modal-overlay').addClass('active');
        }}
    }});
    $('#modal-close, #modal-overlay').on('click', function(e) {{
        if (e.target === this) $('#modal-overlay').removeClass('active');
    }});

    var firstBtn = document.querySelector('.tab-btn');
    if (firstBtn) firstBtn.click();
}});
</script>
</body>
</html>"""


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _gmt_term_names(gmt_path: str | None) -> dict[str, str]:
    """Parse term_id → term_name from the second column of a GMT file.

    The GMT format is:
        term_id <TAB> term_name <TAB> gene1 <TAB> gene2 ...

    Returns an empty dict if the path is None, missing, or unreadable.
    Used as a fallback source of term names when no OBO file is available.
    """
    if not gmt_path or not os.path.isfile(gmt_path):
        return {}
    names: dict[str, str] = {}
    with open(gmt_path) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 2 and fields[0] and fields[1]:
                names[fields[0]] = fields[1]
    return names


def _find_gmt_in_dir(directory: str) -> str | None:
    """Return the path to the first .gmt file found in directory, or None.

    Used to locate the GMT written by run_gsea when --organism is used
    (saved as results_dir/gprofiler_<organism>.gmt).
    """
    if not os.path.isdir(directory):
        return None
    for fname in os.listdir(directory):
        if fname.endswith(".gmt"):
            return os.path.join(directory, fname)
    return None


def _fetch(url: str) -> str:
    ctx = ssl.create_default_context()
    ctx.check_hostname = False
    ctx.verify_mode = ssl.CERT_NONE
    req = urllib.request.Request(
        url, headers={"User-Agent": "Mozilla/5.0 (compatible; ENHYDRA)"}
    )
    with urllib.request.urlopen(req, context=ctx) as r:
        return r.read().decode("utf-8", errors="replace")


def _fetch_cached(url: str, cache_dir: str | None, filename: str) -> str:
    if cache_dir:
        local_path = os.path.join(cache_dir, filename)
        if os.path.isfile(local_path):
            with open(local_path, encoding="utf-8") as fh:
                return fh.read()
    content = _fetch(url)
    if cache_dir:
        os.makedirs(cache_dir, exist_ok=True)
        with open(local_path, "w", encoding="utf-8") as fh:
            fh.write(content)
    return content


def _img_to_base64(img_path: str) -> str:
    with open(img_path, "rb") as fh:
        data = base64.b64encode(fh.read()).decode("utf-8")
    return "data:image/png;base64,%s" % data


def _load_gsea_results(results_dir: str) -> pd.DataFrame | None:
    path = os.path.join(results_dir, "gseapy.gene_set.prerank.report.csv")
    if not os.path.isfile(path):
        logger.warning("GSEA results not found: %s", path)
        return None

    df = pd.read_csv(path)

    if "Tag %" in df.columns:
        def _tag(v):
            try:
                v = str(v).strip()
                if "/" in v:
                    a, b = v.split("/")
                    return round(float(a) / float(b), 4)
                return round(float(v), 4)
            except Exception:
                return None
        df["Tag %"] = df["Tag %"].apply(_tag)

    if "Gene %" in df.columns:
        def _gene(v):
            try:
                return round(float(str(v).strip().rstrip("%")) / 100.0, 4)
            except Exception:
                return None
        df["Gene %"] = df["Gene %"].apply(_gene)

    processed_path = os.path.join(results_dir, "gsea_results_processed.tsv")
    df.to_csv(processed_path, sep="\t", index=False)
    return df


def _build_enrichment_plot_index(results_dir: str) -> dict[str, str]:
    prerank_dir = os.path.join(results_dir, "prerank")
    if not os.path.isdir(prerank_dir):
        return {}
    index = {}
    for filename in os.listdir(prerank_dir):
        if not filename.endswith(".png"):
            continue
        go_id = filename.replace(".png", "").replace("_", ":", 1)
        index[go_id] = _img_to_base64(os.path.join(prerank_dir, filename))
    logger.info("Indexed %d enrichment plot(s).", len(index))
    return index


def _results_table_html(
    df: pd.DataFrame,
    obo_names: dict[str, str],
    plot_index: dict[str, str],
    fdr_threshold: float = 0.25,
    metric: str | None = None,
) -> tuple[str, list[int]]:
    """Build an HTML table from GSEA results.

    Args:
        df:            GSEA results DataFrame.
        obo_names:     GO ID → term name mapping (may already include names
                       merged from the GMT file by the caller).
        plot_index:    GO ID → base64 PNG URI mapping.
        fdr_threshold: FDR threshold for row highlighting.
        metric:        Metric name used to generate a unique table ID and to
                       tag go-link anchors with data-metric so the modal
                       handler can look up the correct enrichment plot when
                       multiple tabs share one page.  None preserves the
                       original single-metric behaviour.
    """
    df = df.copy()

    if "Term" not in df.columns:
        logger.warning("'Term' column not found in GSEA results.")
        return "<p>No results to display.</p>", []

    df["GO Term"] = df["Term"].map(obo_names).fillna(df["Term"])

    for col in ["ES", "NES", "NOM p-val", "FDR q-val", "FWER p-val", "Tag %", "Gene %"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce").apply(
                lambda x: "%.4f" % x if pd.notna(x) else ""
            )

    df["Significant"] = df["FDR q-val"].apply(
        lambda x: "✓" if x != "" and float(x) < fdr_threshold else ""
    )

    col_defs = [
        ("Term",       "GO ID",      "Gene Ontology term identifier."),
        ("GO Term",    "Term name",  "Human-readable name of the GO term."),
        ("NES",        "NES",        "Normalised Enrichment Score. Positive = more conserved, negative = faster evolving."),
        ("NOM p-val",  "p-value",    "Nominal p-value from permutation testing."),
        ("FDR q-val",  "FDR",        "False Discovery Rate q-value. Significant below %.2f." % fdr_threshold),
        ("Tag %",      "Tag %",      "Fraction of gene set genes in the leading edge (0-1)."),
        ("Gene %",     "Gene %",     "Fraction of all ranked genes in the leading edge (0-1)."),
        ("Significant","Sig.",       "Significant at FDR < %.2f." % fdr_threshold),
    ]
    col_defs = [(s, d, t) for s, d, t in col_defs
                if s in df.columns or s in ("GO Term", "Significant")]

    table_df = df[[s for s, _, _ in col_defs]].copy()
    display_names = [d for _, d, _ in col_defs]
    tooltips      = [t for _, _, t in col_defs]
    table_df.columns = display_names

    numeric_display_names = {"NES", "p-value", "FDR", "Tag %", "Gene %"}
    numeric_col_indices   = [i for i, n in enumerate(display_names)
                             if n in numeric_display_names]

    header_cells = "".join(
        '<th>%s <span class="col-tip">?<span class="tip-text">%s</span></span></th>'
        % (name, tip)
        for name, tip in zip(display_names, tooltips)
    )
    filter_cells = "<th></th>" * len(display_names)

    table_id    = ("results-table-%s" % metric) if metric else "results-table"
    metric_attr = (' data-metric="%s"' % metric) if metric else ""

    rows = ""
    for _, row in table_df.iterrows():
        sig_class = ' class="sig-row"' if row.get("Sig.") == "✓" else ""
        go_id     = row.get("GO ID", "")
        cells     = ""
        for col, val in row.items():
            if col == "GO ID" and go_id in plot_index:
                cells += (
                    '<td><a href="#" class="go-link" data-goid="%s"%s>%s</a></td>'
                    % (go_id, metric_attr, val)
                )
            else:
                cells += "<td>%s</td>" % str(val)
        rows += "<tr%s>%s</tr>\n" % (sig_class, cells)

    html = (
        '<table id="%s" class="display compact" style="width:100%%">'
        "<thead><tr>%s</tr><tr class=\"filter-row\">%s</tr></thead>"
        "<tbody>%s</tbody></table>"
    ) % (table_id, header_cells, filter_cells, rows)

    return html, numeric_col_indices


def _plot_section(plots_dir: str, names: list[tuple[str, str]]) -> str:
    html = ""
    for stem, caption in names:
        svg_path = os.path.join(plots_dir, stem + ".svg")
        png_path = os.path.join(plots_dir, stem + ".png")
        if os.path.isfile(svg_path):
            with open(svg_path, encoding="utf-8") as fh:
                svg_content = fh.read()
            html += (
                '<div class="plot-block">'
                '<p class="plot-caption">%s (hover bars for GO definitions)</p>'
                '%s</div>'
            ) % (caption, svg_content)
        elif os.path.isfile(png_path):
            uri = _img_to_base64(png_path)
            html += (
                '<div class="plot-block">'
                '<p class="plot-caption">%s</p>'
                '<img src="%s" alt="%s"/></div>'
            ) % (caption, uri, caption)
    return html


def _resolve_term_names(
    results_dir: str,
    obo_path: str | None,
    gmt_path: str | None,
) -> dict[str, str]:
    """Build the term name lookup used by _results_table_html.

    Merges two sources, with OBO names taking precedence:
    1. GMT second column — always available when a GMT file was used; provides
       term names even when no OBO file has been downloaded.
    2. go-basic.obo — canonical GO term names; used when --obo-cache is set.

    If gmt_path is not provided, the results_dir is scanned for a .gmt file
    written by run_gsea in --organism mode.
    """
    effective_gmt = gmt_path or _find_gmt_in_dir(results_dir)
    gmt_names     = _gmt_term_names(effective_gmt)

    obo_names: dict[str, str] = {}
    if obo_path and os.path.isfile(obo_path):
        obo_names = _parse_obo_names(obo_path)
        logger.info("Loaded %d GO term names from OBO.", len(obo_names))

    # OBO takes precedence: its names are more carefully curated and
    # consistent across tool versions.
    return {**gmt_names, **obo_names}


# ---------------------------------------------------------------------------
# Single-metric report (backward-compatible)
# ---------------------------------------------------------------------------

def build_report(
    results_dir: str,
    plots_dir: str,
    report_path: str,
    obo_path: str | None = None,
    mode: str = "single",
    metric: str = "zscore",
    fdr_threshold: float = 0.25,
    gmt_path: str | None = None,
):
    """Build a self-contained HTML report for a single-metric ENHYDRA run.

    Args:
        results_dir:   Directory containing GSEApy results.
        plots_dir:     Directory containing plot files.
        report_path:   Output path for the HTML file.
        obo_path:      Path to go-basic.obo (optional; used for canonical
                       GO term names).
        mode:          "single" or "differential".
        metric:        Ranking metric label for the report title.
        fdr_threshold: FDR threshold for significance highlighting.
        gmt_path:      Path to the GMT file used for GSEA.  Term names are
                       read from its second column and used when an OBO file
                       is not available.  If None, results_dir is scanned
                       for a .gmt file written by --organism mode.
    """
    logger.info("Building HTML report...")

    term_names = _resolve_term_names(results_dir, obo_path, gmt_path)

    df = _load_gsea_results(results_dir)
    if df is None:
        logger.warning("Cannot build report: no GSEA results found.")
        return

    cache_dir = os.path.dirname(obo_path) if obo_path else None
    logger.info("Fetching DataTables assets (cached after first download)...")
    jquery_js = _fetch_cached(_JQUERY_URL,         cache_dir, "jquery.min.js")
    dt_js     = _fetch_cached(_DATATABLES_JS_URL,  cache_dir, "datatables.min.js")
    dt_css    = _fetch_cached(_DATATABLES_CSS_URL, cache_dir, "datatables.min.css")

    plot_index = _build_enrichment_plot_index(results_dir)
    plot_data_js = "var enrichmentPlots = {%s};" % ",".join(
        '"%s": "%s"' % (go_id, uri) for go_id, uri in plot_index.items()
    )

    if mode == "single":
        plot_names = [
            ("identity_distribution", "Distribution of mean alignment identity"),
            ("gsea_barplot",          "Top enriched gene sets (NES)"),
        ]
        title = "ENHYDRA Single-List Enrichment Report"
    else:
        plot_names = [
            ("identity_scatter",          "Identity comparison between lists"),
            ("differential_distribution", "Differential conservation score distribution"),
            ("identity_distribution",     "Distribution of differential scores"),
            ("gsea_barplot",              "Top differentially enriched gene sets (NES)"),
        ]
        title = "ENHYDRA Differential Enrichment Report"

    plots_html = _plot_section(plots_dir, plot_names)
    table_html, numeric_col_indices = _results_table_html(
        df, term_names, plot_index, fdr_threshold, metric=None
    )

    html = _TEMPLATE.format(
        title=title,
        dt_css=dt_css,
        plots_html=plots_html,
        table_html=table_html,
        jquery_js=jquery_js,
        dt_js=dt_js,
        plot_data_js=plot_data_js,
        numeric_col_indices=numeric_col_indices,
    )

    with open(report_path, "w", encoding="utf-8") as fh:
        fh.write(html)
    logger.info("HTML report written to: %s", report_path)


# ---------------------------------------------------------------------------
# Multi-metric tabbed report
# ---------------------------------------------------------------------------

def build_multi_metric_report(
    metric_data: dict,
    report_path: str,
    obo_path: str | None = None,
    fdr_threshold: float = 0.25,
    mode: str = "single",
    gmt_path: str | None = None,
):
    """Build a self-contained tabbed HTML report covering all ranking metrics.

    Each metric is rendered as a separate tab containing its own plots section
    and an interactive DataTables results table.  DataTables are initialised
    lazily on first tab activation to avoid layout issues with hidden elements.

    Args:
        metric_data:   Dict mapping metric name → {"results_dir": str,
                       "plots_dir": str}.  Tabs appear in insertion order.
        report_path:   Output path for the HTML file.
        obo_path:      Path to go-basic.obo for GO term name lookup (optional).
        fdr_threshold: FDR threshold used for row highlighting.
        mode:          "single" or "differential" — controls plot set and title.
        gmt_path:      Path to the GMT file used for GSEA.  Term names from its
                       second column are used as a fallback when an OBO file is
                       not available.  If None, the first metric's results_dir
                       is scanned for a .gmt file written by --organism mode.
    """
    logger.info(
        "Building multi-metric HTML report (%d metrics)...", len(metric_data)
    )

    # If no explicit GMT path, scan the first metric's results_dir.
    effective_gmt = gmt_path
    if not effective_gmt and metric_data:
        first_results = next(iter(metric_data.values()))["results_dir"]
        effective_gmt = _find_gmt_in_dir(first_results)

    term_names = _resolve_term_names(
        results_dir=next(iter(metric_data.values()))["results_dir"] if metric_data else "",
        obo_path=obo_path,
        gmt_path=effective_gmt,
    )

    cache_dir = os.path.dirname(obo_path) if obo_path else None
    jquery_js = _fetch_cached(_JQUERY_URL,         cache_dir, "jquery.min.js")
    dt_js     = _fetch_cached(_DATATABLES_JS_URL,  cache_dir, "datatables.min.js")
    dt_css    = _fetch_cached(_DATATABLES_CSS_URL, cache_dir, "datatables.min.css")

    title = ("ENHYDRA Multi-Metric Differential Enrichment Report"
             if mode == "differential"
             else "ENHYDRA Multi-Metric Enrichment Report")

    if mode == "single":
        plot_names = [
            ("identity_distribution", "Distribution of mean alignment identity"),
            ("gsea_barplot",          "Top enriched gene sets (NES)"),
        ]
    else:
        plot_names = [
            ("identity_scatter",          "Identity comparison between lists"),
            ("differential_distribution", "Differential conservation score distribution"),
            ("gsea_barplot",              "Top differentially enriched gene sets (NES)"),
        ]

    tab_buttons_parts    = []
    tab_panels_parts     = []
    enrichment_plots_map = {}
    numeric_cols_map     = {}

    first = True
    for metric, paths in metric_data.items():
        label       = METRIC_LABELS.get(metric, metric.capitalize())
        active_cls  = " active" if first else ""
        results_dir = paths["results_dir"]
        plots_dir   = paths["plots_dir"]

        tab_buttons_parts.append(
            '    <button class="tab-btn%s" data-metric="%s" '
            'role="tab" aria-controls="tab-%s">%s</button>'
            % (active_cls, metric, metric, label)
        )

        plot_idx = _build_enrichment_plot_index(results_dir)
        enrichment_plots_map[metric] = plot_idx

        df = _load_gsea_results(results_dir)
        if df is not None:
            tbl_html, num_cols = _results_table_html(
                df, term_names, plot_idx, fdr_threshold, metric=metric
            )
            numeric_cols_map[metric] = num_cols
        else:
            tbl_html = "<p>No GSEA results found for this metric.</p>"
            numeric_cols_map[metric] = []

        plots_html = _plot_section(plots_dir, plot_names)
        desc       = _METRIC_DESCS.get(metric, "")

        tab_panels_parts.append(
            '<div id="tab-{m}" class="tab-panel{ac}" role="tabpanel">\n'
            '  <p class="metric-desc">{desc}</p>\n'
            '  <h3>Plots</h3>\n'
            '  <div class="plot-grid">{plots}</div>\n'
            '  <h3>Enrichment results</h3>\n'
            '  <p>Significant gene sets (FDR&nbsp;&lt;&nbsp;{fdr}) are highlighted '
            'in blue. Click a GO ID to view its enrichment plot. '
            'Use the filter boxes beneath each column header to narrow results.</p>\n'
            '  {tbl}\n'
            '</div>\n'.format(
                m=metric, ac=active_cls, desc=desc,
                plots=plots_html, fdr=fdr_threshold, tbl=tbl_html,
            )
        )
        first = False

    html = _MULTI_TEMPLATE.format(
        title=title,
        dt_css=dt_css,
        tab_buttons="\n".join(tab_buttons_parts),
        tab_panels="\n".join(tab_panels_parts),
        jquery_js=jquery_js,
        dt_js=dt_js,
        enrichment_plots_map=json.dumps(enrichment_plots_map),
        numeric_cols_map=json.dumps(numeric_cols_map),
    )

    with open(report_path, "w", encoding="utf-8") as fh:
        fh.write(html)
    logger.info("Multi-metric HTML report written to: %s", report_path)
