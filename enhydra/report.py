from __future__ import annotations

import os
import base64
import logging
import urllib.request
import ssl
import pandas as pd

logger = logging.getLogger(__name__)

# DataTables CDN URLs — downloaded once and embedded
_DATATABLES_JS_URL  = "https://cdn.datatables.net/1.13.6/js/jquery.dataTables.min.js"
_JQUERY_URL         = "https://code.jquery.com/jquery-3.7.0.min.js"
_DATATABLES_CSS_URL = "https://cdn.datatables.net/1.13.6/css/jquery.dataTables.min.css"


def _fetch(url: str) -> str:
    """Download a URL and return its content as a string."""
    ctx = ssl.create_default_context()
    ctx.check_hostname = False
    ctx.verify_mode = ssl.CERT_NONE
    req = urllib.request.Request(
        url, headers={"User-Agent": "Mozilla/5.0 (compatible; ENHYDRA)"}
    )
    with urllib.request.urlopen(req, context=ctx) as r:
        return r.read().decode("utf-8", errors="replace")


def _fetch_cached(url: str, cache_dir: str | None, filename: str) -> str:
    """Download a URL, caching the result locally if cache_dir is provided."""
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
    """Encode a PNG file as a base64 data URI."""
    with open(img_path, "rb") as fh:
        data = base64.b64encode(fh.read()).decode("utf-8")
    return "data:image/png;base64,%s" % data
    """Encode a PNG file as a base64 data URI."""
    with open(img_path, "rb") as fh:
        data = base64.b64encode(fh.read()).decode("utf-8")
    return "data:image/png;base64,%s" % data


def _parse_obo_names(obo_path: str) -> dict[str, str]:
    """Parse go-basic.obo and return GO ID → term name mapping."""
    names = {}
    current_id = None
    current_name = None
    is_obsolete = False
    with open(obo_path) as fh:
        for line in fh:
            line = line.rstrip()
            if line == "[Term]":
                if current_id and not is_obsolete and current_name:
                    names[current_id] = current_name
                current_id = None
                current_name = None
                is_obsolete = False
            elif line.startswith("id: GO:"):
                current_id = line[4:]
            elif line.startswith("name: "):
                current_name = line[6:]
            elif line == "is_obsolete: true":
                is_obsolete = True
    if current_id and not is_obsolete and current_name:
        names[current_id] = current_name
    return names


def _build_enrichment_plot_index(results_dir: str) -> dict[str, str]:
    """Build a mapping from GO ID to base64-encoded PNG enrichment plot.

    Args:
        results_dir: Directory containing the GSEApy prerank/ subdirectory.

    Returns:
        Dict mapping GO ID (e.g. 'GO:0006412') to base64 PNG data URI.
    """
    prerank_dir = os.path.join(results_dir, "prerank")
    if not os.path.isdir(prerank_dir):
        return {}

    index = {}
    for filename in os.listdir(prerank_dir):
        if not filename.endswith(".png"):
            continue
        # GSEApy names files like GO_0006412.png — convert to GO:0006412
        go_id = filename.replace(".png", "").replace("_", ":", 1)
        png_path = os.path.join(prerank_dir, filename)
        index[go_id] = _img_to_base64(png_path)

    logger.info("Indexed %d enrichment plot(s).", len(index))
    return index
    path = os.path.join(results_dir, "gseapy.gene_set.prerank.report.csv")
    if not os.path.isfile(path):
        logger.warning("GSEA results not found: %s", path)
        return None
    return pd.read_csv(path)


def _results_table_html(
    df: pd.DataFrame,
    obo_names: dict[str, str],
    plot_index: dict[str, str],
) -> str:
    """Build the HTML for the results DataTable."""
    df = df.copy()

    # GSEApy prerank results columns:
    # 'Name' = run name (e.g. 'prerank'), 'Term' = GO ID, 'ES', 'NES', etc.
    if "Term" not in df.columns:
        logger.warning("'Term' column not found in GSEA results.")
        return "<p>No results to display.</p>"

    # Look up term name from OBO using the GO ID in the 'Term' column
    df["GO Term"] = df["Term"].map(obo_names).fillna(df["Term"])

    # Format numeric columns
    for col in ["ES", "NES", "NOM p-val", "FDR q-val", "FWER p-val"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce").apply(
                lambda x: "%.4f" % x if pd.notna(x) else ""
            )

    # Significance flag
    df["Significant"] = df["FDR q-val"].apply(
        lambda x: "✓" if x != "" and float(x) < 0.25 else ""
    )

    display_cols = {
        "Term":        "GO ID",
        "GO Term":     "Term name",
        "NES":         "NES",
        "NOM p-val":   "p-value",
        "FDR q-val":   "FDR",
        "Significant": "Significant",
    }
    table_df = df[[c for c in display_cols if c in df.columns
                   or c in ("GO Term", "Significant")]].copy()
    table_df = table_df.rename(columns=display_cols)

    rows = ""
    for _, row in table_df.iterrows():
        sig_class = ' class="sig-row"' if row.get("Significant") == "✓" else ""
        go_id = row.get("GO ID", "")
        has_plot = go_id in plot_index
        cells = ""
        for col, val in row.items():
            if col == "GO ID" and has_plot:
                cells += '<td><a href="#" class="go-link" data-goid="%s">%s</a></td>' % (go_id, val)
            else:
                cells += "<td>%s</td>" % str(val)
        rows += "<tr%s>%s</tr>\n" % (sig_class, cells)

    headers = "".join("<th>%s</th>" % c for c in table_df.columns)
    return """
    <table id="results-table" class="display compact" style="width:100%%">
        <thead><tr>%s</tr></thead>
        <tbody>%s</tbody>
    </table>""" % (headers, rows)


def _plot_section(plots_dir: str, names: list[tuple[str, str]]) -> str:
    """Build HTML img tags for a list of (filename_stem, caption) pairs."""
    html = ""
    for stem, caption in names:
        png_path = os.path.join(plots_dir, stem + ".png")
        if not os.path.isfile(png_path):
            continue
        uri = _img_to_base64(png_path)
        html += """
        <div class="plot-block">
            <p class="plot-caption">%s</p>
            <img src="%s" alt="%s"/>
        </div>""" % (caption, uri, caption)
    return html


def _html_template(
    title: str,
    jquery_js: str,
    dt_css: str,
    dt_js: str,
    plots_html: str,
    table_html: str,
    mode: str,
    plot_index: dict[str, str],
) -> str:
    # Build JS plot data object
    plot_data_js = "var enrichmentPlots = {%s};" % ",".join(
        '"%s": "%s"' % (go_id, uri)
        for go_id, uri in plot_index.items()
    )
    return """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8"/>
<meta name="viewport" content="width=device-width, initial-scale=1.0"/>
<title>%s</title>
<style>
%s
body { font-family: Arial, sans-serif; margin: 0; padding: 0;
       background: #f5f5f5; color: #222; }
header { background: #1a3a5c; color: white; padding: 24px 40px; }
header h1 { margin: 0; font-size: 1.8em; }
header p  { margin: 4px 0 0; font-size: 0.95em; opacity: 0.85; }
main { max-width: 1200px; margin: 32px auto; padding: 0 24px; }
section { background: white; border-radius: 8px;
          box-shadow: 0 1px 4px rgba(0,0,0,0.1);
          padding: 28px 32px; margin-bottom: 28px; }
h2 { margin-top: 0; color: #1a3a5c; border-bottom: 2px solid #e0e0e0;
     padding-bottom: 8px; }
.plot-block { margin: 20px 0; text-align: center; }
.plot-block img { max-width: 100%%; border: 1px solid #e0e0e0;
                  border-radius: 4px; }
.plot-caption { font-size: 0.9em; color: #555; margin-bottom: 6px; }
.plot-grid { display: grid; grid-template-columns: 1fr 1fr; gap: 20px; }
tr.sig-row { background-color: #eaf3fb !important; font-weight: bold; }
a.go-link { color: #1a3a5c; text-decoration: underline dotted; cursor: pointer; }
/* Modal */
#modal-overlay { display: none; position: fixed; top: 0; left: 0;
                 width: 100%%; height: 100%%; background: rgba(0,0,0,0.6);
                 z-index: 1000; justify-content: center; align-items: center; }
#modal-overlay.active { display: flex; }
#modal-box { background: white; border-radius: 8px; padding: 24px;
             max-width: 700px; width: 90%%; position: relative; }
#modal-title { font-size: 1.1em; font-weight: bold; color: #1a3a5c;
               margin-bottom: 12px; }
#modal-img { width: 100%%; border: 1px solid #e0e0e0; border-radius: 4px; }
#modal-close { position: absolute; top: 12px; right: 16px; font-size: 1.4em;
               cursor: pointer; color: #555; background: none; border: none; }
footer { text-align: center; padding: 20px; font-size: 0.85em; color: #888; }
</style>
</head>
<body>
<header>
  <h1>%s</h1>
  <p>ENHYDRA — Gene Set Enrichment Analysis for evolutionary genomics</p>
</header>

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
%s
  </div>
</section>

<section>
  <h2>Enrichment results</h2>
  <p>Significant gene sets (FDR &lt; 0.25) are highlighted in blue.
     Click a GO ID to view its enrichment plot.</p>
%s
</section>
</main>
<footer>Generated by ENHYDRA</footer>

<script>%s</script>
<script>%s</script>
<script>
%s
$(document).ready(function() {
    $('#results-table').DataTable({
        pageLength: 25,
        order: [[4, 'asc']],
        columnDefs: [{ targets: [2,3,4], type: 'num' }]
    });
    $(document).on('click', '.go-link', function(e) {
        e.preventDefault();
        var goId = $(this).data('goid');
        var uri  = enrichmentPlots[goId];
        if (uri) {
            $('#modal-title').text(goId);
            $('#modal-img').attr('src', uri);
            $('#modal-overlay').addClass('active');
        }
    });
    $('#modal-close, #modal-overlay').on('click', function(e) {
        if (e.target === this) $('#modal-overlay').removeClass('active');
    });
});
</script>
</body>
</html>""" % (
        title, dt_css, title,
        plots_html, table_html,
        jquery_js, dt_js, plot_data_js
    )


def build_report(
    results_dir: str,
    plots_dir: str,
    report_path: str,
    obo_path: str | None = None,
    mode: str = "single",
    metric: str = "zscore",
):
    """Build a self-contained HTML report for ENHYDRA results.

    Args:
        results_dir:  Directory containing GSEApy results CSV.
        plots_dir:    Directory containing PNG plot files.
        report_path:  Output path for the HTML file.
        obo_path:     Path to go-basic.obo for GO term name lookup (optional).
        mode:         'single' or 'differential'.
        metric:       Differential metric used (only relevant in differential mode).
    """
    logger.info("Building HTML report...")

    # Load GO term names
    obo_names = {}
    if obo_path and os.path.isfile(obo_path):
        obo_names = _parse_obo_names(obo_path)
        logger.info("Loaded %d GO term names from OBO.", len(obo_names))
    else:
        logger.warning("OBO file not provided or not found — term names from GMT only.")

    # Load GSEA results
    df = _load_gsea_results(results_dir)
    if df is None:
        logger.warning("Cannot build report: no GSEA results found.")
        return

    # Download JS/CSS (cached in obo_cache dir if available)
    cache_dir = os.path.dirname(obo_path) if obo_path else None
    logger.info("Fetching DataTables assets (cached after first download)...")
    jquery_js = _fetch_cached(_JQUERY_URL,         cache_dir, "jquery.min.js")
    dt_js     = _fetch_cached(_DATATABLES_JS_URL,  cache_dir, "datatables.min.js")
    dt_css    = _fetch_cached(_DATATABLES_CSS_URL, cache_dir, "datatables.min.css")

    # Build plots section
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

    # Build enrichment plot index
    plot_index = _build_enrichment_plot_index(results_dir)

    plots_html = _plot_section(plots_dir, plot_names)
    table_html = _results_table_html(df, obo_names, plot_index)

    html = _html_template(
        title=title,
        jquery_js=jquery_js,
        dt_css=dt_css,
        dt_js=dt_js,
        plots_html=plots_html,
        table_html=table_html,
        mode=mode,
        plot_index=plot_index,
    )

    with open(report_path, "w", encoding="utf-8") as fh:
        fh.write(html)

    logger.info("HTML report written to: %s", report_path)
