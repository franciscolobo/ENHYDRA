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
    with urllib.request.urlopen(url, context=ctx) as r:
        return r.read().decode("utf-8", errors="replace")


def _img_to_base64(img_path: str) -> str:
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


def _load_gsea_results(results_dir: str) -> pd.DataFrame | None:
    path = os.path.join(results_dir, "gseapy.gene_set.prerank.report.csv")
    if not os.path.isfile(path):
        logger.warning("GSEA results not found: %s", path)
        return None
    return pd.read_csv(path)


def _results_table_html(df: pd.DataFrame, obo_names: dict[str, str]) -> str:
    """Build the HTML for the results DataTable."""
    df = df.copy()

    # Add term name from OBO if available
    df["GO Term"] = df["Name"].map(obo_names).fillna(df.get("Term", df["Name"]))

    # Format numeric columns
    for col in ["ES", "NES", "NOM p-val", "FDR q-val", "FWER p-val"]:
        if col in df.columns:
            df[col] = df[col].apply(
                lambda x: "%.4f" % x if pd.notna(x) else ""
            )

    # Significance flag
    df["Significant"] = df["FDR q-val"].apply(
        lambda x: "✓" if pd.notna(x) and float(x) < 0.25 else ""
    )

    display_cols = ["Name", "GO Term", "NES", "NOM p-val", "FDR q-val", "Significant"]
    display_cols = [c for c in display_cols if c in df.columns]
    table_df = df[display_cols].copy()
    table_df.columns = ["GO ID", "Term name", "NES", "p-value", "FDR", "Significant"]

    rows = ""
    for _, row in table_df.iterrows():
        sig_class = ' class="sig-row"' if row["Significant"] == "✓" else ""
        cells = "".join("<td>%s</td>" % str(v) for v in row)
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
) -> str:
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
section { background: white; border-radius: 8px; box-shadow: 0 1px 4px rgba(0,0,0,0.1);
          padding: 28px 32px; margin-bottom: 28px; }
h2 { margin-top: 0; color: #1a3a5c; border-bottom: 2px solid #e0e0e0;
     padding-bottom: 8px; }
.plot-block { margin: 20px 0; text-align: center; }
.plot-block img { max-width: 100%%; border: 1px solid #e0e0e0;
                  border-radius: 4px; }
.plot-caption { font-size: 0.9em; color: #555; margin-bottom: 6px; }
.plot-grid { display: grid; grid-template-columns: 1fr 1fr; gap: 20px; }
tr.sig-row { background-color: #eaf3fb !important; font-weight: bold; }
footer { text-align: center; padding: 20px; font-size: 0.85em; color: #888; }
</style>
</head>
<body>
<header>
  <h1>%s</h1>
  <p>ENHYDRA — Gene Set Enrichment Analysis for evolutionary genomics</p>
</header>
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
     Click column headers to sort. Use the search box to filter.</p>
%s
</section>

</main>
<footer>Generated by ENHYDRA</footer>

<script>%s</script>
<script>%s</script>
<script>
$(document).ready(function() {
    $('#results-table').DataTable({
        pageLength: 25,
        order: [[4, 'asc']],
        columnDefs: [{ targets: [2,3,4], type: 'num' }]
    });
});
</script>
</body>
</html>""" % (title, dt_css, title, plots_html, table_html, jquery_js, dt_js)


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

    # Download JS/CSS
    logger.info("Downloading DataTables assets for embedding...")
    jquery_js = _fetch(_JQUERY_URL)
    dt_js     = _fetch(_DATATABLES_JS_URL)
    dt_css    = _fetch(_DATATABLES_CSS_URL)

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

    plots_html = _plot_section(plots_dir, plot_names)
    table_html = _results_table_html(df, obo_names)

    html = _html_template(
        title=title,
        jquery_js=jquery_js,
        dt_css=dt_css,
        dt_js=dt_js,
        plots_html=plots_html,
        table_html=table_html,
        mode=mode,
    )

    with open(report_path, "w", encoding="utf-8") as fh:
        fh.write(html)

    logger.info("HTML report written to: %s", report_path)
