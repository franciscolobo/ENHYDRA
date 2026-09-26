from __future__ import annotations
from .io import parse_obo_names as _parse_obo_names
from .filtering import display_group_id

import json
import shlex
import os
import html
import base64
import logging
import urllib.request
import ssl
from urllib.parse import quote as _urlquote
import pandas as pd

logger = logging.getLogger(__name__)

_DATATABLES_JS_URL  = "https://cdn.datatables.net/1.13.6/js/jquery.dataTables.min.js"
_JQUERY_URL         = "https://code.jquery.com/jquery-3.7.0.min.js"
_DATATABLES_CSS_URL = "https://cdn.datatables.net/1.13.6/css/jquery.dataTables.min.css"
_XLSX_JS_URL        = "https://cdnjs.cloudflare.com/ajax/libs/xlsx/0.18.5/xlsx.full.min.js"

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

# ---------------------------------------------------------------------------
# Unified report template — one tab shell used for every report, regardless
# of how many ranking metrics were run. A single-metric report simply has
# one metric tab plus the Alignments and Filtering summary tabs; a
# multi-metric report has three metric tabs, a Cross-metric consensus tab,
# plus the same two. This replaces what used to be two separate templates
# (a flat non-tabbed layout for single-metric reports, and a tabbed layout
# for multi-metric reports) — they had been drifting into near-duplicates
# anyway, and the Filtering summary / Alignments tabs need to appear in
# both cases.
# ---------------------------------------------------------------------------

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
h3 {{ color: #2c5282; margin: 24px 0 12px; }}
h4 {{ color: #1a3a5c; margin: 20px 0 6px; }}
h5 {{ color: #2c5282; margin: 14px 0 4px; font-size: 0.95em; }}
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
a.go-link, a.geneset-link, a.leadedge-link, a.detail-link {{
    color: #1a3a5c; text-decoration: underline dotted; cursor: pointer; }}
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
                max-width: 360px; line-height: 1.5; pointer-events: none;
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
#modal-text {{ display: none; max-height: 420px; overflow-y: auto;
               text-align: left; white-space: pre-wrap; word-break: break-word;
               font-family: 'Courier New', monospace; font-size: 12px;
               background: #f7f7f7; border: 1px solid #e0e0e0;
               border-radius: 4px; padding: 10px 14px; margin: 0; }}
#modal-close {{ position: absolute; top: 12px; right: 16px; font-size: 1.4em;
                cursor: pointer; color: #555; background: none; border: none; }}
.export-btn {{ padding: 8px 18px; border: none; border-radius: 4px;
               background: #1a3a5c; color: white; font-size: 13px;
               font-weight: 600; cursor: pointer; margin-bottom: 14px;
               margin-right: 8px; }}
.export-btn:hover {{ background: #12293f; }}
.sig-toggle-btn {{ padding: 8px 18px; border: 2px solid #1a3a5c; border-radius: 4px;
                    background: white; color: #1a3a5c; font-size: 13px;
                    font-weight: 600; cursor: pointer; margin-bottom: 14px;
                    margin-right: 8px; }}
.sig-toggle-btn:hover {{ background: #f0f5fa; }}
.sig-toggle-btn.active {{ background: #1a3a5c; color: white; }}
.funnel-block {{ margin-bottom: 18px; }}
.funnel-row {{ display: flex; align-items: center; flex-wrap: wrap; gap: 6px;
               margin-top: 8px; }}
.funnel-step {{ background: #1a3a5c; color: white; border-radius: 6px;
                padding: 10px 14px; text-align: center; min-width: 92px; }}
.funnel-num {{ font-size: 1.3em; font-weight: 700; }}
.funnel-label {{ font-size: 0.72em; opacity: 0.9; margin-top: 2px; }}
.funnel-arrow {{ font-size: 1.2em; color: #888; padding: 0 4px; }}
.filtering-overlap {{ margin-bottom: 22px; padding: 14px 20px; background: #f0f5fa;
                      border-left: 3px solid #1a3a5c; border-radius: 0 6px 6px 0; }}
.filtering-cols {{ display: flex; gap: 24px; flex-wrap: wrap; }}
.filtering-col {{ flex: 1 1 380px; background: #fafbfc; border: 1px solid #e0e0e0;
                  border-radius: 6px; padding: 16px 20px; min-width: 320px; }}
table.reason-table {{ width: 100%; border-collapse: collapse; margin: 6px 0 16px; }}
table.reason-table th, table.reason-table td {{ border: 1px solid #e0e0e0;
    padding: 5px 10px; font-size: 12.5px; text-align: left; }}
table.reason-table th {{ background: #eef2f6; }}
p.no-drops {{ color: #2f7d3c; font-size: 0.88em; margin: 4px 0 16px; }}
.aln-tree-controls {{ margin-bottom: 14px; }}
.aln-tree-controls input {{ width: 100%; max-width: 480px; padding: 8px 12px;
    font-size: 13px; border: 1px solid #ccc; border-radius: 4px;
    box-sizing: border-box; }}
details.aln-tree-term {{ background: #fafbfc; border: 1px solid #e0e0e0;
    border-radius: 6px; margin-bottom: 8px; padding: 8px 14px; }}
details.aln-tree-term > summary {{ cursor: pointer; font-weight: 600;
    color: #1a3a5c; padding: 4px 0; }}
.aln-tree-row {{ display: flex; gap: 18px; align-items: center; flex-wrap: wrap;
    padding: 6px 4px 6px 20px; border-top: 1px solid #eee; font-size: 13px; }}
.aln-tree-group {{ font-weight: 600; min-width: 160px; }}
.aln-tree-gene {{ color: #555; min-width: 200px; }}
.aln-tree-identity {{ color: #555; min-width: 120px; }}
.aln-tree-diff {{ color: #1a3a5c; font-weight: 600; min-width: 100px; }}
.gene-chip-list {{ display: flex; flex-wrap: wrap; gap: 6px; text-align: left;
    font-family: 'Courier New', monospace; font-size: 12px; }}
.gene-chip {{ background: #eef2f6; border-radius: 4px; padding: 3px 8px;
    display: inline-flex; align-items: center; gap: 5px; }}
.gene-chip-badge {{ background: #1a3a5c; color: white; border-radius: 3px;
    padding: 1px 5px; font-size: 10px; font-weight: 700; text-decoration: none; }}
.gene-chip-badge:hover {{ background: #12293f; }}
a.aln-link {{ color: #1a3a5c; text-decoration: underline dotted; margin-right: 10px; }}
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
    <pre id="modal-text"></pre>
  </div>
</div>
<main>
<section>
  <h2>Results</h2>
  <nav class="tab-nav" role="tablist">
{tab_buttons}
  </nav>
{tab_panels}
</section>
</main>
<footer>Generated by ENHYDRA</footer>
<script>{jquery_js}</script>
<script>{dt_js}</script>
<script>{xlsx_js}</script>
<script>
var enrichmentPlotsMap = {enrichment_plots_map};
var fullGeneSets       = {full_gene_sets_js};
var leadingEdgeMap     = {leading_edge_map_js};
var fullGeneSetsHtml   = {full_gene_sets_html_js};
var leadingEdgeHtmlMap = {leading_edge_html_map_js};
var numericColsMap     = {numeric_cols_map};
var dtInstances        = {{}};
var colFiltersMap      = {{}};
var sigOnlyMap         = {{}};
var sigColIndexMap     = {{}};
// Per-table default initial sort override. Tables not listed here keep the
// historical default of sorting by column index 4 ascending (FDR, for the
// per-metric enrichment tables' fixed column layout). The consensus table
// has a different column shape and is already pre-sorted server-side (by
// number of metrics in agreement, descending), so it opts out of any
// client-side re-sort on load by specifying an empty order array.
var defaultOrderMap    = {{ consensus: [] }};

function getHeaderIndex(tableSelector, name) {{
    var headers = [];
    $(tableSelector).find('thead tr').first().find('th').each(function() {{
        var $clone = $(this).clone();
        $clone.find('.col-tip').remove();
        headers.push($clone.text().trim());
    }});
    return headers.indexOf(name);
}}

$.fn.dataTable.ext.search.push(function(settings, data) {{
    var metric = settings.nTable.id.replace('results-table-', '');
    if (sigOnlyMap[metric]) {{
        var sigCol = sigColIndexMap[metric];
        if (sigCol !== undefined && sigCol !== null &&
            data[sigCol] !== '\u2713') {{
            return false;
        }}
    }}
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
    sigOnlyMap[metric]    = false;
    var numericCols = numericColsMap[metric] || [];
    var order = defaultOrderMap.hasOwnProperty(metric)
        ? defaultOrderMap[metric] : [[4, 'asc']];
    var dt = $('#results-table-' + metric).DataTable({{
        pageLength: 25, orderCellsTop: true, order: order,
        columnDefs: [{{ targets: numericCols, type: 'num' }}],
    }});
    sigColIndexMap[metric] = getHeaderIndex('#results-table-' + metric, 'Sig.');
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
                    var m = val.match(/^([<>=!]=?)\\s*([\\d.eE+\\-]+)$/);
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
function initFilteringSummaryTables() {{
    if (dtInstances['filtering-summary']) return;
    $('#tab-filtering-summary table.species-count-table').each(function() {{
        $(this).DataTable({{ pageLength: 25, order: [] }});
    }});
    dtInstances['filtering-summary'] = true;
}}
function showImageModal(title, uri) {{
    $('#modal-title').text(title);
    $('#modal-text').hide();
    $('#modal-img').attr('src', uri).show();
    $('#modal-overlay').addClass('active');
}}
function showTextModal(title, text) {{
    $('#modal-title').text(title);
    $('#modal-img').hide();
    $('#modal-text').text(text).show();
    $('#modal-overlay').addClass('active');
}}
function showHtmlModal(title, htmlContent) {{
    $('#modal-title').text(title);
    $('#modal-img').hide();
    $('#modal-text').html(htmlContent).show();
    $('#modal-overlay').addClass('active');
}}
function tableToXLSX(tableSelector, numericCols, fullGeneSetsMap, leadEdgeMap, filename, onlyFiltered) {{
    var $table = $(tableSelector);
    if ($table.length === 0) return;

    var headers = [];
    $table.find('thead tr').first().find('th').each(function() {{
        var $clone = $(this).clone();
        $clone.find('.col-tip').remove();
        headers.push($clone.text().trim());
    }});
    var goIdCol     = headers.indexOf('GO ID');
    var fullSetCol  = headers.indexOf('Full gene set');
    var leadEdgeCol = headers.indexOf('Leading edge');

    var dt       = $table.DataTable();
    var selector = onlyFiltered ? {{ search: 'applied' }} : {{}};
    var rows     = dt.rows(selector).nodes();
    var aoa      = [headers];

    $(rows).each(function() {{
        var cells   = $(this).find('td');
        var goId    = goIdCol >= 0 ? $(cells[goIdCol]).text().trim() : null;
        var rowData = [];
        cells.each(function(i) {{
            var text = $(this).text().trim();
            if (i === fullSetCol && goId && fullGeneSetsMap &&
                fullGeneSetsMap[goId] !== undefined) {{
                rowData.push(fullGeneSetsMap[goId].split('\\n').join('; '));
            }} else if (i === leadEdgeCol && goId && leadEdgeMap &&
                       leadEdgeMap[goId] !== undefined) {{
                rowData.push(leadEdgeMap[goId].split('\\n').join('; '));
            }} else if (numericCols.indexOf(i) !== -1) {{
                var v = parseFloat(text);
                rowData.push(isNaN(v) ? text : v);
            }} else {{
                rowData.push(text);
            }}
        }});
        aoa.push(rowData);
    }});

    var ws = XLSX.utils.aoa_to_sheet(aoa);
    var wb = XLSX.utils.book_new();
    XLSX.utils.book_append_sheet(wb, ws, 'Enrichment results');
    XLSX.writeFile(wb, filename);
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
        el.addEventListener('mouseleave', function() {{ svgTip.style.display = 'none'; }});
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
            if (metric === 'filtering-summary') {{
                initFilteringSummaryTables();
            }} else if (numericColsMap.hasOwnProperty(metric)) {{
                initTable(metric);
            }}
        }});
    }});
    var alignFilterInput = document.getElementById('aln-tree-filter');
    if (alignFilterInput) {{
        alignFilterInput.addEventListener('input', function() {{
            var q = this.value.trim().toLowerCase();
            document.querySelectorAll('#aln-tree .aln-tree-term').forEach(function(term) {{
                var termMatch = (term.getAttribute('data-term-search') || '').indexOf(q) !== -1;
                var rows = term.querySelectorAll('.aln-tree-row');
                var anyRowMatch = false;
                rows.forEach(function(row) {{
                    var rowMatch = !q || (row.getAttribute('data-search') || '').indexOf(q) !== -1;
                    row.style.display = rowMatch ? '' : 'none';
                    if (rowMatch) anyRowMatch = true;
                }});
                var show = !q || termMatch || anyRowMatch;
                term.style.display = show ? '' : 'none';
                if (q && show) {{ term.open = true; }}
                if (!q) {{ term.open = false; }}
            }});
        }});
    }}
    $(document).on('click', '.go-link', function(e) {{
        e.preventDefault();
        var goId   = $(this).data('goid');
        var metric = $(this).data('metric');
        var plots  = metric ? enrichmentPlotsMap[metric] : undefined;
        var uri    = plots ? plots[goId] : undefined;
        if (uri) showImageModal(goId, uri);
    }});
    $(document).on('click', '.geneset-link', function(e) {{
        e.preventDefault();
        var goId = $(this).data('goid');
        var htmlContent = fullGeneSetsHtml[goId];
        if (htmlContent !== undefined) {{
            showHtmlModal(goId + ' \u2014 full gene set', htmlContent);
        }}
    }});
    $(document).on('click', '.leadedge-link', function(e) {{
        e.preventDefault();
        var goId   = $(this).data('goid');
        var metric = $(this).data('metric');
        var map    = metric ? leadingEdgeHtmlMap[metric] : undefined;
        var htmlContent = map ? map[goId] : undefined;
        if (htmlContent !== undefined) {{
            showHtmlModal(goId + ' \u2014 leading edge genes (' + metric + ')', htmlContent);
        }}
    }});
    $(document).on('click', '.sig-toggle-btn', function() {{
        var metric = $(this).data('metric');
        sigOnlyMap[metric] = !sigOnlyMap[metric];
        $(this).toggleClass('active', sigOnlyMap[metric]);
        $(this).text(sigOnlyMap[metric] ? 'Showing significant only \u2713' : 'Show only significant');
        if (dtInstances[metric]) dtInstances[metric].draw();
    }});
    $(document).on('click', '.export-xlsx-btn', function() {{
        var metric      = $(this).data('metric');
        var numericCols = numericColsMap[metric] || [];
        var leadMap     = leadingEdgeMap[metric] || {{}};
        tableToXLSX('#results-table-' + metric, numericCols, fullGeneSets, leadMap,
                   'enrichment_results_' + metric + '.xlsx', true);
    }});
    $(document).on('click', '.export-xlsx-all-btn', function() {{
        var metric      = $(this).data('metric');
        var numericCols = numericColsMap[metric] || [];
        var leadMap     = leadingEdgeMap[metric] || {{}};
        tableToXLSX('#results-table-' + metric, numericCols, fullGeneSets, leadMap,
                   'enrichment_results_' + metric + '_all.xlsx', false);
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
# Standalone drop-details pages (one per drop-reasons/skipped-groups TSV)
# ---------------------------------------------------------------------------
# These intentionally do NOT depend on jQuery/DataTables — they are meant to
# stay tiny even when a stage has thousands of dropped rows, and the main
# report should not need to embed or fetch this data at all. Filtering is
# done with plain JS string matching, which is more than sufficient for a
# single-column reason filter plus a free-text search box.

_DETAILS_PAGE_TEMPLATE = """\
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8"/>
<title>{title}</title>
<style>
body {{ font-family: Arial, sans-serif; margin: 24px; color: #222; background: #f5f5f5; }}
h1 {{ color: #1a3a5c; font-size: 1.3em; }}
.controls {{ margin-bottom: 14px; }}
.controls input {{ padding: 6px 10px; font-size: 13px;
    border: 1px solid #ccc; border-radius: 4px; margin-right: 8px; width: 260px; }}
table {{ border-collapse: collapse; width: 100%; background: white; }}
th, td {{ border: 1px solid #e0e0e0; padding: 6px 10px; font-size: 13px; text-align: left; }}
th {{ background: #1a3a5c; color: white; position: sticky; top: 0; }}
tr:nth-child(even) {{ background: #f9f9f9; }}
.count {{ color: #555; font-size: 0.9em; margin-bottom: 10px; }}
</style>
</head>
<body>
<h1>{title}</h1>
<p class="count"><span id="visible-count">{n_rows}</span> of {n_rows} rows shown</p>
<div class="controls">
  <input type="text" id="reason-filter" placeholder="Filter by reason (exact match)..."/>
  <input type="text" id="text-filter" placeholder="Search all columns..."/>
</div>
<table id="details-table">
<thead><tr>{header_html}</tr></thead>
<tbody>
{rows_html}
</tbody>
</table>
<script>
function getQueryParam(name) {{
    var params = new URLSearchParams(window.location.search);
    return params.get(name);
}}
function applyFilters() {{
    var reasonVal = document.getElementById('reason-filter').value.trim().toLowerCase();
    var textVal   = document.getElementById('text-filter').value.trim().toLowerCase();
    var rows      = document.querySelectorAll('#details-table tbody tr');
    var visible   = 0;
    rows.forEach(function(row) {{
        var reasonMatch = !reasonVal || (row.getAttribute('data-reason') || '').toLowerCase() === reasonVal;
        var textMatch   = !textVal || row.textContent.toLowerCase().indexOf(textVal) !== -1;
        var show = reasonMatch && textMatch;
        row.style.display = show ? '' : 'none';
        if (show) visible++;
    }});
    document.getElementById('visible-count').textContent = visible;
}}
document.getElementById('reason-filter').addEventListener('input', applyFilters);
document.getElementById('text-filter').addEventListener('input', applyFilters);
var initialReason = getQueryParam('reason');
if (initialReason) {{
    document.getElementById('reason-filter').value = initialReason;
    applyFilters();
}}
</script>
</body>
</html>"""


def _tsv_to_details_html(tsv_path: str, html_path: str, title: str) -> bool:
    """Render a drop-reasons-style TSV as a standalone, filterable HTML page.

    Returns False (writing nothing) if the source TSV doesn't exist or is
    empty (header-only), so callers can treat 'no page generated' as the
    signal to omit the corresponding 'View list' link entirely.
    """
    if not os.path.isfile(tsv_path):
        return False
    with open(tsv_path) as fh:
        lines = [l.rstrip("\n").split("\t") for l in fh if l.strip()]
    if len(lines) < 2:   # header only, or empty
        return False
    header, rows = lines[0], lines[1:]
    reason_idx = header.index("reason") if "reason" in header else None

    header_html = "".join("<th>%s</th>" % h for h in header)
    row_lines = []
    for row in rows:
        reason_val = row[reason_idx] if reason_idx is not None and reason_idx < len(row) else ""
        cells = "".join("<td>%s</td>" % c for c in row)
        row_lines.append('<tr data-reason="%s">%s</tr>' % (reason_val, cells))

    html_content = _DETAILS_PAGE_TEMPLATE.format(
        title=title, header_html=header_html,
        rows_html="\n".join(row_lines), n_rows=len(rows),
    )
    with open(html_path, "w", encoding="utf-8") as fh:
        fh.write(html_content)
    return True


def _generate_drop_details_pages(
    listdir: str,
    detail_files: dict,
    report_dir: str,
) -> dict[str, str | None]:
    """Write a standalone details page next to each drop-reasons TSV.

    Args:
        listdir:      The list's root directory (dirname of its
                      pipeline_stats.json), used to resolve the relative
                      paths stored in detail_files.
        detail_files: The 'detail_files' dict from pipeline_stats.json,
                      mapping a stage key to a path relative to listdir
                      (or null if that stage's file wasn't produced).
        report_dir:   Directory the report.html itself will be written
                      into, used as the base for the returned link paths
                      (same convention as the enrichment-plot links).

    Returns:
        Dict with the same keys as detail_files, mapped to a path relative
        to report_dir for the generated HTML page — or None if no page was
        generated (source TSV missing/empty, or not applicable — e.g.
        species_counts, which isn't a drop-reasons file and is rendered
        inline instead, is passed through untouched as None here since
        callers handle it separately).
    """
    out: dict[str, str | None] = {}
    for key, rel_tsv in (detail_files or {}).items():
        if key == "group_filter_species_counts":
            continue
        if not rel_tsv:
            out[key] = None
            continue
        tsv_path  = os.path.join(listdir, rel_tsv)
        html_path = os.path.splitext(tsv_path)[0] + ".html"
        title     = key.replace("_", " ").title()
        if _tsv_to_details_html(tsv_path, html_path, title):
            out[key] = os.path.relpath(html_path, report_dir).replace(os.sep, "/")
        else:
            out[key] = None
    return out


# ---------------------------------------------------------------------------
# Filtering summary tab content
# ---------------------------------------------------------------------------

def _load_json(path: str | None) -> dict | None:
    if not path or not os.path.isfile(path):
        return None
    with open(path) as fh:
        return json.load(fh)


def _funnel_html(name: str, stats: dict) -> str:
    steps = [
        ("Input groups",         stats.get("n_input")),
        ("After length filter",  stats.get("n_after_length_filter")),
        ("After group filter",   stats.get("n_after_group_filter")),
    ]
    # The divergence filter is optional and off by default (see
    # stats.aggregate_pipeline_stats()'s doc comment on auto-detection).
    # A run that never enabled it has no meaningful count to show here —
    # inserting a permanently dashed/empty funnel step for every ordinary
    # run would be clutter, not information, so the step is only added
    # when it actually ran for this list.
    if stats.get("n_after_divergence_filter") is not None:
        steps.append((
            "After divergence filter", stats.get("n_after_divergence_filter")
        ))
    steps.append(("Final groups",  stats.get("n_final_groups")))
    steps.append(("Anchor-mapped", stats.get("n_anchor_mapped")))

    parts = []
    for i, (label, val) in enumerate(steps):
        if i > 0:
            parts.append('<div class="funnel-arrow">&rarr;</div>')
        parts.append(
            '<div class="funnel-step"><div class="funnel-num">%s</div>'
            '<div class="funnel-label">%s</div></div>'
            % (val if val is not None else "\u2013", label)
        )
    return (
        '<div class="funnel-block"><h4>%s</h4><div class="funnel-row">%s</div></div>'
        % (name, "".join(parts))
    )


def _reason_table_html(
    reason_counts: dict[str, int],
    detail_html_rel: str | None,
    empty_message: str,
) -> str:
    if not reason_counts:
        return '<p class="no-drops">%s</p>' % empty_message
    rows = ""
    for reason, count in sorted(reason_counts.items(), key=lambda kv: -kv[1]):
        if detail_html_rel:
            link = (
                '<a href="%s?reason=%s" target="_blank" class="detail-link">View list</a>'
                % (detail_html_rel, _urlquote(reason))
            )
        else:
            link = ""
        rows += "<tr><td>%s</td><td>%s</td><td>%s</td></tr>" % (reason, count, link)
    return (
        '<table class="reason-table"><thead><tr><th>Reason</th><th>Count</th>'
        '<th></th></tr></thead><tbody>%s</tbody></table>' % rows
    )


def _species_counts_table_html(tsv_path: str | None, table_id: str) -> str:
    if not tsv_path or not os.path.isfile(tsv_path):
        return "<p>No species count data available.</p>"
    with open(tsv_path) as fh:
        lines = [l.rstrip("\n").split("\t") for l in fh if l.strip()]
    if len(lines) < 2:
        return "<p>No species count data available.</p>"
    header, rows = lines[0], lines[1:]
    header_html = "".join("<th>%s</th>" % h for h in header)
    rows_html   = "".join(
        "<tr>%s</tr>" % "".join("<td>%s</td>" % c for c in row)
        for row in rows
    )
    return (
        '<table id="%s" class="species-count-table display compact" style="width:100%%">'
        '<thead><tr>%s</tr></thead><tbody>%s</tbody></table>'
        % (table_id, header_html, rows_html)
    )


def _divergence_filter_reasons_html(
    ds: dict,
    detail_links: dict[str, str | None],
    heading_tag: str = "h5",
) -> str:
    """Build the two divergence-filter drop-reason subsections, if that
    step ran for this list.

    Both subsections ('sequences removed' and 'groups dropped') link to
    the same drop_reasons.tsv details page — that file carries both
    reasons in one 'reason' column, and the generic '?reason=' query-param
    filtering already used by every other _reason_table_html() call (see
    _generate_drop_details_pages()/_tsv_to_details_html()) narrows the
    details page to just the relevant rows on click, with no special
    handling needed here for this being a shared file.

    Args:
        ds:            The 'drop_summary' dict from pipeline_stats.json.
                      If it has no 'divergence_filter' key, the step did
                      not run for this list and an empty string is
                      returned — no heading, no placeholder text, since
                      showing a permanently-empty divergence filter
                      section on every ordinary run (where this filter is
                      off by default) would be clutter.
        detail_links:  Output of _generate_drop_details_pages(), used to
                      resolve the 'divergence_filter_drop_reasons' link.
        heading_tag:   'h4' in single-list mode, 'h5' in the two-list
                      per-column blocks — matches the surrounding
                      section's own heading level in each caller.

    Returns:
        HTML string, or "" if the divergence filter did not run.
    """
    if "divergence_filter" not in ds:
        return ""
    div_ds = ds["divergence_filter"]
    detail_link = detail_links.get("divergence_filter_drop_reasons")
    return (
        "<%s>Divergence filter \u2014 sequences removed</%s>"
        % (heading_tag, heading_tag)
        + _reason_table_html(
            div_ds.get("sequences_removed", {}),
            detail_link,
            "No sequences were removed by the divergence filter.",
        )
        + "<%s>Divergence filter \u2014 groups dropped</%s>"
        % (heading_tag, heading_tag)
        + _reason_table_html(
            div_ds.get("groups_dropped", {}),
            detail_link,
            "No groups were dropped by the divergence filter.",
        )
    )


def _one_list_filtering_block(
    name: str,
    pipeline_stats_path: str | None,
    report_dir: str,
    table_id_prefix: str,
) -> str:
    stats = _load_json(pipeline_stats_path)
    if stats is None:
        return (
            '<div class="filtering-col"><h4>%s</h4>'
            '<p>No pipeline statistics available.</p></div>' % name
        )
    listdir       = os.path.dirname(os.path.abspath(pipeline_stats_path))
    detail_files  = stats.get("detail_files", {})
    detail_links  = _generate_drop_details_pages(listdir, detail_files, report_dir)
    species_rel   = detail_files.get("group_filter_species_counts")
    species_tsv   = os.path.join(listdir, species_rel) if species_rel else None
    ds            = stats.get("drop_summary", {})

    funnel = _funnel_html(name, stats)
    reasons_html = (
        "<h5>Length filter \u2014 groups skipped entirely</h5>"
        + _reason_table_html(
            ds.get("length_filter", {}).get("groups_skipped", {}),
            detail_links.get("length_filter_skipped_groups"),
            "No groups were skipped entirely at this step.",
        )
        + "<h5>Length filter \u2014 sequences removed as outliers</h5>"
        + _reason_table_html(
            ds.get("length_filter", {}).get("sequences_removed", {}),
            detail_links.get("length_filter_drop_reasons"),
            "No individual sequences were removed as length outliers.",
        )
        + "<h5>Group filter \u2014 groups removed</h5>"
        + _reason_table_html(
            ds.get("group_filter", {}),
            detail_links.get("group_filter_drop_reasons"),
            "No groups were removed at the group filter step.",
        )
        + _divergence_filter_reasons_html(ds, detail_links, heading_tag="h5")
        + "<h5>Tables step \u2014 anomalies</h5>"
        + _reason_table_html(
            ds.get("tables", {}),
            detail_links.get("tables_drop_reasons"),
            "No anomalies at the table generation step.",
        )
    )
    species_html = (
        "<h5>Groups per species (surviving groups)</h5>"
        + _species_counts_table_html(species_tsv, table_id_prefix + "-species-counts")
    )
    return '<div class="filtering-col"><h4>%s</h4>%s%s%s</div>' % (
        name, funnel, reasons_html, species_html,
    )


def _build_filtering_summary_tab_content(
    mode: str,
    label1: str,
    label2: str,
    pipeline_stats_path: str | None,
    pipeline_stats_path1: str | None,
    pipeline_stats_path2: str | None,
    differential_stats_path: str | None,
    report_dir: str,
) -> str:
    """Build the Filtering summary tab body.

    mode == 'differential' corresponds exactly to two-list mode throughout
    this module (see build_report()/build_multi_metric_report() call sites
    in cli.py), so it doubles as the single-list vs two-list switch here —
    no separate flag is needed.
    """
    if mode != "differential":
        stats = _load_json(pipeline_stats_path)
        if stats is None:
            return "<p>No pipeline statistics available for this run.</p>"
        listdir      = os.path.dirname(os.path.abspath(pipeline_stats_path))
        detail_files = stats.get("detail_files", {})
        detail_links = _generate_drop_details_pages(listdir, detail_files, report_dir)
        species_rel  = detail_files.get("group_filter_species_counts")
        species_tsv  = os.path.join(listdir, species_rel) if species_rel else None
        ds           = stats.get("drop_summary", {})

        funnel = _funnel_html("Groups", stats)
        reasons_html = (
            "<h4>Length filter \u2014 groups skipped entirely</h4>"
            + _reason_table_html(
                ds.get("length_filter", {}).get("groups_skipped", {}),
                detail_links.get("length_filter_skipped_groups"),
                "No groups were skipped entirely at this step.",
            )
            + "<h4>Length filter \u2014 sequences removed as outliers</h4>"
            + _reason_table_html(
                ds.get("length_filter", {}).get("sequences_removed", {}),
                detail_links.get("length_filter_drop_reasons"),
                "No individual sequences were removed as length outliers.",
            )
            + "<h4>Group filter \u2014 groups removed</h4>"
            + _reason_table_html(
                ds.get("group_filter", {}),
                detail_links.get("group_filter_drop_reasons"),
                "No groups were removed at the group filter step.",
            )
            + _divergence_filter_reasons_html(ds, detail_links, heading_tag="h4")
            + "<h4>Tables step \u2014 anomalies</h4>"
            + _reason_table_html(
                ds.get("tables", {}),
                detail_links.get("tables_drop_reasons"),
                "No anomalies at the table generation step.",
            )
        )
        species_html = (
            "<h4>Groups per species (surviving groups)</h4>"
            + _species_counts_table_html(species_tsv, "species-counts")
        )
        return '<div class="filtering-summary">%s%s%s</div>' % (
            funnel, reasons_html, species_html,
        )

    # Two-list mode: one column per list, plus an overlap summary on top.
    diff_stats = _load_json(differential_stats_path)
    if diff_stats:
        overlap_html = (
            '<div class="filtering-overlap"><h4>Overlap between lists '
            '(used for differential scoring)</h4>'
            '<div class="funnel-row">'
            '<div class="funnel-step"><div class="funnel-num">%d</div>'
            '<div class="funnel-label">%s groups</div></div>'
            '<div class="funnel-arrow">&cap;</div>'
            '<div class="funnel-step"><div class="funnel-num">%d</div>'
            '<div class="funnel-label">%s groups</div></div>'
            '<div class="funnel-arrow">=</div>'
            '<div class="funnel-step"><div class="funnel-num">%d</div>'
            '<div class="funnel-label">Common groups</div></div>'
            '<div class="funnel-step"><div class="funnel-num">%d</div>'
            '<div class="funnel-label">%s only</div></div>'
            '<div class="funnel-step"><div class="funnel-num">%d</div>'
            '<div class="funnel-label">%s only</div></div>'
            '</div></div>'
            % (diff_stats["n_groups_list1"], diff_stats["list1_name"],
               diff_stats["n_groups_list2"], diff_stats["list2_name"],
               diff_stats["n_common_groups"],
               diff_stats["n_list1_only"], diff_stats["list1_name"],
               diff_stats["n_list2_only"], diff_stats["list2_name"])
        )
    else:
        overlap_html = "<p>No overlap statistics available.</p>"

    col1 = _one_list_filtering_block(label1, pipeline_stats_path1, report_dir, "list1")
    col2 = _one_list_filtering_block(label2, pipeline_stats_path2, report_dir, "list2")
    return (
        '<div class="filtering-summary-two-list">%s'
        '<div class="filtering-cols">%s%s</div></div>'
        % (overlap_html, col1, col2)
    )


# ---------------------------------------------------------------------------
# Run parameters tab
# ---------------------------------------------------------------------------
#
# Presents the reproducibility record written by
# provenance.write_run_parameters() (see cli.py) — the effective,
# post-CLI-override parameter values actually used for this run, not the
# raw config file contents (a CLI flag may have silently overridden a
# config value; this tab shows what actually ran). Explicit (key, label)
# field lists are used rather than a generic dict-key-to-label
# transformation, since the underlying dict shape is small and fixed
# (mirrors provenance.collect_run_parameters()'s own contract), and
# hand-written labels read better than a mechanical snake_case conversion
# (e.g. "FDR threshold" rather than "Fdr threshold").

_RUN_INFO_FIELDS = [
    ("enhydra_version",     "ENHYDRA version"),
    ("timestamp",           "Run timestamp"),
    ("code_config_path",    "Code config file"),
    ("project_config_path", "Project config file"),
    ("outdir",              "Output directory"),
]

_INPUT_FIELDS = [
    ("inputdir",        "Input directory"),
    ("orthofinder_dir", "OrthoFinder directory"),
    ("anchor",          "Anchor species"),
    ("two_list_mode",   "Two-list differential mode"),
]

_LIST_FIELDS = [
    ("path",            "Species list file"),
    ("name",            "Display name"),
    ("n_species",       "Number of species"),
    ("anchor_injected", "Anchor auto-injected into this list"),
]

_FILTERING_FIELDS = [
    ("min_species",          "Minimum species per group"),
    ("min_sequences",        "Minimum sequences per group"),
    ("paralogs",             "Paralog handling"),
    ("length_filter_sd",     "Length filter (SD threshold)"),
    ("divergence_filter_sd", "Divergence filter (SD threshold)"),
    ("trim",                 "Column trimming mode"),
]

_ALIGNMENT_FIELDS = [
    ("aligner",    "Aligner"),
    ("mafft_mode", "MAFFT mode"),
]

_RANKING_FIELDS = [
    ("all_metrics", "All metrics run (--all-metrics)"),
    ("metrics_run", "Metric(s) run"),
]

_GSEA_FIELDS = [
    ("permutations",  "Permutations"),
    ("min_size",      "Minimum gene set size"),
    ("max_size",      "Maximum gene set size"),
    ("seed",          "Random seed"),
    ("fdr_threshold", "FDR threshold"),
    ("top_n",         "Top N gene sets (plots)"),
]

_GENE_SETS_FIELDS = [
    ("gene_sets_path", "Local GMT file"),
    ("organism",       "g:Profiler organism"),
    ("sources",        "g:Profiler data sources"),
]


def get_run_parameter_field_specs() -> dict[str, list[tuple[str, str]]]:
    """Return the (key, label) field specs for each run_parameters.json
    section, grouped by section name.

    Single source of truth for field ordering/labelling, shared by this
    module's own Run parameters tab (_build_run_parameters_tab_content())
    and by compare_runs.py's cross-run parameter diff table, so the two
    presentations of the same underlying provenance record never drift
    out of sync with each other.

    Returns:
        Dict mapping section name to its ordered list of (key, label)
        pairs. 'run_info' fields live at the top level of a
        run_parameters.json record; 'list' fields are nested under
        input.list1 / input.list2 (two-list mode only) rather than under
        a top-level 'list' key — callers extract the relevant sub-dict
        themselves. All other section names match a top-level key in the
        record directly (e.g. record['filtering']).
    """
    return {
        "run_info":  _RUN_INFO_FIELDS,
        "input":     _INPUT_FIELDS,
        "list":      _LIST_FIELDS,
        "filtering": _FILTERING_FIELDS,
        "alignment": _ALIGNMENT_FIELDS,
        "ranking":   _RANKING_FIELDS,
        "gsea":      _GSEA_FIELDS,
        "gene_sets": _GENE_SETS_FIELDS,
    }


def format_provenance_value(value) -> str:
    """Format a single run_parameters.json value for HTML display.

    Public (not module-private) since it is reused outside this module by
    compare_runs.py's cross-run parameter diff table, in addition to this
    module's own _fields_table_html().
    """
    if value is None:
        return '<span style="color:#999;">not set</span>'
    if isinstance(value, bool):
        return "Yes" if value else "No"
    if isinstance(value, list):
        if not value:
            return '<span style="color:#999;">none</span>'
        return html.escape(", ".join(str(v) for v in value))
    return html.escape(str(value))


def _fields_table_html(d: dict, fields: list[tuple[str, str]]) -> str:
    rows = "".join(
        "<tr><td>%s</td><td>%s</td></tr>"
        % (html.escape(label), format_provenance_value(d.get(key)))
        for key, label in fields
    )
    return '<table class="reason-table"><tbody>%s</tbody></table>' % rows


def _build_run_parameters_tab_content(
    run_parameters_path: str | None,
    report_dir: str,
) -> str:
    """Build the Run parameters tab body from run_parameters.json.

    Args:
        run_parameters_path: Path to run_parameters.json, as written by
                             provenance.write_run_parameters(). None or a
                             missing file returns a placeholder message
                             rather than raising — a report may be
                             regenerated (e.g. via --replot on an output
                             directory) from a run predating this feature.
        report_dir:          Directory report.html itself will be written
                             into, used to compute a relative link to the
                             raw JSON file for download.

    Returns:
        HTML string.
    """
    if not run_parameters_path or not os.path.isfile(run_parameters_path):
        return (
            "<p>No run parameters record is available for this run "
            "(run_parameters.json not found \u2014 this report may have "
            "been generated by an older version of ENHYDRA).</p>"
        )

    with open(run_parameters_path) as fh:
        params = json.load(fh)

    rel_json = os.path.relpath(run_parameters_path, report_dir).replace(os.sep, "/")

    cmd = params.get("command_line") or []
    try:
        cmd_str = shlex.join(str(c) for c in cmd)
    except Exception:
        cmd_str = " ".join(str(c) for c in cmd)

    input_section = params.get("input", {}) or {}
    two_list      = bool(input_section.get("two_list_mode"))

    parts = []
    parts.append("<h4>Run info</h4>")
    parts.append(_fields_table_html(params, _RUN_INFO_FIELDS))

    parts.append("<h4>Command line</h4>")
    cmd_display = html.escape(cmd_str) if cmd_str else \
        '<span style="color:#999;">unavailable</span>'
    parts.append(
        '<table class="reason-table"><tbody><tr><td>Invocation</td>'
        '<td style="font-family:\'Courier New\',monospace;">%s</td>'
        "</tr></tbody></table>" % cmd_display
    )

    parts.append("<h4>Input</h4>")
    parts.append(_fields_table_html(input_section, _INPUT_FIELDS))
    if two_list:
        list1 = input_section.get("list1", {}) or {}
        list2 = input_section.get("list2", {}) or {}
        list1_heading = "List 1"
        if list1.get("name"):
            list1_heading += " \u2014 %s" % html.escape(str(list1["name"]))
        list2_heading = "List 2"
        if list2.get("name"):
            list2_heading += " \u2014 %s" % html.escape(str(list2["name"]))
        parts.append("<h5>%s</h5>" % list1_heading)
        parts.append(_fields_table_html(list1, _LIST_FIELDS))
        parts.append("<h5>%s</h5>" % list2_heading)
        parts.append(_fields_table_html(list2, _LIST_FIELDS))

    parts.append("<h4>Filtering</h4>")
    parts.append(_fields_table_html(params.get("filtering", {}) or {}, _FILTERING_FIELDS))

    parts.append("<h4>Alignment</h4>")
    parts.append(_fields_table_html(params.get("alignment", {}) or {}, _ALIGNMENT_FIELDS))

    parts.append("<h4>Ranking</h4>")
    parts.append(_fields_table_html(params.get("ranking", {}) or {}, _RANKING_FIELDS))

    parts.append("<h4>GSEA</h4>")
    parts.append(_fields_table_html(params.get("gsea", {}) or {}, _GSEA_FIELDS))

    parts.append("<h4>Gene sets</h4>")
    parts.append(_fields_table_html(params.get("gene_sets", {}) or {}, _GENE_SETS_FIELDS))

    parts.append(
        '<p style="margin-top:18px;"><a href="%s" target="_blank" '
        'class="detail-link">Download raw run_parameters.json</a></p>' % rel_json
    )

    return "".join(parts)


# ---------------------------------------------------------------------------
# Alignments tab (GO-term-nested alignment page browser)
# ---------------------------------------------------------------------------

def _load_group_to_genes(tables_dir: str) -> dict[str, list[str]]:
    """Load group2anchor.tsv as {group_id: [anchor_gene_id, ...]}.

    Unlike msa_viewer.load_group_anchor() (which keeps only the last
    gene_id per group_id, sufficient for single-line page metadata), this
    preserves every row for a group — a group can map to more than one
    anchor gene under --paralogs all, and GO-term membership needs to be
    checked against all of them, not just whichever row happened to be
    read last.
    """
    path = os.path.join(tables_dir, "group2anchor.tsv")
    result: dict[str, list[str]] = {}
    if not os.path.isfile(path):
        return result
    with open(path) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                continue
            result.setdefault(fields[0], []).append(fields[1])
    return result


def _load_group_identity(tables_dir: str | None) -> dict[str, float]:
    """Load group2mean.tsv as {group_id: mean_identity}."""
    if not tables_dir:
        return {}
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


def _build_alignment_tree_html(
    gmt_gene_sets: dict[str, list[str]],
    tables_dir1: str | None,
    obo_names: dict[str, str],
    alignment_pages1: dict[str, str] | None,
    alignment_pages2: dict[str, str] | None,
    report_dir: str,
    label1: str = "List 1",
    label2: str = "List 2",
    tables_dir2: str | None = None,
) -> str:
    """Build the Alignments tab body: a GO-term-nested tree of alignment pages.

    Membership is determined by list1's (or the single list's)
    group2anchor.tsv — the only source of anchor-gene -> GMT mapping that
    exists — intersected with whichever groups actually received a
    rendered page in alignment_pages1. In two-list mode, alignment_pages1
    is already restricted by cli.py to the group_id column of
    differential_scores.tsv (i.e. bucket 3: groups that actually fed
    GSEA), so no separate intersection check against that file is needed
    here.

    A GO term is only listed if at least one of its genes maps to a group
    with a rendered page; terms are sorted by matched-group count
    descending. If a group has more than one anchor gene under
    --paralogs all and more than one of those genes falls under the same
    term, the group is still listed only once for that term (paired with
    the first gene, in sorted order, that matched) — this is a documented
    simplification to avoid confusing duplicate rows for one group under
    a single term.

    In two-list mode (tables_dir2 given), each row shows both lists'
    mean identity for that group (from each list's own group2mean.tsv)
    plus their signed difference (list1 - list2). Rows within a term are
    sorted by that signed difference descending: groups where list1 is
    more conserved / list2 is more variable (positive diff) appear
    first, sliding down through zero to groups where list2 is more
    conserved / list1 is more variable (negative diff) at the bottom.
    Deliberately signed rather than by |diff| — the direction of the
    difference is the point, not just its magnitude. Rows missing a
    score in either list sort to the end regardless of sign. In
    single-list mode, each row shows a single plain identity value, and
    rows are sorted by group_id (there is nothing to diff against).

    Args:
        gmt_gene_sets:     {term_id: [gene_id, ...]}, as returned by
                          _gmt_gene_sets().
        tables_dir1:       list1's (or the single list's) tables/
                          directory. None disables the tab content
                          (returns an explanatory message) rather than
                          raising.
        obo_names:         {term_id: term_name}, for display alongside
                          each GO ID.
        alignment_pages1:  {group_id: absolute_path}, as returned by
                          build_alignment_pages() for the single list /
                          list1. None or empty disables the tab content.
        alignment_pages2:  {group_id: absolute_path} for list2, or None
                          in single-list mode. When given, each leaf
                          shows a second link labelled with label2,
                          shown only for groups that also have a list2
                          page.
        report_dir:        Directory report.html itself will be written
                          into, used to convert the absolute paths in
                          alignment_pages1/2 into report-relative links
                          (same convention as enrichment plot linking).
        label1:            Display label for the first link/score column
                          (the list name in two-list mode). In
                          single-list mode (alignment_pages2 is None)
                          the link label falls back to plain "View
                          alignment", since naming it after a list would
                          be meaningless with only one list.
        label2:            Display label for the second link/score
                          column (two-list mode only).
        tables_dir2:       list2's tables/ directory, used to load its
                          own group2mean.tsv for the per-row identity
                          comparison. None in single-list mode, or if
                          two-list identity comparison should be skipped
                          (falls back to the single plain-identity
                          display even if alignment_pages2 is given).

    Returns:
        HTML string: a filter input plus the nested <details> tree, or a
        placeholder message if no data is available to build it from.
    """
    if not tables_dir1 or not alignment_pages1:
        return "<p>No alignment pages were generated for this run.</p>"

    group_to_genes = _load_group_to_genes(tables_dir1)
    identity1 = _load_group_identity(tables_dir1)
    identity2 = _load_group_identity(tables_dir2) if tables_dir2 else {}

    gene_to_groups: dict[str, list[str]] = {}
    for gid, genes in group_to_genes.items():
        for gene in genes:
            gene_to_groups.setdefault(gene, []).append(gid)

    rel_pages1 = {
        gid: os.path.relpath(p, report_dir).replace(os.sep, "/")
        for gid, p in alignment_pages1.items()
    }
    rel_pages2 = (
        {gid: os.path.relpath(p, report_dir).replace(os.sep, "/")
         for gid, p in alignment_pages2.items()}
        if alignment_pages2 else {}
    )
    two_list = bool(alignment_pages2) and tables_dir2 is not None

    term_entries: dict[str, list[dict]] = {}
    for term_id, genes in gmt_gene_sets.items():
        seen_groups: set[str] = set()
        entries: list[dict] = []
        for gene in sorted(genes):
            for group_id in gene_to_groups.get(gene, []):
                if group_id in seen_groups or group_id not in rel_pages1:
                    continue
                seen_groups.add(group_id)
                id1 = identity1.get(group_id)
                id2 = identity2.get(group_id) if two_list else None
                diff = (id1 - id2) if (id1 is not None and id2 is not None) else None
                entries.append({
                    "group_id": group_id, "gene_id": gene,
                    "id1": id1, "id2": id2, "diff": diff,
                })
        if entries:
            if two_list:
                entries.sort(key=lambda e: (e["diff"] is None,
                                            -e["diff"] if e["diff"] is not None else 0))
            else:
                entries.sort(key=lambda e: e["group_id"])
            term_entries[term_id] = entries

    if not term_entries:
        return "<p>No GO terms could be matched to rendered alignment pages.</p>"

    sorted_terms = sorted(term_entries.items(), key=lambda kv: -len(kv[1]))

    blocks = []
    for term_id, entries in sorted_terms:
        term_name = obo_names.get(term_id, "")
        summary_label = "%s%s (%d group%s)" % (
            html.escape(term_id),
            " \u2014 %s" % html.escape(term_name) if term_name else "",
            len(entries), "" if len(entries) == 1 else "s",
        )
        rows_html = []
        for e in entries:
            group_id, gene_id = e["group_id"], e["gene_id"]
            link1_label = html.escape(label1) if two_list else "View alignment"
            links = '<a href="%s" target="_blank" class="aln-link">%s</a>' % (
                rel_pages1[group_id], link1_label,
            )
            if two_list and group_id in rel_pages2:
                links += ' <a href="%s" target="_blank" class="aln-link">%s</a>' % (
                    rel_pages2[group_id], html.escape(label2),
                )

            if two_list:
                id1_str = "%.4f" % e["id1"] if e["id1"] is not None else "N/A"
                id2_str = "%.4f" % e["id2"] if e["id2"] is not None else "N/A"
                diff_str = "%+.4f" % e["diff"] if e["diff"] is not None else "N/A"
                score_html = (
                    '<span class="aln-tree-identity">%s: %s</span>'
                    '<span class="aln-tree-identity">%s: %s</span>'
                    '<span class="aln-tree-diff">\u0394: %s</span>'
                    % (html.escape(label1), id1_str,
                       html.escape(label2), id2_str, diff_str)
                )
            else:
                id1_str = "%.4f" % e["id1"] if e["id1"] is not None else "N/A"
                score_html = '<span class="aln-tree-identity">identity: %s</span>' % id1_str

            search_key = html.escape(
                ("%s %s %s %s" % (group_id, gene_id, term_id, term_name)).lower()
            )
            rows_html.append(
                '<div class="aln-tree-row" data-search="%s">'
                '<span class="aln-tree-group">%s</span>'
                '<span class="aln-tree-gene">anchor: %s</span>'
                '%s'
                '<span class="aln-tree-links">%s</span>'
                "</div>"
                % (search_key, html.escape(display_group_id(group_id)),
                   html.escape(gene_id), score_html, links)
            )
        term_search_key = html.escape(("%s %s" % (term_id, term_name)).lower())
        blocks.append(
            '<details class="aln-tree-term" data-term-search="%s">'
            "<summary>%s</summary>%s</details>"
            % (term_search_key, summary_label, "".join(rows_html))
        )

    return (
        '<div class="aln-tree-controls">'
        '<input type="text" id="aln-tree-filter" '
        'placeholder="Filter by GO ID, term name, group ID, or gene ID..."/>'
        "</div>"
        '<div id="aln-tree">%s</div>' % "".join(blocks)
    )


def _build_gene_to_page_links(
    tables_dir1: str | None,
    alignment_pages1: dict[str, str] | None,
    alignment_pages2: dict[str, str] | None,
    report_dir: str,
) -> dict[str, tuple[str | None, str | None]]:
    """Build {gene_id: (list1_relpath_or_None, list2_relpath_or_None)}.

    Used to turn plain gene IDs in the "Full gene set" / "Leading edge"
    modals into links to their pre-rendered alignment page, alongside the
    GO-term tree (Alignments tab) which serves the same underlying data
    but organised by GO term rather than by gene.

    If a gene's group (from list1's group2anchor.tsv) has no rendered
    page in alignment_pages1, that gene is simply absent from the
    returned dict — callers should render it as plain, non-linked text.
    If a gene happens to map to more than one group (possible under
    --paralogs all, though unusual), the first matching group encountered
    is used — a documented simplification, consistent with the same
    choice already made in _build_alignment_tree_html().

    Args:
        tables_dir1:       list1's (or the single list's) tables/
                          directory. None returns an empty dict.
        alignment_pages1:  {group_id: absolute_path} for the single list
                          or list1. None or empty returns an empty dict.
        alignment_pages2:  {group_id: absolute_path} for list2, or None
                          in single-list mode.
        report_dir:        Directory report.html will be written into,
                          used to compute relative link paths.

    Returns:
        Dict mapping gene_id to a (list1_link, list2_link) tuple, where
        each element is a report-relative path string or None if that
        list has no page for this gene's group.
    """
    if not tables_dir1 or not alignment_pages1:
        return {}

    group_to_genes = _load_group_to_genes(tables_dir1)
    rel_pages1 = {
        gid: os.path.relpath(p, report_dir).replace(os.sep, "/")
        for gid, p in alignment_pages1.items()
    }
    rel_pages2 = (
        {gid: os.path.relpath(p, report_dir).replace(os.sep, "/")
         for gid, p in alignment_pages2.items()}
        if alignment_pages2 else {}
    )

    gene_to_page: dict[str, tuple[str | None, str | None]] = {}
    for group_id, genes in group_to_genes.items():
        if group_id not in rel_pages1:
            continue
        link1 = rel_pages1[group_id]
        link2 = rel_pages2.get(group_id)
        for gene in genes:
            if gene not in gene_to_page:
                gene_to_page[gene] = (link1, link2)
    return gene_to_page


def _gene_chip_list_html(
    gene_ids: list[str],
    gene_to_page: dict[str, tuple[str | None, str | None]],
) -> str:
    """Render a list of gene IDs as chips, linked to alignment pages where
    a mapping exists in gene_to_page, plain text otherwise.

    Each linked gene shows small 'L1'/'L2' badge links (rather than full
    list-name labels, which would be too cramped at chip scale) opening
    that list's alignment page in a new tab. A gene with only a list1
    link (single-list mode, or a two-list gene whose group has no list2
    page) shows only the L1 badge.
    """
    chips = []
    for gene in gene_ids:
        link1, link2 = gene_to_page.get(gene, (None, None))
        badges = ""
        if link1:
            badges += (
                '<a href="%s" target="_blank" class="gene-chip-badge" '
                'title="View alignment">L1</a>' % link1
            )
        if link2:
            badges += (
                '<a href="%s" target="_blank" class="gene-chip-badge" '
                'title="View alignment">L2</a>' % link2
            )
        chips.append(
            '<span class="gene-chip">%s%s</span>'
            % (html.escape(gene), (" " + badges) if badges else "")
        )
    return '<div class="gene-chip-list">%s</div>' % "".join(chips)


# ---------------------------------------------------------------------------
# Internal helpers (GMT / GSEA results / enrichment table)
# ---------------------------------------------------------------------------

def _gmt_term_names(gmt_path: str | None) -> dict[str, str]:
    if not gmt_path or not os.path.isfile(gmt_path):
        return {}
    names: dict[str, str] = {}
    with open(gmt_path) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 2 and fields[0] and fields[1]:
                names[fields[0]] = fields[1]
    return names


def _gmt_gene_sets(gmt_path: str | None) -> dict[str, list[str]]:
    """Return term_id -> full list of member gene IDs, read from a GMT file."""
    if not gmt_path or not os.path.isfile(gmt_path):
        return {}
    sets: dict[str, list[str]] = {}
    with open(gmt_path) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 3 and fields[0]:
                sets[fields[0]] = [g for g in fields[2:] if g]
    return sets


def _find_gmt_in_dir(directory: str) -> str | None:
    if not os.path.isdir(directory):
        return None
    for fname in os.listdir(directory):
        if fname.endswith(".gmt"):
            return os.path.join(directory, fname)
    return None


def _resolve_term_names(
    results_dir: str,
    obo_path: str | None,
    gmt_path: str | None,
) -> dict[str, str]:
    effective_gmt = gmt_path or _find_gmt_in_dir(results_dir)
    gmt_names     = _gmt_term_names(effective_gmt)
    if effective_gmt and not gmt_names:
        logger.warning("GMT at '%s' had no parseable term names.", effective_gmt)
    elif not effective_gmt:
        logger.warning("No GMT file found; term names will show as GO IDs.")
    else:
        logger.info("Loaded %d term names from GMT.", len(gmt_names))
    obo_names: dict[str, str] = {}
    if obo_path and os.path.isfile(obo_path):
        obo_names = _parse_obo_names(obo_path)
        logger.info("Loaded %d GO term names from OBO.", len(obo_names))
    return {**gmt_names, **obo_names}


def _normalise_series(scores: pd.Series, metric: str) -> pd.Series:
    if metric == "identity":
        return scores
    elif metric == "zscore":
        return (scores - scores.mean()) / scores.std()
    elif metric == "rank":
        return scores.rank(ascending=True) / len(scores)
    return scores


def _per_term_scores(
    gmt_path: str,
    tables_dir1: str,
    tables_dir2: str | None,
    term_ids: list[str],
    metric: str = "identity",
) -> pd.DataFrame:
    """Compute mean metric scores per GO term.

    Single-list mode (tables_dir2=None) → "Mean score" column.
    Two-list mode → "List 1 score", "List 2 score", "Score diff".
    """
    def _load_anchor2mean(tables_dir: str) -> dict[str, float]:
        path = os.path.join(tables_dir, "anchor2mean.tsv")
        if not os.path.isfile(path):
            return {}
        df = pd.read_csv(path, sep="\t", header=None,
                         names=["gene_id", "score"],
                         dtype={"gene_id": str})
        df["score"] = pd.to_numeric(df["score"], errors="coerce")
        df = df.dropna(subset=["score"]).drop_duplicates("gene_id")
        return dict(_normalise_series(df.set_index("gene_id")["score"], metric))

    def _load_via_group_mapping(tables_dir_scores: str,
                                tables_dir_mapping: str) -> dict[str, float]:
        g2m_path = os.path.join(tables_dir_scores,  "group2mean.tsv")
        g2a_path = os.path.join(tables_dir_mapping, "group2anchor.tsv")
        if not os.path.isfile(g2m_path) or not os.path.isfile(g2a_path):
            return {}
        g2m = pd.read_csv(g2m_path, sep="\t", header=None,
                          names=["group_id", "score"],
                          dtype={"group_id": str})
        g2m["score"] = pd.to_numeric(g2m["score"], errors="coerce")
        g2m = (g2m.dropna(subset=["score"])
                  .drop_duplicates("group_id")
                  .set_index("group_id")["score"])
        g2m = _normalise_series(g2m, metric)

        g2a = pd.read_csv(g2a_path, sep="\t", header=None,
                          names=["group_id", "gene_id"],
                          dtype={"group_id": str, "gene_id": str})
        g2a = (g2a.dropna()
                  .drop_duplicates("group_id")
                  .set_index("group_id")["gene_id"])

        merged = (g2m.rename("score")
                     .reset_index()
                     .merge(g2a.reset_index(), on="group_id", how="inner")
                     .dropna(subset=["gene_id", "score"]))
        return dict(zip(merged["gene_id"], merged["score"]))

    scores1 = _load_anchor2mean(tables_dir1)
    if tables_dir2 is not None:
        scores2 = _load_anchor2mean(tables_dir2)
        if not scores2:
            scores2 = _load_via_group_mapping(tables_dir2, tables_dir1)
    else:
        scores2 = None

    if not scores1:
        logger.warning("No anchor scores loaded from %s.", tables_dir1)
        return pd.DataFrame()

    term_id_set = set(term_ids)
    gmt: dict[str, set[str]] = {}
    with open(gmt_path) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 3 and fields[0] in term_id_set:
                gmt[fields[0]] = set(fields[2:])

    rows = []
    for tid in term_ids:
        genes = gmt.get(tid, set())
        s1    = [scores1[g] for g in genes if g in scores1]
        mean1 = round(sum(s1) / len(s1), 4) if s1 else None
        if scores2 is None:
            rows.append({"Term": tid, "Mean score": mean1})
        else:
            s2    = [scores2[g] for g in genes if g in scores2]
            mean2 = round(sum(s2) / len(s2), 4) if s2 else None
            diff  = (round(mean1 - mean2, 4)
                     if mean1 is not None and mean2 is not None else None)
            rows.append({"Term": tid,
                         "List 1 score": mean1,
                         "List 2 score": mean2,
                         "Score diff":   diff})
    return pd.DataFrame(rows).set_index("Term")


def _fetch(url: str) -> str:
    ctx = ssl.create_default_context()
    ctx.check_hostname = False
    ctx.verify_mode    = ssl.CERT_NONE
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
    df.to_csv(os.path.join(results_dir, "gsea_results_processed.tsv"),
              sep="\t", index=False)
    return df


def _build_enrichment_plot_index(results_dir: str, report_dir: str) -> dict[str, str]:
    """Build a mapping of GO ID -> image path for per-gene-set enrichment plots.

    Records a path relative to the report's own directory rather than
    base64-embedding every PNG, so the browser loads each plot from disk on
    demand and the HTML stays small even with hundreds of significant sets.
    """
    prerank_dir = os.path.join(results_dir, "prerank")
    if not os.path.isdir(prerank_dir):
        return {}
    index = {}
    for filename in os.listdir(prerank_dir):
        if not filename.endswith(".png"):
            continue
        go_id    = filename.replace(".png", "").replace("_", ":", 1)
        abs_path = os.path.abspath(os.path.join(prerank_dir, filename))
        rel_path = os.path.relpath(abs_path, start=report_dir)
        index[go_id] = rel_path.replace(os.sep, "/")
    logger.info(
        "Indexed %d enrichment plot(s) as on-disk links (not embedded).",
        len(index),
    )
    return index


def _augment_with_per_term_scores(
    df: pd.DataFrame,
    gmt_path: str | None,
    tables_dir1: str | None,
    tables_dir2: str | None,
    metric: str,
) -> pd.DataFrame:
    if not (gmt_path and tables_dir1):
        return df
    if not os.path.isfile(gmt_path):
        return df
    try:
        per_term = _per_term_scores(
            gmt_path=gmt_path,
            tables_dir1=tables_dir1,
            tables_dir2=tables_dir2,
            term_ids=df["Term"].tolist(),
            metric=metric,
        )
        if per_term.empty:
            return df
        df = df.set_index("Term").join(per_term).reset_index()
    except Exception:
        import traceback
        logger.warning("Could not compute per-term scores:\n%s",
                       traceback.format_exc())
    return df


_LEAD_GENES_COL_ALIASES = ("lead_genes", "lead genes", "leading_edge_genes")


def _find_lead_genes_col(df: pd.DataFrame) -> str | None:
    for c in df.columns:
        if c.strip().lower() in _LEAD_GENES_COL_ALIASES:
            return c
    return None


def _results_table_html(
    df: pd.DataFrame,
    obo_names: dict[str, str],
    plot_index: dict[str, str],
    fdr_threshold: float = 0.25,
    metric: str | None = None,
    col1_label: str = "List 1",
    col2_label: str = "List 2",
    gene_sets: dict[str, list[str]] | None = None,
    gene_to_page: dict[str, tuple[str | None, str | None]] | None = None,
) -> tuple[str, list[int], dict[str, str], dict[str, str], dict[str, str], dict[str, str]]:
    """Build the enrichment results DataTable HTML for one metric.

    Returns a 6-tuple:
        (table_html, numeric_col_indices,
         full_gene_sets_text, leading_edge_text,
         full_gene_sets_html, leading_edge_html)

    The *_text dicts are unchanged plain newline-joined gene lists, used
    for the Excel export (tableToXLSX() splits these on newline — turning
    them into HTML would corrupt the exported cells). The *_html dicts are
    a parallel rendering of the same gene lists as clickable chips (via
    _gene_chip_list_html()), used only for the on-screen modal display.
    """
    df = df.copy()
    if "Term" not in df.columns:
        logger.warning("'Term' column not found in GSEA results.")
        return "<p>No results to display.</p>", [], {}, {}, {}, {}

    df["GO Term"] = df["Term"].map(obo_names).fillna(df["Term"])

    for col in ["ES", "NES", "NOM p-val", "FDR q-val", "FWER p-val",
                "Tag %", "Gene %", "Mean score",
                "List 1 score", "List 2 score", "Score diff"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce").apply(
                lambda x: "%.4f" % x if pd.notna(x) else ""
            )

    df["Significant"] = df["FDR q-val"].apply(
        lambda x: "✓" if x != "" and float(x) < fdr_threshold else ""
    )

    gene_to_page = gene_to_page or {}
    full_sets_js: dict[str, str] = {}
    leadedge_js:  dict[str, str] = {}
    full_sets_html_js: dict[str, str] = {}
    leadedge_html_js:  dict[str, str] = {}

    if gene_sets:
        def _full_set_cell(term_id):
            genes = gene_sets.get(term_id)
            if not genes:
                return ""
            full_sets_js[term_id] = "\n".join(genes)
            full_sets_html_js[term_id] = _gene_chip_list_html(genes, gene_to_page)
            return "View (%d)" % len(genes)
        df["Full gene set"] = df["Term"].apply(_full_set_cell)

    lead_col = _find_lead_genes_col(df)
    if lead_col:
        def _leadedge_cell(row):
            raw = row[lead_col]
            if pd.isna(raw) or not str(raw).strip():
                return ""
            genes = [g.strip() for g in str(raw).split(";") if g.strip()]
            if not genes:
                return ""
            leadedge_js[row["Term"]] = "\n".join(genes)
            leadedge_html_js[row["Term"]] = _gene_chip_list_html(genes, gene_to_page)
            return "View (%d)" % len(genes)
        df["Leading edge"] = df.apply(_leadedge_cell, axis=1)

    diff_label = "%s \u2212 %s" % (col1_label, col2_label)
    col_defs = [
        ("Term",         "GO ID",
         "Gene Ontology term identifier."),
        ("GO Term",      "Term name",
         "Human-readable name of the GO term."),
        ("Mean score",   "Mean score",
         "Mean metric score for genes in this gene set."),
        ("List 1 score", "%s score" % col1_label,
         "Mean metric score for genes in this set in %s." % col1_label),
        ("List 2 score", "%s score" % col2_label,
         "Mean metric score for genes in this set in %s." % col2_label),
        ("Score diff",   diff_label,
         "Difference in mean score (%s \u2212 %s)." % (col1_label, col2_label)),
        ("NES",          "NES",
         "Normalised Enrichment Score."),
        ("NOM p-val",    "p-value",
         "Nominal p-value from permutation testing."),
        ("FDR q-val",    "FDR",
         "False Discovery Rate q-value. Significant below %.2f." % fdr_threshold),
        ("Tag %",        "Tag %",
         "Fraction of gene set genes in the leading edge (0-1)."),
        ("Gene %",       "Gene %",
         "Fraction of all ranked genes in the leading edge (0-1)."),
        ("Full gene set", "Full gene set",
         "All genes annotated to this gene set. Click to view as a text list."),
        ("Leading edge", "Leading edge",
         "Genes from this set found in the leading edge of the ranked list "
         "(i.e. driving the enrichment score). Click to view as a text list."),
        ("Significant",  "Sig.",
         "Significant at FDR < %.2f." % fdr_threshold),
    ]
    col_defs = [(s, d, t) for s, d, t in col_defs
                if s in df.columns or s in ("GO Term", "Significant")]

    table_df      = df[[s for s, _, _ in col_defs]].copy()
    display_names = [d for _, d, _ in col_defs]
    tooltips      = [t for _, _, t in col_defs]
    table_df.columns = display_names

    numeric_names       = {"NES", "p-value", "FDR", "Tag %", "Gene %",
                           "Mean score", diff_label,
                           "%s score" % col1_label, "%s score" % col2_label}
    numeric_col_indices = [i for i, n in enumerate(display_names)
                           if n in numeric_names]

    header_cells = "".join(
        '<th>%s <span class="col-tip">?<span class="tip-text">%s</span></span></th>'
        % (name, tip)
        for name, tip in zip(display_names, tooltips)
    )
    filter_cells = "<th></th>" * len(display_names)
    table_id     = ("results-table-%s" % metric) if metric else "results-table"
    metric_attr  = (' data-metric="%s"' % metric) if metric else ""

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
            elif col == "Full gene set" and val:
                cells += (
                    '<td><a href="#" class="geneset-link" data-goid="%s">%s</a></td>'
                    % (go_id, val)
                )
            elif col == "Leading edge" and val:
                cells += (
                    '<td><a href="#" class="leadedge-link" data-goid="%s"%s>%s</a></td>'
                    % (go_id, metric_attr, val)
                )
            else:
                cells += "<td>%s</td>" % str(val)
        rows += "<tr%s>%s</tr>\n" % (sig_class, cells)

    html_content = (
        '<table id="%s" class="display compact" style="width:100%%">'
        '<thead><tr>%s</tr><tr class="filter-row">%s</tr></thead>'
        "<tbody>%s</tbody></table>"
    ) % (table_id, header_cells, filter_cells, rows)

    return (html_content, numeric_col_indices,
            full_sets_js, leadedge_js,
            full_sets_html_js, leadedge_html_js)


def build_significance_agreement_table_html(
    column_dfs: dict[str, pd.DataFrame],
    column_labels: dict[str, str],
    obo_names: dict[str, str],
    fdr_threshold: float,
    enrichment_plots_map: dict[str, dict[str, str]] | None = None,
    table_id: str = "results-table-consensus",
    direction_labels: dict[int, str] | None = None,
) -> tuple[str, list[int] | None]:
    """Build a cross-column significance agreement table.

    Originally written specifically for this module's own Cross-metric
    consensus tab (comparing identity/zscore/rank rankings within one
    run), this is now a general-purpose utility: a "column" is any axis
    of comparison whose GSEA output can be expressed as a
    {Term, NES, FDR q-val} DataFrame. That covers both the original use
    (columns = ranking metrics within one run) and a new one
    (columns = different runs, for the same metric — see
    compare_runs.py's cross-run comparison report), without either caller
    needing its own copy of this logic.

    A gene set qualifies for inclusion if it is significant
    (FDR < fdr_threshold) in a *majority* of the columns given — for 3
    columns this means 2 of 3, not requiring unanimous agreement. This
    directly matches the intent of a consensus view: surfacing results
    that hold up under more than one arbitrary choice of comparison axis,
    without demanding unrealistic universal agreement.

    Beyond the raw significance count, each qualifying row's direction is
    also checked: among the columns where it is significant, do they
    agree on the *sign* of NES (positive = more conserved, negative =
    faster evolving)? A gene set significant in two columns with
    opposite-signed NES is not a genuine consensus finding — it is
    contradictory — so this is surfaced explicitly as 'Mixed direction'
    rather than silently lumped in with genuine agreement.

    Deliberately excludes the "Full gene set" / "Leading edge" columns
    present in the per-metric enrichment tables: this view is meant to be
    a compact cross-check, and that level of per-gene detail remains
    available in whichever original per-column report produced each
    column's own data.

    Args:
        column_dfs:  {column_key: DataFrame} with at least 'Term', 'NES',
                    and 'FDR q-val' columns (numeric dtypes, not yet
                    stringified for display). Only columns with usable
                    data are expected to be included; callers are
                    responsible for that filtering.
        column_labels: {column_key: display_label} — e.g. metric display
                    names (METRIC_LABELS) or run display names.
        obo_names:   {term_id: term_name}, for display alongside each GO ID.
        fdr_threshold: Significance cutoff applied per column.
        enrichment_plots_map: {column_key: {go_id: plot_path}}, used to make
                    a qualifying gene set's per-column NES cell clickable
                    (emits a '.go-link' anchor with data-goid and
                    data-metric="<column_key>" attributes; the consuming
                    page's own JS is responsible for wiring that click to
                    an actual plot viewer — see this module's own
                    _TEMPLATE JS for one such wiring, keyed by metric, and
                    compare_runs.py for a second, independent wiring keyed
                    by (table_id, run_name)). None or an empty dict
                    disables clickable NES cells entirely (plain text
                    values only).
        table_id:    HTML id for the rendered <table>, so more than one
                    such table can appear on the same page without id
                    collisions (e.g. one table per metric in
                    compare_runs.py's Significance agreement tab).
        direction_labels: Optional {1: label_for_positive_sign,
                    -1: label_for_negative_sign} overriding the default
                    'Conserved' / 'Faster-evolving' wording, which is
                    only correct when NES sign means "more/less conserved
                    than the genome-wide average" (single-list mode).
                    Two-list/differential callers must override this — a
                    positive NES there means genes are more conserved in
                    one specific list relative to the other (see
                    differential.compute_differential()'s own docstring),
                    not conserved relative to a shared baseline, so both
                    label values need to name a *list*, not describe an
                    absolute rate of evolution — e.g.
                    {1: "Faster-evolving in <list2>",
                     -1: "Faster-evolving in <list1>"}. Values are
                    inserted as-is (not escaped by this function), the
                    same convention already used for column_labels above.
                    Defaults to {1: "Conserved", -1: "Faster-evolving"}
                    if not given.

    Returns:
        Tuple of (html, numeric_col_indices). If no gene set reaches
        majority significance, html is an explanatory <p> message and
        numeric_col_indices is None — callers should treat a None second
        element as "do not register this table with the DataTables JS",
        since no <table> element is present in that case.
    """
    columns   = list(column_dfs.keys())
    n_columns = len(columns)
    majority  = n_columns // 2 + 1
    enrichment_plots_map = enrichment_plots_map or {}
    _direction_labels = direction_labels or {1: "Conserved", -1: "Faster-evolving"}

    term_data: dict[str, dict[str, tuple[float, float]]] = {}
    for col, df in column_dfs.items():
        for _, row in df.iterrows():
            term = row.get("Term")
            nes, fdr = row.get("NES"), row.get("FDR q-val")
            if term is None or pd.isna(nes) or pd.isna(fdr):
                continue
            term_data.setdefault(term, {})[col] = (float(nes), float(fdr))

    rows_data: list[tuple[str, dict, int, str]] = []
    for term, per_col in term_data.items():
        sig_cols = [
            c for c in columns
            if c in per_col and per_col[c][1] < fdr_threshold
        ]
        n_sig = len(sig_cols)
        if n_sig < majority:
            continue
        signs = {(1 if per_col[c][0] > 0 else -1) for c in sig_cols}
        if len(signs) == 1:
            direction = _direction_labels[next(iter(signs))]
        else:
            direction = "Mixed direction"
        rows_data.append((term, per_col, n_sig, direction))

    if not rows_data:
        return (
            "<p>No gene sets reached majority significance "
            "(significant in \u2265 %d of %d columns at "
            "FDR &lt; %.2f).</p>"
            % (majority, n_columns, fdr_threshold)
        ), None

    # Full agreement first, then by term ID for a stable secondary order.
    rows_data.sort(key=lambda r: (-r[2], r[0]))

    headers  = ["GO ID", "Term name"]
    tooltips = [
        "Gene Ontology term identifier.",
        "Human-readable name of the GO term.",
    ]
    for c in columns:
        label = column_labels.get(c, str(c))
        headers += ["%s NES" % label, "%s FDR" % label]
        tooltips += [
            "Normalised Enrichment Score for %s." % label,
            "FDR q-value for %s (bold if significant, FDR < %.2f)."
            % (label, fdr_threshold),
        ]
    headers += ["Significant in", "Direction"]
    tooltips += [
        "Number of columns (out of %d) in which this gene set was "
        "significant." % n_columns,
        "Whether the direction of enrichment (%s vs. %s) agrees across "
        "the columns in which this gene set was significant."
        % (_direction_labels[1], _direction_labels[-1]),
    ]

    numeric_col_indices = [
        i for i, h in enumerate(headers)
        if h.endswith("NES") or h.endswith("FDR")
    ]

    header_cells = "".join(
        '<th>%s <span class="col-tip">?<span class="tip-text">%s</span></span></th>'
        % (h, t) for h, t in zip(headers, tooltips)
    )
    filter_cells = "<th></th>" * len(headers)

    body_rows = []
    for term, per_col, n_sig, direction in rows_data:
        term_name = obo_names.get(term, term)
        row_cells = "<td>%s</td><td>%s</td>" % (
            html.escape(term), html.escape(term_name),
        )
        for c in columns:
            if c in per_col:
                nes, fdr = per_col[c]
                is_sig  = fdr < fdr_threshold
                nes_str = "%.4f" % nes
                fdr_str = "%.4f" % fdr
                if enrichment_plots_map.get(c, {}).get(term):
                    nes_cell = (
                        '<a href="#" class="go-link" data-goid="%s" '
                        'data-metric="%s">%s</a>' % (term, c, nes_str)
                    )
                else:
                    nes_cell = nes_str
                fdr_cell = ("<strong>%s</strong>" % fdr_str) if is_sig else fdr_str
                row_cells += "<td>%s</td><td>%s</td>" % (nes_cell, fdr_cell)
            else:
                row_cells += "<td>\u2014</td><td>\u2014</td>"
        row_cells += "<td>%d / %d</td><td>%s</td>" % (n_sig, n_columns, direction)
        row_class = ' class="sig-row"' if n_sig == n_columns else ""
        body_rows.append("<tr%s>%s</tr>" % (row_class, row_cells))

    table_html = (
        '<table id="%s" class="display compact" style="width:100%%">'
        '<thead><tr>%s</tr><tr class="filter-row">%s</tr></thead>'
        "<tbody>%s</tbody></table>"
    ) % (table_id, header_cells, filter_cells, "".join(body_rows))

    return table_html, numeric_col_indices


def _plot_section(plots_dir: str, names: list[tuple[str, str]]) -> str:
    html_content = ""
    for stem, caption in names:
        svg_path = os.path.join(plots_dir, stem + ".svg")
        png_path = os.path.join(plots_dir, stem + ".png")
        if os.path.isfile(svg_path):
            with open(svg_path, encoding="utf-8") as fh:
                svg_content = fh.read()
            html_content += (
                '<div class="plot-block">'
                '<p class="plot-caption">%s (hover for details)</p>'
                '%s</div>'
            ) % (caption, svg_content)
        elif os.path.isfile(png_path):
            uri = _img_to_base64(png_path)
            html_content += (
                '<div class="plot-block">'
                '<p class="plot-caption">%s</p>'
                '<img src="%s" alt="%s"/></div>'
            ) % (caption, uri, caption)
    return html_content


# ---------------------------------------------------------------------------
# Shared report builder
# ---------------------------------------------------------------------------

def _build_report_impl(
    metric_data: dict,
    report_path: str,
    obo_path: str | None,
    fdr_threshold: float,
    mode: str,
    gmt_path: str | None,
    tables_dir1: str | None,
    tables_dir2: str | None,
    label1: str,
    label2: str,
    pipeline_stats_path: str | None,
    pipeline_stats_path1: str | None,
    pipeline_stats_path2: str | None,
    differential_stats_path: str | None,
    alignment_pages1: dict[str, str] | None = None,
    alignment_pages2: dict[str, str] | None = None,
    run_parameters_path: str | None = None,
):
    logger.info("Building HTML report (%d metric tab(s))...", len(metric_data))
    first_results = next(iter(metric_data.values()))["results_dir"] if metric_data else ""
    effective_gmt = gmt_path or _find_gmt_in_dir(first_results)
    term_names    = _resolve_term_names(first_results, obo_path, effective_gmt)
    gene_sets_all = _gmt_gene_sets(effective_gmt)
    cache_dir = os.path.dirname(obo_path) if obo_path else None
    jquery_js = _fetch_cached(_JQUERY_URL,         cache_dir, "jquery.min.js")
    dt_js     = _fetch_cached(_DATATABLES_JS_URL,  cache_dir, "datatables.min.js")
    dt_css    = _fetch_cached(_DATATABLES_CSS_URL, cache_dir, "datatables.min.css")
    xlsx_js   = _fetch_cached(_XLSX_JS_URL,        cache_dir, "xlsx.full.min.js")

    single_metric = len(metric_data) == 1
    if mode == "differential":
        title = ("ENHYDRA Differential Enrichment Report" if single_metric
                  else "ENHYDRA Multi-Metric Differential Enrichment Report")
        plot_names = [
            ("identity_scatter",          "Identity comparison between lists"),
            ("differential_distribution", "Differential conservation score distribution"),
            ("gsea_barplot",              "Top differentially enriched gene sets (NES)"),
        ]
    else:
        title = ("ENHYDRA Single-List Enrichment Report" if single_metric
                  else "ENHYDRA Multi-Metric Enrichment Report")
        plot_names = [
            ("identity_distribution", "Distribution of mean alignment identity"),
            ("gsea_barplot",          "Top enriched gene sets (NES)"),
        ]

    report_dir = os.path.dirname(os.path.abspath(report_path))

    # Built once (not per metric) — feeds the "Full gene set" / "Leading
    # edge" modal chip links below. Independent of which metric is being
    # rendered, since it only depends on group2anchor.tsv + the rendered
    # alignment page maps.
    gene_to_page = _build_gene_to_page_links(
        tables_dir1, alignment_pages1, alignment_pages2, report_dir,
    )

    tab_buttons_parts    = []
    tab_panels_parts     = []
    enrichment_plots_map = {}
    numeric_cols_map      = {}
    full_gene_sets_accum  = {}
    leading_edge_map      = {}
    full_gene_sets_html_accum = {}
    leading_edge_html_map     = {}
    # Raw (still-numeric) per-metric GSEA result frames, captured purely to
    # feed the cross-metric consensus tab below — kept separate from the
    # `df` variable inside the loop, which gets progressively mutated
    # in-place-by-rebinding (augmented with per-term scores, then stringified
    # inside _results_table_html()'s own internal copy) for that metric's
    # own results table.
    metric_raw_dfs: dict[str, pd.DataFrame] = {}
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
        plot_idx = _build_enrichment_plot_index(results_dir, report_dir)
        enrichment_plots_map[metric] = plot_idx
        df = _load_gsea_results(results_dir)
        if df is not None:
            metric_raw_dfs[metric] = df[["Term", "NES", "FDR q-val"]].copy()
            df = _augment_with_per_term_scores(
                df, effective_gmt, tables_dir1, tables_dir2, metric
            )
            (tbl_html, num_cols, full_sets_js, leadedge_js,
             full_sets_html_js, leadedge_html_js) = _results_table_html(
                df, term_names, plot_idx, fdr_threshold,
                metric=metric, col1_label=label1, col2_label=label2,
                gene_sets=gene_sets_all, gene_to_page=gene_to_page,
            )
            numeric_cols_map[metric] = num_cols
            full_gene_sets_accum.update(full_sets_js)
            leading_edge_map[metric] = leadedge_js
            full_gene_sets_html_accum.update(full_sets_html_js)
            leading_edge_html_map[metric] = leadedge_html_js
        else:
            tbl_html = "<p>No GSEA results found for this metric.</p>"
            numeric_cols_map[metric] = []
            leading_edge_map[metric] = {}
            leading_edge_html_map[metric] = {}
        plots_html = _plot_section(plots_dir, plot_names)
        desc       = "" if single_metric else _METRIC_DESCS.get(metric, "")
        desc_html  = ('<p class="metric-desc">%s</p>' % desc) if desc else ""
        tab_panels_parts.append(
            '<div id="tab-{m}" class="tab-panel{ac}" role="tabpanel">\n'
            '  {desc_html}\n'
            '  <h3>Plots</h3>\n'
            '  <div class="plot-grid">{plots}</div>\n'
            '  <h3>Enrichment results</h3>\n'
            '  <p>Significant gene sets (FDR&nbsp;&lt;&nbsp;{fdr}) highlighted '
            'in blue. Click a GO ID to view its enrichment plot, or "View" '
            'under Full&nbsp;gene&nbsp;set / Leading&nbsp;edge to see the '
            'gene lists, with links to alignment pages where available.</p>\n'
            '  <p>'
            '<button class="sig-toggle-btn" data-metric="{m}">Show only significant</button> '
            '<button class="export-btn export-xlsx-btn" data-metric="{m}">'
            '&#8681; Export table to Excel (.xlsx)</button> '
            '<button class="export-btn export-xlsx-all-btn" data-metric="{m}">'
            '&#8681; Export all rows (ignore filters)</button>'
            '</p>\n'
            '  {tbl}\n'
            '</div>\n'.format(
                m=metric, ac=active_cls, desc_html=desc_html,
                plots=plots_html, fdr=fdr_threshold, tbl=tbl_html,
            )
        )
        first = False

    # Cross-metric consensus tab — only meaningful with 2+ metrics that
    # each produced usable GSEA results (i.e. an --all-metrics run; a
    # single-metric report has nothing to be consensus about). Surfaces
    # gene sets significant in a majority of the ranking metrics run, so a
    # result can be trusted as more than an artifact of one particular
    # ranking choice (identity vs. zscore vs. rank). See
    # _build_consensus_table_html() for the full significance/direction
    # rules.
    if len(metric_raw_dfs) >= 2:
        column_labels = {
            m: METRIC_LABELS.get(m, m.capitalize()) for m in metric_raw_dfs
        }
        # In two-list/differential mode, NES sign reflects which list a
        # gene set is relatively more conserved in (see
        # differential.compute_differential()'s own docstring: a positive
        # differential score means higher conservation in list 1), not
        # "more/less conserved than the genome-wide average" the way
        # single-list mode's NES sign does. Reusing single-list mode's
        # "Conserved" / "Faster-evolving" wording here would misrepresent
        # what a positive or negative NES actually means for a
        # differential run — both direction values need to name a
        # *list*, not describe an absolute rate of evolution. Sign
        # convention matches plotting.plot_differential_distribution()'s
        # own histogram, which labels scores >= 0 as "More variable in
        # <list2>".
        if mode == "differential":
            esc_label1, esc_label2 = html.escape(label1), html.escape(label2)
            direction_labels = {
                1:  "Faster-evolving in %s" % esc_label2,
                -1: "Faster-evolving in %s" % esc_label1,
            }
            direction_desc = (
                "\u201cDirection\u201d reports, for gene sets significant "
                "under a majority of metrics, whether %s or %s appears "
                "faster-evolving for that gene set (i.e. which list shows "
                "lower relative conservation) \u2014 \u201cMixed "
                "direction\u201d means the metrics disagree on which "
                "list, and the result should be treated with caution "
                "rather than as consensus support."
                % (esc_label1, esc_label2)
            )
        else:
            direction_labels = None
            direction_desc = (
                "\u201cDirection\u201d reports whether the sign of "
                "enrichment (conserved vs. faster-evolving) agrees "
                "across the metrics in which the gene set was "
                "significant \u2014 \u201cMixed direction\u201d means "
                "the metrics disagree and the result should be treated "
                "with caution rather than as consensus support."
            )
        consensus_html, consensus_numeric_cols = build_significance_agreement_table_html(
            column_dfs=metric_raw_dfs,
            column_labels=column_labels,
            obo_names=term_names,
            fdr_threshold=fdr_threshold,
            enrichment_plots_map=enrichment_plots_map,
            table_id="results-table-consensus",
            direction_labels=direction_labels,
        )
        tab_buttons_parts.append(
            '    <button class="tab-btn" data-metric="consensus" '
            'role="tab" aria-controls="tab-consensus">Cross-metric consensus</button>'
        )
        if consensus_numeric_cols is not None:
            numeric_cols_map["consensus"] = consensus_numeric_cols
        n_metrics_run = len(metric_raw_dfs)
        majority_n    = n_metrics_run // 2 + 1
        tab_panels_parts.append(
            '<div id="tab-consensus" class="tab-panel" role="tabpanel">\n'
            '  <p class="metric-desc">Gene sets reaching significance '
            '(FDR&nbsp;&lt;&nbsp;{fdr}) in at least {maj} of the {n} ranking '
            'metrics run. Rows highlighted in blue are significant across '
            'all {n} metrics. {ddesc} Full per-gene detail for any of '
            'these gene sets remains available in that metric\u2019s own '
            'tab.</p>\n'
            '  {tbl}\n'
            '</div>\n'.format(
                fdr=fdr_threshold, maj=majority_n, n=n_metrics_run,
                ddesc=direction_desc, tbl=consensus_html,
            )
        )

    # Alignments tab — appended after the metric tabs (and the consensus
    # tab, if present), before Filtering summary. Built from list1's (or
    # the single list's) group2anchor.tsv plus whichever groups actually
    # received a rendered alignment page; see _build_alignment_tree_html()
    # for the full membership rules.
    tab_buttons_parts.append(
        '    <button class="tab-btn" data-metric="alignments" '
        'role="tab" aria-controls="tab-alignments">Alignments</button>'
    )
    alignment_tree_html = _build_alignment_tree_html(
        gmt_gene_sets=gene_sets_all,
        tables_dir1=tables_dir1,
        obo_names=term_names,
        alignment_pages1=alignment_pages1,
        alignment_pages2=alignment_pages2,
        report_dir=report_dir,
        label1=label1,
        label2=label2,
        tables_dir2=tables_dir2,
    )
    tab_panels_parts.append(
        '<div id="tab-alignments" class="tab-panel" role="tabpanel">\n'
        '  <p class="metric-desc">Alignments for every orthogroup that fed '
        'GSEA, nested by GO term. Columns shown hatched/dimmed, if any, '
        'were removed by trimAl before identity estimation.</p>\n'
        '  %s\n'
        '</div>\n' % alignment_tree_html
    )

    # Filtering summary tab — appended after Alignments.
    tab_buttons_parts.append(
        '    <button class="tab-btn" data-metric="filtering-summary" '
        'role="tab" aria-controls="tab-filtering-summary">Filtering summary</button>'
    )
    filtering_html = _build_filtering_summary_tab_content(
        mode=mode, label1=label1, label2=label2,
        pipeline_stats_path=pipeline_stats_path,
        pipeline_stats_path1=pipeline_stats_path1,
        pipeline_stats_path2=pipeline_stats_path2,
        differential_stats_path=differential_stats_path,
        report_dir=report_dir,
    )
    tab_panels_parts.append(
        '<div id="tab-filtering-summary" class="tab-panel" role="tabpanel">\n'
        '  <p class="metric-desc">Counts of input groups and sequences '
        'retained or excluded at each pipeline stage, and why. "View list" '
        'links open a separate page listing the affected group or sequence '
        'IDs.</p>\n'
        '  %s\n'
        '</div>\n' % filtering_html
    )

    # Run parameters tab — always appended last, after Filtering summary.
    # Presents the reproducibility record written by
    # provenance.write_run_parameters() (see cli.py): the effective
    # parameters used to produce this exact report. Never the
    # default-active tab, matching Filtering summary's own placement, and
    # requires no DataTables/JS registration since it has no <table>
    # element of its own (aside from the plain, unsorted key/value tables
    # rendered inline).
    tab_buttons_parts.append(
        '    <button class="tab-btn" data-metric="run-parameters" '
        'role="tab" aria-controls="tab-run-parameters">Run parameters</button>'
    )
    run_parameters_html = _build_run_parameters_tab_content(
        run_parameters_path=run_parameters_path,
        report_dir=report_dir,
    )
    tab_panels_parts.append(
        '<div id="tab-run-parameters" class="tab-panel" role="tabpanel">\n'
        '  <p class="metric-desc">Effective parameters used to produce this '
        'report, captured after any command-line overrides were applied to '
        'the configuration files \u2014 for reproducibility.</p>\n'
        '  %s\n'
        '</div>\n' % run_parameters_html
    )

    html_content = _TEMPLATE.format(
        title=title, dt_css=dt_css,
        tab_buttons="\n".join(tab_buttons_parts),
        tab_panels="\n".join(tab_panels_parts),
        jquery_js=jquery_js, dt_js=dt_js, xlsx_js=xlsx_js,
        enrichment_plots_map=json.dumps(enrichment_plots_map),
        full_gene_sets_js=json.dumps(full_gene_sets_accum),
        leading_edge_map_js=json.dumps(leading_edge_map),
        full_gene_sets_html_js=json.dumps(full_gene_sets_html_accum),
        leading_edge_html_map_js=json.dumps(leading_edge_html_map),
        numeric_cols_map=json.dumps(numeric_cols_map),
    )
    with open(report_path, "w", encoding="utf-8") as fh:
        fh.write(html_content)
    logger.info("HTML report written to: %s", report_path)


# ---------------------------------------------------------------------------
# Public API
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
    tables_dir1: str | None = None,
    tables_dir2: str | None = None,
    label1: str = "List 1",
    label2: str = "List 2",
    pipeline_stats_path: str | None = None,
    pipeline_stats_path1: str | None = None,
    pipeline_stats_path2: str | None = None,
    differential_stats_path: str | None = None,
    alignment_pages1: dict[str, str] | None = None,
    alignment_pages2: dict[str, str] | None = None,
    run_parameters_path: str | None = None,
):
    """Build a single-metric HTML report (one enrichment tab + Alignments +
    Filtering summary + Run parameters).

    Thin wrapper around _build_report_impl() with a one-entry metric_data
    dict, so single-metric and multi-metric reports share one implementation
    and one visual shell. With only one metric, _build_report_impl() never
    builds a Cross-metric consensus tab (that requires 2+ metrics), so
    single-metric reports are unaffected by that feature.

    Args:
        run_parameters_path: Path to run_parameters.json (see
                             provenance.write_run_parameters()), used to
                             populate the Run parameters tab. None omits
                             the record and shows a placeholder message
                             instead of raising.
    """
    metric_data = {metric: {"results_dir": results_dir, "plots_dir": plots_dir}}
    _build_report_impl(
        metric_data=metric_data, report_path=report_path, obo_path=obo_path,
        fdr_threshold=fdr_threshold, mode=mode, gmt_path=gmt_path,
        tables_dir1=tables_dir1, tables_dir2=tables_dir2,
        label1=label1, label2=label2,
        pipeline_stats_path=pipeline_stats_path,
        pipeline_stats_path1=pipeline_stats_path1,
        pipeline_stats_path2=pipeline_stats_path2,
        differential_stats_path=differential_stats_path,
        alignment_pages1=alignment_pages1,
        alignment_pages2=alignment_pages2,
        run_parameters_path=run_parameters_path,
    )


def build_multi_metric_report(
    metric_data: dict,
    report_path: str,
    obo_path: str | None = None,
    fdr_threshold: float = 0.25,
    mode: str = "single",
    gmt_path: str | None = None,
    tables_dir1: str | None = None,
    tables_dir2: str | None = None,
    label1: str = "List 1",
    label2: str = "List 2",
    pipeline_stats_path: str | None = None,
    pipeline_stats_path1: str | None = None,
    pipeline_stats_path2: str | None = None,
    differential_stats_path: str | None = None,
    alignment_pages1: dict[str, str] | None = None,
    alignment_pages2: dict[str, str] | None = None,
    run_parameters_path: str | None = None,
):
    """Build a multi-metric (identity/zscore/rank) tabbed HTML report,
    including the Cross-metric consensus, Alignments, Filtering summary,
    and Run parameters tabs appended after the per-metric tabs.

    Args:
        run_parameters_path: Path to run_parameters.json (see
                             provenance.write_run_parameters()), used to
                             populate the Run parameters tab. None omits
                             the record and shows a placeholder message
                             instead of raising.
    """
    _build_report_impl(
        metric_data=metric_data, report_path=report_path, obo_path=obo_path,
        fdr_threshold=fdr_threshold, mode=mode, gmt_path=gmt_path,
        tables_dir1=tables_dir1, tables_dir2=tables_dir2,
        label1=label1, label2=label2,
        pipeline_stats_path=pipeline_stats_path,
        pipeline_stats_path1=pipeline_stats_path1,
        pipeline_stats_path2=pipeline_stats_path2,
        differential_stats_path=differential_stats_path,
        alignment_pages1=alignment_pages1,
        alignment_pages2=alignment_pages2,
        run_parameters_path=run_parameters_path,
    )
