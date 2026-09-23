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
a.go-link, a.geneset-link, a.leadedge-link {{
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
  <h2>Plots</h2>
  <div class="plot-grid">
{plots_html}
  </div>
</section>
<section>
  <h2>Enrichment results</h2>
  <p>Significant gene sets are highlighted in blue.
     Click a GO ID to view its enrichment plot, or "View" under
     Full&nbsp;gene&nbsp;set / Leading&nbsp;edge to see the gene lists as text.
     Use the filter boxes below each column header to filter by that column.</p>
  <p>
    <button id="sig-toggle-btn" class="sig-toggle-btn">Show only significant</button>
    <button id="export-xlsx-btn" class="export-btn">&#8681; Export table to Excel (.xlsx)</button>
    <button id="export-xlsx-all-btn" class="export-btn">&#8681; Export all rows (ignore filters)</button>
  </p>
{table_html}
</section>
</main>
<footer>Generated by ENHYDRA</footer>
<script>{jquery_js}</script>
<script>{dt_js}</script>
<script>{xlsx_js}</script>
<script>
{plot_data_js}
{gene_data_js}

function getHeaderIndex(tableSelector, name) {{
    var headers = [];
    $(tableSelector).find('thead tr').first().find('th').each(function() {{
        var $clone = $(this).clone();
        $clone.find('.col-tip').remove();
        headers.push($clone.text().trim());
    }});
    return headers.indexOf(name);
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
    var tip = document.getElementById('svg-tooltip');
    document.querySelectorAll('[data-tip]').forEach(function(el) {{
        el.addEventListener('mousemove', function(e) {{
            tip.innerHTML = this.getAttribute('data-tip');
            tip.style.display = 'block';
            tip.style.left = (e.clientX + 15) + 'px';
            tip.style.top  = (e.clientY + 15) + 'px';
        }});
        el.addEventListener('mouseleave', function() {{ tip.style.display = 'none'; }});
    }});
    var numericCols = {numeric_col_indices};
    var colFilters  = {{}};
    var sigOnly     = false;
    var sigColIndex = null;
    $.fn.dataTable.ext.search.push(function(settings, data) {{
        if (sigOnly && sigColIndex !== null && data[sigColIndex] !== '\u2713') {{
            return false;
        }}
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
        pageLength: 25, orderCellsTop: true, order: [[4, 'asc']],
        columnDefs: [{{ targets: numericCols, type: 'num' }}]
    }});
    sigColIndex = getHeaderIndex('#results-table', 'Sig.');
    $('#results-table thead tr.filter-row th').each(function(i) {{
        var isNumeric   = numericCols.indexOf(i) !== -1;
        var placeholder = isNumeric ? "e.g. < 0.05" : "Filter...";
        var input       = $('<input type="text" placeholder="' + placeholder + '"/>');
        $(this).html(input);
        input.on('keyup change', function() {{
            var val = this.value.trim();
            if (val === '') {{
                delete colFilters[i];
            }} else if (isNumeric) {{
                var m = val.match(/^([<>=!]=?)\\s*([\\d.eE+\\-]+)$/);
                if (m) colFilters[i] = {{ op: m[1], num: parseFloat(m[2]) }};
                else   delete colFilters[i];
            }} else {{
                colFilters[i] = {{ text: val.toLowerCase() }};
            }}
            table.draw();
        }});
    }});
    $('#sig-toggle-btn').on('click', function() {{
        sigOnly = !sigOnly;
        $(this).toggleClass('active', sigOnly);
        $(this).text(sigOnly ? 'Showing significant only \u2713' : 'Show only significant');
        table.draw();
    }});
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
    $(document).on('click', '.go-link', function(e) {{
        e.preventDefault();
        var goId = $(this).data('goid');
        var uri  = enrichmentPlots[goId];
        if (uri) showImageModal(goId, uri);
    }});
    $(document).on('click', '.geneset-link', function(e) {{
        e.preventDefault();
        var goId = $(this).data('goid');
        var text = fullGeneSets[goId];
        if (text !== undefined) showTextModal(goId + ' \u2014 full gene set', text);
    }});
    $(document).on('click', '.leadedge-link', function(e) {{
        e.preventDefault();
        var goId = $(this).data('goid');
        var text = leadingEdge[goId];
        if (text !== undefined) showTextModal(goId + ' \u2014 leading edge genes', text);
    }});
    $('#modal-close, #modal-overlay').on('click', function(e) {{
        if (e.target === this) $('#modal-overlay').removeClass('active');
    }});
    $(document).on('click', '#export-xlsx-btn', function() {{
        tableToXLSX('#results-table', numericCols, fullGeneSets, leadingEdge,
                   'enrichment_results.xlsx', true);
    }});
    $(document).on('click', '#export-xlsx-all-btn', function() {{
        tableToXLSX('#results-table', numericCols, fullGeneSets, leadingEdge,
                   'enrichment_results_all.xlsx', false);
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
a.go-link, a.geneset-link, a.leadedge-link {{
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
<script>{xlsx_js}</script>
<script>
var enrichmentPlotsMap = {enrichment_plots_map};
var fullGeneSets       = {full_gene_sets_js};
var leadingEdgeMap     = {leading_edge_map_js};
var numericColsMap     = {numeric_cols_map};
var dtInstances        = {{}};
var colFiltersMap      = {{}};
var sigOnlyMap         = {{}};
var sigColIndexMap     = {{}};

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
    var dt = $('#results-table-' + metric).DataTable({{
        pageLength: 25, orderCellsTop: true, order: [[4, 'asc']],
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
            initTable(metric);
        }});
    }});
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
        var text = fullGeneSets[goId];
        if (text !== undefined) showTextModal(goId + ' \u2014 full gene set', text);
    }});
    $(document).on('click', '.leadedge-link', function(e) {{
        e.preventDefault();
        var goId   = $(this).data('goid');
        var metric = $(this).data('metric');
        var map    = metric ? leadingEdgeMap[metric] : undefined;
        var text   = map ? map[goId] : undefined;
        if (text !== undefined) {{
            showTextModal(goId + ' \u2014 leading edge genes (' + metric + ')', text);
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
# Internal helpers
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
    """Return term_id -> full list of member gene IDs, read from a GMT file.

    Used to power the "Full gene set" link column: the complete membership
    of each gene set, independent of ranking metric.
    """
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

    Gene IDs are always read as str to prevent pandas inferring numeric
    columns as int/float, which would break matching against GMT identifiers.
    """
    def _load_anchor2mean(tables_dir: str) -> dict[str, float]:
        path = os.path.join(tables_dir, "anchor2mean.tsv")
        if not os.path.isfile(path):
            return {}
        df = pd.read_csv(path, sep="\t", header=None,
                         names=["gene_id", "score"],
                         dtype={"gene_id": str})          # always str
        df["score"] = pd.to_numeric(df["score"], errors="coerce")
        df = df.dropna(subset=["score"]).drop_duplicates("gene_id")
        return dict(_normalise_series(df.set_index("gene_id")["score"], metric))

    def _load_via_group_mapping(tables_dir_scores: str,
                                tables_dir_mapping: str) -> dict[str, float]:
        """Fallback for lists without anchor sequences."""
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
                          dtype={"group_id": str, "gene_id": str})  # always str
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

    Rather than reading and base64-encoding every PNG (which used to embed
    the full binary content of every plot directly into the HTML — the
    dominant cause of multi-hundred-MB reports on large analyses such as
    vertebrate genomes), this now records a path *relative to the report's
    own directory* on disk. The browser then loads each plot on demand,
    directly from the enrichment/prerank/ folder, only when the user clicks
    a GO ID. The images are never copied or embedded — the report simply
    keeps a pointer to where they already live.

    Note: this means the resulting report.html is no longer a single
    self-contained file. It must stay alongside the output directory
    structure (specifically the relevant enrichment/prerank/ subfolder) for
    the modal images to resolve.

    Args:
        results_dir: Directory passed to run_gsea() for this metric/list
                     (i.e. the parent of the 'prerank' subfolder).
        report_dir:  Directory where the report.html file itself will be
                     written, used as the base for the relative paths.

    Returns:
        Dict mapping GO ID -> path (relative to report_dir) of its PNG.
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
        # Use forward slashes regardless of platform, since this path is
        # used as a URL/src attribute inside the HTML, not a filesystem call.
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
) -> tuple[str, list[int], dict[str, str], dict[str, str]]:
    """Build the results table HTML.

    In addition to the rendered table, this returns two lookup dicts used
    to power the "Full gene set" and "Leading edge" link columns (and, in
    turn, the "Export to Excel" button, which reads gene lists from these
    same dicts rather than from the truncated "View (N)" cell text):
      - full_sets_js: GO ID -> newline-joined full gene set membership
                      (from the GMT; identical regardless of metric).
      - leadedge_js:  GO ID -> newline-joined leading-edge gene list for
                      this specific ranking metric (from GSEApy's own
                      Lead_genes column).
    Neither gene list is written into the table cell itself — only a
    "View (N)" link is, and the actual text is injected into the page as a
    small JS lookup object, shown in a modal on click (or exported to xlsx
    on export). The "Sig." column (rendered as a checkmark, '\u2713') is
    also what the "Show only significant" toggle filters on client-side.
    """
    df = df.copy()
    if "Term" not in df.columns:
        logger.warning("'Term' column not found in GSEA results.")
        return "<p>No results to display.</p>", [], {}, {}

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

    # --- Full gene set / leading edge link columns --------------------
    full_sets_js: dict[str, str] = {}
    leadedge_js:  dict[str, str] = {}

    if gene_sets:
        def _full_set_cell(term_id):
            genes = gene_sets.get(term_id)
            if not genes:
                return ""
            full_sets_js[term_id] = "\n".join(genes)
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

    html = (
        '<table id="%s" class="display compact" style="width:100%%">'
        '<thead><tr>%s</tr><tr class="filter-row">%s</tr></thead>'
        "<tbody>%s</tbody></table>"
    ) % (table_id, header_cells, filter_cells, rows)

    return html, numeric_col_indices, full_sets_js, leadedge_js


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
                '<p class="plot-caption">%s (hover for details)</p>'
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


# ---------------------------------------------------------------------------
# Single-metric report
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
):
    logger.info("Building HTML report...")
    effective_gmt = gmt_path or _find_gmt_in_dir(results_dir)
    term_names    = _resolve_term_names(results_dir, obo_path, effective_gmt)
    gene_sets     = _gmt_gene_sets(effective_gmt)
    df = _load_gsea_results(results_dir)
    if df is None:
        logger.warning("Cannot build report: no GSEA results found.")
        return
    df = _augment_with_per_term_scores(
        df, effective_gmt, tables_dir1, tables_dir2, metric
    )
    cache_dir    = os.path.dirname(obo_path) if obo_path else None
    jquery_js    = _fetch_cached(_JQUERY_URL,         cache_dir, "jquery.min.js")
    dt_js        = _fetch_cached(_DATATABLES_JS_URL,  cache_dir, "datatables.min.js")
    dt_css       = _fetch_cached(_DATATABLES_CSS_URL, cache_dir, "datatables.min.css")
    xlsx_js      = _fetch_cached(_XLSX_JS_URL,        cache_dir, "xlsx.full.min.js")

    # Plot links are resolved relative to the directory the report itself
    # will be written into, so images are read from disk on click rather
    # than being embedded as base64 blobs in the HTML.
    report_dir   = os.path.dirname(os.path.abspath(report_path))
    plot_index   = _build_enrichment_plot_index(results_dir, report_dir)
    plot_data_js = "var enrichmentPlots = {%s};" % ",".join(
        '"%s": "%s"' % (go_id, path) for go_id, path in plot_index.items()
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
    table_html, numeric_col_indices, full_sets_js, leadedge_js = _results_table_html(
        df, term_names, plot_index, fdr_threshold,
        metric=None, col1_label=label1, col2_label=label2,
        gene_sets=gene_sets,
    )
    # Full gene sets are metric-independent (straight from the GMT); leading
    # edge genes are specific to this metric's ranking. Both are injected as
    # plain JS lookup objects — never written directly into table cells —
    # so the table (and the resulting Excel export) stays fast to build even
    # when gene sets contain hundreds of genes.
    gene_data_js = "var fullGeneSets = %s;\nvar leadingEdge = %s;" % (
        json.dumps(full_sets_js), json.dumps(leadedge_js)
    )
    html = _TEMPLATE.format(
        title=title, dt_css=dt_css, plots_html=plots_html,
        table_html=table_html, jquery_js=jquery_js, dt_js=dt_js,
        xlsx_js=xlsx_js,
        plot_data_js=plot_data_js, gene_data_js=gene_data_js,
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
    tables_dir1: str | None = None,
    tables_dir2: str | None = None,
    label1: str = "List 1",
    label2: str = "List 2",
):
    logger.info("Building multi-metric HTML report (%d metrics)...", len(metric_data))
    first_results = next(iter(metric_data.values()))["results_dir"] if metric_data else ""
    effective_gmt = gmt_path or _find_gmt_in_dir(first_results)
    term_names    = _resolve_term_names(first_results, obo_path, effective_gmt)
    gene_sets_all = _gmt_gene_sets(effective_gmt)
    cache_dir = os.path.dirname(obo_path) if obo_path else None
    jquery_js = _fetch_cached(_JQUERY_URL,         cache_dir, "jquery.min.js")
    dt_js     = _fetch_cached(_DATATABLES_JS_URL,  cache_dir, "datatables.min.js")
    dt_css    = _fetch_cached(_DATATABLES_CSS_URL, cache_dir, "datatables.min.css")
    xlsx_js   = _fetch_cached(_XLSX_JS_URL,        cache_dir, "xlsx.full.min.js")
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

    # Plot links are resolved relative to the directory the report itself
    # will be written into (see build_report for the rationale).
    report_dir = os.path.dirname(os.path.abspath(report_path))

    tab_buttons_parts    = []
    tab_panels_parts     = []
    enrichment_plots_map = {}
    numeric_cols_map     = {}
    # Full gene set membership is the same across metrics (same GMT), so it
    # is accumulated into one flat lookup shared by all tabs. Leading edge
    # genes differ per metric (different ranking -> different leading edge),
    # so that lookup stays keyed by metric, mirroring enrichment_plots_map.
    full_gene_sets_accum: dict[str, str] = {}
    leading_edge_map: dict[str, dict[str, str]] = {}
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
            df = _augment_with_per_term_scores(
                df, effective_gmt, tables_dir1, tables_dir2, metric
            )
            tbl_html, num_cols, full_sets_js, leadedge_js = _results_table_html(
                df, term_names, plot_idx, fdr_threshold,
                metric=metric, col1_label=label1, col2_label=label2,
                gene_sets=gene_sets_all,
            )
            numeric_cols_map[metric] = num_cols
            full_gene_sets_accum.update(full_sets_js)
            leading_edge_map[metric] = leadedge_js
        else:
            tbl_html = "<p>No GSEA results found for this metric.</p>"
            numeric_cols_map[metric] = []
            leading_edge_map[metric] = {}
        plots_html = _plot_section(plots_dir, plot_names)
        desc       = _METRIC_DESCS.get(metric, "")
        tab_panels_parts.append(
            '<div id="tab-{m}" class="tab-panel{ac}" role="tabpanel">\n'
            '  <p class="metric-desc">{desc}</p>\n'
            '  <h3>Plots</h3>\n'
            '  <div class="plot-grid">{plots}</div>\n'
            '  <h3>Enrichment results</h3>\n'
            '  <p>Significant gene sets (FDR&nbsp;&lt;&nbsp;{fdr}) highlighted '
            'in blue. Click a GO ID to view its enrichment plot, or "View" '
            'under Full&nbsp;gene&nbsp;set / Leading&nbsp;edge to see the '
            'gene lists as text.</p>\n'
            '  <p>'
            '<button class="sig-toggle-btn" data-metric="{m}">Show only significant</button> '
            '<button class="export-btn export-xlsx-btn" data-metric="{m}">'
            '&#8681; Export table to Excel (.xlsx)</button> '
            '<button class="export-btn export-xlsx-all-btn" data-metric="{m}">'
            '&#8681; Export all rows (ignore filters)</button>'
            '</p>\n'
            '  {tbl}\n'
            '</div>\n'.format(
                m=metric, ac=active_cls, desc=desc,
                plots=plots_html, fdr=fdr_threshold, tbl=tbl_html,
            )
        )
        first = False
    html = _MULTI_TEMPLATE.format(
        title=title, dt_css=dt_css,
        tab_buttons="\n".join(tab_buttons_parts),
        tab_panels="\n".join(tab_panels_parts),
        jquery_js=jquery_js, dt_js=dt_js, xlsx_js=xlsx_js,
        enrichment_plots_map=json.dumps(enrichment_plots_map),
        full_gene_sets_js=json.dumps(full_gene_sets_accum),
        leading_edge_map_js=json.dumps(leading_edge_map),
        numeric_cols_map=json.dumps(numeric_cols_map),
    )
    with open(report_path, "w", encoding="utf-8") as fh:
        fh.write(html)
    logger.info("Multi-metric HTML report written to: %s", report_path)
