#!/usr/bin/env python3
"""Builds the static Rejuvenation Atlas site from pipeline results.

Unlike the Streamlit companion app (app/streamlit_app.py), this produces a
single self-contained, server-free HTML page: open it directly in a browser
(file://), no Python runtime required to *view* it. Only Plotly.js is loaded
from a local file next to index.html, so the whole output folder works fully
offline.

    python site/build.py --results results/ --out site/public

Data, CSS, and JS are inlined into index.html; PNG figures are embedded as
base64 data URIs; CSVs referenced by "Download" links are copied alongside.
"""
from __future__ import annotations

import argparse
import base64
import html
import json
import re
import shutil
from pathlib import Path

import pandas as pd

from theme_tokens import BLUE, DIVERGING, templates_payload

HERE = Path(__file__).parent


def load_csv(path: Path) -> pd.DataFrame | None:
    return pd.read_csv(path) if path.exists() else None


def image_data_uri(path: Path) -> str | None:
    if not path.exists():
        return None
    encoded = base64.b64encode(path.read_bytes()).decode("ascii")
    return f"data:image/png;base64,{encoded}"


def df_to_table_html(df: pd.DataFrame, float_format: str = "{:.3f}") -> str:
    cols = list(df.columns)
    thead = "".join(f"<th>{html.escape(str(c))}</th>" for c in cols)
    rows = []
    for _, row in df.iterrows():
        cells = []
        for c in cols:
            v = row[c]
            if isinstance(v, float):
                cells.append(f"<td>{float_format.format(v) if pd.notna(v) else ''}</td>")
            else:
                cells.append(f"<td>{html.escape(str(v))}</td>")
        rows.append(f"<tr>{''.join(cells)}</tr>")
    return (
        '<div class="table-scroll"><table class="data-table">'
        f"<thead><tr>{thead}</tr></thead><tbody>{''.join(rows)}</tbody></table></div>"
    )


def figure_caption(number: int, text: str) -> str:
    return f'<div class="figcaption"><b>Figure {number}.</b> {html.escape(text)}</div>'


def table_caption(number: int, text: str) -> str:
    return f'<div class="tablecaption"><b>Table {number}.</b> {html.escape(text)}</div>'


def download_link(rel_href: str, label: str) -> str:
    return f'<a class="download-link" href="{rel_href}" download>⬇ {html.escape(label)}</a>'


def empty_state(name: str) -> str:
    return (
        '<div class="empty-state">No results found for '
        f"<code>{html.escape(name)}</code> yet. Run the pipeline first: "
        f"<code>nextflow run main.nf -profile docker --module {html.escape(name)}</code></div>"
    )


class Module:
    """Accumulates one module's HTML fragment + chart JSON payload."""

    def __init__(self, mid: str, title: str, badge: str, lede: str):
        self.mid = mid
        self.title = title
        self.badge = badge
        self.lede = lede
        self.body_html: list[str] = []
        self.charts: list[dict] = []
        self.variants: dict[str, dict] = {}

    def nav_html(self, subtitle: str) -> str:
        return (
            f'<button class="nav-item" data-module="{self.mid}">'
            f'<span class="nav-num">{self.mid}</span>'
            f'<span><span>{html.escape(self.title)}</span>'
            f'<span class="nav-sub">{html.escape(subtitle)}</span></span>'
            "</button>"
        )

    def section_html(self) -> str:
        return (
            f'<section class="module-section" data-module="{self.mid}">'
            f'<span class="module-badge">{html.escape(self.badge)}</span>'
            f'<h1 class="module-title">{html.escape(self.title)}</h1>'
            f'<p class="module-lede">{html.escape(self.lede)}</p>'
            + "".join(self.body_html)
            + "</section>"
        )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", default="results", type=Path)
    parser.add_argument("--out", default=str(HERE / "public"), type=Path)
    args = parser.parse_args()

    results: Path = args.results
    out: Path = args.out
    if out.exists():
        shutil.rmtree(out)
    (out / "data").mkdir(parents=True)
    (out / "assets").mkdir(parents=True)

    modules: dict[str, Module] = {}
    nav_html_parts: list[str] = []

    # ---------------------------------------------------------------
    # Module 1 — epigenetic clock
    # ---------------------------------------------------------------
    m = Module("1", "Epigenetic age across transient reprogramming",
               "GSE165179 · Gill et al. 2022, eLife",
               "Does maturation-phase transient reprogramming (MPTR) reverse the DNA-methylation age of human fibroblasts?")
    df = load_csv(results / "module1" / "dnam_age_timecourse.csv")
    stats_path = results / "module1" / "delta_age_stats.json"
    if df is None:
        m.body_html.append(empty_state("rejuvenation_clock"))
    else:
        chart = {
            "id": "chart-m1-line",
            "data": [{"type": "scatter", "mode": "lines+markers", "name": "DNAm age",
                      "x": df["day"].tolist(), "y": df["dnam_age"].tolist(),
                      "line": {"width": 2}, "marker": {"size": 8}}],
            "layout": {"xaxis": {"title": {"text": "MPTR day"}},
                       "yaxis": {"title": {"text": "DNAm age (years)"}}},
            "color": BLUE,
        }
        m.charts.append(chart)
        m.body_html.append(
            '<div class="card"><div class="card-head"><span class="card-title">DNAm age vs. MPTR day</span></div>'
            f'<div id="{chart["id"]}" class="chart-canvas"></div>'
            + figure_caption(1, "GSE165179 — DNAm age falls below baseline after maturation phase "
                                 "transient reprogramming (MPTR).")
            + "</div>"
        )
        if stats_path.exists():
            stats = json.loads(stats_path.read_text())
            ci = stats.get("ci95", [0, 0])
            m.body_html.append(
                '<div class="stat-row">'
                f'<div class="stat-tile"><div class="stat-label">Δ age ({html.escape(stats.get("contrast", "peak vs day 0"))})</div>'
                f'<div class="stat-value">{stats.get("delta_years", 0):.1f} y</div></div>'
                f'<div class="stat-tile"><div class="stat-label">95% CI</div>'
                f'<div class="stat-value">[{ci[0]:.1f}, {ci[1]:.1f}]</div></div>'
                f'<div class="stat-tile"><div class="stat-label">Cohen’s d</div>'
                f'<div class="stat-value">{stats.get("cohens_dz", 0):.2f}</div></div>'
                "</div>"
            )
    modules["1"] = m

    # ---------------------------------------------------------------
    # Module 2 — reprogramming trajectory
    # ---------------------------------------------------------------
    m = Module("2", "Fate probabilities from optimal transport",
               "GSE122662 · Schiebinger et al. 2019, Cell",
               "Which day-0 fibroblasts are on a trajectory that completes reprogramming to iPSC?")
    df = load_csv(results / "module2" / "fate_probabilities.csv")
    if df is None:
        m.body_html.append(empty_state("trajectory"))
    else:
        chart = {
            "id": "chart-m2-hist",
            "data": [{"type": "histogram", "name": "cells", "x": df["ipsc_fate_prob"].tolist(), "nbinsx": 50}],
            "layout": {"xaxis": {"title": {"text": "P(reach iPSC fate)"}},
                       "yaxis": {"title": {"text": "count"}}},
            "color": BLUE,
        }
        m.charts.append(chart)
        m.body_html.append(
            '<div class="card"><div class="card-head"><span class="card-title">Day-0 fate-probability distribution</span></div>'
            f'<div id="{chart["id"]}" class="chart-canvas"></div>'
            + figure_caption(2, "Day-0 cells — probability of successful reprogramming under the "
                                 "fitted optimal-transport coupling.")
            + "</div>"
        )
        img = image_data_uri(results / "module2" / "fate_analysis.png")
        if img:
            m.body_html.append(
                f'<div class="card"><img src="{img}" alt="Fate analysis">'
                + figure_caption(3, "UMAP of reprogramming endpoints and Waddington-OT fate probabilities "
                                     "(Schiebinger et al. 2019, rebuilt).")
                + "</div>"
            )
    modules["2"] = m

    # ---------------------------------------------------------------
    # Module 3 — cross-study integration
    # ---------------------------------------------------------------
    m = Module("3", "Does partial reprogramming move aged tissue toward a young state?",
               "GSE190983 × GSE149590 · Browder et al. 2022, Nat Aging",
               "Browder 2022 bulk RNA-seq (6 tissues, 7-month cyclic OSK) projected onto the Tabula Muris "
               "Senis aging axis. Negative Δ = 4F shifted young.")
    df = load_csv(results / "module3" / "tissue_scores.csv")
    if df is None:
        m.body_html.append(empty_state("integration"))
    else:
        csv_name = "module3_tissue_scores.csv"
        shutil.copy(results / "module3" / "tissue_scores.csv", out / "data" / csv_name)
        m.body_html.append(
            '<div class="card"><div class="card-head"><span class="card-title">Per-tissue aging-axis scores</span>'
            + download_link(f"data/{csv_name}", "tissue_scores.csv") + "</div>"
            + df_to_table_html(df)
            + table_caption(1, "Per-tissue aging-axis scores, 4F vs. control, with 95% CI, Cohen's d "
                                "and BH-adjusted p-value.")
            + "</div>"
        )
        img = image_data_uri(results / "module3" / "aging_score.png")
        if img:
            m.body_html.append(
                f'<div class="card"><img src="{img}" alt="Aging score">'
                + figure_caption(4, "Tissue-level aging-score shift, 4F vs. control.")
                + "</div>"
            )
    modules["3"] = m

    # ---------------------------------------------------------------
    # Module 4 — spatial aging (age-group variants of one heatmap)
    # ---------------------------------------------------------------
    m = Module("4", "Spatial neighborhood changes with age (MERFISH)",
               "GSE207848 · Allen et al. 2023, Cell",
               "Cell-type neighborhood enrichment z-scores across the aging mouse brain, by age cohort.")
    mod4_dir = results / "module4"
    files = sorted(mod4_dir.glob("nhood_enrichment_*.csv")) if mod4_dir.exists() else []
    files.sort(key=lambda f: int(re.match(r"\d+", f.stem.replace("nhood_enrichment_", "")).group() or 0))
    if not files:
        m.body_html.append(empty_state("spatial"))
    else:
        ages = [f.stem.replace("nhood_enrichment_", "") for f in files]
        chart_id = "chart-m4-heatmap"
        for f, age in zip(files, ages):
            z = pd.read_csv(f, index_col=0)
            zmax = float(z.abs().max().max())
            variant = {
                "id": chart_id,
                "data": [{"type": "heatmap", "z": z.values.tolist(), "x": list(z.columns), "y": list(z.index),
                          "zmin": -zmax, "zmax": zmax}],
                "layout": {
                    "margin": {"l": 24, "r": 24, "t": 24, "b": 24},
                    "xaxis": {"automargin": True, "tickangle": -40},
                    "yaxis": {"automargin": True},
                },
                "colorscale": DIVERGING,
            }
            m.variants[age] = variant
        default_age = ages[0]
        m.charts.append(m.variants[default_age])
        options = "".join(f'<option value="{html.escape(a)}">{html.escape(a)}</option>' for a in ages)
        m.body_html.append(
            '<div class="card"><div class="card-head"><span class="card-title">Neighborhood enrichment z-scores</span>'
            f'<select class="control variant-select" data-module="4" data-chart="{chart_id}">{options}</select></div>'
            f'<div id="{chart_id}" class="chart-canvas" style="min-height:480px"></div>'
            + figure_caption(5, "Cell-type neighborhood enrichment z-scores. Positive values indicate "
                                 "co-localization above chance.")
            + "</div>"
        )
    modules["4"] = m

    # ---------------------------------------------------------------
    # Module 5 — SASP proteomics
    # ---------------------------------------------------------------
    m = Module("5", "Core SASP from proteomics",
               "PXD013721 · Basisty et al. 2020, PLoS Biology",
               "Proteins reproducibly induced across senescence inducers in the SASP Atlas.")
    df = load_csv(results / "module5" / "core_sasp.csv")
    if df is None:
        m.body_html.append(empty_state("sasp"))
    else:
        m.body_html.append(
            '<div class="stat-row"><div class="stat-tile"><div class="stat-label">Core SASP proteins</div>'
            f'<div class="stat-value">{len(df)}</div></div></div>'
        )
        chips = "".join(f'<span class="chip">{html.escape(str(v))}</span>' for v in df.iloc[:, 0].tolist())
        csv_name = "module5_core_sasp.csv"
        shutil.copy(results / "module5" / "core_sasp.csv", out / "data" / csv_name)
        m.body_html.append(
            '<div class="card"><div class="card-head"><span class="card-title">Protein list</span>'
            + download_link(f"data/{csv_name}", "core_sasp.csv") + "</div>"
            + f'<div class="chip-box">{chips}</div>'
            + table_caption(2, "Proteins reproducibly induced across senescence inducers "
                                "(Basisty et al. 2020 SASP Atlas).")
            + "</div>"
        )
        img = image_data_uri(results / "module5" / "sasp_per_inducer.png")
        if img:
            m.body_html.append(
                f'<div class="card"><img src="{img}" alt="SASP per inducer">'
                + figure_caption(6, "Core SASP protein abundance by senescence inducer.")
                + "</div>"
            )
    modules["5"] = m

    # ---------------------------------------------------------------
    # Assemble page
    # ---------------------------------------------------------------
    subtitles = {
        "1": "Gill 2022", "2": "Schiebinger 2019", "3": "Browder 2022 × TMS",
        "4": "Allen 2023", "5": "Basisty 2020",
    }
    for mid, mod in modules.items():
        nav_html_parts.append(mod.nav_html(subtitles[mid]))

    payload = {
        "templates": templates_payload(),
        "modules": {
            mid: {"charts": mod.charts, "variants": mod.variants}
            for mid, mod in modules.items()
        },
    }

    css = (HERE / "static" / "style.css").read_text()
    js = (HERE / "static" / "app.js").read_text()

    plotly_src = Path(pd.__file__).parent.parent / "plotly" / "package_data" / "plotly.min.js"
    if plotly_src.exists():
        shutil.copy(plotly_src, out / "assets" / "plotly.min.js")
        plotly_tag = '<script src="assets/plotly.min.js"></script>'
    else:
        plotly_tag = '<script src="https://cdn.plot.ly/plotly-2.32.0.min.js"></script>'

    sections_html = "".join(mod.section_html() for mod in modules.values())

    html_doc = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Cell Rejuvenation Atlas — Explorer</title>
<meta name="description" content="Static, offline-capable companion to the Cell Rejuvenation Atlas reproducible reanalysis pipeline.">
<style>
{css}
</style>
</head>
<body>
<header class="top-bar">
  <div class="brand"><span class="brand-mark">CRA</span>Cell Rejuvenation Atlas<small>&nbsp;Explorer</small></div>
  <div class="top-bar-actions">
    <a class="btn" href="https://github.com/lynchaos/cell-rejuvenation-atlas" target="_blank" rel="noopener">Repository</a>
    <button class="btn" id="theme-toggle" type="button">● Dark</button>
  </div>
</header>
<div class="layout-body">
  <nav class="sidebar">
    <div class="sidebar-label">Module</div>
    {''.join(nav_html_parts)}
    <div class="sidebar-footer">Figures regenerated from versioned code — never pasted by hand.</div>
  </nav>
  <main class="content">
    {sections_html}
  </main>
</div>
<footer class="site-footer">
  Reproducible reanalyses of peer-reviewed aging &amp; reprogramming data. See
  <a href="https://github.com/lynchaos/cell-rejuvenation-atlas/blob/main/DATASETS.md">DATASETS.md</a>
  for full citations and accessions. Code licensed AGPL-3.0; datasets remain the property of their
  original studies.
</footer>
<script type="application/json" id="site-data">{json.dumps(payload)}</script>
{plotly_tag}
<script>
{js}
</script>
</body>
</html>
"""
    (out / "index.html").write_text(html_doc)
    size_kb = (out / "index.html").stat().st_size / 1024
    print(f"Wrote {out / 'index.html'} ({size_kb:.0f} KB)")
    print(f"Open it directly in a browser, or serve: python -m http.server -d {out} 8877")


if __name__ == "__main__":
    main()
