"""Shared academic visual theme for the Rejuvenation Explorer app.

Registers a Plotly template ("academic") and provides the page CSS plus
small HTML-rendering helpers, so figures and tables read like journal
supplementary-material panels rather than a default dashboard.
"""
from __future__ import annotations

import html as _html

import pandas as pd
import plotly.graph_objects as go
import plotly.io as pio

# -- Palette (validated 8-hue categorical order; blue/red diverging pair) ----
BLUE = "#2a78d6"
GREEN = "#008300"
MAGENTA = "#e87ba4"
YELLOW = "#eda100"
AQUA = "#1baf7a"
ORANGE = "#eb6834"
VIOLET = "#4a3aa7"
RED = "#e34948"

CATEGORICAL = [BLUE, GREEN, MAGENTA, YELLOW, AQUA, ORANGE, VIOLET, RED]
DIVERGING = [[0.0, BLUE], [0.5, "#f0efec"], [1.0, RED]]

SURFACE = "#fcfcfb"
PAGE_PLANE = "#f9f9f7"
INK_PRIMARY = "#0b0b0b"
INK_SECONDARY = "#52514e"
INK_MUTED = "#898781"
GRIDLINE = "#e1e0d9"
BASELINE = "#c3c2b7"
BORDER = "rgba(11, 11, 11, 0.10)"

_FONT = "Inter, -apple-system, 'Segoe UI', system-ui, sans-serif"


def _build_template() -> go.layout.Template:
    axis = dict(
        gridcolor=GRIDLINE,
        zerolinecolor=BASELINE,
        linecolor=BASELINE,
        showline=True,
        ticks="outside",
        tickcolor=GRIDLINE,
        tickfont=dict(color=INK_MUTED, size=12),
        title=dict(font=dict(color=INK_SECONDARY, size=12)),
    )
    template = go.layout.Template()
    template.layout = go.Layout(
        font=dict(family=_FONT, size=13, color=INK_SECONDARY),
        title=dict(font=dict(family=_FONT, size=15, color=INK_PRIMARY)),
        paper_bgcolor=SURFACE,
        plot_bgcolor=SURFACE,
        colorway=CATEGORICAL,
        margin=dict(l=56, r=24, t=24, b=48),
        xaxis=axis,
        yaxis=axis,
        legend=dict(bgcolor="rgba(0,0,0,0)", bordercolor=BORDER, borderwidth=1,
                    font=dict(size=12, color=INK_SECONDARY)),
        hoverlabel=dict(bgcolor=SURFACE, bordercolor=BORDER,
                         font=dict(family=_FONT, size=12, color=INK_PRIMARY)),
        hovermode="closest",
    )
    template.data.scatter = [go.Scatter(marker=dict(size=8, line=dict(width=0)), line=dict(width=2))]
    template.data.histogram = [go.Histogram(marker=dict(color=BLUE, line=dict(color=SURFACE, width=1)))]
    return template


pio.templates["academic"] = _build_template()
pio.templates.default = "academic"


PAGE_CSS = """
<style>
@import url('https://fonts.googleapis.com/css2?family=Source+Serif+4:opsz,wght@8..60,400;8..60,600;8..60,700&family=Inter:wght@400;500;600&family=JetBrains+Mono:wght@400;500&display=swap');

html, body, [class*="css"] { font-family: 'Inter', -apple-system, 'Segoe UI', sans-serif; }

.block-container { max-width: 980px; padding-top: 2.5rem; }

h1, h2, h3 { font-family: 'Source Serif 4', Georgia, serif !important; color: #0b0b0b; letter-spacing: -0.01em; }
h1 { font-weight: 700; border-bottom: 1px solid #e1e0d9; padding-bottom: 0.6rem; }
h2 { font-weight: 600; font-size: 1.3rem !important; margin-top: 0.25rem; }

[data-testid="stCaptionContainer"] p { font-style: italic; color: #898781; }

.acad-figcaption {
    font-family: 'Source Serif 4', Georgia, serif;
    font-style: italic;
    font-size: 0.92rem;
    color: #52514e;
    border-top: 1px solid #e1e0d9;
    margin-top: -0.6rem;
    padding-top: 0.6rem;
}
.acad-figcaption b { font-style: normal; color: #0b0b0b; }

[data-testid="stMetric"] {
    background: #fcfcfb;
    border: 1px solid rgba(11, 11, 11, 0.10);
    border-radius: 6px;
    padding: 0.9rem 1rem 0.7rem;
}
[data-testid="stMetricValue"] { font-variant-numeric: tabular-nums; color: #0b0b0b; }
[data-testid="stMetricLabel"] { color: #52514e; }

section[data-testid="stSidebar"] {
    background: #f9f9f7;
    border-right: 1px solid #e1e0d9;
}
section[data-testid="stSidebar"] label p { font-family: 'Inter', sans-serif; font-size: 0.92rem; }

table.acad-table {
    width: 100%;
    border-collapse: collapse;
    font-family: 'Inter', sans-serif;
    font-size: 0.85rem;
    font-variant-numeric: tabular-nums;
}
table.acad-table th {
    text-align: left;
    color: #52514e;
    font-weight: 600;
    border-bottom: 1.5px solid #c3c2b7;
    padding: 0.4rem 0.6rem;
    white-space: nowrap;
}
table.acad-table td {
    padding: 0.35rem 0.6rem;
    border-bottom: 1px solid #e1e0d9;
    color: #0b0b0b;
    white-space: nowrap;
}
table.acad-table tr:last-child td { border-bottom: none; }
.acad-table-scroll { overflow-x: auto; border: 1px solid #e1e0d9; border-radius: 6px; }

.acad-chip-box {
    max-height: 220px;
    overflow-y: auto;
    border: 1px solid #e1e0d9;
    border-radius: 6px;
    padding: 0.7rem;
    background: #fcfcfb;
    line-height: 2.1;
}
.acad-chip {
    display: inline-block;
    font-family: 'JetBrains Mono', ui-monospace, monospace;
    font-size: 0.78rem;
    color: #0b0b0b;
    background: #f9f9f7;
    border: 1px solid #e1e0d9;
    border-radius: 4px;
    padding: 0.12rem 0.45rem;
    margin: 0 0.28rem 0.28rem 0;
}
</style>
"""


def inject() -> None:
    """Apply the academic CSS theme to the current Streamlit page."""
    import streamlit as st

    st.markdown(PAGE_CSS, unsafe_allow_html=True)


def figure_caption(number: int, text: str) -> str:
    return f'<div class="acad-figcaption"><b>Figure {number}.</b> {text}</div>'


def table_caption(number: int, text: str) -> str:
    return f'<div class="acad-figcaption"><b>Table {number}.</b> {text}</div>'


def df_to_html(df: pd.DataFrame, float_format: str = "{:.3f}") -> str:
    """Render a DataFrame as a styled HTML table (no pyarrow involved)."""
    formatted = df.copy()
    for col in formatted.select_dtypes(include=["float", "float64", "float32"]).columns:
        formatted[col] = formatted[col].map(lambda v: float_format.format(v) if pd.notna(v) else "")
    table = formatted.to_html(index=False, border=0, classes="acad-table", escape=True)
    return f'<div class="acad-table-scroll">{table}</div>'


def chip_list_html(items: list[str]) -> str:
    chips = "".join(f'<span class="acad-chip">{_html.escape(str(item))}</span>' for item in items)
    return f'<div class="acad-chip-box">{chips}</div>'
