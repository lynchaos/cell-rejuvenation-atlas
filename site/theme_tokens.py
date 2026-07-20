"""Design tokens for the static Explorer (site/build.py).

One validated categorical/diverging palette, one Plotly template builder,
kept separate from the HTML/CSS assembly in build.py.
"""
from __future__ import annotations

FONT_SANS = "Inter, -apple-system, 'Segoe UI', system-ui, sans-serif"

LIGHT = {
    "surface": "#fcfcfb", "grid": "#e1e0d9", "baseline": "#c3c2b7",
    "ink": "#0b0b0b", "ink2": "#52514e", "ink3": "#898781",
    "series": ["#2a78d6", "#008300", "#e87ba4", "#eda100", "#1baf7a", "#eb6834", "#4a3aa7", "#e34948"],
}
DARK = {
    "surface": "#1a1a19", "grid": "#2c2c2a", "baseline": "#383835",
    "ink": "#ffffff", "ink2": "#c3c2b7", "ink3": "#898781",
    "series": ["#3987e5", "#008300", "#d55181", "#c98500", "#199e70", "#d95926", "#9085e9", "#e66767"],
}
DIVERGING = {
    "light": [[0, "#2a78d6"], [0.5, "#f0efec"], [1, "#e34948"]],
    "dark": [[0, "#3987e5"], [0.5, "#383835"], [1, "#e66767"]],
}
BLUE = {"light": LIGHT["series"][0], "dark": DARK["series"][0]}


def axis_layout(colors: dict) -> dict:
    return {
        "gridcolor": colors["grid"], "zerolinecolor": colors["baseline"], "linecolor": colors["baseline"],
        "showline": True, "ticks": "outside", "tickcolor": colors["grid"],
        "tickfont": {"color": colors["ink3"], "size": 12},
        "title": {"font": {"color": colors["ink2"], "size": 12}},
    }


def build_template(colors: dict) -> dict:
    ax = axis_layout(colors)
    return {
        "layout": {
            "paper_bgcolor": colors["surface"], "plot_bgcolor": colors["surface"],
            "font": {"family": FONT_SANS, "size": 13, "color": colors["ink2"]},
            "colorway": colors["series"],
            "margin": {"l": 56, "r": 24, "t": 24, "b": 48},
            "xaxis": ax, "yaxis": ax,
            "legend": {"bgcolor": "rgba(0,0,0,0)", "bordercolor": "rgba(128,128,128,0.25)",
                       "borderwidth": 1, "font": {"size": 12, "color": colors["ink2"]}},
            "hoverlabel": {"bgcolor": colors["surface"], "bordercolor": "rgba(128,128,128,0.3)",
                           "font": {"family": FONT_SANS, "size": 12, "color": colors["ink"]}},
            "hovermode": "closest",
        }
    }


def templates_payload() -> dict:
    return {"light": build_template(LIGHT), "dark": build_template(DARK)}
