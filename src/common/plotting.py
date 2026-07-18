"""Consistent figure style across all modules."""
from __future__ import annotations

import matplotlib.pyplot as plt

PALETTE = {
    "control": "#8c8c8c",
    "aged": "#d1495b",
    "reprogrammed": "#2e86ab",
    "rejuvenated": "#1b998b",
    "young": "#5bc0eb",
}


def use_style() -> None:
    plt.rcParams.update(
        {
            "figure.dpi": 120,
            "savefig.dpi": 300,
            "savefig.bbox": "tight",
            "font.size": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": False,
            "legend.frameon": False,
        }
    )
