"""Fate analysis: couple consecutive days and score iPSC fate probability.

Usage:
    python -m src.module2_reprogramming_trajectory.fate_analysis \
        --adata results/module2/adata.h5ad --outdir results/module2
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

from src.common.plotting import PALETTE, use_style
from src.module2_reprogramming_trajectory.ot_transport import (
    estimate_growth_rates,
    fate_probability,
    transport_map,
)

# Canonical markers: pluripotency (success) vs stromal diversion
IPSC_MARKERS = ["Nanog", "Pou5f1", "Sox2", "Esrrb", "Dppa5a"]
STROMAL_MARKERS = ["Col1a1", "Col1a2", "Dcn", "Lum", "Sparc"]


def build_couplings(
    adata: ad.AnnData, day_column: str, max_cells_per_day: int = 2000, seed: int = 0
) -> tuple[list[np.ndarray], list[float]]:
    """Compute OT couplings between consecutive days in PCA space."""
    rng = np.random.default_rng(seed)
    days = np.sort(adata.obs[day_column].unique())
    couplings: list[np.ndarray] = []
    for d0, d1 in zip(days[:-1], days[1:]):
        src = adata[adata.obs[day_column] == d0]
        tgt = adata[adata.obs[day_column] == d1]
        if src.n_obs > max_cells_per_day:
            src = src[rng.choice(src.n_obs, max_cells_per_day, replace=False)]
        if tgt.n_obs > max_cells_per_day:
            tgt = tgt[rng.choice(tgt.n_obs, max_cells_per_day, replace=False)]
        g = estimate_growth_rates(src.n_obs, tgt.n_obs, dt=float(d1 - d0))
        couplings.append(
            transport_map(src.obsm["X_pca"], tgt.obsm["X_pca"], growth=g, dt=float(d1 - d0))
        )
    return couplings, [float(d) for d in days]


def score_fates(adata: ad.AnnData) -> pd.Series:
    """Classify cells by marker programs: iPSC, stromal, other."""
    present_ipsc = [g for g in IPSC_MARKERS if g in adata.var_names]
    present_stromal = [g for g in STROMAL_MARKERS if g in adata.var_names]
    if not present_ipsc or not present_stromal:
        raise ValueError("Marker genes missing from matrix; check gene naming.")
    sc.tl.score_genes(adata, present_ipsc, score_name="ipsc_score")
    sc.tl.score_genes(adata, present_stromal, score_name="stromal_score")
    fate = pd.Series("other", index=adata.obs_names)
    fate[adata.obs["ipsc_score"] > adata.obs["ipsc_score"].quantile(0.75)] = "iPSC"
    fate[
        (adata.obs["stromal_score"] > adata.obs["stromal_score"].quantile(0.75))
        & (fate == "other")
    ] = "stromal"
    return fate


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--adata", required=True)
    ap.add_argument("--outdir", default="results/module2")
    ap.add_argument("--day-column", default="day")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    use_style()

    adata = ad.read_h5ad(args.adata)
    adata.obs["fate"] = score_fates(adata)
    couplings, days = build_couplings(adata, args.day_column)

    last_day = days[-1]
    terminal = adata.obs[adata.obs[args.day_column] == last_day]
    ipsc_mask = (terminal["fate"] == "iPSC").to_numpy()
    # Fate probabilities for cells at the first day (subsampled to match coupling size)
    n0 = couplings[0].shape[0]
    probs = fate_probability(couplings, ipsc_mask)
    pd.DataFrame({"cell_index": np.arange(n0), "ipsc_fate_prob": probs}).to_csv(
        outdir / "fate_probabilities.csv", index=False
    )

    # ---- Figures ----
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.4))
    sc.pl.umap(adata, color="fate", ax=axes[0], show=False,
               palette={"iPSC": PALETTE["rejuvenated"], "stromal": PALETTE["aged"],
                        "other": PALETTE["control"]})
    axes[0].set_title("Reprogramming endpoints")
    axes[1].hist(probs, bins=50, color=PALETTE["reprogrammed"])
    axes[1].set_xlabel("P(reach iPSC fate) — day 0 cells")
    axes[1].set_ylabel("cells")
    axes[1].set_title("Waddington-OT fate probabilities")
    fig.suptitle("Schiebinger et al. 2019 (GSE122662) — rebuilt with WOT-style transport")
    fig.tight_layout()
    fig.savefig(outdir / "fate_analysis.pdf")
    fig.savefig(outdir / "fate_analysis.png")

    (outdir / "summary.json").write_text(
        json.dumps(
            {
                "days": days,
                "n_couplings": len(couplings),
                "terminal_fate_counts": terminal["fate"].value_counts().to_dict(),
                "mean_ipsc_fate_prob_day0": float(np.mean(probs)),
            },
            indent=2,
        )
    )
    print(f"[done] results in {outdir}")


if __name__ == "__main__":
    main()
