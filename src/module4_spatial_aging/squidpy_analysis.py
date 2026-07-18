"""Spatial analysis of MERFISH aging-brain data with Squidpy.

Computes:
  * spatial kNN graphs per section
  * neighborhood enrichment (do cell types co-localize more with age?)
  * Moran's I for aging-associated markers (spatial autocorrelation)

Usage:
    python -m src.module4_spatial_aging.squidpy_analysis \
        --input data/module4 --outdir results/module4 \
        --celltype-key cell_type --age-key age
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
import squidpy as sq

from src.common.plotting import use_style

AGING_MARKERS = ["C4b", "Cst3", "B2m", "Gfap", "Aqp1", "Il33", "Vim"]


def load_merfish(indir: Path) -> ad.AnnData:
    h5ad = sorted(indir.rglob("*.h5ad"))
    if h5ad:
        return ad.read_h5ad(h5ad[0])
    counts = sorted(indir.rglob("*counts*.csv*"))
    coords = sorted(indir.rglob("*coord*.csv*"))
    if counts and coords:
        a = sc.read_csv(counts[0]).T
        xy = pd.read_csv(coords[0], index_col=0)
        a.obsm["spatial"] = xy.loc[a.obs_names].to_numpy()
        return a
    raise FileNotFoundError(f"No MERFISH matrix found under {indir}")


def neighborhood_enrichment_by_age(
    adata: ad.AnnData, celltype_key: str, age_key: str, spatial_key: str = "spatial"
) -> dict[str, pd.DataFrame]:
    """sq.gr.nhood_enrichment per age group; returns z-score tables."""
    out: dict[str, pd.DataFrame] = {}
    for age in adata.obs[age_key].unique():
        sub = adata[adata.obs[age_key] == age].copy()
        if sub.n_obs < 500:
            continue
        sq.gr.spatial_neighbors(sub, coord_type="generic", spatial_key=spatial_key)
        sq.gr.nhood_enrichment(sub, cluster_key=celltype_key)
        out[str(age)] = pd.DataFrame(
            sub.uns[f"{celltype_key}_nhood_enrichment"]["zscore"],
            index=sub.obs[celltype_key].cat.categories,
            columns=sub.obs[celltype_key].cat.categories,
        )
    return out


def morans_i(adata: ad.AnnData, genes: list[str]) -> pd.Series:
    present = [g for g in genes if g in adata.var_names]
    if not present:
        return pd.Series(dtype=float)
    sq.gr.spatial_autocorr(adata, mode="moran", genes=present)
    return adata.uns["moranI"]["I"]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", default="data/module4")
    ap.add_argument("--outdir", default="results/module4")
    ap.add_argument("--celltype-key", default="cell_type")
    ap.add_argument("--age-key", default="age")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    use_style()

    adata = load_merfish(Path(args.input))
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)

    enrich = neighborhood_enrichment_by_age(adata, args.celltype_key, args.age_key)
    summary: dict = {"ages_analyzed": list(enrich), "n_cells": int(adata.n_obs)}

    n = max(len(enrich), 1)
    fig, axes = plt.subplots(1, n, figsize=(4.2 * n, 3.8), squeeze=False)
    for ax, (age, z) in zip(axes[0], enrich.items()):
        im = ax.imshow(z, cmap="RdBu_r", vmin=-np.nanmax(np.abs(z)), vmax=np.nanmax(np.abs(z)))
        ax.set_xticks(range(z.shape[1]), z.columns, rotation=90, fontsize=6)
        ax.set_yticks(range(z.shape[0]), z.index, fontsize=6)
        ax.set_title(f"neighborhood enrichment (z) — {age}")
        fig.colorbar(im, ax=ax, shrink=0.7)
        z.to_csv(outdir / f"nhood_enrichment_{age}.csv")
    fig.suptitle("Allen et al. 2023 (GSE207848) — MERFISH aging brain")
    fig.tight_layout()
    fig.savefig(outdir / "nhood_enrichment.pdf")
    fig.savefig(outdir / "nhood_enrichment.png")

    mi = morans_i(adata, AGING_MARKERS)
    mi.to_csv(outdir / "morans_i_aging_markers.csv")
    summary["morans_i"] = mi.to_dict() if len(mi) else {}
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
