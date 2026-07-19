"""Spatial analysis of MERFISH aging-brain data with Squidpy.

Computes, per age group:
  * neighborhood enrichment (do cell types co-localize more with age?)
  * Moran's I for aging-associated markers (spatial autocorrelation)

All spatial graphs are built WITHIN a (donor, slice) section and the resulting
z-scores / Moran's I values are averaged across sections of the same age.
Pooling cells across sections would create false spatial neighbors between
cells that were never on the same tissue slice.

Usage:
    python -m src.module4_spatial_aging.squidpy_analysis \
        --input data/module4 --outdir results/module4 \
        --celltype-key cell_type --age-key age
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq

from src.common.plotting import use_style

AGING_MARKERS = ["C4b", "Cst3", "B2m", "Gfap", "Aqp1", "Il33", "Vim"]
GROUP_KEYS = ("donor_id", "slice")
MIN_CELLS_PER_SECTION = 500


def load_merfish(indir: Path) -> ad.AnnData:
    preferred = sorted(indir.rglob("allen_merfish_brain.h5ad"))
    h5ad = preferred or sorted(indir.rglob("*.h5ad"))
    if h5ad:
        adata = ad.read_h5ad(h5ad[0])
        _ensure_spatial_key(adata)
        return adata
    counts = sorted(indir.rglob("*counts*.csv*"))
    coords = sorted(indir.rglob("*coord*.csv*"))
    if counts and coords:
        a = sc.read_csv(counts[0]).T
        xy = pd.read_csv(coords[0], index_col=0)
        a.obsm["spatial"] = xy.loc[a.obs_names].to_numpy()
        return a
    raise FileNotFoundError(f"No MERFISH matrix found under {indir}")


def _ensure_spatial_key(adata: ad.AnnData) -> None:
    """Map whatever coordinate array the object ships to obsm['spatial']."""
    if "spatial" in adata.obsm:
        return
    for key in list(adata.obsm.keys()):
        if "spatial" in key.lower() or "coord" in key.lower():
            adata.obsm["spatial"] = np.asarray(adata.obsm[key])
            return
    raise KeyError(f"No spatial coordinates found in .obsm (keys: {list(adata.obsm.keys())})")


def _age_order(values: pd.Series) -> list[str]:
    """Order age labels by their numeric prefix: 4wk < 24wk < 90wk."""
    def key(v: str) -> tuple[int, str]:
        m = re.match(r"\s*(\d+)", v)
        return (int(m.group(1)) if m else 0, v)

    return sorted({str(v) for v in values}, key=key)


def _iter_sections(
    adata: ad.AnnData,
    age: str,
    age_key: str,
    group_keys: tuple[str, ...] = GROUP_KEYS,
    min_cells: int = MIN_CELLS_PER_SECTION,
):
    """Yield (section_name, section_adata) for one age group, per (donor, slice)."""
    sub = adata[adata.obs[age_key].astype(str) == age]
    keys = [k for k in group_keys if k in sub.obs.columns]
    if not keys:
        if sub.n_obs >= min_cells:
            yield "all", sub.copy()
        return
    for name, idx in sub.obs.groupby(list(keys), observed=True).groups.items():
        if len(idx) >= min_cells:
            yield name, sub[sub.obs_names.isin(idx)].copy()


def neighborhood_enrichment_by_age(
    adata: ad.AnnData,
    celltype_key: str,
    age_key: str,
    group_keys: tuple[str, ...] = GROUP_KEYS,
) -> dict[str, pd.DataFrame]:
    """Mean sq.gr.nhood_enrichment z-scores per age, averaged over sections."""
    cats = pd.Categorical(adata.obs[celltype_key]).categories
    out: dict[str, pd.DataFrame] = {}
    for age in _age_order(adata.obs[age_key]):
        zs = []
        for _, section in _iter_sections(adata, age, age_key, group_keys):
            section.obs[celltype_key] = section.obs[celltype_key].astype("category")
            sq.gr.spatial_neighbors(section, coord_type="generic", spatial_key="spatial")
            sq.gr.nhood_enrichment(section, cluster_key=celltype_key)
            z = pd.DataFrame(
                section.uns[f"{celltype_key}_nhood_enrichment"]["zscore"],
                index=section.obs[celltype_key].cat.categories,
                columns=section.obs[celltype_key].cat.categories,
            )
            zs.append(z.reindex(index=cats, columns=cats).to_numpy(dtype=float))
        if zs:
            out[age] = pd.DataFrame(np.nanmean(np.stack(zs), axis=0), index=cats, columns=cats)
    return out


def morans_i_by_age(
    adata: ad.AnnData,
    genes: list[str],
    age_key: str,
    group_keys: tuple[str, ...] = GROUP_KEYS,
) -> pd.DataFrame:
    """Moran's I per gene per age, averaged over sections (genes x ages)."""
    present = [g for g in genes if g in adata.var_names]
    if not present:
        return pd.DataFrame()
    cols: dict[str, pd.Series] = {}
    for age in _age_order(adata.obs[age_key]):
        vals = []
        for _, section in _iter_sections(adata, age, age_key, group_keys):
            sq.gr.spatial_neighbors(section, coord_type="generic", spatial_key="spatial")
            sq.gr.spatial_autocorr(section, mode="moran", genes=present)
            vals.append(section.uns["moranI"]["I"].reindex(present))
        if vals:
            cols[age] = pd.concat(vals, axis=1).mean(axis=1)
    return pd.DataFrame(cols)


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
    n_sections = {
        age: sum(1 for _ in _iter_sections(adata, age, args.age_key))
        for age in _age_order(adata.obs[args.age_key])
    }
    summary: dict = {
        "ages_analyzed": list(enrich),
        "n_cells": int(adata.n_obs),
        "n_sections_per_age": n_sections,
    }

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

    mi = morans_i_by_age(adata, AGING_MARKERS, args.age_key)
    if not mi.empty:
        mi.to_csv(outdir / "morans_i_aging_markers.csv")
        summary["morans_i"] = mi.to_dict()
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
