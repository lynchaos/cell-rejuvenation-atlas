"""Standard scanpy preprocessing for the reprogramming time course.

Usage:
    python -m src.module2_reprogramming_trajectory.preprocess \
        --input data/module2 --out results/module2/adata.h5ad
"""
from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import scanpy as sc


def load_matrix(indir: Path) -> ad.AnnData:
    """Load 10x-style mtx / h5 / csv counts from the download directory."""
    h5 = sorted(indir.rglob("*.h5"))
    mtx = sorted(indir.rglob("*.mtx*"))
    csv = sorted(indir.rglob("*counts*.csv*"))
    if h5:
        return sc.read_10x_h5(h5[0])
    if mtx:
        a = sc.read_mtx(mtx[0]).T
        return a
    if csv:
        return sc.read_csv(csv[0]).T
    raise FileNotFoundError(f"No recognized count matrix under {indir}")


def preprocess(adata: ad.AnnData, n_hvg: int = 2000, n_pcs: int = 50) -> ad.AnnData:
    """QC -> normalize -> HVG -> PCA -> neighbors -> UMAP."""
    adata.var_names_make_unique()
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
    adata = adata[adata.obs["pct_counts_mt"] < 20].copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg, flavor="seurat")
    sc.pp.scale(adata, max_value=10)
    sc.pp.pca(adata, n_comps=n_pcs, use_highly_variable=True)
    sc.pp.neighbors(adata, n_pcs=min(n_pcs, 30))
    sc.tl.umap(adata)
    return adata


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", default="data/module2")
    ap.add_argument("--out", default="results/module2/adata.h5ad")
    ap.add_argument("--day-column", default="day", help="obs column with time point")
    args = ap.parse_args()

    adata = load_matrix(Path(args.input))
    adata = preprocess(adata)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    adata.write(args.out)
    print(f"[done] {adata.n_obs} cells x {adata.n_vars} genes -> {args.out}")


if __name__ == "__main__":
    main()
