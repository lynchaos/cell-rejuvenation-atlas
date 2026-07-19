"""Standard scanpy preprocessing for the reprogramming time course.

Usage:
    python -m src.module2_reprogramming_trajectory.preprocess \
        --input data/module2 --out results/module2/adata.h5ad
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

from src.module2_reprogramming_trajectory.fate_analysis import IPSC_MARKERS, STROMAL_MARKERS


def parse_sample_name(name: str) -> tuple[float, str]:
    """'GSM2836276_D9-1-2i' -> (9.0, '2i'); 'GSM2836269_D2-2' -> (2.0, 'serum');
    'GSM3195648_D0_Dox_C1_gene_bc_mat' -> (0.0, 'dox')."""
    m = re.search(r"_D(\d+(?:\.\d+)?)", name)
    day = float(m.group(1)) if m else np.nan
    low = name.lower()
    cond = "dox" if "dox" in low else ("2i" if "2i" in low else "serum")
    return day, cond


def _read_triplet(mtx: Path, indir: Path) -> ad.AnnData:
    """Read one GEO 10x triplet (matrix + barcodes + genes) into AnnData."""
    prefix = mtx.name.split(".matrix.mtx")[0].split(".mtx")[0]
    a = sc.read_mtx(mtx).T  # genes x cells -> cells x genes
    barcodes = next(indir.rglob(f"{prefix}.barcodes.tsv*"), None)
    genes = next(indir.rglob(f"{prefix}.genes.tsv*"), None)
    if barcodes is not None:
        a.obs_names = pd.read_csv(barcodes, header=None)[0].to_numpy()
    if genes is not None:
        g = pd.read_csv(genes, header=None, sep="\t")
        a.var_names = g[1].to_numpy() if g.shape[1] > 1 else g[0].to_numpy()
    return _annotate(a, prefix)


def _annotate(a: ad.AnnData, name: str) -> ad.AnnData:
    day, cond = parse_sample_name(name)
    a.obs["day"] = day
    a.obs["condition"] = cond
    a.obs["sample"] = name
    return a


def load_matrix(indir: Path) -> ad.AnnData:
    """Load the Schiebinger reprogramming time course from the download dir.

    GSE122662_RAW.tar bundles TWO experiments: the Dox-inducible time course
    (D0..D18 every 12 h, *_gene_bc_mat.h5 — the canonical WOT dataset) and a
    smaller serum/2i media comparison (10x triplets). They are never mixed:
    the Dox h5 series takes precedence; triplets/csv are fallbacks.
    """
    h5 = sorted(indir.rglob("*_gene_bc_mat.h5"))
    mtx = sorted(indir.rglob("*.matrix.mtx*")) or sorted(indir.rglob("*.mtx*"))
    csv = sorted(indir.rglob("*counts*.csv*"))
    if h5:
        adatas = [_annotate(sc.read_10x_h5(f), f.name.split("_gene_bc_mat")[0]) for f in h5]
    elif mtx:
        adatas = [_read_triplet(m, indir) for m in mtx]
    elif csv:
        return sc.read_csv(csv[0]).T
    else:
        raise FileNotFoundError(f"No recognized count matrix under {indir}")
    for a in adatas:
        a.var_names_make_unique()
    combined = ad.concat(adatas, join="inner") if len(adatas) > 1 else adatas[0]
    combined = combined[~combined.obs["day"].isna()].copy()
    combined.obs_names_make_unique()  # 10x barcodes repeat across samples
    return combined


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
    # Subset before scale(): scaling densifies the matrix (250k x 20k would
    # not fit in memory); PCA uses highly_variable genes anyway. Fate marker
    # genes are kept explicitly — fate_analysis.score_fates needs them.
    keep = adata.var["highly_variable"] | adata.var_names.isin(IPSC_MARKERS + STROMAL_MARKERS)
    adata = adata[:, keep].copy()
    sc.pp.scale(adata, max_value=10)
    sc.pp.pca(adata, n_comps=n_pcs)
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
