"""Transparent baseline to the deep model: per-modality PCA + CCA consensus.

A lightweight MOFA-style factor view. Where scVI is powerful but opaque, this
module provides the interpretable cross-check reviewers ask for: do the deep
latent factors and the linear consensus axes tell the same story?

Usage:
    python -m src.module3_multiomics_integration.factor_model \
        --modality-a results/module3/tms_pca.csv \
        --modality-b results/module3/browder_pca.csv --outdir results/module3
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.cross_decomposition import CCA
from sklearn.decomposition import PCA


def modality_pcs(expr: pd.DataFrame, n_pcs: int = 30) -> pd.DataFrame:
    """PCA of a log-normalized expression matrix (cells x genes)."""
    x = expr.to_numpy(dtype=float)
    if np.nanmin(x) >= 0:  # count-like data: variance-stabilize
        x = np.log1p(x)
    x = (x - x.mean(axis=0)) / (x.std(axis=0) + 1e-8)
    x = np.nan_to_num(x)
    pcs = PCA(n_components=n_pcs, random_state=0).fit_transform(x)
    return pd.DataFrame(pcs, index=expr.index, columns=[f"PC{i+1}" for i in range(n_pcs)])


def consensus_axes(pcs_a: pd.DataFrame, pcs_b: pd.DataFrame, n_components: int = 10) -> pd.DataFrame:
    """CCA between two modalities' PC spaces over shared cells -> shared axes."""
    common = pcs_a.index.intersection(pcs_b.index)
    if len(common) < n_components * 5:
        raise ValueError(f"Too few shared cells ({len(common)}) for CCA.")
    cca = CCA(n_components=n_components)
    xa, xb = cca.fit_transform(pcs_a.loc[common], pcs_b.loc[common])
    corrs = [float(np.corrcoef(xa[:, i], xb[:, i])[0, 1]) for i in range(n_components)]
    return pd.DataFrame({"axis": range(1, n_components + 1), "canonical_corr": corrs})


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--modality-a", required=True, help="CSV: cells x genes (modality A)")
    ap.add_argument("--modality-b", required=True, help="CSV: cells x genes (modality B)")
    ap.add_argument("--outdir", default="results/module3")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    a = pd.read_csv(args.modality_a, index_col=0)
    b = pd.read_csv(args.modality_b, index_col=0)
    axes = consensus_axes(modality_pcs(a), modality_pcs(b))
    axes.to_csv(outdir / "cca_consensus_axes.csv", index=False)
    print(axes.to_string(index=False))


if __name__ == "__main__":
    main()
