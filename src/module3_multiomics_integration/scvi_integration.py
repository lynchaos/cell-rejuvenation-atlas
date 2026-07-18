"""scVI integration of Tabula Muris Senis with Browder et al. reprogramming data.

scVI (Lopez et al. 2018, Nature Methods) fits a deep generative model (VAE)
to raw counts with batch as a covariate — the generative-modeling component
of this project. Heavy dependencies (torch/scvi-tools) are imported lazily so
the pure scoring functions stay light and testable.

Usage:
    python -m src.module3_multiomics_integration.scvi_integration \
        --tms data/module3/tabula_muris_senis/tms.h5ad \
        --browder data/module3/browder_long7m/browder.h5ad \
        --outdir results/module3
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from src.common.plotting import PALETTE, use_style
from src.module3_multiomics_integration.rejuvenation_score import (
    aging_signature,
    latent_shift,
    project_onto_signature,
)


def integrate(tms: ad.AnnData, browder: ad.AnnData, max_epochs: int = 400) -> ad.AnnData:
    """Train scVI on the concatenated object; return object with latent + UMAP."""
    import scvi  # lazy: heavy dependency
    import scanpy as sc

    common_genes = tms.var_names.intersection(browder.var_names)
    tms = tms[:, common_genes].copy()
    browder = browder[:, common_genes].copy()
    tms.obs["study"] = "tabula_muris_senis"
    browder.obs["study"] = "browder_partial_reprogramming"
    adata = ad.concat([tms, browder], label="dataset", join="inner")

    scvi.model.SCVI.setup_anndata(adata, batch_key="study")
    model = scvi.model.SCVI(adata, n_latent=30, gene_likelihood="nb")
    model.train(max_epochs=max_epochs, early_stopping=True)
    adata.obsm["X_scVI"] = model.get_latent_representation()

    sc.pp.neighbors(adata, use_rep="X_scVI")
    sc.tl.umap(adata)
    return adata


def transfer_labels(adata: ad.AnnData, label_key: str, ref_study: str = "tabula_muris_senis") -> None:
    """kNN label transfer from the annotated atlas to reprogramming cells."""
    from sklearn.neighbors import KNeighborsClassifier

    is_ref = adata.obs["study"] == ref_study
    x_ref, x_qry = adata.obsm["X_scVI"][is_ref], adata.obsm["X_scVI"][~is_ref]
    knn = KNeighborsClassifier(n_neighbors=15, weights="distance").fit(
        x_ref, adata.obs.loc[is_ref, label_key]
    )
    adata.obs.loc[~is_ref, f"transferred_{label_key}"] = knn.predict(x_qry)
    adata.obs.loc[~is_ref, f"transfer_confidence"] = knn.predict_proba(x_qry).max(axis=1)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tms", required=True, help="Tabula Muris Senis h5ad")
    ap.add_argument("--browder", required=True, help="Browder et al. h5ad")
    ap.add_argument("--outdir", default="results/module3")
    ap.add_argument("--age-key", default="age")
    ap.add_argument("--celltype-key", default="cell_ontology_class")
    ap.add_argument("--max-epochs", type=int, default=400)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    use_style()

    tms = ad.read_h5ad(args.tms)
    browder = ad.read_h5ad(args.browder)
    adata = integrate(tms, browder, max_epochs=args.max_epochs)
    transfer_labels(adata, args.celltype_key)

    # ---- Rejuvenation scoring per cell type ----
    results = []
    expr = pd.DataFrame(
        adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X,
        index=adata.obs_names, columns=adata.var_names,
    )
    for ct in adata.obs[args.celltype_key].dropna().unique():
        mask = adata.obs[args.celltype_key] == ct
        sub = expr.loc[mask]
        groups = adata.obs.loc[mask, args.age_key]
        if groups.nunique() < 2 or mask.sum() < 50:
            continue
        sig = aging_signature(sub, groups)
        adata.obs.loc[mask, "aging_score"] = project_onto_signature(sub, sig)
        results.append({"cell_type": ct, "n_cells": int(mask.sum())})
    pd.DataFrame(results).to_csv(outdir / "scored_cell_types.csv", index=False)

    # ---- Latent shift: old+OSK relative to old control ----
    latent = pd.DataFrame(adata.obsm["X_scVI"], index=adata.obs_names)
    browder_obs = adata.obs[adata.obs["study"] == "browder_partial_reprogramming"]
    condition_cols = [c for c in browder_obs.columns if "treat" in c.lower() or "osk" in c.lower()]
    summary = {"n_cells_integrated": int(adata.n_obs), "cell_types_scored": len(results)}
    if condition_cols:
        cond = condition_cols[0]
        treated = browder_obs[browder_obs[cond].astype(str).str.contains("OSK|osk", regex=True)]
        control = browder_obs[~browder_obs.index.isin(treated.index)]
        if len(treated) and len(control):
            shift = latent_shift(latent.loc[control.index], latent.loc[treated.index])
            summary["mean_latent_shift_norm_osk_vs_control"] = float(np.linalg.norm(shift))

    # ---- Figure ----
    import scanpy as sc

    fig, axes = plt.subplots(1, 3, figsize=(13, 3.6))
    sc.pl.umap(adata, color="study", ax=axes[0], show=False)
    sc.pl.umap(adata, color=args.celltype_key, ax=axes[1], show=False, legend_loc="on data")
    if "aging_score" in adata.obs:
        sc.pl.umap(adata, color="aging_score", ax=axes[2], show=False, cmap="RdBu_r")
    fig.suptitle("scVI integration: Tabula Muris Senis x Browder 2022 (GSE149590 x GSE190986)")
    fig.tight_layout()
    fig.savefig(outdir / "scvi_integration.pdf")
    fig.savefig(outdir / "scvi_integration.png")

    adata.write(outdir / "integrated.h5ad")
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
