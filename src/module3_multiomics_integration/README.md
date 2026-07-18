# Module 3 — Cross-study integration: aging atlas × in-vivo partial reprogramming

## Biology

Does **in-vivo OSK partial reprogramming** push aged cells toward a young state, cell type by cell type? We integrate two peer-reviewed resources — the **Tabula Muris Senis** aging atlas (2020, *Nature*) and **Browder et al.** (2022, *Nature Aging*) in-vivo partial reprogramming data — and ask, in a shared latent space, whether reprogrammed cells from old animals move toward young-animal states.

## Data

| Resource | Accession |
|---|---|
| Browder et al. 2022 (super-series) | GSE190986 |
| Browder subseries (1-, 7-, 10-month treatments) | GSE190984, GSE190983, GSE190985 |
| Tabula Muris Senis | GSE149590 |

## What the code does

1. `download_data.py` — fetch processed matrices for both studies.
2. `scvi_integration.py` — train **scVI** (deep generative model: VAE over counts with batch as a covariate) on the concatenated datasets; harmonize cell-type labels via kNN label transfer; produce a shared UMAP.
3. `rejuvenation_score.py` — per cell type: an *aging signature* (young vs. old DE in Tabula Muris Senis) projected onto reprogrammed cells; tests whether OSK shifts cells along the young direction (scGen-style latent arithmetic included).
4. `factor_model.py` — a lightweight MOFA-style consensus: per-modality PCA followed by CCA to expose shared axes of variation (documented as the transparent baseline against the deep model).

## Reproduced vs. published

Expected: partial reprogramming moves aged cells toward young states in a subset of cell types (e.g., kidney mesangial/adipose lineages per Browder et al.), with effect sizes heterogenous across lineages. All deviations logged in `reports/module3_report.md`.

## Run

```bash
python -m src.module3_multiomics_integration.download_data --outdir data/module3
python -m src.module3_multiomics_integration.scvi_integration \
    --tms data/module3/tabula_muris_senis.h5ad \
    --browder data/module3/browder.h5ad --outdir results/module3
```
