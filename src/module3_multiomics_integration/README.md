# Module 3 — Cross-study integration: aging atlas × in-vivo partial reprogramming

## Biology

Does **in-vivo OSK partial reprogramming** push aged tissue toward a young state? We combine two peer-reviewed resources — the **Tabula Muris Senis** aging atlas (2020, *Nature*) and **Browder et al.** (2022, *Nature Aging*) in-vivo partial reprogramming data — and ask whether tissue from OSK-treated old animals looks transcriptomically younger than matched controls.

## Data

| Resource | Accession |
|---|---|
| Browder et al. 2022, 7-month cohort — **bulk** RNA-seq, 6 tissues | GSE190983 (`GSE190983_count.tsv.gz`) |
| Tabula Muris Senis — annotated FACS object | GSE149590 (`GSM4505405` h5ad) |

## What the code does

1. `download_data.py` — fetch the Browder count matrix and the annotated TMS object.
2. `aging_signature_score.py` — derive a transcriptomic **aging axis** from TMS (mean log-CPM old − young per gene), map Browder Ensembl IDs to symbols (mygene.info, cached), and project each bulk sample onto the axis (`rejuvenation_score.project_onto_signature`). Per tissue: Welch t-test 4F vs. control, bootstrap CI, Cohen's d, BH FDR. Outputs `tissue_scores.csv`, `sample_scores.csv`, `aging_score.png`, `summary.json`.

Note: Browder GSE190983 is **bulk** RNA-seq (genes × samples), so single-cell integration with scVI is not applicable; the signature-projection design is the statistically valid way to combine the two studies. `rejuvenation_score.py` (signature algebra) and `factor_model.py` (PCA+CCA consensus baseline) remain as shared, unit-tested helpers.

## Reproduced vs. published

Expected: negative Δ (4F younger than control) in the tissues where Browder et al. report age-reversal, with heterogeneity across tissues. All deviations logged in `reports/module3_report.md`.

## Run

```bash
python -m src.module3_multiomics_integration.download_data --outdir data/module3
python -m src.module3_multiomics_integration.aging_signature_score \
    --tms data/module3/tabula_muris_senis/tms.h5ad \
    --browder data/module3/browder_long7m/GSE190983_count.tsv.gz \
    --outdir results/module3
```
