# Module 5 — The senescence-associated secretory phenotype (proteomics)

## Biology

Basisty et al. (2020, *PLoS Biology*) built the **SASP Atlas**: DIA mass-spectrometry of secretomes from senescent human cells across inducers (irradiation, RAS, drugs) and cell types. This module reanalyzes the protein-level data and connects it back to the transcriptomic modules: which SASP factors are robust at the protein level, and do they form a cross-modality aging signature?

## Data

* ProteomeXchange: **PXD013721** (DIA-NN-processed protein matrices via the authors' SASP Atlas portal / PRIDE)

## What the code does

1. `download_data.py` — fetch the processed protein matrix (PRIDE/portal).
2. `sasp_analysis.py` — log2FC of senescent vs. control per inducer; core SASP = proteins upregulated across >= 2 inducers (FDR < 0.05); compares against the published core SASP; enrichment via gseapy (secretory/ECM terms); correlation with the transcriptomic aging signature from module 3.

## Run

```bash
python -m src.module5_proteomics_sasp.download_data --outdir data/module5
python -m src.module5_proteomics_sasp.sasp_analysis --input data/module5 --outdir results/module5
```
