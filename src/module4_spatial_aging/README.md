# Module 4 — Spatial aging of the mouse brain (MERFISH)

## Biology

Allen et al. (2023, *Cell*) used MERFISH to map how cell types and their **spatial organization** change in the aging mouse brain — e.g., activation of glial states with stereotyped spatial distributions. This module rebuilds the spatial-neighborhood analysis.

## Data

* GEO: **GSE207848** (MERFISH count matrices + spatial coordinates + cell metadata)

## What the code does

1. `download_data.py` — fetch processed MERFISH matrices.
2. `squidpy_analysis.py` — build spatial kNN graphs, compute neighborhood enrichment, Moran's I spatial autocorrelation of aging markers, and spatially-variable gene detection per age group.

## Run

```bash
python -m src.module4_spatial_aging.download_data --outdir data/module4
python -m src.module4_spatial_aging.squidpy_analysis --input data/module4 --outdir results/module4
```
