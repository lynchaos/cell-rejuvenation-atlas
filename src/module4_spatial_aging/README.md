# Module 4 — Spatial aging of the mouse brain (MERFISH)

## Biology

Allen et al. (2023, *Cell*) used MERFISH to map how cell types and their **spatial organization** change in the aging mouse brain — e.g., activation of glial states with stereotyped spatial distributions. This module rebuilds the spatial-neighborhood analysis.

## Data

* Annotated object (378,918 cells × 374 genes, with spatial coordinates, donor, slice and age labels) via **CELLxGENE** collection 31937775-0602-4e52-a799-b6acdd2bac2e ("BrainAgingSpatialAtlas_MERFISH"). GEO **GSE207848** ships only raw per-run archives without unified cell typing, so the curated object is used.

## What the code does

1. `download_data.py` — fetch the annotated h5ad from CELLxGENE.
2. `squidpy_analysis.py` — spatial kNN graphs, neighborhood enrichment and Moran's I for aging markers, per age group. **All graphs are built within a (donor, slice) section** and statistics averaged across sections — pooling cells across slices would create false spatial neighbors.

## Run

```bash
python -m src.module4_spatial_aging.download_data --outdir data/module4
python -m src.module4_spatial_aging.squidpy_analysis --input data/module4 --outdir results/module4
```
