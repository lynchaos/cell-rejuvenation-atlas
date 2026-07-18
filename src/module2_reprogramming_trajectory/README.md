# Module 2 — Reprogramming trajectories and optimal transport

## Biology

Schiebinger et al. (2019, *Cell*) profiled ~250k cells over a 39-point time course of mouse embryonic fibroblast (MEF) → iPSC reprogramming and used **Waddington-OT** to show that reprogramming is not a single path: most cells divert into a stromal/trophic state, and only a narrow trajectory reaches pluripotency. This module rebuilds that analysis and adds a fate-probability layer.

## Data

* GEO: **GSE122662** (processed count matrices + cell metadata)

## What the code does

1. `download_data.py` — fetch processed matrices from GEO.
2. `preprocess.py` — scanpy workflow: QC, normalization, HVGs, PCA, kNN graph, UMAP; day annotations from metadata.
3. `ot_transport.py` — a compact WOT-style engine: PCA-space cost matrices, entropic optimal-transport maps between consecutive days (POT/sinkhorn), growth-rate weighting, and ancestor/descendant distributions.
4. `fate_analysis.py` — classify endpoints (iPSC vs. stromal diversion), compute fate probabilities by pushing day-0 mass through the transport maps, and rank driver genes along the successful trajectory.

## Reproduced vs. published

Expected landmarks: bifurcation into iPSC and stromal branches after day ~8; MET wave; proliferation as the dominant growth signal. Deviations are logged in `reports/module2_report.md`.

## Run

```bash
python -m src.module2_reprogramming_trajectory.download_data --outdir data/module2
python -m src.module2_reprogramming_trajectory.preprocess --input data/module2 --out results/module2/adata.h5ad
python -m src.module2_reprogramming_trajectory.fate_analysis --adata results/module2/adata.h5ad --outdir results/module2
```
