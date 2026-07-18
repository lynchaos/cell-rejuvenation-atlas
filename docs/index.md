# Cell Rejuvenation Atlas — Documentation

Welcome. This site documents the design decisions behind each module.

## Philosophy

1. **Real, peer-reviewed data only.** Every analysis starts from a published study (see `DATASETS.md`).
2. **Reproduce, then extend.** Each module first rebuilds a published result, then adds something the paper didn't do (cross-study integration, fate probabilities, interactive exploration).
3. **Engineering is a scientific value.** Containers, CI, checksums, and tests are not polish — they are what makes a result trustworthy.

## Module guides

| Module | Question | Key methods |
|---|---|---|
| 1 · Rejuvenation clock | Does MPTR reverse DNAm age? | elastic-net clocks, Horvath transform, paired bootstrap |
| 2 · Trajectories | Which cells complete reprogramming? | scanpy, entropic optimal transport (WOT-style), fate probabilities |
| 3 · Integration | Does OSK move aged cells toward young states? | scVI (VAE), kNN label transfer, signature projection, PCA+CCA baseline |
| 4 · Spatial aging | How does brain architecture change with age? | Squidpy, neighborhood enrichment, Moran's I |
| 5 · SASP proteomics | What is the protein-level core SASP? | DIA-MS reanalysis, per-inducer DE, multi-inducer core set |

## Operations

* [AWS setup](aws_setup.md) — one-time cloud setup for Batch runs
* `nextflow run main.nf -profile test,docker` — smoke test in minutes
* `streamlit run app/streamlit_app.py` — interactive explorer
