# Datasets

All data reanalyzed here are public and peer-reviewed. Accessions were verified against NCBI GEO / ProteomeXchange at project creation.

| Module | Study | Journal (year) | Data type | Repository / accession |
|---|---|---|---|---|
| 1 | Gill et al., "Multi-omic rejuvenation of human cells by maturation phase transient reprogramming" | eLife (2022) | Human fibroblast RNA-seq + Illumina EPIC methylation arrays, transient (MPTR) and Sendai reprogramming time courses | GEO: **GSE165180** (super-series); RNA-seq **GSE165176**, **GSE165177**; methylation **GSE165178**, **GSE165179** |
| 1 | Horvath, "DNA methylation age of human tissues and cell types" | Genome Biology (2013) | Methylation clock training data (multiple tissues) | GEO (see module README for the series list used) |
| 2 | Schiebinger et al., "Optimal-transport analysis of single-cell gene expression identifies developmental trajectories in reprogramming" | Cell (2019) | ~250k-cell scRNA-seq time course of MEF → iPSC reprogramming (39 time points) | GEO: **GSE122662** |
| 3 | Browder et al., "In vivo partial reprogramming alters manifestation of age-associated molecular changes in aging mice" | Nature Aging (2022) | snRNA-seq / scRNA-seq of mouse tissues, OSK partial reprogramming, multiple treatment durations | GEO: **GSE190986** (super-series); **GSE190983**, **GSE190984**, **GSE190985** |
| 3 | Tabula Muris Senis Consortium, "A single-cell transcriptomic atlas characterizes ageing tissues in the mouse" | Nature (2020) | scRNA-seq (droplet + FACS) across mouse tissues and ages | GEO: **GSE149590** |
| 4 | Allen et al., "Molecular and spatial signatures of mouse brain aging at single-cell resolution" | Cell (2023) | MERFISH spatial transcriptomics, mouse brain, young vs. aged | GEO: **GSE207848** |
| 5 | Basisty et al., "A proteomic atlas of senescence-associated secretomes for aging biomarker development" | PLoS Biology (2020) | DIA LC-MS/MS proteomics of senescent-cell secretomes (SASP Atlas) | ProteomeXchange: **PXD013721** |

## Retrieval policy

* `src/module*/download_data.py` scripts fetch processed matrices directly (GEO supplementary files / cellxgene / PRIDE) — no manual browser steps.
* Raw sequencing reads (SRA) are **optional** and only needed for the from-fastq Nextflow demonstration path (`--from_fastq true`).
* Every download records SHA-256 checksums into `reports/data_manifest.json` for provenance.

## Licensing note

Code in this repository is MIT-licensed. The datasets remain the property of their respective studies; use them under the terms of the hosting repositories and cite the original publications.
