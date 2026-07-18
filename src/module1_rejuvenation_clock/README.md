# Module 1 — Rejuvenation benchmark: epigenetic clocks on transient-reprogramming data

## Biology

Gill et al. (2022, *eLife*) showed that **maturation-phase transient reprogramming (MPTR)** — exposing human fibroblasts to OSKMLN reprogramming factors for ~13 days, then withdrawing them — restores a youthful transcriptome and **reverses epigenetic age by ~30 years** without loss of cell identity. This module independently re-derives that result from the raw public data.

## Data

| Content | Accession |
|---|---|
| Super-series | GSE165180 |
| EPIC methylation arrays, transient time course | GSE165179 |
| EPIC methylation arrays, Sendai time course | GSE165178 |
| RNA-seq, transient time course | GSE165177 |
| RNA-seq, Sendai time course | GSE165176 |

## What the code does

1. `download_data.py` — pulls processed methylation matrices and sample metadata from GEO programmatically (checksums recorded).
2. `clock.py` — epigenetic-clock machinery: Horvath-style age transformation, elastic-net clock training with nested cross-validation, DNAm-age prediction with bootstrap CIs. Accepts published clock coefficients (CSV: `probe,coef`) or trains a clock on user-supplied training data.
3. `analyze_gill.py` — applies the clock across the MPTR time course; paired bootstrap test of DNAm-age change at day 13 vs. day 0; effect size (Cohen's dz); transcriptomic concordance check via fibroblast-identity markers from the RNA-seq arm.

## Reproduced vs. published

The headline figure (`results/module1/dnam_age_timecourse.pdf`) should show DNAm age rising during early phases, then falling below baseline by ~day 13 of MPTR — the signature result of Gill et al. Any deviation is documented in `reports/module1_report.md` with hypotheses (clock choice, array batch, probe overlap).

## Run

```bash
python -m src.module1_rejuvenation_clock.download_data --outdir data/module1
python -m src.module1_rejuvenation_clock.analyze_gill \
    --beta data/module1/ --metadata data/module1/metadata.csv \
    --coef path/to/clock_coefficients.csv --outdir results/module1
```

Or via Nextflow: `nextflow run main.nf -profile docker --module rejuvenation_clock`
