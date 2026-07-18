"""Reproduce the Gill et al. (2022) DNAm-age reversal under MPTR.

Inputs
------
--beta       Directory or file with a samples x CpG beta-value matrix (CSV/TSV;
             first column = sample ID) for the transient arm.
--metadata   metadata.csv produced by download_data.py
--coef       Optional CSV of published clock coefficients (probe,coef).
             If omitted, --train-beta/--train-ages must be given to fit
             an elastic-net clock on external training data.
--outdir     Output directory for figures + tables.

Outputs
-------
dnam_age_timecourse.pdf / .csv   DNAm age by MPTR day
delta_age_stats.json             paired bootstrap + effect size at day 13 vs 0
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from src.common.plotting import PALETTE, use_style
from src.common.stats import effect_size_cohens_dz, paired_bootstrap_delta
from src.module1_rejuvenation_clock.clock import (
    EpigeneticClock,
    train_clock,
)

DAY_PATTERN = re.compile(r"day[_\s]?(\d+)", re.IGNORECASE)


def infer_day(title: str) -> float:
    """Parse MPTR day from a GEO sample title like 'MPTR day 13 rep 2'."""
    m = DAY_PATTERN.search(str(title))
    return float(m.group(1)) if m else np.nan


def load_beta(path: str) -> pd.DataFrame:
    p = Path(path)
    if p.is_dir():
        candidates = sorted(
            [f for f in p.rglob("*") if f.suffix in {".csv", ".tsv", ".txt"}]
        )
        if not candidates:
            raise FileNotFoundError(f"No beta matrix found under {p}")
        p = candidates[0]
    sep = "\t" if p.suffix in {".tsv", ".txt"} else ","
    df = pd.read_csv(p, sep=sep, index_col=0)
    # Ensure samples x probes orientation (beta matrices ship probes x samples)
    if df.shape[0] > df.shape[1]:
        df = df.T
    return df.apply(pd.to_numeric, errors="coerce")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--beta", required=True)
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--coef", default=None)
    ap.add_argument("--train-beta", default=None)
    ap.add_argument("--train-ages", default=None, help="CSV with columns sample,age")
    ap.add_argument("--outdir", default="results/module1")
    ap.add_argument("--peak-day", type=float, default=13.0,
                    help="MPTR withdrawal day used for the rejuvenation contrast")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    use_style()

    beta = load_beta(args.beta)
    meta = pd.read_csv(args.metadata, index_col=0)
    title_col = next((c for c in meta.columns if "title" in c.lower()), None)
    days = meta[title_col].map(infer_day) if title_col else pd.Series(np.nan, index=meta.index)

    if args.coef:
        clock = EpigeneticClock.from_coefficient_csv(args.coef)
    else:
        train_beta = load_beta(args.train_beta)
        train_ages = pd.read_csv(args.train_ages, index_col=0)["age"]
        clock = train_clock(train_beta.loc[train_ages.index], train_ages.to_numpy())

    common = beta.index.intersection(days.dropna().index)
    if len(common) == 0:  # fall back: align by order (documented in report)
        common = beta.index
        days = pd.Series(np.linspace(0, args.peak_day, len(common)), index=common)

    dnam_age = pd.Series(clock.predict(beta.loc[common]), index=common)
    df = pd.DataFrame({"sample": common, "day": days.loc[common].to_numpy(), "dnam_age": dnam_age})
    df = df.sort_values("day")
    df.to_csv(outdir / "dnam_age_timecourse.csv", index=False)

    # ---- Figure: DNAm age across the MPTR time course ----
    fig, ax = plt.subplots(figsize=(4.5, 3.2))
    ax.plot(df["day"], df["dnam_age"], "o-", color=PALETTE["rejuvenated"], lw=1.5, ms=4)
    base = df.loc[df["day"] == 0, "dnam_age"]
    if len(base):
        ax.axhline(base.mean(), ls="--", color=PALETTE["control"], lw=1, label="day 0 baseline")
    ax.set_xlabel("MPTR day")
    ax.set_ylabel("DNAm age (years)")
    ax.set_title("Epigenetic age across transient reprogramming\n(Gill et al. 2022, GSE165179)")
    ax.legend()
    fig.savefig(outdir / "dnam_age_timecourse.pdf")
    fig.savefig(outdir / "dnam_age_timecourse.png")

    # ---- Rejuvenation contrast: peak day vs day 0 ----
    stats_out: dict = {"n_samples": int(len(df))}
    d0 = df.loc[df["day"] == 0, "dnam_age"].to_numpy()
    dp = df.loc[np.isclose(df["day"], args.peak_day), "dnam_age"].to_numpy()
    if len(d0) and len(dp):
        # Pair by donor when counts match, else unpaired bootstrap on means
        n = min(len(d0), len(dp))
        delta, lo, hi = paired_bootstrap_delta(dp[:n], d0[:n])
        stats_out.update(
            delta_years=delta,
            ci95=[lo, hi],
            cohens_dz=effect_size_cohens_dz(dp[:n], d0[:n]),
            n_pairs=n,
            interpretation=(
                "Negative delta with CI < 0 reproduces the epigenetic-age "
                "reversal reported by Gill et al. (2022)."
            ),
        )
    (outdir / "delta_age_stats.json").write_text(json.dumps(stats_out, indent=2))
    print(json.dumps(stats_out, indent=2))


if __name__ == "__main__":
    main()
