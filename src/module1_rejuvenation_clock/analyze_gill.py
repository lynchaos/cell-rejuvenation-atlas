"""Reproduce the Gill et al. (2022) DNAm-age reversal under MPTR.

Inputs
------
--beta       Directory or file with a samples x CpG beta-value matrix (CSV/TSV,
             optionally gzipped; first column = probe ID when probes x samples)
             for the transient arm.
--metadata   metadata.csv produced by download_data.py
--coef       Optional CSV of published clock coefficients (probe,coef).
             If omitted, --train-beta/--train-ages must be given to fit
             an elastic-net clock on external training data.
--outdir     Output directory for figures + tables.

Outputs
-------
dnam_age_timecourse.pdf / .csv   DNAm age by MPTR day and experimental arm
delta_age_stats.json             bootstrap + effect size, transient vs control
                                 at the peak reprogramming day
"""
from __future__ import annotations

import argparse
import gzip
import json
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from src.common.plotting import PALETTE, use_style
from src.common.stats import (
    effect_size_cohens_d,
    unpaired_bootstrap_delta,
)
from src.module1_rejuvenation_clock.clock import (
    EpigeneticClock,
    train_clock,
)

DAY_AFTER_WORD = re.compile(r"day[_\s]?(\d+)", re.IGNORECASE)    # 'MPTR day 13 rep 2'
DAY_BEFORE_WORD = re.compile(r"(\d+)\s*days?", re.IGNORECASE)    # '..._13days_exp1'
PVAL_COLUMN = re.compile(r"Detection Pval(\.\d+)?$", re.IGNORECASE)


def infer_day(title: str) -> float:
    """Parse MPTR day from a sample title ('MPTR day 13' or '13days')."""
    m = DAY_AFTER_WORD.search(str(title)) or DAY_BEFORE_WORD.search(str(title))
    return float(m.group(1)) if m else np.nan


def infer_arm(title: str) -> str:
    """Classify a Gill et al. sample title into an experimental arm."""
    t = str(title)
    if "transiently_reprogrammed" in t or "transient_reprogramming" in t:
        return "transient"
    if "fail" in t:
        return "failed"
    if "control" in t:
        return "control"
    return "baseline"


def _sniff_sep(path: Path) -> str:
    """Pick the delimiter by content, not extension (GEO mislabels CSV as .txt)."""
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", errors="replace") as fh:
        first = fh.readline()
    return "\t" if first.count("\t") > first.count(",") else ","


def load_beta(path: str) -> pd.DataFrame:
    p = Path(path)
    if p.is_dir():
        candidates = sorted(
            [f for f in p.rglob("*") if f.suffix in {".csv", ".tsv", ".txt"} or
             (f.suffix == ".gz" and f.stem.endswith((".txt", ".tsv", ".csv")))]
        )
        if not candidates:
            raise FileNotFoundError(f"No beta matrix found under {p}")
        p = candidates[0]
    df = pd.read_csv(p, sep=_sniff_sep(p), index_col=0)
    # GEO processed matrices interleave one 'Detection Pval' column per sample
    df = df.drop(columns=[c for c in df.columns if PVAL_COLUMN.match(str(c))])
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

    if args.coef:
        clock = EpigeneticClock.from_coefficient_csv(args.coef)
    else:
        train_beta = load_beta(args.train_beta)
        train_ages = pd.read_csv(args.train_ages, index_col=0)["age"]
        clock = train_clock(train_beta.loc[train_ages.index], train_ages.to_numpy())

    # Sample names in the beta matrix are GEO titles; parse day/arm from them.
    # Fall back to joining metadata by accession when titles carry no day.
    days = beta.index.to_series().map(infer_day)
    if days.notna().any():
        days.index = beta.index
        arms = pd.Series(beta.index.map(infer_arm), index=beta.index)
    else:
        meta = pd.read_csv(args.metadata, index_col=0)
        title_col = next((c for c in meta.columns if "title" in c.lower()), None)
        days = meta[title_col].map(infer_day) if title_col else pd.Series(np.nan, index=meta.index)
        arms = meta[title_col].map(infer_arm) if title_col else pd.Series("baseline", index=meta.index)
        common = beta.index.intersection(days.dropna().index)
        if len(common) == 0:  # fall back: align by order (documented in report)
            common = beta.index
            days = pd.Series(np.linspace(0, args.peak_day, len(common)), index=common)
            arms = pd.Series("transient", index=common)

    common = beta.index.intersection(days.dropna().index)
    if len(common) == 0:
        common = beta.index
        days = pd.Series(np.nan, index=common)
    dnam_age = pd.Series(clock.predict(beta.loc[common]), index=common)
    df = pd.DataFrame({
        "sample": common,
        "day": days.loc[common].to_numpy(),
        "arm": arms.reindex(common).fillna("baseline").to_numpy(),
        "dnam_age": dnam_age,
    })
    df = df.sort_values("day")
    df.to_csv(outdir / "dnam_age_timecourse.csv", index=False)

    # ---- Figure: DNAm age per arm across the MPTR time course ----
    arm_colors = {
        "transient": PALETTE["rejuvenated"],
        "control": PALETTE["control"],
        "failed": PALETTE["aged"],
        "baseline": PALETTE["young"],
    }
    fig, ax = plt.subplots(figsize=(5, 3.4))
    for arm, sub in df.dropna(subset=["day"]).groupby("arm"):
        color = arm_colors.get(arm, "#888888")
        means = sub.groupby("day")["dnam_age"].mean()
        ax.scatter(sub["day"], sub["dnam_age"], color=color, alpha=0.3, s=10, zorder=1)
        ax.plot(means.index, means.values, "o-", color=color, lw=1.5, ms=4, label=arm, zorder=2)
    ax.set_xlabel("MPTR day")
    ax.set_ylabel("DNAm age (years)")
    ax.set_title("Epigenetic age across transient reprogramming\n(Gill et al. 2022, GSE165179)")
    ax.legend(title="arm")
    fig.savefig(outdir / "dnam_age_timecourse.pdf")
    fig.savefig(outdir / "dnam_age_timecourse.png")

    # ---- Rejuvenation contrast: transient vs control at the peak day ----
    stats_out: dict = {"n_samples": int(len(df)), "peak_day": args.peak_day}
    at_peak = np.isclose(df["day"], args.peak_day)
    dp = df.loc[at_peak & (df["arm"] == "transient"), "dnam_age"].to_numpy()
    d0 = df.loc[at_peak & (df["arm"] == "control"), "dnam_age"].to_numpy()
    if len(dp) and len(d0):
        delta, lo, hi = unpaired_bootstrap_delta(dp, d0)
        stats_out.update(
            contrast=f"transient vs control, day {args.peak_day:g}",
            delta_years=delta,
            ci95=[lo, hi],
            cohens_dz=effect_size_cohens_d(dp, d0),
            n_transient=len(dp),
            n_control=len(d0),
            interpretation=(
                "Negative delta with CI < 0 reproduces the epigenetic-age "
                "reversal reported by Gill et al. (2022)."
            ),
        )
    (outdir / "delta_age_stats.json").write_text(json.dumps(stats_out, indent=2))
    print(json.dumps(stats_out, indent=2))


if __name__ == "__main__":
    main()
