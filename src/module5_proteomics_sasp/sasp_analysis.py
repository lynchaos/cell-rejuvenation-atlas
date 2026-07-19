"""Reanalysis of the SASP Atlas proteomics (Basisty et al. 2020).

Core logic is pure pandas/scipy (unit-testable):
  * per-inducer differential abundance (senescent vs control), Welch t-test
  * core SASP = up in >= 2 inducers at FDR < 0.05 and log2FC > 0.5
  * enrichment against secretory/ECM gene sets via gseapy (optional)

The public processed data (PLoS S1 Table / MassIVE MSV000083750) is
group-level: one protein per row with the publication's own SEN/CTL
statistics per experiment sheet. main() uses those directly when the
workbook is present; the per-sample t-test path (--design) remains for
synthetic/test inputs.

Usage:
    python -m src.module5_proteomics_sasp.sasp_analysis \
        --input data/module5 --outdir results/module5
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
from src.common.stats import benjamini_hochberg, two_group_de

# Published exemplar SASP factors (Basisty et al. 2020; Coppé et al. 2010)
PUBLISHED_SASP = ["GDF15", "STC1", "MMP1", "MMP3", "SERPINE1", "IL6", "CXCL8", "TGFB1"]

# Sheets of the S1 Table workbook forming the paper's fibroblast soluble-SASP
# triad, against which the core SASP is defined
CORE_SASP_SHEETS = ["IR SASP", "RAS SASP", "ATV SASP"]


def load_protein_matrix(indir: Path) -> pd.DataFrame:
    """samples x proteins intensity matrix (first column = sample id)."""
    skip = ("listing", "design", "metadata")
    files = sorted(
        f
        for f in indir.rglob("*")
        if f.suffix in {".csv", ".tsv", ".txt"} and not any(s in f.name.lower() for s in skip)
    )
    if not files:
        raise FileNotFoundError(f"No protein matrix under {indir}")
    preferred = [f for f in files if any(k in f.name.lower() for k in ("protein", "matrix", "sasp"))]
    p = (preferred or files)[0]
    sep = "\t" if p.suffix in {".tsv", ".txt"} else ","
    df = pd.read_csv(p, sep=sep, index_col=0)
    if df.shape[0] > df.shape[1]:
        df = df.T
    return df.apply(pd.to_numeric, errors="coerce")


def load_sasp_workbook(indir: Path) -> dict[str, pd.DataFrame]:
    """Read the Basisty S1 Table workbook into per-experiment DE tables.

    The public processed data is group-level: one protein per row with the
    publication's own SEN/CTL statistics. Each sheet is normalized to the
    differential_abundance() output shape (index=gene; log2fc, p, padj) so
    core_sasp() applies unchanged.
    """
    files = sorted(indir.rglob("*.xlsx"))
    if not files:
        raise FileNotFoundError(f"No SASP workbook under {indir}")
    sheets = pd.read_excel(files[0], sheet_name=None, skiprows=2)
    per_experiment = {}
    for sheet, raw in sheets.items():
        cols = {str(c).lower().replace(" ", "").replace("-", ""): c for c in raw.columns}
        gene_col = next((c for k, c in cols.items() if k.startswith("gene")), None)
        lfc_col = next((c for k, c in cols.items() if k.startswith("avglog")), None)
        p_col = next((c for k, c in cols.items() if k in ("pvalue", "pval")), None)
        q_col = next((c for k, c in cols.items() if k in ("qvalue", "qval")), None)
        if gene_col is None or lfc_col is None or q_col is None:
            continue
        de = pd.DataFrame({
            "log2fc": pd.to_numeric(raw[lfc_col], errors="coerce"),
            "p": pd.to_numeric(raw[p_col], errors="coerce") if p_col else np.nan,
            "padj": pd.to_numeric(raw[q_col], errors="coerce"),
        })
        de.index = raw[gene_col].astype(str).str.split().str[0]
        de = de.dropna(subset=["padj"])
        per_experiment[str(sheet)] = de[~de.index.duplicated(keep="first")]
    return per_experiment


def differential_abundance(
    matrix: pd.DataFrame, design: pd.Series
) -> pd.DataFrame:
    """Per-protein Welch t-test; log2FC = treatment vs reference.

    The group whose label contains 'control' (case-insensitive) is the
    reference; otherwise the alphabetically first group is the reference.
    """
    groups = sorted(design.unique(), key=lambda g: "control" not in str(g).lower())
    if len(groups) != 2:
        raise ValueError("design must have exactly 2 groups")
    a = matrix.loc[design[design == groups[0]].index].dropna(axis=1, how="all")
    b = matrix.loc[design[design == groups[1]].index]
    t, p = two_group_de(b[a.columns].to_numpy(), a.to_numpy())
    out = pd.DataFrame(
        {
            "protein": a.columns,
            "log2fc": np.log2((b[a.columns].mean() + 1e-9) / (a.mean() + 1e-9)),
            "t": t,
            "p": p,
        }
    ).set_index("protein")
    out["padj"] = benjamini_hochberg(out["p"].fillna(1).to_numpy())
    return out.sort_values("padj")


def core_sasp(per_inducer: dict[str, pd.DataFrame], min_inducers: int = 2,
              padj: float = 0.05, lfc: float = 0.5) -> pd.Index:
    """Proteins significantly upregulated in >= min_inducers contrasts."""
    hits = []
    for name, de in per_inducer.items():
        up = de[(de["padj"] < padj) & (de["log2fc"] > lfc)].index
        hits.append(pd.Series(1, index=up, name=name))
    if not hits:
        return pd.Index([])
    counts = pd.concat(hits, axis=1).fillna(0).sum(axis=1)
    return counts[counts >= min_inducers].index


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", default="data/module5")
    ap.add_argument("--outdir", default="results/module5")
    ap.add_argument("--design", default=None,
                    help="CSV: sample,inducer,condition (condition in {senescent,control})")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    use_style()

    indir = Path(args.input)
    summary: dict = {}

    if list(indir.rglob("*.xlsx")):
        # Real-data path: published group-level SEN/CTL statistics per sheet.
        per_inducer = load_sasp_workbook(indir)
        summary["n_experiments"] = len(per_inducer)
        summary["proteins_per_experiment"] = {k: int(len(v)) for k, v in per_inducer.items()}
    elif args.design:
        # Per-sample path (synthetic/test inputs): Welch t-test per inducer.
        matrix = load_protein_matrix(indir)
        summary.update(n_samples=int(matrix.shape[0]), n_proteins=int(matrix.shape[1]))
        design = pd.read_csv(args.design, index_col=0)
        per_inducer = {}
        for inducer in design["inducer"].unique():
            sub = design[design["inducer"] == inducer]
            per_inducer[inducer] = differential_abundance(
                matrix.loc[sub.index], sub["condition"]
            )
    else:
        summary["note"] = "Provide the S1 Table workbook or --design to run differential analysis."
        (outdir / "summary.json").write_text(json.dumps(summary, indent=2))
        print(json.dumps(summary, indent=2))
        return

    for name, de in per_inducer.items():
        safe = re.sub(r"[^A-Za-z0-9]+", "_", name).strip("_")
        de.to_csv(outdir / f"de_{safe}.csv")
    core_keys = [s for s in CORE_SASP_SHEETS if s in per_inducer]
    core_pool = {k: per_inducer[k] for k in core_keys} if core_keys else per_inducer
    core = core_sasp(core_pool)
    pd.Series(core, name="protein").to_csv(outdir / "core_sasp.csv", index=False)
    summary["core_sasp_sheets"] = core_keys or sorted(per_inducer)
    summary["core_sasp_size"] = len(core)
    summary["published_sasp_recovered"] = sorted(set(PUBLISHED_SASP) & set(core))

    # UpSet-style bar: overlap of significant-up sets across inducers
    fig, ax = plt.subplots(figsize=(4.6, 3.2))
    sizes = {k: int(((v["padj"] < 0.05) & (v["log2fc"] > 0.5)).sum()) for k, v in per_inducer.items()}
    ax.bar(sizes.keys(), sizes.values(), color=PALETTE["aged"])
    ax.set_ylabel("proteins up in senescence (FDR<0.05)")
    ax.set_title("SASP Atlas (PXD013721) — per-inducer secretomes")
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    fig.savefig(outdir / "sasp_per_inducer.pdf")
    fig.savefig(outdir / "sasp_per_inducer.png")

    (outdir / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
