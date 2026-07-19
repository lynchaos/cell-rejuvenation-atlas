"""Project Browder 2022 bulk RNA-seq onto a Tabula Muris Senis aging axis.

GSE190983 (Browder in-vivo OSK, 7-month cohort) is BULK RNA-seq — a genes x
samples count matrix over six tissues — so single-cell integration against
the Tabula Muris Senis (TMS) atlas (e.g. scVI) is not valid. Instead we
derive a transcriptomic aging axis from TMS FACS (mean log-CPM old − young
per gene) and score each Browder sample by projection onto that axis. If
long-term cyclic OSK reverts age-associated expression, 4F samples should
score BELOW their tissue-matched controls.

Usage:
    python -m src.module3_multiomics_integration.aging_signature_score \
        --tms data/module3/tabula_muris_senis/tms.h5ad \
        --browder data/module3/browder_long7m/GSE190983_count.tsv.gz \
        --outdir results/module3
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from src.common.stats import (
    benjamini_hochberg,
    effect_size_cohens_d,
    two_group_de,
    unpaired_bootstrap_delta,
)
from src.module3_multiomics_integration.rejuvenation_score import project_onto_signature

MYGENE_URL = "https://mygene.info/v3/query"
YOUNG_AGES = ("3m",)
OLD_AGES = ("18m", "21m", "24m")


def tms_aging_signature(tms_path: str | Path) -> pd.Series:
    """Old − young mean log-CPM per gene over the whole TMS FACS atlas."""
    import scanpy as sc  # heavy import: keep out of module scope for tests

    adata = sc.read_h5ad(tms_path)
    sc.pp.normalize_total(adata, target_sum=1e6)
    sc.pp.log1p(adata)
    age = adata.obs["age"].astype(str)
    x = adata.X

    def col_mean(mask: np.ndarray) -> np.ndarray:
        return np.asarray(x[mask].mean(axis=0)).ravel()

    old = col_mean(age.isin(OLD_AGES).to_numpy())
    young = col_mean(age.isin(YOUNG_AGES).to_numpy())
    return pd.Series(old - young, index=adata.var_names.astype(str), name="aging_signature")


def parse_browder_sample(name: str) -> tuple[str, str]:
    """'Skeletal_Muscle_4F_m3' -> ('Skeletal_Muscle', '4F')."""
    parts = str(name).split("_")
    for i, tok in enumerate(parts):
        if tok.lower() in ("4f", "control"):
            return "_".join(parts[:i]), "4F" if tok.lower() == "4f" else "Control"
    raise ValueError(f"Cannot parse Browder sample name: {name!r}")


def ensembl_to_symbol(ensembl_ids: list[str], cache_path: str | Path) -> dict[str, str]:
    """ENSMUSG -> gene symbol via mygene.info (batch POST), cached as JSON."""
    cache_path = Path(cache_path)
    if cache_path.exists():
        return json.loads(cache_path.read_text())
    import requests  # heavy import: keep out of module scope for tests

    ids = list(dict.fromkeys(ensembl_ids))
    mapping: dict[str, str] = {}
    for i in range(0, len(ids), 1000):
        batch = ids[i : i + 1000]
        r = requests.post(
            MYGENE_URL,
            data={
                "q": ",".join(batch),
                "scopes": "ensembl.gene",
                "fields": "symbol",
                "species": "mouse",
            },
            timeout=120,
        )
        r.raise_for_status()
        for hit in r.json():
            if hit.get("notfound"):
                continue
            symbol = hit.get("symbol")
            if symbol:
                mapping[hit["query"]] = symbol
    cache_path.write_text(json.dumps(mapping, indent=0, sort_keys=True))
    return mapping


def collapse_by_symbol(counts: pd.DataFrame, id_to_symbol: dict[str, str]) -> pd.DataFrame:
    """Map Ensembl-indexed counts (genes x samples) to symbols, summing duplicates."""
    genes = counts.index.astype(str).str.replace(r"\.\d+$", "", regex=True)  # strip version
    symbols = pd.Series(genes, index=counts.index).map(id_to_symbol)
    mapped = counts.loc[symbols.notna()]
    mapped.index = symbols[symbols.notna()].to_numpy()
    return mapped.groupby(level=0).sum()


def normalize_log_cpm(counts: pd.DataFrame) -> pd.DataFrame:
    """log1p(CPM) of a genes x samples count matrix."""
    lib = counts.sum(axis=0)
    return np.log1p(counts / lib * 1e6)


def tissue_stats(scores: pd.Series, meta: pd.DataFrame) -> pd.DataFrame:
    """Per tissue: 4F vs Control on projected aging scores (Welch t + effect size)."""
    rows = []
    for tissue, sub in meta.groupby("tissue"):
        s_4f = scores.reindex(sub.index[sub["condition"] == "4F"]).dropna()
        s_ct = scores.reindex(sub.index[sub["condition"] == "Control"]).dropna()
        row = {
            "tissue": tissue,
            "n_control": int(s_ct.size),
            "n_4f": int(s_4f.size),
            "mean_control": float(s_ct.mean()),
            "mean_4f": float(s_4f.mean()),
        }
        if s_4f.size >= 2 and s_ct.size >= 2:
            t, p = two_group_de(s_4f.to_numpy(), s_ct.to_numpy())
            delta, lo, hi = unpaired_bootstrap_delta(s_4f.to_numpy(), s_ct.to_numpy())
            row.update(
                delta=delta,
                ci95_lo=lo,
                ci95_hi=hi,
                pvalue=float(np.atleast_1d(p)[0]),
                statistic=float(np.atleast_1d(t)[0]),
                cohens_d=effect_size_cohens_d(s_4f.to_numpy(), s_ct.to_numpy()),
            )
        else:
            row.update(delta=np.nan, ci95_lo=np.nan, ci95_hi=np.nan,
                       pvalue=np.nan, statistic=np.nan, cohens_d=np.nan)
        rows.append(row)
    out = pd.DataFrame(rows)
    out["padj"] = benjamini_hochberg(out["pvalue"].fillna(1.0).to_numpy())
    return out.sort_values("delta").reset_index(drop=True)


def plot_scores(meta: pd.DataFrame, stats: pd.DataFrame, out_png: Path) -> None:
    """Strip + mean plot of aging scores per tissue, 4F vs Control."""
    import matplotlib.pyplot as plt

    from src.common.plotting import use_style

    use_style()
    order = stats["tissue"].tolist()
    fig, ax = plt.subplots(figsize=(8.5, 4.5))
    rng = np.random.default_rng(0)
    styles = {"Control": ("C7", "o", "control"), "4F": ("C3", "s", "4F (OSK 7m)")}
    for pos, tissue in enumerate(order):
        sub = meta[meta["tissue"] == tissue]
        for cond, (color, marker, label) in styles.items():
            vals = sub.loc[sub["condition"] == cond, "aging_score"].to_numpy()
            xoff = -0.15 if cond == "Control" else 0.15
            ax.scatter(pos + xoff + rng.uniform(-0.05, 0.05, vals.size), vals,
                       c=color, marker=marker, s=28, zorder=3,
                       label=label if pos == 0 else None)
            if vals.size:
                ax.hlines(vals.mean(), pos + xoff - 0.12, pos + xoff + 0.12,
                          colors="black", linewidth=1.4, zorder=4)
        p = stats.loc[stats["tissue"] == tissue, "padj"].iloc[0]
        if pd.notna(p) and p < 0.05:
            ax.text(pos, ax.get_ylim()[1], "*", ha="center", va="bottom", fontsize=14)
    ax.set_xticks(range(len(order)), order, rotation=30, ha="right")
    ax.set_ylabel("TMS aging-axis score (log-CPM projection)")
    ax.set_title("Browder 2022 (GSE190983): 4F vs control projected on the TMS aging axis")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    fig.savefig(out_png.with_suffix(".pdf"))
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tms", required=True, help="TMS FACS h5ad (annotated, raw counts)")
    ap.add_argument("--browder", required=True, help="GSE190983_count.tsv.gz (genes x samples)")
    ap.add_argument("--outdir", default="results/module3")
    args = ap.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    signature = tms_aging_signature(args.tms)
    signature.to_csv(outdir / "tms_aging_signature.csv")

    counts = pd.read_csv(args.browder, sep="\t", index_col=0)
    mapping = ensembl_to_symbol(list(counts.index.astype(str)),
                                outdir / "ensembl_to_symbol.json")
    symbols = collapse_by_symbol(counts, mapping)
    expr = normalize_log_cpm(symbols).T  # samples x genes
    expr.index = expr.index.astype(str)

    scores = project_onto_signature(expr, signature)
    meta = pd.DataFrame([parse_browder_sample(s) for s in scores.index],
                        columns=["tissue", "condition"], index=scores.index)
    meta["aging_score"] = scores
    meta.reset_index(names="sample").to_csv(outdir / "sample_scores.csv", index=False)

    stats = tissue_stats(scores, meta)
    stats.to_csv(outdir / "tissue_scores.csv", index=False)
    plot_scores(meta, stats, outdir / "aging_score.png")

    summary = {
        "n_samples": int(meta.shape[0]),
        "n_mapped_genes": int(symbols.shape[0]),
        "n_signature_genes_shared": int((signature.reindex(expr.columns).notna()).sum()),
        "tissues_4f_significantly_younger": stats.loc[
            (stats["padj"] < 0.05) & (stats["delta"] < 0), "tissue"
        ].tolist(),
        "interpretation": "delta = mean(4F) - mean(control) on the TMS aging axis; "
                          "negative = 4F shifted toward the young state",
    }
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2))
    print(stats.to_string(index=False))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
