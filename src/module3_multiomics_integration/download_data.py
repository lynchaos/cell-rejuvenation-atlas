"""Download Browder 2022 + Tabula Muris Senis processed matrices.

Usage:
    python -m src.module3_multiomics_integration.download_data --outdir data/module3
"""
from __future__ import annotations

import argparse
from pathlib import Path

from src.common.geo import download_supplementary, list_supplementary

SERIES = {
    "browder_long7m": "GSE190983",
    "browder_short1m": "GSE190984",
    "browder_long10m": "GSE190985",
    "tabula_muris_senis": "GSE149590",
}
MANIFEST = "reports/data_manifest.json"
PATTERNS = [r"\.h5", r"\.h5ad", r"\.rds", r"\.mtx", r"counts", r"metadata", r"\.csv", r"\.tar"]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--outdir", default="data/module3")
    args = ap.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    for name, gse in SERIES.items():
        dest = outdir / name
        dest.mkdir(exist_ok=True)
        available = list_supplementary(gse)
        print(f"[{gse}] available:\n  " + "\n  ".join(available))
        (dest / f"{gse}_supplementary_listing.txt").write_text("\n".join(available))
        try:
            files = download_supplementary(gse, dest, patterns=PATTERNS, manifest_path=MANIFEST)
            print(f"  downloaded: {[f.name for f in files]}")
        except FileNotFoundError as e:
            print(f"  warning: {e}")


if __name__ == "__main__":
    main()
