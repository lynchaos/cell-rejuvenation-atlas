"""Download Allen et al. (2023) MERFISH aging-brain data (GSE207848)."""
from __future__ import annotations

import argparse
from pathlib import Path

from src.common.geo import download_supplementary, list_supplementary

GSE = "GSE207848"
MANIFEST = "reports/data_manifest.json"
PATTERNS = [r"\.h5ad", r"\.h5", r"counts", r"metadata", r"coordinates", r"\.csv", r"\.tar"]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--outdir", default="data/module4")
    args = ap.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    available = list_supplementary(GSE)
    print(f"[{GSE}] available:\n  " + "\n  ".join(available))
    (outdir / f"{GSE}_supplementary_listing.txt").write_text("\n".join(available))
    files = download_supplementary(GSE, outdir, patterns=PATTERNS, manifest_path=MANIFEST)
    print(f"[done] downloaded: {[f.name for f in files]}")


if __name__ == "__main__":
    main()
