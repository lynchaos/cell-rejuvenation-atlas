"""Download the annotated Allen et al. (2023) MERFISH aging-brain object.

GEO GSE207848 ships only raw per-run archives (no unified cell typing), so we
fetch the curated, cell-type-annotated object the authors published on
CELLxGENE (collection 31937775-0602-4e52-a799-b6acdd2bac2e, dataset
"BrainAgingSpatialAtlas_MERFISH"): 378,918 cells x 374 genes with spatial
coordinates, donor, slice and age annotations.

Usage:
    python -m src.module4_spatial_aging.download_data --outdir data/module4
"""
from __future__ import annotations

import argparse
import urllib.request
from pathlib import Path

from src.common.geo import record_checksum

CELLXGENE_H5AD_URL = (
    "https://datasets.cellxgene.cziscience.com/"
    "c93d78c2-ee17-4504-8d1c-17cf093ad7b5.h5ad"
)
MANIFEST = "reports/data_manifest.json"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--outdir", default="data/module4")
    args = ap.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    dest = outdir / "allen_merfish_brain.h5ad"
    if not dest.exists():
        urllib.request.urlretrieve(CELLXGENE_H5AD_URL, dest)
    record_checksum(MANIFEST, dest, source_url=CELLXGENE_H5AD_URL, study="GSE207848")
    print(f"[done] {dest}")


if __name__ == "__main__":
    main()
