"""Download Browder 2022 bulk RNA-seq + Tabula Muris Senis reference.

GSE190983 (Browder in vivo OSK, 7-month) is BULK RNA-seq: genes x samples
(tissue_condition_mouseID), not single-cell. The reference atlas is the
annotated TMS FACS object from GSE149590 (sample-level file; the series
RAW tar bundles the same object).

Usage:
    python -m src.module3_multiomics_integration.download_data --outdir data/module3
"""
from __future__ import annotations

import argparse
import gzip
import shutil
import urllib.request
from pathlib import Path

from src.common.geo import record_checksum

TMS_H5AD_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4505nnn/GSM4505405/"
    "suppl/GSM4505405_tabula-muris-senis-facs-official-raw-obj.h5ad.gz"
)
BROWDER_COUNTS_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE190nnn/GSE190983/"
    "suppl/GSE190983_count.tsv.gz"
)
MANIFEST = "reports/data_manifest.json"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--outdir", default="data/module3")
    args = ap.parse_args()
    outdir = Path(args.outdir)
    (outdir / "tabula_muris_senis").mkdir(parents=True, exist_ok=True)
    (outdir / "browder_long7m").mkdir(parents=True, exist_ok=True)

    tms_gz = outdir / "tabula_muris_senis" / "tms.h5ad.gz"
    tms_out = tms_gz.with_suffix("")  # strip .gz -> tms.h5ad
    if not tms_out.exists():
        if not tms_gz.exists():
            urllib.request.urlretrieve(TMS_H5AD_URL, tms_gz)
        with gzip.open(tms_gz, "rb") as src, open(tms_out, "wb") as dst:
            shutil.copyfileobj(src, dst)
    record_checksum(MANIFEST, tms_out, source_url=TMS_H5AD_URL, study="GSE149590")

    browder = outdir / "browder_long7m" / "GSE190983_count.tsv.gz"
    if not browder.exists():
        urllib.request.urlretrieve(BROWDER_COUNTS_URL, browder)
    record_checksum(MANIFEST, browder, source_url=BROWDER_COUNTS_URL, study="GSE190983")
    print(f"[done] {tms_out}\n[done] {browder}")


if __name__ == "__main__":
    main()
