"""Download the SASP Atlas processed proteomics (Basisty et al. 2020).

The study's ProteomeXchange submission (PXD013721) was removed from PRIDE;
the processed quantification survives as the paper's S1 Table workbook
(group-level per-protein SEN/CTL statistics, one sheet per experiment),
served over HTTPS by PLoS and mirrored on MassIVE (MSV000083750).

Usage:
    python -m src.module5_proteomics_sasp.download_data --outdir data/module5
"""
from __future__ import annotations

import argparse
import urllib.request
from pathlib import Path

from src.common.geo import record_checksum

# PLoS Biology S1 Table (DOI-stable); MassIVE quant/ mirror of the same workbook
S1_TABLE_URLS = [
    "https://journals.plos.org/plosbiology/article/file"
    "?id=10.1371/journal.pbio.3000599.s007&type=supplementary",
    "https://massive.ucsd.edu/ProteoSAFe/DownloadResultFile"
    "?file=f.MSV000083750/quant/Supplemental_Table_1.xlsx&forceDownload=true",
]
S1_TABLE_NAME = "pbio.3000599.s007.xlsx"
MANIFEST = "reports/data_manifest.json"


def download_first_available(urls: list[str], dest: Path) -> str:
    """Try each URL in order; return the one that succeeded."""
    errors = []
    for url in urls:
        try:
            urllib.request.urlretrieve(url, dest)
            return url
        except Exception as e:  # noqa: BLE001 - report all mirrors then fail
            errors.append(f"{url}: {e}")
    raise ConnectionError("all mirrors failed:\n  " + "\n  ".join(errors))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--outdir", default="data/module5")
    args = ap.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    dest = outdir / S1_TABLE_NAME
    if dest.exists():
        url = S1_TABLE_URLS[0]
        print(f"[skip] {dest} already downloaded")
    else:
        url = download_first_available(S1_TABLE_URLS, dest)
        print(f"[done] {S1_TABLE_NAME} <- {url}")
    record_checksum(MANIFEST, dest, source_url=url, study="PXD013721")


if __name__ == "__main__":
    main()
