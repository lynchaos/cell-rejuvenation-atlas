"""Download Gill et al. (2022) processed data from GEO.

Usage:
    python -m src.module1_rejuvenation_clock.download_data --outdir data/module1
"""
from __future__ import annotations

import argparse
import gzip
import shutil
import tarfile
import tempfile
from pathlib import Path

import pandas as pd

from src.common.geo import download_supplementary, list_supplementary, record_checksum

METHYLATION_SERIES = {"transient": "GSE165179", "sendai": "GSE165178"}
RNASEQ_SERIES = {"transient": "GSE165177", "sendai": "GSE165176"}
MANIFEST = "reports/data_manifest.json"


def fetch_series_matrix(gse: str, outdir: Path) -> Path:
    """Fetch the GEO series-matrix tarball (always contains sample metadata)."""
    import urllib.request

    num = gse.replace("GSE", "")
    url = (
        f"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE{num[:-3]}nnn/"
        f"{gse}/matrix/{gse}_series_matrix.txt.gz"
    )
    dest = outdir / f"{gse}_series_matrix.txt.gz"
    if not dest.exists():
        urllib.request.urlretrieve(url, dest)
    record_checksum(MANIFEST, dest, source_url=url, study=gse)
    return dest


def parse_series_matrix_metadata(path: Path) -> pd.DataFrame:
    """Extract per-sample metadata rows (!Sample_*) from a series matrix file."""
    rows: dict[str, list[str]] = {}
    with gzip.open(path, "rt", errors="replace") as fh:
        for line in fh:
            if not line.startswith("!Sample_"):
                if line.startswith('"ID_REF"'):
                    break
                continue
            key, *vals = line.rstrip("\n").split("\t")
            rows[key] = [v.strip('"') for v in vals]
    df = pd.DataFrame(rows)
    if "!Sample_geo_accession" in df.columns:
        df = df.set_index("!Sample_geo_accession")
    return df


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--outdir", default="data/module1")
    args = ap.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    metas = []
    for arm, gse in list(METHYLATION_SERIES.items()) + list(RNASEQ_SERIES.items()):
        print(f"[download] {gse} ({arm})")
        matrix = fetch_series_matrix(gse, outdir)
        meta = parse_series_matrix_metadata(matrix)
        meta["arm"] = arm
        meta["series"] = gse
        metas.append(meta)
        try:
            files = download_supplementary(
                gse, outdir / gse, patterns=[r"\.txt", r"\.csv", r"\.tsv", r"\.tar"],
                manifest_path=MANIFEST,
            )
            print(f"  supplementary: {[f.name for f in files]}")
        except FileNotFoundError as e:
            print(f"  warning: {e}")
        available = list_supplementary(gse)
        (outdir / f"{gse}_supplementary_listing.txt").write_text("\n".join(available))

    pd.concat(metas).to_csv(outdir / "metadata.csv")
    print(f"[done] metadata -> {outdir / 'metadata.csv'}")


if __name__ == "__main__":
    main()
