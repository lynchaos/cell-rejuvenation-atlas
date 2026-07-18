"""Download the SASP Atlas processed proteomics (PXD013721).

The processed DIA protein matrix is distributed with Basisty et al. 2020
(PLoS Biology) via PRIDE/ProteomeXchange. PRIDE exposes files over HTTPS:
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/01/PXD013721/
"""
from __future__ import annotations

import argparse
import re
import urllib.request
from pathlib import Path

from src.common.geo import record_checksum

PRIDE_URL = "https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/01/PXD013721/"
MANIFEST = "reports/data_manifest.json"


def list_pride_files(url: str = PRIDE_URL) -> list[str]:
    html = urllib.request.urlopen(url, timeout=60).read().decode()
    return re.findall(r'href="([^"?/][^"]*)"', html)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--outdir", default="data/module5")
    ap.add_argument("--patterns", nargs="*", default=[r"\.csv", r"\.tsv", r"\.txt", r"\.xlsx"])
    args = ap.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    files = list_pride_files()
    print("[PXD013721] available:\n  " + "\n  ".join(files))
    (outdir / "PXD013721_listing.txt").write_text("\n".join(files))

    selected = [f for f in files if any(re.search(p, f) for p in args.patterns)]
    if not selected:
        raise FileNotFoundError(f"No PRIDE files matched {args.patterns}")
    for name in selected:
        dest = outdir / name
        if not dest.exists():
            urllib.request.urlretrieve(PRIDE_URL + name, dest)
        record_checksum(MANIFEST, dest, source_url=PRIDE_URL + name, study="PXD013721")
        print(f"  downloaded {name}")


if __name__ == "__main__":
    main()
