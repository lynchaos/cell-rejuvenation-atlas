"""Download Schiebinger et al. (2019) processed reprogramming data (GSE122662).

Usage:
    python -m src.module2_reprogramming_trajectory.download_data --outdir data/module2
"""
from __future__ import annotations

import argparse
import tarfile
from pathlib import Path

from src.common.geo import download_supplementary, list_supplementary

GSE = "GSE122662"
MANIFEST = "reports/data_manifest.json"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--outdir", default="data/module2")
    ap.add_argument(
        "--patterns",
        nargs="*",
        default=[r"metadata", r"counts", r"matrix", r"\.h5", r"\.mtx", r"\.csv", r"\.txt", r"\.tar"],
        help="Regexes selecting supplementary files (default: processed matrices).",
    )
    args = ap.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    available = list_supplementary(GSE)
    print(f"[{GSE}] available supplementary files:\n  " + "\n  ".join(available))
    (outdir / f"{GSE}_supplementary_listing.txt").write_text("\n".join(available))

    files = download_supplementary(GSE, outdir, patterns=args.patterns, manifest_path=MANIFEST)
    print(f"[done] downloaded: {[f.name for f in files]}")

    # GEO RAW bundles per-sample 10x triplets in a single tar; unpack in place
    for f in files:
        if f.suffix == ".tar":
            with tarfile.open(f) as tf:
                try:
                    tf.extractall(outdir, filter="data")
                except TypeError:  # python < 3.11.4 lacks the filter argument
                    tf.extractall(outdir)
            print(f"  extracted {f.name}")


if __name__ == "__main__":
    main()
