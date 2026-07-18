"""Programmatic retrieval of GEO series supplementary files.

No browser steps: GEO exposes series supplementary files over HTTPS at
https://ftp.ncbi.nlm.nih.gov/geo/series/GSEnnnnnn/GSEXXXXXX/suppl/
"""
from __future__ import annotations

import hashlib
import json
import re
import urllib.request
from pathlib import Path

GEO_SUPPL = "https://ftp.ncbi.nlm.nih.gov/geo/series/{bucket}/{gse}/suppl/"


def _bucket(gse: str) -> str:
    """GEO groups series FTP dirs as GSE123nnn for GSE123456."""
    num = gse.upper().replace("GSE", "")
    return f"GSE{num[:-3]}nnn"


def list_supplementary(gse: str) -> list[str]:
    """Return file names available in a GEO series' supplementary directory."""
    url = GEO_SUPPL.format(bucket=_bucket(gse), gse=gse.upper())
    html = urllib.request.urlopen(url, timeout=60).read().decode()
    return re.findall(r'href="([^"?/][^"]*)"', html)


def download_supplementary(
    gse: str,
    outdir: str | Path,
    patterns: list[str] | None = None,
    manifest_path: str | Path | None = None,
) -> list[Path]:
    """Download supplementary files whose names match any regex in `patterns`.

    Records SHA-256 checksums and source URLs into the data manifest.
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    names = list_supplementary(gse)
    if patterns:
        names = [n for n in names if any(re.search(p, n) for p in patterns)]
    if not names:
        raise FileNotFoundError(
            f"No supplementary files of {gse} matched {patterns}. "
            f"Available: {list_supplementary(gse)}"
        )
    downloaded: list[Path] = []
    base = GEO_SUPPL.format(bucket=_bucket(gse), gse=gse.upper())
    for name in names:
        dest = outdir / name
        if not dest.exists():
            urllib.request.urlretrieve(base + name, dest)
        downloaded.append(dest)
        if manifest_path:
            record_checksum(manifest_path, dest, source_url=base + name, study=gse)
    return downloaded


def sha256(path: str | Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def record_checksum(manifest_path: str | Path, path: str | Path, **meta) -> None:
    """Append a file's checksum + provenance metadata to the JSON manifest."""
    manifest_path = Path(manifest_path)
    manifest: dict = {}
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text())
    manifest[str(path)] = {"sha256": sha256(path), **meta}
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True))
