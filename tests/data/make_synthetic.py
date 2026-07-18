"""Generate the tiny synthetic dataset used by CI (`data/test/`).

Mimics each module's *input shapes* so pipeline processes can smoke-test
end-to-end in minutes without downloading real data.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

OUT = Path(__file__).resolve().parents[2] / "data" / "test"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(42)

    # Module 1: samples x CpG beta matrix + metadata + ages
    n, p = 24, 2000
    beta = pd.DataFrame(rng.uniform(0, 1, (n, p)), columns=[f"cg{i:06d}" for i in range(p)])
    beta.index = [f"GSM{i:07d}" for i in range(n)]
    beta.to_csv(OUT / "module1_beta.csv")
    pd.DataFrame(
        {"!Sample_title": [f"MPTR day {d} rep {r}" for d in [0, 7, 10, 13] * 6 for r in [1]],
         "arm": "transient"},
        index=beta.index,
    ).to_csv(OUT / "module1_metadata.csv")
    # chronological ages for coef-less clock training in CI
    ages = pd.DataFrame({"age": np.linspace(30, 70, n)}, index=beta.index)
    ages.to_csv(OUT / "module1_train_ages.csv")

    # Module 5: samples x proteins + design
    prot = pd.DataFrame(rng.normal(100, 15, (16, 300)), columns=[f"PROT{i}" for i in range(300)])
    prot.index = [f"S{i}" for i in range(16)]
    senescent_rows = [4, 5, 6, 7, 12, 13, 14, 15]  # IR + RAS senescent
    prot.iloc[senescent_rows, :10] += 150  # planted SASP proteins in both inducers
    prot.to_csv(OUT / "module5_proteins.csv")
    pd.DataFrame(
        {"inducer": ["IR"] * 8 + ["RAS"] * 8,
         "condition": ["control"] * 4 + ["senescent"] * 4 + ["control"] * 4 + ["senescent"] * 4},
        index=prot.index,
    ).to_csv(OUT / "module5_design.csv")

    print(f"synthetic test data -> {OUT}")


if __name__ == "__main__":
    main()
