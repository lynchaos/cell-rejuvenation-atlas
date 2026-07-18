"""Statistical helpers shared across modules."""
from __future__ import annotations

import numpy as np
from scipy import stats


def benjamini_hochberg(pvals: np.ndarray) -> np.ndarray:
    """FDR-adjusted p-values (Benjamini-Hochberg)."""
    p = np.asarray(pvals, dtype=float)
    n = p.size
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    adj = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(adj, 0, 1)
    return out


def paired_bootstrap_delta(
    treated: np.ndarray,
    control: np.ndarray,
    n_boot: int = 10_000,
    seed: int = 0,
) -> tuple[float, float, float]:
    """Bootstrap CI for the mean paired difference (treated - control).

    Returns (point estimate, lower 95% CI, upper 95% CI).
    """
    treated = np.asarray(treated, float)
    control = np.asarray(control, float)
    if treated.shape != control.shape:
        raise ValueError("Paired arrays must have identical shapes.")
    d = treated - control
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, d.size, size=(n_boot, d.size))
    boots = d[idx].mean(axis=1)
    return float(d.mean()), float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))


def effect_size_cohens_dz(treated: np.ndarray, control: np.ndarray) -> float:
    """Standardized effect size for paired designs (dz)."""
    d = np.asarray(treated, float) - np.asarray(control, float)
    sd = d.std(ddof=1)
    return float(d.mean() / sd) if sd > 0 else np.nan


def two_group_de(
    x_treated: np.ndarray,
    x_control: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Welch t-test per feature (columns). Returns (statistic, raw p-values)."""
    t, p = stats.ttest_ind(x_treated, x_control, axis=0, equal_var=False, nan_policy="omit")
    return np.asarray(t), np.asarray(p)
