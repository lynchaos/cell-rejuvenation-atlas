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


def unpaired_bootstrap_delta(
    treated: np.ndarray,
    control: np.ndarray,
    n_boot: int = 10_000,
    seed: int = 0,
) -> tuple[float, float, float]:
    """Bootstrap CI for the difference in means (treated - control), unpaired.

    Returns (point estimate, lower 95% CI, upper 95% CI).
    """
    treated = np.asarray(treated, float)
    control = np.asarray(control, float)
    rng = np.random.default_rng(seed)
    ti = rng.integers(0, treated.size, size=(n_boot, treated.size))
    ci = rng.integers(0, control.size, size=(n_boot, control.size))
    boots = treated[ti].mean(axis=1) - control[ci].mean(axis=1)
    delta = float(treated.mean() - control.mean())
    return delta, float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))


def effect_size_cohens_d(treated: np.ndarray, control: np.ndarray) -> float:
    """Standardized effect size for two independent groups (Cohen's d, pooled SD)."""
    a, b = np.asarray(treated, float), np.asarray(control, float)
    na, nb = a.size, b.size
    if na < 2 or nb < 2:
        return np.nan
    pooled = ((na - 1) * a.var(ddof=1) + (nb - 1) * b.var(ddof=1)) / (na + nb - 2)
    return float((a.mean() - b.mean()) / np.sqrt(pooled)) if pooled > 0 else np.nan


def two_group_de(
    x_treated: np.ndarray,
    x_control: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Welch t-test per feature (columns). Returns (statistic, raw p-values)."""
    t, p = stats.ttest_ind(x_treated, x_control, axis=0, equal_var=False, nan_policy="omit")
    return np.asarray(t), np.asarray(p)
