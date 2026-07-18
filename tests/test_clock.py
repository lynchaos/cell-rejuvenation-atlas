"""Unit tests for the epigenetic-clock machinery (synthetic data)."""
import numpy as np
import pandas as pd
import pytest

from src.common.stats import benjamini_hochberg, paired_bootstrap_delta
from src.module1_rejuvenation_clock.clock import (
    EpigeneticClock,
    horvath_transform,
    inverse_horvath_transform,
    train_clock,
)


def test_horvath_transform_roundtrip():
    ages = np.array([0.5, 5, 20, 40, 80], dtype=float)
    assert np.allclose(inverse_horvath_transform(horvath_transform(ages)), ages)


def _synthetic_methylation(n=200, n_probes=500, seed=1):
    rng = np.random.default_rng(seed)
    ages = rng.uniform(20, 80, n)
    beta = rng.uniform(0.05, 0.95, size=(n, n_probes))
    # plant 40 age-correlated probes
    signal = np.zeros((n, 40))
    for j in range(40):
        coef = rng.choice([-1, 1]) * rng.uniform(0.002, 0.005)
        signal[:, j] = 0.5 + coef * (ages - 50) + rng.normal(0, 0.02, n)
    beta[:, :40] = np.clip(signal, 0.01, 0.99)
    cols = [f"cg{i:06d}" for i in range(n_probes)]
    return pd.DataFrame(beta, columns=cols), ages


def test_trained_clock_recovers_age():
    beta, ages = _synthetic_methylation()
    clock = train_clock(beta, ages, cv_folds=3)
    pred = clock.predict(beta)
    corr = np.corrcoef(pred, ages)[0, 1]
    assert corr > 0.9, f"trained clock should track age (r={corr:.2f})"


def test_clock_prediction_handles_missing_probes():
    beta, ages = _synthetic_methylation()
    clock = train_clock(beta, ages, cv_folds=3)
    subset = beta[beta.columns[:100]]
    pred = clock.predict(subset)  # should not raise; rescales weights
    assert len(pred) == len(beta)


def test_clock_rejects_empty_overlap():
    clock = EpigeneticClock(probes=["cgX"], coef=np.array([1.0]), intercept=0.0)
    with pytest.raises(ValueError):
        clock.predict(pd.DataFrame({"cgY": [0.5]}))


def test_paired_bootstrap_detects_rejuvenation():
    rng = np.random.default_rng(0)
    day0 = rng.normal(60, 2, 12)
    day13 = day0 - rng.normal(25, 3, 12)  # rejuvenation by ~25y
    delta, lo, hi = paired_bootstrap_delta(day13, day0, n_boot=2000)
    assert delta < -15 and hi < 0


def test_benjamini_hochberg_monotone():
    p = np.array([0.001, 0.01, 0.04, 0.2, 0.9])
    q = benjamini_hochberg(p)
    assert np.all(np.diff(q[np.argsort(p)]) >= -1e-12)
    assert q[0] < 0.01
