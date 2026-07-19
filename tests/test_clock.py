"""Unit tests for the epigenetic-clock machinery (synthetic data)."""
import gzip

import numpy as np
import pandas as pd
import pytest

from src.common.stats import (
    benjamini_hochberg,
    effect_size_cohens_d,
    paired_bootstrap_delta,
    unpaired_bootstrap_delta,
)
from src.module1_rejuvenation_clock.analyze_gill import (
    infer_arm,
    infer_day,
    load_beta,
)
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


def test_infer_day_title_formats():
    assert infer_day("MPTR day 13 rep 2") == 13.0
    assert infer_day("O2_transiently_reprogrammed_13days_exp1") == 13.0
    assert np.isnan(infer_day("O1 Fib"))


def test_infer_arm_gill_titles():
    assert infer_arm("O3_transiently_reprogrammed_15days_exp1") == "transient"
    assert infer_arm("O2_transient_reprogramming_intermediate_17days_exp2") == "transient"
    assert infer_arm("O1_failed_to_transiently_reprogram_15days_exp1") == "failed"
    assert infer_arm("O2_failing_to_transiently_reprogram_intermediate_17days_exp2") == "failed"
    assert infer_arm("O1_negative_control_15days_exp1") == "control"
    assert infer_arm("O1 Fib") == "baseline"


def test_load_beta_sniffs_delimiter_and_drops_pval_columns(tmp_path):
    # GEO ships comma-separated content under a .txt.gz name, with a
    # 'Detection Pval' column interleaved after every sample.
    with gzip.open(tmp_path / "matrix.txt.gz", "wt") as fh:
        fh.write(
            "ID_REF,s1,Detection Pval,s2,Detection Pval\n"
            "cg1,0.5,0.001,0.6,0.002\n"
            "cg2,0.7,0.0,0.8,0.0\n"
            "cg3,0.9,0.0,0.1,0.0\n"
        )
    df = load_beta(str(tmp_path))
    assert df.shape == (2, 3)  # samples x probes after transpose
    assert set(df.columns) == {"cg1", "cg2", "cg3"}
    assert set(df.index) == {"s1", "s2"}


def test_unpaired_bootstrap_detects_rejuvenation():
    rng = np.random.default_rng(0)
    control = rng.normal(58, 3, 9)
    transient = control - rng.normal(20, 3, 9)
    delta, lo, hi = unpaired_bootstrap_delta(transient, control, n_boot=2000)
    assert delta < -15 and hi < 0
    assert effect_size_cohens_d(transient, control) < -2
