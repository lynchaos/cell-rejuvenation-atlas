"""Unit tests for the WOT-style transport engine (synthetic populations)."""
import numpy as np

from src.module2_reprogramming_trajectory.ot_transport import (
    estimate_growth_rates,
    fate_probability,
    push_forward,
    transport_map,
)


def _two_clusters(seed=0):
    rng = np.random.default_rng(seed)
    src = rng.normal(0, 0.1, size=(50, 3))
    tgt = np.vstack([rng.normal(1, 0.1, size=(25, 3)), rng.normal(-1, 0.1, size=(25, 3))])
    return src, tgt


def test_transport_map_rows_match_source_mass():
    src, tgt = _two_clusters()
    g = transport_map(src, tgt, growth=0.0, dt=1.0)
    assert g.shape == (50, 50)
    np.testing.assert_allclose(g.sum(axis=1).sum(), 1.0, rtol=0.1)


def test_growth_increases_source_mass():
    src, tgt = _two_clusters()
    g0 = transport_map(src, tgt, growth=0.0)
    g1 = transport_map(src, tgt, growth=1.0, dt=1.0)
    assert g1.sum() > g0.sum()


def test_estimate_growth_rates():
    assert np.isclose(estimate_growth_rates(100, 200, 1.0), np.log(2))


def test_fate_probability_prefers_nearby_fate():
    rng = np.random.default_rng(3)
    # day0: two groups; day1: two groups aligned by position
    d0 = np.vstack([rng.normal(0, 0.05, (30, 2)), rng.normal(5, 0.05, (30, 2))])
    d1 = np.vstack([rng.normal(0, 0.05, (30, 2)), rng.normal(5, 0.05, (30, 2))])
    coupling = transport_map(d0, d1, epsilon=0.01)
    mask = np.zeros(60, dtype=bool)
    mask[:30] = True  # fate = the cluster near position 0
    probs = fate_probability([coupling], mask)
    assert probs[:30].mean() > probs[30:].mean()


def test_push_forward_conserves_mass():
    src, tgt = _two_clusters()
    g = transport_map(src, tgt)
    out = push_forward(np.ones(50), [g])
    assert np.isclose(out.sum(), 1.0)
