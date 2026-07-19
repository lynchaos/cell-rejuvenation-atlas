"""Unit tests for integration scoring and SASP core logic (synthetic data)."""
import numpy as np
import pandas as pd
import pytest

from src.module3_multiomics_integration.aging_signature_score import (
    collapse_by_symbol,
    normalize_log_cpm,
    parse_browder_sample,
    tissue_stats,
)
from src.module3_multiomics_integration.factor_model import consensus_axes, modality_pcs
from src.module3_multiomics_integration.rejuvenation_score import (
    aging_signature,
    latent_shift,
    project_onto_signature,
    signature_consistency,
)
from src.module5_proteomics_sasp.sasp_analysis import core_sasp, differential_abundance


def _aging_data(seed=0, n_genes=300):
    rng = np.random.default_rng(seed)
    young = rng.normal(10, 1, size=(40, n_genes))
    old = young.mean(axis=0) + rng.normal(0, 0.3, size=(40, n_genes)) + rng.normal(0, 1, size=(40, n_genes))
    # plant aging genes up in old
    old[:, :30] += 2.0
    expr = pd.DataFrame(
        np.vstack([young, old]), columns=[f"G{i}" for i in range(n_genes)],
        index=[f"y{i}" for i in range(40)] + [f"o{i}" for i in range(40)],
    )
    group = pd.Series(["young"] * 40 + ["old"] * 40, index=expr.index)
    return expr, group


def test_aging_signature_and_projection():
    expr, group = _aging_data()
    sig = aging_signature(expr, group)
    assert sig.iloc[:30].mean() > sig.iloc[30:].mean()  # planted genes on top
    scores = project_onto_signature(expr, sig)
    assert scores.loc[[i for i in expr.index if i.startswith("o")]].mean() > \
           scores.loc[[i for i in expr.index if i.startswith("y")]].mean()


def test_signature_consistency_identical():
    expr, group = _aging_data()
    sig = aging_signature(expr, group)
    assert signature_consistency(sig, sig) > 0.99


def test_latent_shift_direction():
    zc = np.zeros((10, 5))
    zt = np.ones((10, 5))
    assert np.allclose(latent_shift(zc, zt), 1.0)


def test_factor_model_recovers_shared_axis():
    rng = np.random.default_rng(0)
    shared = rng.normal(size=100)
    a = pd.DataFrame({"g1": shared + rng.normal(0, 0.1, 100),
                      "g2": rng.normal(0, 1, 100)})
    b = pd.DataFrame({"p1": shared + rng.normal(0, 0.1, 100),
                      "p2": rng.normal(0, 1, 100)})
    axes = consensus_axes(modality_pcs(a, n_pcs=2), modality_pcs(b, n_pcs=2), n_components=2)
    assert axes["canonical_corr"].iloc[0] > 0.8


def test_core_sasp_requires_two_inducers():
    idx = pd.Index(["P1", "P2", "P3"])
    de1 = pd.DataFrame({"padj": [0.01, 0.5, 0.01], "log2fc": [1.0, 0.1, 2.0]}, index=idx)
    de2 = pd.DataFrame({"padj": [0.01, 0.01, 0.5], "log2fc": [1.5, 1.0, 0.2]}, index=idx)
    core = core_sasp({"IR": de1, "RAS": de2})
    assert list(core) == ["P1"]


def test_differential_abundance_finds_planted_effect():
    rng = np.random.default_rng(1)
    ctrl = rng.normal(100, 10, size=(8, 50))
    sen = rng.normal(100, 10, size=(8, 50))
    sen[:, :5] += 100  # strongly up proteins
    m = pd.DataFrame(np.vstack([ctrl, sen]), columns=[f"P{i}" for i in range(50)])
    design = pd.Series(["control"] * 8 + ["senescent"] * 8, index=m.index)
    de = differential_abundance(m, design)
    assert set(de.head(5).index) == {f"P{i}" for i in range(5)}


def test_parse_browder_sample():
    assert parse_browder_sample("Skeletal_Muscle_4F_m3") == ("Skeletal_Muscle", "4F")
    assert parse_browder_sample("Lung_Control_m12") == ("Lung", "Control")
    assert parse_browder_sample("Skin_4F_m1") == ("Skin", "4F")
    with pytest.raises(ValueError):
        parse_browder_sample("no_condition_token")


def test_collapse_by_symbol_strips_versions_and_sums():
    counts = pd.DataFrame(
        {"s1": [1, 2, 4, 8], "s2": [0, 1, 0, 2]},
        index=["ENSMUSG00000000001.4", "ENSMUSG00000000003", "ENSMUSG00000000001.4", "ENSMUSG99999999999"],
    )
    mapping = {"ENSMUSG00000000001": "Nans", "ENSMUSG00000000003": "Oct4"}
    out = collapse_by_symbol(counts, mapping)
    assert set(out.index) == {"Nans", "Oct4"}  # unmapped id dropped
    assert out.loc["Nans", "s1"] == 5  # duplicated symbol summed


def test_normalize_log_cpm():
    counts = pd.DataFrame({"a": [10, 0], "b": [5, 5]}, index=["g1", "g2"])
    out = normalize_log_cpm(counts)
    assert np.isclose(out.loc["g1", "a"], np.log1p(1e6))
    assert np.isclose(out.loc["g1", "b"], np.log1p(0.5e6))


def test_tissue_stats_detects_younger_4f():
    rng = np.random.default_rng(2)
    meta = pd.DataFrame(
        {"tissue": ["Lung"] * 8,
         "condition": ["Control"] * 4 + ["4F"] * 4},
        index=[f"s{i}" for i in range(8)],
    )
    scores = pd.Series(
        np.concatenate([rng.normal(10, 0.5, 4), rng.normal(8, 0.5, 4)]),
        index=meta.index,
    )
    stats = tissue_stats(scores, meta)
    row = stats.iloc[0]
    assert row["delta"] < 0
    assert row["pvalue"] < 0.05
    assert row["cohens_d"] < 0
