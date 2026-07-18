"""Signature-based rejuvenation scoring (pure numpy/pandas — unit-testable).

The idea: derive an aging axis per cell type from an independent atlas
(young vs. old differential expression), then project treatment cells onto
that axis. A negative shift for reprogrammed cells = movement toward the
young state.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def aging_signature(
    expr: pd.DataFrame, group: pd.Series, young_label: str = "young", old_label: str = "old"
) -> pd.Series:
    """Mean log-expression difference old - young per gene (the 'aging axis')."""
    young = expr.loc[group[group == young_label].index]
    old = expr.loc[group[group == old_label].index]
    return old.mean(axis=0) - young.mean(axis=0)


def project_onto_signature(expr: pd.DataFrame, signature: pd.Series) -> pd.Series:
    """Score each cell by its projection on the aging axis (cosine-normalized)."""
    common = expr.columns.intersection(signature.index)
    if len(common) < 10:
        raise ValueError(f"Too few shared genes for projection ({len(common)}).")
    x = expr[common].to_numpy(float)
    s = signature[common].to_numpy(float)
    s_norm = np.linalg.norm(s)
    if s_norm == 0:
        raise ValueError("Degenerate signature (zero vector).")
    scores = (x @ s) / s_norm
    return pd.Series(scores, index=expr.index, name="aging_score")


def signature_consistency(sig_a: pd.Series, sig_b: pd.Series) -> float:
    """Spearman correlation between two aging signatures (cross-study check)."""
    from scipy.stats import spearmanr

    common = sig_a.index.intersection(sig_b.index)
    if len(common) < 10:
        return np.nan
    return float(spearmanr(sig_a[common], sig_b[common]).statistic)


def latent_shift(z_control: np.ndarray, z_treated: np.ndarray) -> np.ndarray:
    """Mean latent displacement vector control -> treated (scGen-style arithmetic)."""
    return np.asarray(z_treated, float).mean(axis=0) - np.asarray(z_control, float).mean(axis=0)
