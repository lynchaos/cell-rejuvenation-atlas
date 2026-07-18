"""A compact Waddington-OT-style transport engine.

Implements the core of Schiebinger et al. (2019): consecutive time points are
coupled by entropic optimal transport in PCA space, with cell mass reweighted
by estimated growth (death/birth) rates. Pure functions over numpy arrays so
the math is unit-testable on small synthetic populations.
"""
from __future__ import annotations

import numpy as np
import ot  # POT: Python Optimal Transport
from scipy.spatial.distance import cdist


def estimate_growth_rates(n_source: int, n_target: int, dt: float) -> float:
    """Net per-day growth g such that n_target ~ n_source * exp(g * dt)."""
    if dt <= 0:
        raise ValueError("dt must be positive")
    return np.log(max(n_target, 1) / max(n_source, 1)) / dt


def transport_map(
    x_source: np.ndarray,
    x_target: np.ndarray,
    growth: float = 0.0,
    dt: float = 1.0,
    epsilon: float | None = None,
) -> np.ndarray:
    """Entropic OT coupling between two time points.

    Parameters
    ----------
    x_source, x_target : cells x PCs arrays for the two time points.
    growth : per-day net growth rate; source mass is scaled by exp(g*dt).
    epsilon : entropic regularization; default = 0.05 * median pairwise cost.

    Returns
    -------
    (n_source, n_target) coupling matrix whose rows sum to the source mass.
    """
    cost = cdist(x_source, x_target, metric="sqeuclidean")
    cost /= max(cost.max(), 1e-12)
    if epsilon is None:
        epsilon = 0.05 * float(np.median(cost))
    a = np.full(x_source.shape[0], np.exp(growth * dt) / x_source.shape[0])
    b = np.full(x_target.shape[0], 1.0 / x_target.shape[0])
    return ot.bregman.sinkhorn(a, b, cost, reg=epsilon)


def push_forward(mass: np.ndarray, couplings: list[np.ndarray]) -> np.ndarray:
    """Propagate an initial mass distribution through a chain of couplings.

    Each coupling row-sums to the source mass; we normalize rows to get
    transition probabilities and multiply through the chain.
    """
    p = np.asarray(mass, dtype=float)
    p = p / p.sum()
    for gamma in couplings:
        row_sums = gamma.sum(axis=1, keepdims=True)
        transition = np.divide(gamma, row_sums, where=row_sums > 0)
        p = p @ transition
    return p


def fate_probability(
    couplings: list[np.ndarray], terminal_fate_mask: np.ndarray
) -> np.ndarray:
    """Probability that each initial cell ends in a given terminal fate.

    terminal_fate_mask : boolean vector over cells at the final time point.
    Returns one fate probability per cell at the *first* time point.
    """
    n0 = couplings[0].shape[0]
    eye = np.eye(n0)
    mask = np.asarray(terminal_fate_mask, dtype=float)
    probs = np.empty(n0)
    for i in range(n0):
        traj = push_forward(eye[i], couplings)
        probs[i] = traj @ mask
    return probs
