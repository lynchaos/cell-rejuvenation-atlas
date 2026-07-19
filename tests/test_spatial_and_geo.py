"""Unit tests for spatial helpers and GEO utilities."""
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from src.common.geo import _bucket
from src.module4_spatial_aging.spatial_analysis import (
    AGING_MARKERS,
    _age_order,
    morans_i,
    neighborhood_enrichment_z,
    spatial_connectivities,
)


def test_geo_bucket():
    assert _bucket("GSE165179") == "GSE165nnn"
    assert _bucket("gse122662") == "GSE122nnn"


def test_aging_markers_are_canonical():
    # Guard against typos drifting from the Allen et al. marker panel
    assert {"C4b", "Gfap", "B2m"}.issubset(set(AGING_MARKERS))


def test_age_order_numeric_prefix():
    assert _age_order(pd.Series(["90wk", "4wk", "24wk"])) == ["4wk", "24wk", "90wk"]


def _path_graph(n: int) -> sp.csr_matrix:
    """Adjacency of a 1-D path: i connected to i-1 and i+1."""
    rows = list(range(n - 1)) + list(range(1, n))
    cols = list(range(1, n)) + list(range(n - 1))
    return sp.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(n, n))


def test_morans_i_sign_matches_clustering():
    graph = _path_graph(6)
    clustered = np.array([1.0, 1.0, 1.0, 0.0, 0.0, 0.0])
    alternating = np.array([1.0, 0.0, 1.0, 0.0, 1.0, 0.0])
    assert morans_i(graph, clustered) > 0
    assert morans_i(graph, alternating) < 0


def test_morans_i_constant_input_is_nan():
    assert np.isnan(morans_i(_path_graph(4), np.ones(4)))


def test_spatial_connectivities_symmetric_no_self_loops():
    rng = np.random.default_rng(0)
    graph = spatial_connectivities(rng.normal(size=(50, 2)), n_neigh=6)
    assert (graph != graph.T).nnz == 0  # symmetric
    assert graph.diagonal().sum() == 0  # no self loops
    assert graph.data.max() == 1  # binary


def test_nhood_enrichment_segregated_clusters():
    rng = np.random.default_rng(1)
    a = rng.normal(0, 1, size=(60, 2))
    b = rng.normal(100, 1, size=(60, 2))
    coords = np.vstack([a, b])
    labels = pd.Series(["A"] * 60 + ["B"] * 60)
    z = neighborhood_enrichment_z(spatial_connectivities(coords), labels, n_perms=200, seed=0)
    assert z.loc["A", "A"] > 0 and z.loc["B", "B"] > 0  # within-cluster contacts enriched
    assert z.loc["A", "B"] < 0  # between-cluster contacts depleted
