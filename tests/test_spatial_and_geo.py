"""Unit tests for spatial helpers and GEO utilities."""
import numpy as np
import pandas as pd
import pytest

from src.common.geo import _bucket
from src.module4_spatial_aging.squidpy_analysis import AGING_MARKERS, _age_order


def test_geo_bucket():
    assert _bucket("GSE165179") == "GSE165nnn"
    assert _bucket("gse122662") == "GSE122nnn"


def test_aging_markers_are_canonical():
    # Guard against typos drifting from the Allen et al. marker panel
    assert {"C4b", "Gfap", "B2m"}.issubset(set(AGING_MARKERS))


def test_age_order_numeric_prefix():
    assert _age_order(pd.Series(["90wk", "4wk", "24wk"])) == ["4wk", "24wk", "90wk"]
