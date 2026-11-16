"""The shared density-image raster used by the density figures and the GUI."""

from __future__ import annotations

import numpy as np
import pytest

from cryorole.core.display_policy import visible_value_raster


def test_each_pixel_keeps_the_value_a_scatter_would_show():
    x = np.array([0.1, 0.15, 0.9, 0.5, 2.0])  # the last point lies outside the range
    y = np.array([0.1, 0.12, 0.9, 0.5, 0.5])
    values = np.array([1.0, 3.0, 2.0, np.nan, 9.0])  # NaN is never drawn
    ascending = visible_value_raster(x, y, values, x_range=(0, 1), y_range=(0, 1), width=4, height=2)
    assert ascending.shape == (2, 4)
    assert ascending[0, 0] == 3.0  # two points share the cell: the larger value is on top
    assert ascending[1, 3] == 2.0  # row 0 is at y_range[0], so y=0.9 is in the last row
    assert np.isnan(ascending).sum() == 6
    descending = visible_value_raster(
        x, y, values, x_range=(0, 1), y_range=(0, 1), width=4, height=2, order="descending"
    )
    assert descending[0, 0] == 1.0


def test_none_order_behaves_like_ascending_and_bad_order_is_rejected():
    x = y = np.array([0.5, 0.5])
    values = np.array([4.0, 1.0])
    grid = visible_value_raster(x, y, values, x_range=(0, 1), y_range=(0, 1), width=1, height=1, order="none")
    assert grid[0, 0] == 4.0
    with pytest.raises(ValueError):
        visible_value_raster(x, y, values, x_range=(0, 1), y_range=(0, 1), width=1, height=1, order="random")
