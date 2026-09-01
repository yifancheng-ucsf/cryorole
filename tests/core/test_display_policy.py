from __future__ import annotations

import numpy as np

from cryorole.core.display_policy import (
    EULER_PROJECTIONS,
    LEGACY_DISPLAY_STYLE,
    resolve_color_scale,
    resolve_display_indices,
    sort_display_indices,
)


def test_threshold_and_top_fraction_are_stable_and_source_ordered() -> None:
    values = np.asarray([1.0, 2.0, 2.0, 3.0])

    assert resolve_display_indices(values, threshold=2.0).tolist() == [1, 2, 3]
    assert resolve_display_indices(values, top_fraction=0.5).tolist() == [1, 2, 3]


def test_top_fraction_retains_cutoff_ties() -> None:
    values = np.asarray([1.0, 2.0, 3.0, 3.0, 4.0])

    assert resolve_display_indices(values, top_fraction=0.40).tolist() == [2, 3, 4]
    assert sort_display_indices(
        np.asarray([0, 1, 2, 3]),
        values,
        order="ascending",
    ).tolist() == [0, 1, 2, 3]


def test_legacy_style_projection_order_and_tail_jump_color_scale() -> None:
    assert tuple(item[0] for item in EULER_PROJECTIONS) == (
        "alpha_beta",
        "beta_gamma",
        "alpha_gamma",
    )
    assert LEGACY_DISPLAY_STYLE["color_map"] == "rainbow_r"
    assert LEGACY_DISPLAY_STYLE["point_size"] == 1.0
    assert LEGACY_DISPLAY_STYLE["point_alpha"] == 1.0
    assert LEGACY_DISPLAY_STYLE["aspect"] == "equal"
    assert LEGACY_DISPLAY_STYLE["colorbar_position"] == "bottom"

    scale = resolve_color_scale(
        np.asarray([1.0, 2.0, 100.0]),
        display_outlier_mask=np.asarray([False, False, True]),
        visual_style="legacy",
        color_vmin=None,
        color_vmax=None,
        display_threshold=None,
    )
    assert scale.vmax == 2.0
    assert scale.vmax_source == "tail_jump_display_outliers"
