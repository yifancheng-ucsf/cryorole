from __future__ import annotations

from pathlib import Path

import pytest
from PIL import Image

from cryorole.animation.compositor import (
    compose_frame_sequences,
    resolve_composite_layout,
    validate_dual_structure_horizontal_crop,
    validate_png_sequence,
)


def _write_frames(
    directory: Path,
    *,
    count: int,
    size: tuple[int, int],
    color: tuple[int, int, int],
) -> tuple[Path, ...]:
    directory.mkdir(parents=True, exist_ok=True)
    paths = []
    for index in range(count):
        path = directory / f"frame_{index:06d}.png"
        Image.new("RGB", size, color).save(path)
        paths.append(path)
    return tuple(paths)


def test_compositor_is_deterministic_contained_and_preserves_sources(tmp_path):
    landscape = _write_frames(
        tmp_path / "landscape",
        count=2,
        size=(300, 100),
        color=(255, 0, 0),
    )
    structure = _write_frames(
        tmp_path / "structure",
        count=2,
        size=(100, 100),
        color=(0, 0, 255),
    )
    before = {path: path.read_bytes() for path in landscape + structure}

    result = compose_frame_sequences(
        landscape_dir=tmp_path / "landscape",
        structure_dir=tmp_path / "structure",
        output_dir=tmp_path / "composite",
        expected_count=2,
        canvas_size=(1920, 1080),
        layout_name="side_by_side",
        background_color="#ffffff",
    )

    assert result.frame_count == 2
    assert result.dimensions == (1920, 1080)
    assert [path.name for path in result.paths] == [
        "frame_000000.png",
        "frame_000001.png",
    ]
    assert result.layout.resize_filter == "lanczos"
    assert result.layout.landscape_destination[2] / result.layout.landscape_destination[3] == 3
    assert result.layout.structure_destination[2] == result.layout.structure_destination[3]
    assert all(path.read_bytes() == before[path] for path in landscape + structure)

    with Image.open(result.paths[0]) as image:
        assert image.mode == "RGB"
        assert image.size == (1920, 1080)
        lx, ly, lw, lh = result.layout.landscape_destination
        sx, sy, sw, sh = result.layout.structure_destination
        assert image.getpixel((lx + lw // 2, ly + lh // 2)) == (255, 0, 0)
        assert image.getpixel((sx + sw // 2, sy + sh // 2)) == (0, 0, 255)
        assert image.getpixel((0, 0)) == (255, 255, 255)


def test_stacked_layout_is_structure_first_fixed_and_uses_64_36_regions():
    layout = resolve_composite_layout(
        canvas_size=(1920, 1080),
        landscape_dimensions=(1800, 600),
        structure_dimensions=(900, 900),
        layout_name="stacked",
        background_color="#ffffff",
    )

    assert layout.name == "stacked"
    assert layout.layout_version == "4"
    assert layout.structure_fraction == pytest.approx(0.64)
    assert layout.landscape_fraction == pytest.approx(0.36)
    assert layout.structure_region[0] == layout.landscape_region[0] == 24
    assert layout.structure_region[2] == layout.landscape_region[2] == 1872
    assert layout.structure_region[1] < layout.landscape_region[1]
    assert (
        layout.structure_region[1]
        + layout.structure_region[3]
        + layout.gap_pixels
        == layout.landscape_region[1]
    )
    usable_height = 1080 - 2 * layout.padding_pixels - layout.gap_pixels
    assert layout.structure_region[3] / usable_height == pytest.approx(0.64, abs=0.001)
    assert layout.landscape_region[3] / usable_height == pytest.approx(0.36, abs=0.001)


def test_side_by_side_layout_remains_supported():
    layout = resolve_composite_layout(
        canvas_size=(1920, 1080),
        landscape_dimensions=(1800, 600),
        structure_dimensions=(900, 900),
        layout_name="side_by_side",
        background_color="#ffffff",
    )
    assert layout.name == "side_by_side"
    assert layout.landscape_region[0] < layout.structure_region[0]
    assert layout.structure_fraction is None
    assert layout.landscape_fraction is None


def test_dual_stacked_layout_has_equal_fixed_structure_regions_and_small_gap(tmp_path):
    layout = resolve_composite_layout(
        canvas_size=(1920, 1080),
        landscape_dimensions=(1800, 600),
        structure_dimensions=(900, 900),
        secondary_structure_dimensions=(900, 900),
        layout_name="stacked",
        background_color="#ffffff",
    )
    primary = layout.primary_structure_region
    secondary = layout.secondary_structure_region
    assert layout.structure_view_count == 2
    assert primary is not None and secondary is not None
    assert primary[2] == secondary[2]
    assert secondary[0] - (primary[0] + primary[2]) == layout.dual_view_gap_pixels
    assert layout.dual_view_gap_pixels <= 1920 * 0.015
    assert layout.primary_structure_destination is not None
    assert layout.secondary_structure_destination is not None

    _write_frames(tmp_path / "landscape", count=2, size=(300, 100), color=(255, 0, 0))
    _write_frames(tmp_path / "primary", count=2, size=(100, 100), color=(0, 0, 255))
    _write_frames(tmp_path / "secondary", count=2, size=(100, 100), color=(0, 255, 0))
    result = compose_frame_sequences(
        landscape_dir=tmp_path / "landscape",
        structure_dir=tmp_path / "primary",
        secondary_structure_dir=tmp_path / "secondary",
        output_dir=tmp_path / "composite",
        expected_count=2,
        canvas_size=(1920, 1080),
        layout_name="stacked",
        background_color="#ffffff",
    )
    assert result.secondary_structure_validation is not None
    with Image.open(result.paths[0]) as image:
        px, py, pw, ph = result.layout.primary_structure_destination
        sx, sy, sw, sh = result.layout.secondary_structure_destination
        assert image.getpixel((px + pw // 2, py + ph // 2)) == (0, 0, 255)
        assert image.getpixel((sx + sw // 2, sy + sh // 2)) == (0, 255, 0)


def test_dual_structure_view_rejects_side_by_side():
    with pytest.raises(ValueError, match="dual structure view.*stacked"):
        resolve_composite_layout(
            canvas_size=(1920, 1080),
            landscape_dimensions=(1800, 600),
            structure_dimensions=(900, 900),
            secondary_structure_dimensions=(900, 900),
            layout_name="side_by_side",
            background_color="#ffffff",
        )


@pytest.mark.parametrize(
    "fraction",
    [-0.01, 0.5, 1.0, float("nan"), float("inf"), float("-inf")],
)
def test_dual_structure_horizontal_crop_rejects_invalid_fraction(fraction):
    with pytest.raises(ValueError, match="finite.*0 <=.*< 0.5"):
        validate_dual_structure_horizontal_crop(
            fraction,
            has_secondary=True,
            layout_name="stacked",
        )


def test_nonzero_dual_structure_crop_requires_dual_stacked_layout():
    with pytest.raises(ValueError, match="only.*dual-view.*stacked"):
        resolve_composite_layout(
            canvas_size=(1920, 1080),
            landscape_dimensions=(1800, 600),
            structure_dimensions=(1200, 840),
            layout_name="stacked",
            background_color="#ffffff",
            dual_structure_horizontal_crop=0.1,
        )
    with pytest.raises(ValueError, match="dual structure view.*stacked"):
        resolve_composite_layout(
            canvas_size=(1920, 1080),
            landscape_dimensions=(1800, 600),
            structure_dimensions=(1200, 840),
            secondary_structure_dimensions=(1200, 840),
            layout_name="side_by_side",
            background_color="#ffffff",
            dual_structure_horizontal_crop=0.1,
        )


def test_zero_crop_preserves_existing_dual_geometry():
    common = {
        "canvas_size": (1920, 1080),
        "landscape_dimensions": (1800, 600),
        "structure_dimensions": (1200, 840),
        "secondary_structure_dimensions": (1200, 840),
        "layout_name": "stacked",
        "background_color": "#ffffff",
    }
    existing = resolve_composite_layout(**common)
    explicit_zero = resolve_composite_layout(
        **common,
        dual_structure_horizontal_crop=0,
    )
    assert explicit_zero == existing


def test_ten_percent_dual_crop_is_fixed_and_aligned_inward():
    layout = resolve_composite_layout(
        canvas_size=(1920, 1080),
        landscape_dimensions=(1800, 600),
        structure_dimensions=(1200, 840),
        secondary_structure_dimensions=(1200, 840),
        layout_name="stacked",
        background_color="#ffffff",
        dual_structure_horizontal_crop=0.1,
    )

    assert layout.dual_structure_horizontal_crop_fraction == pytest.approx(0.1)
    assert layout.primary_structure_source_crop == (120, 0, 1080, 840)
    assert layout.secondary_structure_source_crop == (120, 0, 1080, 840)
    assert layout.primary_structure_cropped_dimensions == (960, 840)
    assert layout.secondary_structure_cropped_dimensions == (960, 840)
    assert layout.primary_structure_alignment == "right"
    assert layout.secondary_structure_alignment == "left"
    primary_region = layout.primary_structure_region
    secondary_region = layout.secondary_structure_region
    primary_destination = layout.primary_structure_destination
    secondary_destination = layout.secondary_structure_destination
    assert primary_region is not None and secondary_region is not None
    assert primary_destination is not None and secondary_destination is not None
    assert primary_destination[0] + primary_destination[2] == (
        primary_region[0] + primary_region[2]
    )
    assert secondary_destination[0] == secondary_region[0]
    assert secondary_destination[0] - (
        primary_destination[0] + primary_destination[2]
    ) == layout.dual_view_gap_pixels


def test_dual_crop_reuses_geometry_and_does_not_modify_sources(tmp_path):
    landscape = _write_frames(
        tmp_path / "landscape",
        count=2,
        size=(300, 100),
        color=(255, 0, 0),
    )
    primary = _write_frames(
        tmp_path / "primary",
        count=2,
        size=(120, 84),
        color=(0, 0, 255),
    )
    secondary = _write_frames(
        tmp_path / "secondary",
        count=2,
        size=(120, 84),
        color=(0, 255, 0),
    )
    sources = landscape + primary + secondary
    before = {path: path.read_bytes() for path in sources}

    result = compose_frame_sequences(
        landscape_dir=tmp_path / "landscape",
        structure_dir=tmp_path / "primary",
        secondary_structure_dir=tmp_path / "secondary",
        output_dir=tmp_path / "composite",
        expected_count=2,
        canvas_size=(1920, 1080),
        layout_name="stacked",
        background_color="#ffffff",
        dual_structure_horizontal_crop=0.1,
    )

    assert result.layout.primary_structure_source_crop == (12, 0, 108, 84)
    assert result.layout.secondary_structure_source_crop == (12, 0, 108, 84)
    assert all(path.read_bytes() == before[path] for path in sources)
    output_sizes = []
    for path in result.paths:
        with Image.open(path) as image:
            output_sizes.append(image.size)
    assert len(set(output_sizes)) == 1


def test_compositor_rejects_missing_invalid_and_inconsistent_frames(tmp_path):
    _write_frames(tmp_path / "landscape", count=2, size=(30, 10), color=(1, 2, 3))
    _write_frames(tmp_path / "structure", count=1, size=(10, 10), color=(4, 5, 6))
    with pytest.raises(RuntimeError, match="structure frame indices"):
        compose_frame_sequences(
            landscape_dir=tmp_path / "landscape",
            structure_dir=tmp_path / "structure",
            output_dir=tmp_path / "composite",
            expected_count=2,
            canvas_size=(1920, 1080),
            layout_name="side_by_side",
            background_color="white",
        )

    invalid = tmp_path / "invalid"
    invalid.mkdir()
    (invalid / "frame_000000.png").write_bytes(b"not-png")
    with pytest.raises(RuntimeError, match="valid PNG"):
        validate_png_sequence(invalid, expected_count=1, label="test")

    inconsistent = tmp_path / "inconsistent"
    _write_frames(inconsistent, count=2, size=(10, 10), color=(1, 1, 1))
    Image.new("RGB", (11, 10), "black").save(inconsistent / "frame_000001.png")
    with pytest.raises(RuntimeError, match="dimensions are inconsistent"):
        validate_png_sequence(inconsistent, expected_count=2, label="test")


@pytest.mark.parametrize("canvas_size", [(0, 1080), (1921, 1080), (1920, 1079)])
def test_compositor_requires_positive_even_yuv420_canvas(tmp_path, canvas_size):
    _write_frames(tmp_path / "landscape", count=1, size=(30, 10), color=(1, 2, 3))
    _write_frames(tmp_path / "structure", count=1, size=(10, 10), color=(4, 5, 6))
    with pytest.raises(ValueError, match="positive even"):
        compose_frame_sequences(
            landscape_dir=tmp_path / "landscape",
            structure_dir=tmp_path / "structure",
            output_dir=tmp_path / "composite",
            expected_count=1,
            canvas_size=canvas_size,
            layout_name="side_by_side",
            background_color="white",
        )
