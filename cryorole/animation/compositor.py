"""Deterministic Phase 4 composition of validated animation frame sequences."""

from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path

from PIL import Image, ImageColor, UnidentifiedImageError


COMPOSITE_LAYOUT_VERSION = "4"
STACKED_STRUCTURE_FRACTION = 0.64
STACKED_LANDSCAPE_FRACTION = 0.36
DUAL_VIEW_GAP_FRACTION = 0.0125


@dataclass(frozen=True)
class PNGSequenceValidation:
    """Validated contiguous PNG frame sequence."""

    label: str
    frame_count: int
    dimensions: tuple[int, int]
    paths: tuple[Path, ...]


@dataclass(frozen=True)
class CompositeLayout:
    """Auditable fixed geometry using ``(x, y, width, height)`` rectangles."""

    name: str
    canvas_size: tuple[int, int]
    landscape_region: tuple[int, int, int, int]
    structure_region: tuple[int, int, int, int]
    landscape_destination: tuple[int, int, int, int]
    structure_destination: tuple[int, int, int, int]
    background_color: str
    primary_structure_region: tuple[int, int, int, int] | None = None
    secondary_structure_region: tuple[int, int, int, int] | None = None
    primary_structure_destination: tuple[int, int, int, int] | None = None
    secondary_structure_destination: tuple[int, int, int, int] | None = None
    primary_structure_source_crop: tuple[int, int, int, int] | None = None
    secondary_structure_source_crop: tuple[int, int, int, int] | None = None
    primary_structure_cropped_dimensions: tuple[int, int] | None = None
    secondary_structure_cropped_dimensions: tuple[int, int] | None = None
    primary_structure_alignment: str = "center"
    secondary_structure_alignment: str | None = None
    dual_structure_horizontal_crop_fraction: float = 0.0
    structure_view_count: int = 1
    dual_view_gap_pixels: int = 0
    layout_version: str = COMPOSITE_LAYOUT_VERSION
    structure_fraction: float | None = None
    landscape_fraction: float | None = None
    resize_filter: str = "lanczos"
    aspect_policy: str = "contain_letterbox_no_crop"
    padding_pixels: int = 24
    gap_pixels: int = 24


@dataclass(frozen=True)
class CompositeResult:
    """Composite frames plus validated input and output geometry."""

    frame_count: int
    dimensions: tuple[int, int]
    paths: tuple[Path, ...]
    layout: CompositeLayout
    landscape_validation: PNGSequenceValidation
    structure_validation: PNGSequenceValidation
    secondary_structure_validation: PNGSequenceValidation | None = None


def validate_png_sequence(
    directory: str | Path,
    *,
    expected_count: int,
    label: str,
) -> PNGSequenceValidation:
    """Validate exact contiguous, non-empty, fixed-size PNG frames."""

    if expected_count <= 0:
        raise ValueError("expected frame count must be positive")
    root = Path(directory)
    paths = tuple(sorted(root.glob("frame_*.png"))) if root.is_dir() else ()
    expected_names = tuple(f"frame_{index:06d}.png" for index in range(expected_count))
    names = tuple(path.name for path in paths)
    if names != expected_names:
        raise RuntimeError(
            f"{label} frame indices are not contiguous or do not match "
            f"trajectory count {expected_count}: {names}"
        )
    dimensions: set[tuple[int, int]] = set()
    for path in paths:
        if path.stat().st_size <= 0:
            raise RuntimeError(f"{label} frame is empty: {path}")
        try:
            with Image.open(path) as image:
                if image.format != "PNG":
                    raise RuntimeError(f"{label} frame is not a valid PNG: {path}")
                dimensions.add(tuple(map(int, image.size)))
                image.verify()
        except (OSError, UnidentifiedImageError) as exc:
            raise RuntimeError(f"{label} frame is not a valid PNG: {path}") from exc
    if len(dimensions) != 1:
        raise RuntimeError(f"{label} frame dimensions are inconsistent")
    dimensions_value = next(iter(dimensions))
    if dimensions_value[0] <= 0 or dimensions_value[1] <= 0:
        raise RuntimeError(f"{label} frame dimensions must be positive")
    return PNGSequenceValidation(
        label=label,
        frame_count=len(paths),
        dimensions=dimensions_value,
        paths=paths,
    )


def compose_frame_sequences(
    *,
    landscape_dir: str | Path,
    structure_dir: str | Path,
    secondary_structure_dir: str | Path | None = None,
    output_dir: str | Path,
    expected_count: int,
    canvas_size: tuple[int, int],
    layout_name: str,
    background_color: str,
    dual_structure_horizontal_crop: float = 0.0,
) -> CompositeResult:
    """Compose validated frames with one fixed, optional dual-view source crop."""

    canvas_width, canvas_height = _validate_canvas_size(canvas_size)
    if layout_name not in {"stacked", "side_by_side"}:
        raise ValueError("layout_name must be 'stacked' or 'side_by_side'")
    try:
        background_rgb = ImageColor.getrgb(background_color)
    except ValueError as exc:
        raise ValueError(f"Invalid composite background color: {background_color!r}") from exc
    if len(background_rgb) not in {3, 4}:
        raise ValueError("Composite background color must resolve to RGB or RGBA")
    background_rgb = tuple(background_rgb[:3])

    landscape = validate_png_sequence(
        landscape_dir,
        expected_count=expected_count,
        label="landscape",
    )
    structure = validate_png_sequence(
        structure_dir,
        expected_count=expected_count,
        label="primary structure" if secondary_structure_dir is not None else "structure",
    )
    secondary_structure = (
        validate_png_sequence(
            secondary_structure_dir,
            expected_count=expected_count,
            label="secondary structure",
        )
        if secondary_structure_dir is not None
        else None
    )
    layout = resolve_composite_layout(
        canvas_size=(canvas_width, canvas_height),
        landscape_dimensions=landscape.dimensions,
        structure_dimensions=structure.dimensions,
        secondary_structure_dimensions=(
            secondary_structure.dimensions
            if secondary_structure is not None
            else None
        ),
        layout_name=layout_name,
        background_color=background_color,
        dual_structure_horizontal_crop=dual_structure_horizontal_crop,
    )
    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    outputs: list[Path] = []
    for index, (landscape_path, structure_path) in enumerate(
        zip(landscape.paths, structure.paths)
    ):
        canvas = Image.new("RGB", (canvas_width, canvas_height), background_rgb)
        _paste_contained(canvas, landscape_path, layout.landscape_destination)
        _paste_contained(
            canvas,
            structure_path,
            layout.primary_structure_destination or layout.structure_destination,
            source_crop=layout.primary_structure_source_crop,
        )
        if secondary_structure is not None:
            _paste_contained(
                canvas,
                secondary_structure.paths[index],
                layout.secondary_structure_destination,
                source_crop=layout.secondary_structure_source_crop,
            )
        output = destination / f"frame_{index:06d}.png"
        canvas.save(output, format="PNG")
        outputs.append(output)
    composite = validate_png_sequence(
        destination,
        expected_count=expected_count,
        label="composite",
    )
    if composite.dimensions != (canvas_width, canvas_height):
        raise RuntimeError(
            "Composite frame dimensions do not match requested canvas: "
            f"{composite.dimensions} != {(canvas_width, canvas_height)}"
        )
    return CompositeResult(
        frame_count=composite.frame_count,
        dimensions=composite.dimensions,
        paths=composite.paths,
        layout=layout,
        landscape_validation=landscape,
        structure_validation=structure,
        secondary_structure_validation=secondary_structure,
    )


def resolve_composite_layout(
    *,
    canvas_size: tuple[int, int],
    landscape_dimensions: tuple[int, int],
    structure_dimensions: tuple[int, int],
    secondary_structure_dimensions: tuple[int, int] | None = None,
    layout_name: str,
    background_color: str,
    padding_pixels: int = 24,
    gap_pixels: int = 24,
    dual_structure_horizontal_crop: float = 0.0,
) -> CompositeLayout:
    """Resolve fixed source crops and aspect-preserving destination rectangles."""

    width, height = _validate_canvas_size(canvas_size)
    if layout_name not in {"stacked", "side_by_side"}:
        raise ValueError("layout_name must be 'stacked' or 'side_by_side'")
    if secondary_structure_dimensions is not None and layout_name != "stacked":
        raise ValueError("dual structure view is supported only with stacked layout")
    crop_fraction = validate_dual_structure_horizontal_crop(
        dual_structure_horizontal_crop,
        has_secondary=secondary_structure_dimensions is not None,
        layout_name=layout_name,
    )
    if padding_pixels < 0 or gap_pixels < 0:
        raise ValueError("Composite padding and gap must be non-negative")
    if layout_name == "stacked":
        inner_width = width - 2 * padding_pixels
        inner_height = height - 2 * padding_pixels - gap_pixels
        if inner_width < 1 or inner_height < 2:
            raise ValueError("Composite canvas is too small for configured padding and gap")
        structure_fraction = STACKED_STRUCTURE_FRACTION
        landscape_fraction = STACKED_LANDSCAPE_FRACTION
        structure_height = min(
            inner_height - 1,
            max(1, int(round(inner_height * structure_fraction))),
        )
        landscape_height = inner_height - structure_height
        structure_region = (
            padding_pixels,
            padding_pixels,
            inner_width,
            structure_height,
        )
        landscape_region = (
            padding_pixels,
            padding_pixels + structure_height + gap_pixels,
            inner_width,
            landscape_height,
        )
    else:
        inner_width = width - 2 * padding_pixels - gap_pixels
        inner_height = height - 2 * padding_pixels
        if inner_width < 2 or inner_height < 1:
            raise ValueError("Composite canvas is too small for configured padding and gap")
        structure_fraction = None
        landscape_fraction = None
        structure_width = min(inner_height, inner_width // 2)
        landscape_width = inner_width - structure_width
        landscape_region = (
            padding_pixels,
            padding_pixels,
            landscape_width,
            inner_height,
        )
        structure_region = (
            padding_pixels + landscape_width + gap_pixels,
            padding_pixels,
            structure_width,
            inner_height,
        )
    if secondary_structure_dimensions is None:
        primary_structure_region = structure_region
        secondary_structure_region = None
        dual_view_gap_pixels = 0
    else:
        dual_view_gap_pixels = _dual_view_gap(width, structure_region[2])
        view_width = (structure_region[2] - dual_view_gap_pixels) // 2
        primary_structure_region = (
            structure_region[0],
            structure_region[1],
            view_width,
            structure_region[3],
        )
        secondary_structure_region = (
            structure_region[0] + view_width + dual_view_gap_pixels,
            structure_region[1],
            view_width,
            structure_region[3],
        )
    primary_crop, primary_cropped_dimensions = _horizontal_crop_rectangle(
        structure_dimensions,
        crop_fraction,
    )
    secondary_crop_and_dimensions = (
        _horizontal_crop_rectangle(
            secondary_structure_dimensions,
            crop_fraction,
        )
        if secondary_structure_dimensions is not None
        else None
    )
    primary_alignment = "right" if crop_fraction > 0 else "center"
    secondary_alignment = (
        "left" if crop_fraction > 0 else "center"
        if secondary_structure_dimensions is not None
        else None
    )
    primary_destination = _contain_rectangle(
        primary_cropped_dimensions,
        primary_structure_region,
        horizontal_alignment=primary_alignment,
    )
    secondary_destination = (
        _contain_rectangle(
            secondary_crop_and_dimensions[1],
            secondary_structure_region,
            horizontal_alignment=secondary_alignment,
        )
        if secondary_crop_and_dimensions is not None
        else None
    )
    return CompositeLayout(
        name=layout_name,
        canvas_size=(width, height),
        landscape_region=landscape_region,
        structure_region=structure_region,
        landscape_destination=_contain_rectangle(
            landscape_dimensions,
            landscape_region,
        ),
        structure_destination=primary_destination,
        background_color=background_color,
        primary_structure_region=primary_structure_region,
        secondary_structure_region=secondary_structure_region,
        primary_structure_destination=primary_destination,
        secondary_structure_destination=secondary_destination,
        primary_structure_source_crop=primary_crop,
        secondary_structure_source_crop=(
            secondary_crop_and_dimensions[0]
            if secondary_crop_and_dimensions is not None
            else None
        ),
        primary_structure_cropped_dimensions=primary_cropped_dimensions,
        secondary_structure_cropped_dimensions=(
            secondary_crop_and_dimensions[1]
            if secondary_crop_and_dimensions is not None
            else None
        ),
        primary_structure_alignment=primary_alignment,
        secondary_structure_alignment=secondary_alignment,
        dual_structure_horizontal_crop_fraction=crop_fraction,
        structure_view_count=2 if secondary_structure_dimensions is not None else 1,
        dual_view_gap_pixels=dual_view_gap_pixels,
        aspect_policy=(
            "fixed_source_crop_then_contain_no_stretch"
            if crop_fraction > 0
            else "contain_letterbox_no_crop"
        ),
        structure_fraction=structure_fraction,
        landscape_fraction=landscape_fraction,
        padding_pixels=padding_pixels,
        gap_pixels=gap_pixels,
    )


def validate_dual_structure_horizontal_crop(
    fraction: float,
    *,
    has_secondary: bool,
    layout_name: str,
) -> float:
    """Validate the presentation-only fixed crop policy."""

    try:
        value = float(fraction)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            "dual structure horizontal crop must be finite with 0 <= FRACTION < 0.5"
        ) from exc
    if not math.isfinite(value) or not 0 <= value < 0.5:
        raise ValueError(
            "dual structure horizontal crop must be finite with 0 <= FRACTION < 0.5"
        )
    if value > 0 and (not has_secondary or layout_name != "stacked"):
        raise ValueError(
            "non-zero dual structure horizontal crop is allowed only for "
            "dual-view stacked layout"
        )
    return value


def _horizontal_crop_rectangle(
    source_dimensions: tuple[int, int],
    fraction: float,
) -> tuple[tuple[int, int, int, int], tuple[int, int]]:
    """Return a Pillow crop rectangle using deterministic half-up rounding."""

    source_width, source_height = map(int, source_dimensions)
    if source_width <= 0 or source_height <= 0:
        raise ValueError("Source frame dimensions must be positive")
    crop_pixels = int(math.floor(source_width * fraction + 0.5))
    cropped_width = source_width - 2 * crop_pixels
    if cropped_width < 1:
        raise ValueError("dual structure horizontal crop leaves no source pixels")
    return (
        (crop_pixels, 0, source_width - crop_pixels, source_height),
        (cropped_width, source_height),
    )


def _dual_view_gap(canvas_width: int, region_width: int) -> int:
    """Resolve a small gap while keeping both dual-view regions exactly equal."""

    maximum = max(1, int(canvas_width * 0.015))
    gap = min(maximum, max(1, int(round(canvas_width * DUAL_VIEW_GAP_FRACTION))))
    if (region_width - gap) % 2:
        gap = gap + 1 if gap < maximum else gap - 1
    if gap < 0 or region_width - gap < 2:
        raise ValueError("Composite structure region is too small for dual view")
    return gap


def _validate_canvas_size(canvas_size: tuple[int, int]) -> tuple[int, int]:
    raw_values = tuple(canvas_size)
    if len(raw_values) != 2 or any(
        isinstance(value, bool) or not isinstance(value, int)
        for value in raw_values
    ):
        raise ValueError("Composite canvas dimensions must be positive even integers")
    values = tuple(raw_values)
    if any(value <= 0 or value % 2 for value in values):
        raise ValueError("Composite canvas dimensions must be positive even integers")
    return values


def _contain_rectangle(
    source_dimensions: tuple[int, int],
    region: tuple[int, int, int, int],
    *,
    horizontal_alignment: str = "center",
) -> tuple[int, int, int, int]:
    source_width, source_height = map(int, source_dimensions)
    region_x, region_y, region_width, region_height = region
    if source_width <= 0 or source_height <= 0:
        raise ValueError("Source frame dimensions must be positive")
    scale = min(region_width / source_width, region_height / source_height)
    width = max(1, int(source_width * scale))
    height = max(1, int(source_height * scale))
    if horizontal_alignment == "left":
        x = region_x
    elif horizontal_alignment == "right":
        x = region_x + region_width - width
    elif horizontal_alignment == "center":
        x = region_x + (region_width - width) // 2
    else:
        raise ValueError("horizontal_alignment must be left, center, or right")
    y = region_y + (region_height - height) // 2
    return x, y, width, height


def _paste_contained(
    canvas: Image.Image,
    source_path: Path,
    destination: tuple[int, int, int, int],
    *,
    source_crop: tuple[int, int, int, int] | None = None,
) -> None:
    x, y, width, height = destination
    with Image.open(source_path) as source:
        cropped = source.crop(source_crop) if source_crop is not None else source
        converted = cropped.convert("RGB")
        resized = converted.resize((width, height), resample=Image.Resampling.LANCZOS)
        canvas.paste(resized, (x, y))
