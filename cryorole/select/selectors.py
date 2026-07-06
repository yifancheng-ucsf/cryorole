"""Policy-driven selectors for orientation Landscapes."""

from __future__ import annotations

import math

import numpy as np
import pandas as pd
from scipy.spatial.transform import Rotation

from cryorole.core.euler_conventions import resolve_euler_convention
from cryorole.models.landscape import Landscape
from cryorole.models.policies import SelectionPolicy
from cryorole.models.selection import Selection


def select_particles(
    landscape: Landscape,
    *,
    policy: SelectionPolicy | None = None,
    source_metadata=None,
) -> Selection:
    """Select particle keys from a Landscape without modifying it."""

    policy = policy or SelectionPolicy()
    policy = _normalize_euler_policy(policy)
    _validate_policy(policy)
    if policy.selection_mode == "top_fraction_by_density":
        selected_keys, artifact_metadata = _select_top_fraction_by_density(landscape, policy)
        selection_basis = f"density:{policy.density_support_field}"
        center_evaluated = None
        resolved_evaluation_space = None
        resolved_range_source = None
    elif policy.selection_mode == "threshold_by_density":
        selected_keys = _select_threshold_by_density(landscape, policy)
        artifact_metadata = _density_artifact_metadata(
            len(landscape.data),
            policy=policy,
            applied=False,
        )
        selection_basis = f"density:{policy.density_support_field}"
        center_evaluated = None
        resolved_evaluation_space = None
        resolved_range_source = None
    elif policy.selection_mode == "radius_around_center":
        center_evaluated = _evaluate_center_input(policy)
        coordinate_column, resolved_evaluation_space = _coordinate_column_for_space(
            landscape,
            policy.evaluation_space,
        )
        selected_keys = _select_radius_around_center(
            landscape,
            policy,
            coordinate_column=coordinate_column,
            center_evaluated=center_evaluated,
        )
        artifact_metadata = _density_artifact_metadata(
            len(landscape.data),
            policy=policy,
            applied=False,
        )
        selection_basis = f"coordinates_{resolved_evaluation_space}"
        resolved_range_source = None
    elif policy.selection_mode == "range_by_coordinates":
        coordinate_column, resolved_range_source = _coordinate_column_for_space(
            landscape,
            policy.range_coordinate_source,
        )
        selected_keys = _select_range_by_coordinates(
            landscape,
            policy,
            coordinate_column=coordinate_column,
        )
        if policy.range_representation == "euler":
            selection_basis = f"display_coordinates:{resolved_range_source}:euler"
        else:
            selection_basis = f"coordinates_{resolved_range_source}:rotvec"
        center_evaluated = None
        resolved_evaluation_space = None
        artifact_metadata = _density_artifact_metadata(
            len(landscape.data),
            policy=policy,
            applied=False,
        )
    elif policy.selection_mode == "random":
        selected_keys = _select_random_fraction(landscape, policy)
        artifact_metadata = _density_artifact_metadata(
            len(landscape.data),
            policy=policy,
            applied=False,
        )
        selection_basis = "random_fraction"
        center_evaluated = None
        resolved_evaluation_space = None
        resolved_range_source = None
        metadata_metadata = _empty_metadata_metadata()
    elif policy.selection_mode == "metadata_value":
        selected_keys, metadata_metadata = _select_by_source_metadata(
            landscape,
            policy,
            source_metadata=source_metadata,
        )
        artifact_metadata = _density_artifact_metadata(
            len(landscape.data),
            policy=policy,
            applied=False,
        )
        selection_basis = f"source_metadata:{policy.metadata_domain}:{policy.metadata_column}"
        center_evaluated = None
        resolved_evaluation_space = None
        resolved_range_source = None
    else:
        raise ValueError(f"Unsupported selection_mode: {policy.selection_mode}")

    if policy.selection_mode != "metadata_value":
        metadata_metadata = _empty_metadata_metadata()

    return Selection(
        selection_id=policy.selection_id or _default_selection_id(policy),
        parent_landscape_id=policy.parent_landscape_id,
        parent_landscape_metadata=dict(policy.parent_landscape_metadata),
        selection_mode=policy.selection_mode,
        selection_basis=selection_basis,
        metric=(
            "coordinate_box"
            if policy.selection_mode == "range_by_coordinates"
            else (
                "random_without_replacement"
                if policy.selection_mode == "random"
                else (
                    "metadata_value_match"
                    if policy.selection_mode == "metadata_value"
                    else policy.metric
                )
            )
        ),
        center_input=policy.center_input,
        center_input_representation=(
            policy.center_input_representation
            if policy.selection_mode == "radius_around_center"
            else None
        ),
        center_input_space=(
            policy.center_input_space if policy.selection_mode == "radius_around_center" else None
        ),
        center_euler_sequence=(
            policy.center_euler_sequence
            if policy.selection_mode == "radius_around_center"
            and policy.center_input_representation == "euler"
            else None
        ),
        center_euler_convention=(
            policy.center_euler_convention
            if policy.selection_mode == "radius_around_center"
            and policy.center_input_representation == "euler"
            else None
        ),
        center_scipy_euler_sequence=(
            policy.center_scipy_euler_sequence
            if policy.selection_mode == "radius_around_center"
            and policy.center_input_representation == "euler"
            else None
        ),
        center_degrees=(
            policy.center_degrees
            if policy.selection_mode == "radius_around_center"
            and policy.center_input_representation == "euler"
            else None
        ),
        center_evaluated=center_evaluated,
        evaluation_space=(
            policy.evaluation_space if policy.selection_mode == "radius_around_center" else None
        ),
        resolved_evaluation_space=resolved_evaluation_space,
        radius=policy.radius,
        radius_unit=policy.radius_unit,
        threshold=policy.threshold if policy.selection_mode == "threshold_by_density" else None,
        sld_min=policy.sld_min if policy.selection_mode == "threshold_by_density" else None,
        sld_max=policy.sld_max if policy.selection_mode == "threshold_by_density" else None,
        threshold_operator=(
            policy.threshold_operator if policy.selection_mode == "threshold_by_density" else None
        ),
        density_support_field=policy.density_support_field,
        density_artifact_policy=policy.density_artifact_policy,
        density_artifact_flag_field=artifact_metadata["flag_field"],
        density_artifact_candidate_count_before=artifact_metadata["candidate_count_before"],
        density_artifact_candidate_count_after=artifact_metadata["candidate_count_after"],
        density_artifact_excluded_count=artifact_metadata["excluded_count"],
        top_fraction=(
            policy.top_fraction if policy.selection_mode == "top_fraction_by_density" else None
        ),
        random_fraction=(
            policy.random_fraction if policy.selection_mode == "random" else None
        ),
        random_seed=policy.random_seed if policy.selection_mode == "random" else None,
        random_candidate_count=(
            len(landscape.data) if policy.selection_mode == "random" else None
        ),
        metadata_domain=metadata_metadata["metadata_domain"],
        metadata_source_file=metadata_metadata["metadata_source_file"],
        metadata_column=metadata_metadata["metadata_column"],
        metadata_values=metadata_metadata["metadata_values"],
        metadata_source_row_id_field=metadata_metadata["metadata_source_row_id_field"],
        metadata_candidate_count=metadata_metadata["metadata_candidate_count"],
        metadata_missing_count=metadata_metadata["metadata_missing_count"],
        metadata_invalid_count=metadata_metadata["metadata_invalid_count"],
        tie_break_rule=(
            policy.tie_break_rule if policy.selection_mode == "top_fraction_by_density" else None
        ),
        range_coordinate_source=(
            policy.range_coordinate_source
            if policy.selection_mode == "range_by_coordinates"
            else None
        ),
        resolved_range_coordinate_source=resolved_range_source,
        range_representation=(
            policy.range_representation if policy.selection_mode == "range_by_coordinates" else None
        ),
        range_euler_sequence=(
            policy.range_euler_sequence
            if policy.selection_mode == "range_by_coordinates"
            and policy.range_representation == "euler"
            else None
        ),
        range_euler_convention=(
            policy.range_euler_convention
            if policy.selection_mode == "range_by_coordinates"
            and policy.range_representation == "euler"
            else None
        ),
        range_scipy_euler_sequence=(
            policy.range_scipy_euler_sequence
            if policy.selection_mode == "range_by_coordinates"
            and policy.range_representation == "euler"
            else None
        ),
        range_degrees=(
            policy.range_degrees
            if policy.selection_mode == "range_by_coordinates"
            and policy.range_representation == "euler"
            else None
        ),
        range_bounds=(
            dict(policy.range_bounds)
            if policy.selection_mode == "range_by_coordinates"
            else None
        ),
        selected_particle_keys=selected_keys,
        selected_count=len(selected_keys),
        total_count=len(landscape.data),
        active_policy=policy,
    )


def _normalize_euler_policy(policy: SelectionPolicy) -> SelectionPolicy:
    center = _resolve_selection_euler_policy(
        policy.center_euler_convention,
        scipy_euler_sequence=policy.center_scipy_euler_sequence,
        legacy_euler_sequence=policy.center_euler_sequence,
    )
    range_ = _resolve_selection_euler_policy(
        policy.range_euler_convention,
        scipy_euler_sequence=policy.range_scipy_euler_sequence,
        legacy_euler_sequence=policy.range_euler_sequence,
    )
    return SelectionPolicy(
        selection_mode=policy.selection_mode,
        density_support_field=policy.density_support_field,
        density_artifact_policy=policy.density_artifact_policy,
        top_fraction=policy.top_fraction,
        random_fraction=policy.random_fraction,
        random_seed=policy.random_seed,
        metadata_domain=policy.metadata_domain,
        metadata_column=policy.metadata_column,
        metadata_values=tuple(policy.metadata_values),
        metadata_source_file=policy.metadata_source_file,
        metadata_source_row_id_field=policy.metadata_source_row_id_field,
        split_by_metadata=policy.split_by_metadata,
        threshold=policy.threshold,
        sld_min=policy.sld_min,
        sld_max=policy.sld_max,
        threshold_operator=policy.threshold_operator,
        tie_break_rule=policy.tie_break_rule,
        center_input=policy.center_input,
        center_input_representation=policy.center_input_representation,
        center_input_space=policy.center_input_space,
        center_euler_sequence=center.scipy_euler_sequence,
        center_euler_convention=center.euler_convention,
        center_scipy_euler_sequence=center.scipy_euler_sequence,
        center_degrees=policy.center_degrees,
        evaluation_space=policy.evaluation_space,
        metric=policy.metric,
        radius=policy.radius,
        radius_unit=policy.radius_unit,
        range_coordinate_source=policy.range_coordinate_source,
        range_representation=policy.range_representation,
        range_euler_sequence=range_.scipy_euler_sequence,
        range_euler_convention=range_.euler_convention,
        range_scipy_euler_sequence=range_.scipy_euler_sequence,
        range_degrees=policy.range_degrees,
        range_bounds=policy.range_bounds,
        parent_landscape_id=policy.parent_landscape_id,
        parent_landscape_metadata=policy.parent_landscape_metadata,
        selection_id=policy.selection_id,
    )


def _resolve_selection_euler_policy(
    euler_convention: str,
    *,
    scipy_euler_sequence: str,
    legacy_euler_sequence: str,
):
    if legacy_euler_sequence != scipy_euler_sequence:
        return resolve_euler_convention(
            scipy_euler_sequence=legacy_euler_sequence,
            source="selection_policy_legacy_sequence",
        )
    return resolve_euler_convention(
        euler_convention,
        scipy_euler_sequence=scipy_euler_sequence,
        source="selection_policy",
    )


def _validate_policy(policy: SelectionPolicy) -> None:
    if policy.selection_mode not in {
        "top_fraction_by_density",
        "threshold_by_density",
        "radius_around_center",
        "range_by_coordinates",
        "random",
        "metadata_value",
        "metadata_group",
    }:
        raise ValueError(f"Unsupported selection_mode: {policy.selection_mode}")
    if policy.evaluation_space not in {"analysis", "canonical", "canonical_if_available"}:
        raise ValueError(f"Unsupported selection evaluation_space: {policy.evaluation_space}")
    if policy.center_input_space not in {"analysis_display", "canonical_display", "evaluation"}:
        raise ValueError(f"Unsupported selection center_input_space: {policy.center_input_space}")
    _resolve_selection_euler_policy(
        policy.center_euler_convention,
        scipy_euler_sequence=policy.center_scipy_euler_sequence,
        legacy_euler_sequence=policy.center_euler_sequence,
    )
    _resolve_selection_euler_policy(
        policy.range_euler_convention,
        scipy_euler_sequence=policy.range_scipy_euler_sequence,
        legacy_euler_sequence=policy.range_euler_sequence,
    )
    if policy.selection_mode != "range_by_coordinates" and policy.metric not in {
        "so3_geodesic",
        "rotvec_euclidean",
        "euclidean",
    }:
        raise ValueError(f"Unsupported selection metric: {policy.metric}")
    if policy.density_support_field == "sld_display":
        raise ValueError("sld_display is display-only and cannot drive selection")
    if policy.density_artifact_policy not in {"include_all", "exclude_display_outliers"}:
        raise ValueError(
            "density_artifact_policy must be include_all or exclude_display_outliers"
        )
    if policy.radius is not None and (not np.isfinite(policy.radius) or policy.radius <= 0.0):
        raise ValueError("Selection radius must be a positive finite value")
    if policy.selection_mode == "top_fraction_by_density":
        if policy.threshold is not None:
            raise ValueError("top_fraction_by_density selection must not use threshold")
        if not 0.0 < policy.top_fraction <= 1.0:
            raise ValueError("Selection top_fraction must be > 0 and <= 1")
        if policy.tie_break_rule != "stable_input_order":
            raise ValueError(f"Unsupported selection tie_break_rule: {policy.tie_break_rule}")
    if policy.selection_mode == "random":
        if policy.random_fraction is None:
            raise ValueError("random selection requires random_fraction")
        if not 0.0 < policy.random_fraction <= 1.0:
            raise ValueError("Selection random_fraction must be > 0 and <= 1")
        if policy.random_seed is not None and isinstance(policy.random_seed, bool):
            raise ValueError("Selection random_seed must be an integer or None")
    if policy.selection_mode in {"metadata_value", "metadata_group"}:
        if policy.metadata_domain not in {"ref", "mov"}:
            raise ValueError("metadata selection requires metadata_domain ref or mov")
        if not policy.metadata_column:
            raise ValueError("metadata selection requires metadata_column")
        expected_row_field = f"{policy.metadata_domain}_source_row_id"
        if (
            policy.metadata_source_row_id_field is not None
            and policy.metadata_source_row_id_field != expected_row_field
        ):
            raise ValueError(
                "metadata_source_row_id_field must match metadata_domain "
                f"({expected_row_field})"
            )
        if policy.selection_mode == "metadata_value" and not policy.metadata_values:
            raise ValueError("metadata_value selection requires metadata_values")
        if policy.selection_mode == "metadata_group" and not policy.split_by_metadata:
            raise ValueError("metadata_group selection requires split_by_metadata")
    if policy.selection_mode == "threshold_by_density":
        if policy.threshold is None and policy.sld_min is None and policy.sld_max is None:
            raise ValueError(
                "threshold_by_density selection requires threshold, sld_min, or sld_max"
            )
        for field_name, value in (
            ("threshold", policy.threshold),
            ("sld_min", policy.sld_min),
            ("sld_max", policy.sld_max),
        ):
            if value is not None and (not np.isfinite(value) or value < 0.0):
                raise ValueError(f"Selection {field_name} must be finite non-negative")
        if (
            policy.sld_min is not None
            and policy.sld_max is not None
            and policy.sld_min > policy.sld_max
        ):
            raise ValueError("Selection sld_min must be <= sld_max")
        if policy.threshold_operator != ">=":
            raise ValueError(f"Unsupported selection threshold_operator: {policy.threshold_operator}")
    if policy.selection_mode == "radius_around_center":
        if policy.center_input is None:
            raise ValueError("radius_around_center selection requires center_input")
        if policy.radius is None:
            raise ValueError("radius_around_center selection requires radius")
        _evaluate_center_input(policy)
    if policy.selection_mode == "range_by_coordinates":
        _validate_range_policy(policy)


def _select_top_fraction_by_density(
    landscape: Landscape,
    policy: SelectionPolicy,
) -> tuple[tuple[object, ...], dict[str, int | str | None]]:
    data = landscape.data
    if policy.density_support_field not in data.columns:
        raise ValueError(f"Landscape missing density support field: {policy.density_support_field}")
    support_values = np.asarray(data[policy.density_support_field], dtype=float)
    if not np.isfinite(support_values).all():
        raise ValueError(f"{policy.density_support_field} contains non-finite values")

    candidate_indices, artifact_metadata = _density_artifact_candidate_indices(data, policy)
    if candidate_indices.size == 0:
        return (), artifact_metadata
    n_select = max(1, int(math.ceil(candidate_indices.size * policy.top_fraction)))
    n_select = min(n_select, candidate_indices.size)
    row_indices = candidate_indices
    ranked_within_candidates = np.lexsort((row_indices, -support_values[candidate_indices]))
    ranked_indices = candidate_indices[ranked_within_candidates]
    selected_indices = set(ranked_indices[:n_select])
    return tuple(
        particle_key
        for row_index, particle_key in enumerate(data["particle_key"])
        if row_index in selected_indices
    ), artifact_metadata


def _density_artifact_candidate_indices(
    data,
    policy: SelectionPolicy,
) -> tuple[np.ndarray, dict[str, int | str | None]]:
    all_indices = np.arange(len(data), dtype=int)
    if policy.density_artifact_policy == "include_all":
        return all_indices, _density_artifact_metadata(len(data), policy=policy, applied=False)
    flag_field = "sld_display_is_outlier"
    if flag_field not in data.columns:
        raise ValueError(
            "exclude_display_outliers requires landscape field sld_display_is_outlier"
        )
    clipped = np.asarray(data[flag_field], dtype=bool)
    candidate_indices = np.where(~clipped)[0]
    return candidate_indices, {
        "flag_field": flag_field,
        "candidate_count_before": int(len(data)),
        "candidate_count_after": int(candidate_indices.size),
        "excluded_count": int(clipped.sum()),
    }


def _density_artifact_metadata(
    n_rows: int,
    *,
    policy: SelectionPolicy,
    applied: bool,
) -> dict[str, int | str | None]:
    if policy.density_artifact_policy == "exclude_display_outliers" or applied:
        return {
            "flag_field": "sld_display_is_outlier",
            "candidate_count_before": int(n_rows),
            "candidate_count_after": int(n_rows),
            "excluded_count": 0,
        }
    return {
        "flag_field": None,
        "candidate_count_before": int(n_rows),
        "candidate_count_after": int(n_rows),
        "excluded_count": 0,
    }


def _select_threshold_by_density(
    landscape: Landscape,
    policy: SelectionPolicy,
) -> tuple[object, ...]:
    data = landscape.data
    if policy.density_support_field not in data.columns:
        raise ValueError(f"Landscape missing density support field: {policy.density_support_field}")
    support_values = np.asarray(data[policy.density_support_field], dtype=float)
    if not np.isfinite(support_values).all():
        raise ValueError(f"{policy.density_support_field} contains non-finite values")
    selected = np.ones(len(data), dtype=bool)
    lower = policy.sld_min if policy.sld_min is not None else policy.threshold
    upper = policy.sld_max
    if lower is not None:
        selected &= support_values >= float(lower)
    if upper is not None:
        selected &= support_values <= float(upper)
    return tuple(data.loc[selected, "particle_key"])


def _select_random_fraction(
    landscape: Landscape,
    policy: SelectionPolicy,
) -> tuple[object, ...]:
    data = landscape.data
    n_rows = len(data)
    if n_rows == 0:
        return ()
    n_select = max(1, int(math.ceil(n_rows * float(policy.random_fraction))))
    n_select = min(n_select, n_rows)
    rng = np.random.default_rng(policy.random_seed)
    selected_indices = set(rng.choice(n_rows, size=n_select, replace=False).tolist())
    return tuple(
        particle_key
        for row_index, particle_key in enumerate(data["particle_key"])
        if row_index in selected_indices
    )


def _select_by_source_metadata(
    landscape: Landscape,
    policy: SelectionPolicy,
    *,
    source_metadata,
) -> tuple[tuple[object, ...], dict[str, object]]:
    if source_metadata is None:
        raise ValueError("metadata selection requires source_metadata")
    source_row_id_field = policy.metadata_source_row_id_field or (
        f"{policy.metadata_domain}_source_row_id"
    )
    if source_row_id_field not in landscape.data.columns:
        raise ValueError(f"Landscape missing metadata source-row field: {source_row_id_field}")
    values_by_row_id = _metadata_values_by_source_row_id(
        source_metadata,
        metadata_column=str(policy.metadata_column),
    )
    requested_values = tuple(_metadata_value_key(value) for value in policy.metadata_values)
    requested_set = set(requested_values)
    selected_keys: list[object] = []
    missing_count = 0
    invalid_count = 0
    for particle_key, source_row_id in zip(
        landscape.data["particle_key"],
        landscape.data[source_row_id_field],
    ):
        row_id = _coerce_source_row_id(source_row_id)
        if row_id is None or row_id < 0:
            invalid_count += 1
            missing_count += 1
            continue
        if row_id not in values_by_row_id:
            missing_count += 1
            continue
        metadata_value = values_by_row_id[row_id]
        if metadata_value is None:
            missing_count += 1
            continue
        if metadata_value in requested_set:
            selected_keys.append(particle_key)
    metadata = {
        "metadata_domain": policy.metadata_domain,
        "metadata_source_file": policy.metadata_source_file,
        "metadata_column": policy.metadata_column,
        "metadata_values": requested_values,
        "metadata_source_row_id_field": source_row_id_field,
        "metadata_candidate_count": int(len(landscape.data)),
        "metadata_missing_count": int(missing_count),
        "metadata_invalid_count": int(invalid_count),
    }
    return tuple(selected_keys), metadata


def _metadata_values_by_source_row_id(source_metadata, *, metadata_column: str) -> dict[int, str | None]:
    if isinstance(source_metadata, pd.DataFrame):
        column = _resolve_metadata_column_name(source_metadata.columns, metadata_column)
        return {
            int(index): _metadata_value_key(value)
            for index, value in source_metadata[column].items()
        }
    if isinstance(source_metadata, dict):
        return {
            int(key): _metadata_value_key(value)
            for key, value in source_metadata.items()
        }
    raise ValueError("source_metadata must be a pandas DataFrame or row-id mapping")


def _resolve_metadata_column_name(columns, requested: str) -> str:
    if requested in columns:
        return requested
    prefixed = requested if requested.startswith("_") else f"_{requested}"
    if prefixed in columns:
        return prefixed
    unprefixed = requested[1:] if requested.startswith("_") else requested
    if unprefixed in columns:
        return unprefixed
    raise ValueError(f"Source metadata missing column: {requested}")


def _coerce_source_row_id(value) -> int | None:
    if isinstance(value, (bool, np.bool_)):
        return None
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return None
    if not np.isfinite(numeric) or not numeric.is_integer():
        return None
    return int(numeric)


def _metadata_value_key(value) -> str | None:
    if value is None or pd.isna(value):
        return None
    text = str(value).strip()
    try:
        numeric = float(text)
    except ValueError:
        return text
    if np.isfinite(numeric) and numeric.is_integer():
        return str(int(numeric))
    return text


def _empty_metadata_metadata() -> dict[str, object]:
    return {
        "metadata_domain": None,
        "metadata_source_file": None,
        "metadata_column": None,
        "metadata_values": (),
        "metadata_source_row_id_field": None,
        "metadata_candidate_count": None,
        "metadata_missing_count": None,
        "metadata_invalid_count": None,
    }


def _select_radius_around_center(
    landscape: Landscape,
    policy: SelectionPolicy,
    *,
    coordinate_column,
    center_evaluated: tuple[float, float, float],
) -> tuple[object, ...]:
    coordinates = np.vstack([np.asarray(coords, dtype=float) for coords in coordinate_column])
    center = np.asarray(center_evaluated, dtype=float)
    if policy.metric == "so3_geodesic":
        distances = _so3_geodesic_distances(coordinates, center)
    else:
        distances = np.linalg.norm(coordinates - center, axis=1)
    selected = distances <= _radius_threshold_in_radians(policy)
    return tuple(landscape.data.loc[selected, "particle_key"])


def _radius_threshold_in_radians(policy: SelectionPolicy) -> float:
    radius = float(policy.radius)
    if policy.radius_unit == "degrees":
        return math.radians(radius)
    return radius


def _so3_geodesic_distances(
    rotvec_coordinates: np.ndarray,
    center_rotvec: np.ndarray,
) -> np.ndarray:
    """Return SO(3) geodesic distances in radians from center to each rotvec."""

    center_rotation = Rotation.from_rotvec(center_rotvec)
    point_rotations = Rotation.from_rotvec(rotvec_coordinates)
    relative = center_rotation.inv() * point_rotations
    return relative.magnitude()


def _select_range_by_coordinates(
    landscape: Landscape,
    policy: SelectionPolicy,
    *,
    coordinate_column,
) -> tuple[object, ...]:
    rotvec_coordinates = np.vstack([np.asarray(coords, dtype=float) for coords in coordinate_column])
    if policy.range_representation == "euler":
        coordinates = Rotation.from_rotvec(rotvec_coordinates).as_euler(
            policy.range_scipy_euler_sequence,
            degrees=policy.range_degrees,
        )
        axis_names = ("alpha", "beta", "gamma")
    elif policy.range_representation == "rotvec":
        coordinates = rotvec_coordinates
        axis_names = ("x", "y", "z")
    else:
        raise ValueError(f"Unsupported range_representation: {policy.range_representation}")

    selected = np.ones(len(landscape.data), dtype=bool)
    for axis_index, axis_name in enumerate(axis_names):
        bounds = policy.range_bounds.get(axis_name)
        if bounds is None:
            continue
        selected &= _range_mask(
            coordinates[:, axis_index],
            bounds,
            wraparound=policy.range_representation == "euler",
        )
    return tuple(landscape.data.loc[selected, "particle_key"])


def _coordinate_column_for_space(landscape: Landscape, evaluation_space: str):
    if evaluation_space == "analysis":
        return landscape.data["coordinates_analysis"], "analysis"
    if evaluation_space == "canonical":
        if "coordinates_canonical" not in landscape.data.columns:
            raise ValueError("canonical coordinate selection requires coordinates_canonical")
        return landscape.data["coordinates_canonical"], "canonical"
    if evaluation_space == "canonical_if_available":
        if "coordinates_canonical" in landscape.data.columns:
            return landscape.data["coordinates_canonical"], "canonical"
        return landscape.data["coordinates_analysis"], "analysis"
    raise ValueError(f"Unsupported selection evaluation_space: {evaluation_space}")


def _evaluate_center_input(policy: SelectionPolicy) -> tuple[float, float, float]:
    center_input = np.asarray(policy.center_input, dtype=float)
    if policy.center_input_representation == "rotvec":
        if center_input.shape != (3,) or not np.isfinite(center_input).all():
            raise ValueError("rotvec center_input must be a finite length-3 vector")
        return tuple(float(value) for value in center_input)
    if policy.center_input_representation == "euler":
        if center_input.shape != (3,) or not np.isfinite(center_input).all():
            raise ValueError("euler center_input must be a finite length-3 vector")
        rotvec = Rotation.from_euler(
            policy.center_scipy_euler_sequence,
            center_input,
            degrees=policy.center_degrees,
        ).as_rotvec()
        return tuple(float(value) for value in rotvec)
    raise ValueError(
        f"Unsupported center_input_representation: {policy.center_input_representation}"
    )


def _validate_range_policy(policy: SelectionPolicy) -> None:
    if policy.range_coordinate_source not in {"analysis", "canonical", "canonical_if_available"}:
        raise ValueError(
            f"Unsupported selection range_coordinate_source: {policy.range_coordinate_source}"
        )
    if policy.range_representation not in {"euler", "rotvec"}:
        raise ValueError(f"Unsupported range_representation: {policy.range_representation}")
    allowed_axes = (
        {"alpha", "beta", "gamma"} if policy.range_representation == "euler" else {"x", "y", "z"}
    )
    unknown_axes = set(policy.range_bounds) - allowed_axes
    if unknown_axes:
        raise ValueError(f"Unsupported range axis names: {sorted(unknown_axes)}")
    has_constrained_axis = False
    for axis_name, bounds in policy.range_bounds.items():
        if bounds is None:
            continue
        if len(bounds) != 2:
            raise ValueError(f"range_bounds for {axis_name} must be a (lower, upper) pair")
        lower, upper = bounds
        if lower is not None or upper is not None:
            has_constrained_axis = True
        if lower is not None and not np.isfinite(lower):
            raise ValueError(f"range lower bound for {axis_name} must be finite or None")
        if upper is not None and not np.isfinite(upper):
            raise ValueError(f"range upper bound for {axis_name} must be finite or None")
        if (
            policy.range_representation == "rotvec"
            and lower is not None
            and upper is not None
            and lower > upper
        ):
            raise ValueError("rotvec range bounds must have lower <= upper")
    if not has_constrained_axis:
        raise ValueError("range_by_coordinates requires at least one constrained axis")


def _range_mask(
    values: np.ndarray,
    bounds: tuple[float | None, float | None],
    *,
    wraparound: bool,
) -> np.ndarray:
    lower, upper = bounds
    if lower is None and upper is None:
        return np.ones(values.shape, dtype=bool)
    if lower is None:
        return values <= float(upper)
    if upper is None:
        return values >= float(lower)
    if wraparound and lower > upper:
        return (values >= float(lower)) | (values <= float(upper))
    return (values >= float(lower)) & (values <= float(upper))


def _default_selection_id(policy: SelectionPolicy) -> str:
    if policy.selection_mode == "top_fraction_by_density":
        artifact_suffix = (
            ""
            if policy.density_artifact_policy == "include_all"
            else f":{policy.density_artifact_policy}"
        )
        return (
            f"selection:{policy.selection_mode}:{policy.density_support_field}:"
            f"{policy.top_fraction:g}{artifact_suffix}"
        )
    if policy.selection_mode == "threshold_by_density":
        return (
            f"selection:{policy.selection_mode}:{policy.density_support_field}:"
            f"{policy.threshold_operator}{policy.threshold:g}"
        )
    if policy.selection_mode == "range_by_coordinates":
        return (
            f"selection:{policy.selection_mode}:{policy.range_coordinate_source}:"
            f"{policy.range_representation}"
        )
    if policy.selection_mode == "random":
        seed = "none" if policy.random_seed is None else str(policy.random_seed)
        return f"selection:random:{policy.random_fraction:g}:seed{seed}"
    if policy.selection_mode == "metadata_value":
        values = "_".join(str(value) for value in policy.metadata_values) or "value"
        return (
            f"selection:metadata:{policy.metadata_domain}:"
            f"{policy.metadata_column}:{values}"
        )
    return f"selection:{policy.selection_mode}:{policy.evaluation_space}:{policy.radius:g}"
