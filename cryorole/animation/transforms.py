"""Pure numerical transforms for offline animation.

Matrices in this module act on column vectors. Canonical landscape RVs remain
row vectors at the artifact boundary: ``canonical_rv = raw_rv @ C``.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.spatial.transform import Rotation

from cryorole.canonicalize.transforms import validate_canonical_transform


def canonical_rotation_to_raw(
    canonical_rotation: Rotation,
    canonical_transform: np.ndarray,
) -> Rotation:
    """Map a canonical waypoint rotation back to the physical raw RO space.

    The canonical frame is a linear coordinate reparameterization of the
    rotation vector, not a direct substitution of canonical Euler values into
    raw Euler space.
    """

    transform = np.asarray(canonical_transform, dtype=float)
    validate_canonical_transform(transform)
    canonical_rv = np.asarray(canonical_rotation.as_rotvec(), dtype=float)
    raw_rv = canonical_rv @ transform.T
    if not np.isfinite(raw_rv).all():
        raise ValueError("canonical-to-raw mapping produced non-finite values")
    return Rotation.from_rotvec(raw_rv)


def ro_target_to_active_delta(
    target_ro: Rotation,
    baseline_ro: Rotation,
) -> Rotation:
    """Return the active density delta for a target passive RO assertion.

    ``delta_passive = RO_target @ inverse(RO_baseline)``, followed by
    ``delta_active = inverse(delta_passive)``.
    """

    delta_passive = target_ro * baseline_ro.inv()
    delta_active = delta_passive.inv()
    _validate_rotation_matrix(delta_active.as_matrix(), "active RO delta")
    return delta_active


def conjugate_active_delta(
    delta_active_raw: np.ndarray,
    raw_to_scene: np.ndarray,
) -> np.ndarray:
    """Change a raw-frame active delta into the ChimeraX scene frame.

    With column vectors ``x_scene = S @ x_raw``, the scene-space active
    rotation is ``S @ delta_active_raw @ inverse(S)``.
    """

    delta = _validate_rotation_matrix(delta_active_raw, "raw active delta")
    scene = _validate_rotation_matrix(raw_to_scene, "raw_to_scene")
    result = scene @ delta @ scene.T
    return _validate_rotation_matrix(result, "scene active delta")


def pivoted_scene_delta(
    delta_active_scene: np.ndarray,
    pivot_scene: np.ndarray,
) -> np.ndarray:
    """Return a 4x4 scene-space delta rotating around ``pivot_scene``.

    ``x_target = pivot + R @ (x_baseline - pivot)``. The result is intended
    for left-composition with a saved local-to-scene baseline transform.
    """

    rotation = _validate_rotation_matrix(delta_active_scene, "scene active delta")
    pivot = np.asarray(pivot_scene, dtype=float)
    if pivot.shape != (3,):
        raise ValueError(f"pivot must have shape (3,), got {pivot.shape}")
    if not np.isfinite(pivot).all():
        raise ValueError("pivot contains non-finite values")
    transform = np.eye(4, dtype=float)
    transform[:3, :3] = rotation
    transform[:3, 3] = pivot - rotation @ pivot
    if not np.isfinite(transform).all():
        raise ValueError("pivot transform contains non-finite values")
    return transform


def validate_scene_basis(matrix: np.ndarray, *, name: str = "raw_to_scene") -> np.ndarray:
    """Validate a physical right-handed orthonormal frame transform."""

    return _validate_rotation_matrix(matrix, name)


def parse_ro_assertion(
    value: str,
    *,
    scipy_euler_sequence: str,
    first_waypoint_raw: Rotation,
) -> tuple[Rotation, dict[str, object]]:
    """Parse the required saved-session baseline RO assertion."""

    text = str(value).strip()
    if text == "identity":
        return Rotation.identity(), {"kind": "identity", "value": text}
    if text == "first-waypoint":
        return first_waypoint_raw, {"kind": "first-waypoint", "value": text}
    if text.startswith("ea:"):
        coordinates = _parse_triplet(text[3:], "baseline ea")
        return (
            Rotation.from_euler(scipy_euler_sequence, coordinates, degrees=True),
            {"kind": "ea", "value": coordinates.tolist(), "units": "degrees"},
        )
    if text.startswith("rv:"):
        coordinates = _parse_triplet(text[3:], "baseline rv")
        return (
            Rotation.from_rotvec(coordinates),
            {"kind": "rv", "value": coordinates.tolist(), "units": "radians"},
        )
    raise ValueError(
        "baseline_ro must be identity, first-waypoint, ea:A,B,G, or rv:X,Y,Z"
    )


def load_explicit_scene_basis(path: str | Path) -> np.ndarray:
    """Load an audited raw-to-scene rotation matrix from JSON, NPY, or NPZ."""

    source = Path(path)
    if not source.is_file():
        raise ValueError(f"Map-frame transform does not exist: {source}")
    suffix = source.suffix.lower()
    if suffix == ".json":
        with source.open("r", encoding="utf-8") as handle:
            payload = json.load(handle)
        if isinstance(payload, dict):
            payload = payload.get("raw_to_scene", payload.get("transform"))
        matrix = np.asarray(payload, dtype=float)
    elif suffix == ".npy":
        matrix = np.asarray(np.load(source, allow_pickle=False), dtype=float)
    elif suffix == ".npz":
        with np.load(source, allow_pickle=False) as archive:
            key = "raw_to_scene" if "raw_to_scene" in archive else "transform"
            if key not in archive:
                raise ValueError("Map-frame NPZ requires raw_to_scene or transform")
            matrix = np.asarray(archive[key], dtype=float)
    else:
        raise ValueError("Map-frame transform must be JSON, NPY, or NPZ")
    return validate_scene_basis(matrix)


def _parse_triplet(value: str, name: str) -> np.ndarray:
    parts = value.split(",")
    if len(parts) != 3:
        raise ValueError(f"{name} must contain exactly three comma-separated values")
    try:
        result = np.asarray([float(part.strip()) for part in parts], dtype=float)
    except ValueError as exc:
        raise ValueError(f"{name} values must be numeric") from exc
    if not np.isfinite(result).all():
        raise ValueError(f"{name} values must be finite")
    return result


def _validate_rotation_matrix(matrix: np.ndarray, name: str) -> np.ndarray:
    value = np.asarray(matrix, dtype=float)
    if value.shape != (3, 3):
        raise ValueError(f"{name} must have shape (3, 3), got {value.shape}")
    if not np.isfinite(value).all():
        raise ValueError(f"{name} contains non-finite values")
    if not np.allclose(value.T @ value, np.eye(3), atol=1e-8):
        raise ValueError(f"{name} must be orthonormal")
    if not np.isclose(np.linalg.det(value), 1.0, atol=1e-8):
        raise ValueError(f"{name} must be right-handed with determinant +1")
    return value
