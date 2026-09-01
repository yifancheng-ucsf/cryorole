"""Create visualization-compatible landscapes after a controlled RO rotation."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
from typing import Any, Sequence

import numpy as np
from scipy.spatial.transform import Rotation

from cryorole.canonicalize.transforms import validate_canonical_transform
from cryorole.core.euler_conventions import (
    RAW_EULER_ANGLE_COLUMNS,
    resolve_euler_convention,
)
from cryorole.export import (
    read_landscape,
    resolve_landscape_path,
    write_canonical_landscape_csv,
    write_json_artifact,
    write_landscape_npz,
    write_raw_landscape_csv,
)
from cryorole.models.landscape import Landscape


ROTATION_OPERATION = "moving_domain_local_right_rotation"
RO_TRANSFORM = "R_ro_shifted = R_ro @ G_raw"
CANONICAL_TRANSFORM_DIRECTION = "canonical_rv = raw_rv @ canonical_transform"


def build_parser() -> argparse.ArgumentParser:
    """Build the standalone rotated-landscape diagnostic parser."""

    parser = argparse.ArgumentParser(
        prog="rotate_landscape",
        description=(
            "Apply one extrinsic-ZYX rotation to every RO in a raw or canonical "
            "landscape and write a bundle readable by 'cryorole visualize'."
        ),
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Parent cryoROLE run directory or raw/canonical landscape CSV/NPZ.",
    )
    parser.add_argument(
        "--space",
        choices=("raw", "canonical"),
        default="raw",
        help="Coordinate frame in which the input rotation is defined. Default: raw.",
    )
    parser.add_argument(
        "--rotation-euler",
        nargs=3,
        required=True,
        type=float,
        metavar=("ALPHA", "BETA", "GAMMA"),
        help="Extrinsic fixed-axis ZYX rotation in degrees.",
    )
    parser.add_argument("--output-dir", required=True)
    parser.add_argument(
        "--canonical-id",
        default="default",
        help="Canonical landscape/frame id for input resolution and output. Default: default.",
    )
    parser.add_argument(
        "--canonical-frame",
        help=(
            "Optional canonical_frame.json/.npz. Required only when canonical "
            "coordinates do not carry or determine their transform."
        ),
    )
    parser.add_argument("--overwrite", action="store_true")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Run the standalone rotated-landscape diagnostic."""

    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        rotate_landscape_bundle(
            args.input,
            args.output_dir,
            rotation_euler_deg=args.rotation_euler,
            space=args.space,
            canonical_id=args.canonical_id,
            canonical_frame=args.canonical_frame,
            overwrite=args.overwrite,
        )
    except (ValueError, FileExistsError, IsADirectoryError) as exc:
        parser.exit(2, f"rotate_landscape: error: {exc}\n")
    print(str(Path(args.output_dir).resolve()))
    return 0


def rotate_landscape_bundle(
    input_path: str | Path,
    output_dir: str | Path,
    *,
    rotation_euler_deg: Sequence[float],
    space: str = "raw",
    canonical_id: str = "default",
    canonical_frame: str | Path | None = None,
    overwrite: bool = False,
) -> dict[str, Any]:
    """Rotate a landscape and write a derived bundle consumable by visualize.

    The same rotation is right-multiplied onto every physical RO. SLD and
    particle/source-row provenance are inherited unchanged.
    """

    if space not in {"raw", "canonical"}:
        raise ValueError("space must be 'raw' or 'canonical'")
    if not canonical_id or not str(canonical_id).strip():
        raise ValueError("canonical_id must not be empty")
    angles = np.asarray(rotation_euler_deg, dtype=float)
    if angles.shape != (3,) or not np.isfinite(angles).all():
        raise ValueError("rotation_euler_deg must contain three finite values")

    requested_input = Path(input_path)
    if not requested_input.exists():
        raise ValueError(f"Input landscape or run directory does not exist: {requested_input}")
    source_path = resolve_landscape_path(
        requested_input,
        space=space,
        canonical_id=canonical_id,
    )
    source_path = source_path.resolve()
    parent_landscape = read_landscape(source_path)
    raw_coordinates = _stack_coordinates(
        parent_landscape,
        "coordinates_analysis",
    )
    if len(raw_coordinates) == 0:
        raise ValueError("Input landscape must contain at least one row")
    if not np.isfinite(raw_coordinates).all():
        raise ValueError("Input landscape raw rotation vectors contain non-finite values")

    canonical_coordinates = _optional_coordinates(
        parent_landscape,
        "coordinates_canonical",
    )
    if canonical_coordinates is not None and not np.isfinite(canonical_coordinates).all():
        raise ValueError("Input landscape canonical rotation vectors contain non-finite values")
    if space == "canonical" and canonical_coordinates is None:
        raise ValueError(
            "canonical rotation requires canonical coordinates in the input landscape"
        )

    transform, transform_source = _resolve_canonical_transform(
        requested_input=requested_input,
        parent_landscape=parent_landscape,
        raw_coordinates=raw_coordinates,
        canonical_coordinates=canonical_coordinates,
        canonical_id=canonical_id,
        explicit_frame=canonical_frame,
    )
    if space == "canonical" and transform is None:
        raise ValueError(
            "canonical rotation requires a canonical transform; provide a standard "
            "canonical CSV/NPZ, a parent run bundle, or --canonical-frame"
        )

    euler = resolve_euler_convention("extrinsic_zyx", source="rotate_landscape_fixed")
    requested_increment = Rotation.from_euler(
        euler.scipy_euler_sequence,
        angles,
        degrees=True,
    ).as_matrix()
    if space == "canonical":
        assert transform is not None
        raw_increment = transform @ requested_increment @ transform.T
    else:
        raw_increment = requested_increment

    source_rotations = Rotation.from_rotvec(raw_coordinates).as_matrix()
    rotated_raw = Rotation.from_matrix(source_rotations @ raw_increment).as_rotvec()
    rotated_canonical = rotated_raw @ transform if transform is not None else None
    derived_landscape = _derived_landscape(
        parent_landscape,
        rotated_raw=rotated_raw,
        rotated_canonical=rotated_canonical,
        canonical_transform=transform,
        rotation_policy={
            "operation": ROTATION_OPERATION,
            "ro_transform": RO_TRANSFORM,
            "requested_space": space,
            "rotation_euler_deg": angles.tolist(),
            "euler_convention": euler.euler_convention,
            "scipy_euler_sequence": euler.scipy_euler_sequence,
        },
    )

    output_path = Path(output_dir).resolve()
    _validate_output_target(
        requested_input=requested_input.resolve(),
        source_path=source_path,
        output_dir=output_path,
        overwrite=overwrite,
    )
    _remove_stale_canonical_artifacts(
        output_path,
        canonical_id=canonical_id,
        canonical_output_written=transform is not None,
        overwrite=overwrite,
    )

    raw_dir = output_path / "data"
    raw_dir.mkdir(parents=True, exist_ok=True)
    raw_npz = write_landscape_npz(
        derived_landscape,
        raw_dir / "raw_landscape.npz",
        overwrite=overwrite,
        artifact_type="rotated_raw_landscape",
    )
    raw_csv = write_raw_landscape_csv(
        derived_landscape,
        raw_dir / "raw_landscape.csv",
        overwrite=overwrite,
        euler_sequence=euler.scipy_euler_sequence,
    )

    canonical_outputs: dict[str, Path] = {}
    if transform is not None:
        canonical_dir = output_path / "canonical" / canonical_id
        canonical_dir.mkdir(parents=True, exist_ok=True)
        canonical_npz = write_landscape_npz(
            derived_landscape,
            canonical_dir / "canonical_landscape.npz",
            overwrite=overwrite,
            artifact_type="rotated_canonical_landscape",
        )
        canonical_csv = write_canonical_landscape_csv(
            derived_landscape,
            canonical_dir / "canonical_landscape.csv",
            overwrite=overwrite,
            euler_sequence=euler.scipy_euler_sequence,
        )
        frame_npz = canonical_dir / "canonical_frame.npz"
        if frame_npz.exists() and not overwrite:
            raise FileExistsError(f"Output already exists: {frame_npz}")
        np.savez(
            frame_npz,
            canonical_transform=np.asarray(transform, dtype=float),
        )
        frame_payload = {
            "artifact_type": "inherited_canonical_frame",
            "schema_version": "1",
            "canonical_transform": transform,
            "transform_direction": CANONICAL_TRANSFORM_DIRECTION,
            "source": transform_source,
            "source_landscape_path": source_path,
            "canonical_frame_npz": frame_npz,
        }
        frame_json = write_json_artifact(
            frame_payload,
            canonical_dir / "canonical_frame.json",
            overwrite=overwrite,
        )
        canonical_outputs = {
            "canonical_landscape_npz": canonical_npz,
            "canonical_landscape_csv": canonical_csv,
            "canonical_frame_npz": frame_npz,
            "canonical_frame_json": frame_json,
        }

    parent_metadata = _read_parent_metadata(requested_input, source_path)
    requested_sld_metric, resolved_sld_metric = _inherited_sld_metrics(parent_metadata)
    now = datetime.now(timezone.utc).isoformat()
    output_paths: dict[str, Any] = {
        "raw_landscape_npz": raw_npz,
        "raw_landscape_csv": raw_csv,
        **canonical_outputs,
    }
    rotation_report_path = output_path / "rotation_report.json"
    run_summary_path = output_path / "run_summary.json"
    manifest_path = output_path / "run_manifest.json"
    output_paths.update(
        {
            "rotation_report": rotation_report_path,
            "run_summary": run_summary_path,
            "run_manifest": manifest_path,
        }
    )
    euler_metadata = euler.metadata(euler_angle_columns=RAW_EULER_ANGLE_COLUMNS)
    report: dict[str, Any] = {
        "artifact_type": "rotated_landscape_report",
        "schema_version": "1",
        "created_at_utc": now,
        "operation": ROTATION_OPERATION,
        "ro_transform": RO_TRANSFORM,
        "source_landscape_path": source_path,
        "source_row_count": len(parent_landscape.data),
        "requested_space": space,
        "rotation_euler_deg": {
            "alpha": float(angles[0]),
            "beta": float(angles[1]),
            "gamma": float(angles[2]),
        },
        "euler_convention": euler.euler_convention,
        "scipy_euler_sequence": euler.scipy_euler_sequence,
        "rotation_matrix_requested_space": requested_increment,
        "rotation_matrix_raw": raw_increment,
        "canonical_id": canonical_id if transform is not None else None,
        "canonical_transform": transform,
        "canonical_transform_direction": (
            CANONICAL_TRANSFORM_DIRECTION if transform is not None else None
        ),
        "canonical_transform_source": transform_source,
        "density_source": "parent_landscape",
        "sld_recomputed": False,
        "requested_sld_metric": requested_sld_metric,
        "resolved_sld_metric": resolved_sld_metric,
        "source_rows_or_metadata_modified": False,
        "warnings": _density_inheritance_warnings(resolved_sld_metric),
        "output_paths": output_paths,
        "visualize_commands": _visualize_commands(
            output_path,
            canonical_id=canonical_id if transform is not None else None,
        ),
    }
    summary = {
        "artifact_type": "rotated_landscape_summary",
        "schema_version": "1",
        "created_at_utc": now,
        "derived_bundle": True,
        "source_landscape_path": source_path,
        "row_count": len(parent_landscape.data),
        "operation": ROTATION_OPERATION,
        "requested_space": space,
        "density_source": "parent_landscape",
        "sld_recomputed": False,
        "requested_sld_metric": requested_sld_metric,
        "resolved_sld_metric": resolved_sld_metric,
        **euler_metadata,
        "output_artifacts": output_paths,
    }
    manifest = {
        "artifact_type": "rotated_landscape_manifest",
        "schema_version": "1",
        "created_at_utc": now,
        "workflow_name": "rotate_landscape",
        "derived_bundle": True,
        "source_landscape_path": source_path,
        "row_counts": {
            "source_landscape": len(parent_landscape.data),
            "rotated_landscape": len(derived_landscape.data),
        },
        "active_policies": {
            "rotation_policy": {
                "operation": ROTATION_OPERATION,
                "ro_transform": RO_TRANSFORM,
                "requested_space": space,
                "rotation_euler_deg": angles,
                "euler_convention": euler.euler_convention,
                "scipy_euler_sequence": euler.scipy_euler_sequence,
            },
            "density_policy": {
                "density_source": "parent_landscape",
                "sld_recomputed": False,
                "requested_sld_metric": requested_sld_metric,
                "resolved_sld_metric": resolved_sld_metric,
            },
        },
        "results": {
            "canonical_transform_source": transform_source,
            "euler_metadata": euler_metadata,
        },
        "output_artifacts": output_paths,
    }
    write_json_artifact(report, rotation_report_path, overwrite=overwrite)
    write_json_artifact(summary, run_summary_path, overwrite=overwrite)
    write_json_artifact(manifest, manifest_path, overwrite=overwrite)
    return report


def _derived_landscape(
    parent: Landscape,
    *,
    rotated_raw: np.ndarray,
    rotated_canonical: np.ndarray | None,
    canonical_transform: np.ndarray | None,
    rotation_policy: dict[str, Any],
) -> Landscape:
    data = parent.data.copy(deep=True)
    data["coordinates_analysis"] = [row.copy() for row in rotated_raw]
    data["coordinates_display"] = [row.copy() for row in rotated_raw]
    if rotated_canonical is None:
        if "coordinates_canonical" in data.columns:
            data = data.drop(columns=["coordinates_canonical"])
    else:
        data["coordinates_canonical"] = [row.copy() for row in rotated_canonical]
    active_policies = dict(parent.active_policies or {})
    active_policies["rotation_policy"] = rotation_policy
    return Landscape(
        data=data,
        canonical_transform=canonical_transform,
        active_policies=active_policies,
    )


def _stack_coordinates(landscape: Landscape, column: str) -> np.ndarray:
    values = [np.asarray(value, dtype=float) for value in landscape.data[column]]
    if not values:
        return np.empty((0, 3), dtype=float)
    coordinates = np.vstack(values)
    if coordinates.shape != (len(values), 3):
        raise ValueError(f"{column} must have shape (n, 3)")
    return coordinates


def _optional_coordinates(landscape: Landscape, column: str) -> np.ndarray | None:
    if column not in landscape.data.columns:
        return None
    return _stack_coordinates(landscape, column)


def _resolve_canonical_transform(
    *,
    requested_input: Path,
    parent_landscape: Landscape,
    raw_coordinates: np.ndarray,
    canonical_coordinates: np.ndarray | None,
    canonical_id: str,
    explicit_frame: str | Path | None,
) -> tuple[np.ndarray | None, str | None]:
    if explicit_frame is not None:
        transform = _read_canonical_frame(Path(explicit_frame))
        if (
            parent_landscape.canonical_transform is not None
            and not np.allclose(
                transform,
                parent_landscape.canonical_transform,
                atol=1e-8,
            )
        ):
            raise ValueError(
                "--canonical-frame disagrees with the transform stored in the landscape"
            )
        return transform, "explicit_canonical_frame"

    if parent_landscape.canonical_transform is not None:
        transform = np.asarray(parent_landscape.canonical_transform, dtype=float)
        validate_canonical_transform(transform)
        return transform, "landscape_artifact"

    if requested_input.is_dir():
        frame = _find_parent_canonical_frame(requested_input, canonical_id)
        if frame is not None:
            return _read_canonical_frame(frame), "parent_run_canonical_frame"

    if canonical_coordinates is not None:
        return (
            _infer_canonical_transform(raw_coordinates, canonical_coordinates),
            "inferred_from_landscape_coordinates",
        )
    return None, None


def _find_parent_canonical_frame(run_dir: Path, canonical_id: str) -> Path | None:
    base = run_dir / "canonical" / canonical_id
    for candidate in (base / "canonical_frame.npz", base / "canonical_frame.json"):
        if candidate.exists():
            return candidate
    return None


def _read_canonical_frame(path: Path) -> np.ndarray:
    if not path.exists():
        raise ValueError(f"Canonical frame does not exist: {path}")
    if path.is_dir():
        for candidate in (path / "canonical_frame.npz", path / "canonical_frame.json"):
            if candidate.exists():
                return _read_canonical_frame(candidate)
        raise ValueError(f"No canonical_frame.json/.npz found under: {path}")
    if path.suffix.lower() == ".npz":
        with np.load(path, allow_pickle=False) as payload:
            if "canonical_transform" not in payload:
                raise ValueError(f"Canonical frame NPZ missing canonical_transform: {path}")
            transform = np.asarray(payload["canonical_transform"], dtype=float)
    elif path.suffix.lower() == ".json":
        payload = json.loads(path.read_text(encoding="utf-8"))
        if "canonical_transform" in payload:
            transform = np.asarray(payload["canonical_transform"], dtype=float)
        elif payload.get("canonical_frame_npz"):
            nested = Path(str(payload["canonical_frame_npz"]))
            if not nested.is_absolute():
                nested = path.parent / nested
            return _read_canonical_frame(nested)
        else:
            raise ValueError(f"Canonical frame JSON missing canonical_transform: {path}")
    else:
        raise ValueError("Canonical frame must be a .json or .npz artifact")
    validate_canonical_transform(transform)
    return transform


def _infer_canonical_transform(
    raw_coordinates: np.ndarray,
    canonical_coordinates: np.ndarray,
) -> np.ndarray:
    if raw_coordinates.shape != canonical_coordinates.shape:
        raise ValueError("raw and canonical coordinate arrays must have matching shape")
    least_squares, _, rank, _ = np.linalg.lstsq(
        raw_coordinates,
        canonical_coordinates,
        rcond=None,
    )
    if rank < 3:
        raise ValueError(
            "Cannot infer a unique canonical transform from rank-deficient coordinates; "
            "provide --canonical-frame"
        )
    u, _, vt = np.linalg.svd(least_squares)
    transform = u @ vt
    if np.linalg.det(transform) < 0.0:
        u[:, -1] *= -1.0
        transform = u @ vt
    validate_canonical_transform(transform)
    if not np.allclose(
        raw_coordinates @ transform,
        canonical_coordinates,
        rtol=1e-7,
        atol=1e-9,
    ):
        raise ValueError(
            "Canonical coordinates are inconsistent with one proper canonical transform; "
            "provide the matching --canonical-frame"
        )
    return transform


def _validate_output_target(
    *,
    requested_input: Path,
    source_path: Path,
    output_dir: Path,
    overwrite: bool,
) -> None:
    if output_dir == requested_input or output_dir == source_path:
        raise ValueError("Output directory must differ from the input landscape/run directory")
    if _is_within(source_path, output_dir):
        raise ValueError("Output directory contains the input landscape; choose a separate output")
    if output_dir.exists() and not output_dir.is_dir():
        raise FileExistsError(f"Output path exists and is not a directory: {output_dir}")
    if output_dir.exists() and any(output_dir.iterdir()) and not overwrite:
        raise FileExistsError(
            f"Output directory already exists and is not empty: {output_dir}; "
            "use --overwrite to replace generated artifacts"
        )


def _is_within(path: Path, directory: Path) -> bool:
    try:
        path.relative_to(directory)
    except ValueError:
        return False
    return True


def _remove_stale_canonical_artifacts(
    output_dir: Path,
    *,
    canonical_id: str,
    canonical_output_written: bool,
    overwrite: bool,
) -> None:
    if not overwrite or canonical_output_written:
        return
    canonical_dir = output_dir / "canonical" / canonical_id
    for filename in (
        "canonical_landscape.npz",
        "canonical_landscape.csv",
        "canonical_frame.npz",
        "canonical_frame.json",
    ):
        path = canonical_dir / filename
        if path.is_file():
            path.unlink()


def _read_parent_metadata(requested_input: Path, source_path: Path) -> list[dict[str, Any]]:
    candidates: list[Path] = []
    if requested_input.is_dir():
        candidates.append(requested_input)
    for ancestor in source_path.parents:
        if ancestor not in candidates:
            candidates.append(ancestor)
        if len(candidates) >= 5:
            break
    payloads: list[dict[str, Any]] = []
    for directory in candidates:
        for filename in ("run_summary.json", "run_manifest.json"):
            path = directory / filename
            if path.exists():
                payload = json.loads(path.read_text(encoding="utf-8"))
                if isinstance(payload, dict):
                    payloads.append(payload)
        if payloads:
            break
    return payloads


def _inherited_sld_metrics(
    payloads: Sequence[dict[str, Any]],
) -> tuple[str | None, str | None]:
    requested: str | None = None
    resolved: str | None = None
    for payload in payloads:
        requested_value = payload.get("requested_sld_metric")
        resolved_value = payload.get("resolved_sld_metric")
        if isinstance(requested_value, str):
            requested = requested_value
        if isinstance(resolved_value, str):
            resolved = resolved_value
        active_policies = payload.get("active_policies")
        if isinstance(active_policies, dict):
            density = active_policies.get("density_policy")
            if isinstance(density, dict):
                metric = density.get("sld_metric")
                if isinstance(metric, str):
                    requested = requested or metric
                    resolved = resolved or metric
    return requested, resolved


def _visualize_commands(
    output_dir: Path,
    *,
    canonical_id: str | None,
) -> list[str]:
    commands = [
        f'cryorole visualize --run-dir "{output_dir}" --space raw --representation both'
    ]
    if canonical_id is not None:
        commands.append(
            "cryorole visualize "
            f'--run-dir "{output_dir}" --space canonical '
            f'--canonical-id "{canonical_id}" --representation both'
        )
    return commands


def _density_inheritance_warnings(resolved_sld_metric: str | None) -> list[str]:
    warning = (
        "SLD fields are inherited from the parent landscape so colors remain "
        "particle-comparable; density was not recomputed after rotation."
    )
    if resolved_sld_metric == "rotvec_euclidean":
        warning += (
            " RV-Euclidean distances are not generally invariant under this "
            "SO(3) composition, so inherited SLD is parent-density provenance."
        )
    elif resolved_sld_metric == "so3_geodesic":
        warning += (
            " Pairwise SO(3) geodesic distances are invariant under the common "
            "right rotation."
        )
    return [warning]


if __name__ == "__main__":
    raise SystemExit(main())
