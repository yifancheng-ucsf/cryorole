"""Scientific selection workflow and artifact service."""

from __future__ import annotations

import secrets

from dataclasses import dataclass, fields, replace
from datetime import datetime, timezone
import json
from pathlib import Path
from typing import Any, Sequence

import numpy as np
import pandas as pd
from scipy.spatial.transform import Rotation

from cryorole.core.density import compute_sld_display_values, compute_sld_values
from cryorole.core.euler_conventions import DEFAULT_EULER_CONVENTION, LEGACY_MISSING_EULER_CONVENTION_SOURCE, RAW_EULER_ANGLE_COLUMNS, resolve_euler_convention
from cryorole.export import (
    landscape_from_arrays,
    read_landscape,
    read_landscape_metadata,
    read_landscape_npz_arrays,
    write_json_artifact,
    write_landscape_npz,
)
from cryorole.export.serialization import to_json_safe
from cryorole.io.readers import read_relion_star
from cryorole.io.readers.cs_reader import read_cryosparc_cs_column
from cryorole.provenance.source_identity import verify_source_identity
from cryorole.select.metadata import CsMetadataColumn, CS_METADATA_MAX_GROUPS
from cryorole.models.landscape import Landscape
from cryorole.models.landscape_arrays import LandscapeArrays
from cryorole.models.policies import DensityPolicy, SelectionPolicy
from cryorole.run_bundle import validate_completed_run_bundle
from cryorole.select.artifacts import SelectedRowProvenance, SelectionArtifactRequest, write_selection_artifact
from cryorole.select.selectors import (
    select_particles, _metadata_value_key as _cli_metadata_value_key,
    _coerce_source_row_id as _coerce_cli_source_row_id, _resolve_metadata_column_name,
)
from cryorole.workflows.input_policy import resolve_source_type
from cryorole.logs import warn_user


@dataclass(frozen=True)
class SelectRequest:
    run_dir: str
    selection_id: str | None = None
    space: str = "raw"
    canonical_id: str = "default"
    selection_mode: str = "radius"
    center: tuple[float, float, float] | None = None
    center_representation: str = "euler"
    radius: float | None = None
    radius_rad: float | None = None
    metric: str = "so3"
    sld_min: float | None = None
    sld_max: float | None = None
    range_bound: tuple[Any, ...] = ()
    fraction: float | None = None
    seed: int | None = None
    metadata_domain: str | None = None
    metadata_column: str | None = None
    metadata_value: str | None = None
    split_by_value: bool = False
    write_selected_landscape: bool = False
    recompute_sld: bool = False
    overwrite: bool = False
    euler_convention: str | None = None
    explicit_options: tuple[str, ...] = ()
    # How omitted --run-dir / --canonical-id were filled in (cryorole.workflow.resolve).
    resolved_by: dict[str, dict[str, str]] | None = None
    # "user" or "generated": random mode always records the seed it used.
    seed_source: str | None = None

    @classmethod
    def from_namespace(cls, namespace: Any) -> "SelectRequest":
        values = vars(namespace)
        return cls(**{field.name: values[field.name] for field in fields(cls) if field.name in values})


@dataclass(frozen=True)
class SelectResult:
    output_dir: Path
    selection_dirs: tuple[Path, ...]
    selected_counts: tuple[int, ...] = ()


@dataclass
class _ArraySelectionLandscape:
    """Minimal compatibility view backed by compact NPZ arrays."""

    data: pd.DataFrame
    canonical_transform: np.ndarray | None = None
    active_policies: dict[str, object] | None = None
    canonicalization_report: object | None = None


def _validate_select_request(request: SelectRequest) -> None:
    """Validate mode intent before reading inputs or creating artifacts."""
    if not request.selection_id or not request.selection_id.strip():
        raise ValueError("Choose a name for this selection: add --selection-id region_01. No default name is assigned.")
    if request.selection_id in {".", ".."} or any(
        char in "/\\:" or ord(char) < 32 for char in request.selection_id
    ):
        raise ValueError("--selection-id must be a name, not a path; for example region_01")
    modes = {
        "radius": ("center", "center_representation", "radius", "radius_rad", "metric"),
        "threshold": ("sld_min", "sld_max"),
        "range": ("range_bound",),
        "random": ("fraction", "seed"),
        "metadata": ("metadata_domain", "metadata_column", "metadata_value", "split_by_value"),
    }
    mode = request.selection_mode
    if mode not in modes:
        raise ValueError(f"Unsupported --mode {mode!r}; choose from {', '.join(modes)}")
    defaults = {"center_representation": "euler", "metric": "so3", "range_bound": (), "split_by_value": False}
    for owner, names in modes.items():
        for name in names:
            value = getattr(request, name)
            supplied = name in request.explicit_options or (
                bool(value) if name == "range_bound" else value != defaults.get(name)
            )
            if supplied and owner != mode:
                option = "--" + name.replace("_", "-")
                raise ValueError(f"{option} applies only to --mode {owner}; current mode is {mode!r}. Remove it or change --mode.")
    if request.space != "canonical" and (
        "canonical_id" in request.explicit_options
        or request.canonical_id not in (None, "default")
    ):
        raise ValueError("--canonical-id requires --space canonical")
    required = {
        "radius": [(request.center is not None, "--center A B C"),
                   (request.radius is not None or request.radius_rad is not None, "--radius DEG or --radius-rad RAD")],
        "threshold": [(request.sld_min is not None or request.sld_max is not None, "--sld-min VALUE, --sld-max VALUE, or both")],
        "range": [(bool(request.range_bound), "--range-bound AXIS:LOWER:UPPER")],
        "random": [(request.fraction is not None, "--fraction F")],
        "metadata": [(bool(request.metadata_domain), "--metadata-domain ref|mov"),
                     (bool(request.metadata_column), "--metadata-column COLUMN"),
                     (bool(request.metadata_value) or request.split_by_value, "--metadata-value VALUES or --split-by-value")],
    }
    for present, usage in required[mode]:
        if not present:
            raise ValueError(f"Mode {mode!r} requires {usage}. See cryorole select --help for a complete example.")
    _resolve_radius_args(request)
    _internal_selection_mode(request)
    if request.fraction is not None and (not np.isfinite(request.fraction) or not 0 < request.fraction <= 1):
        raise ValueError("--fraction must be finite and in (0, 1]")
    for name in ("sld_min", "sld_max"):
        value = getattr(request, name)
        if value is not None and (not np.isfinite(value) or value < 0):
            raise ValueError(f"--{name.replace('_', '-')} must be finite and non-negative")
    if request.sld_min is not None and request.sld_max is not None and request.sld_min > request.sld_max:
        raise ValueError("--sld-min must be <= --sld-max")
    if request.seed is not None and request.seed < 0:
        raise ValueError("--seed must be a non-negative integer")
    if request.recompute_sld and not request.write_selected_landscape:
        raise ValueError("--recompute-sld requires --write-selected-landscape")


def create_selection(request: SelectRequest) -> SelectResult:
    """Load a saved landscape, select particles, and export selection artifacts."""

    _validate_select_request(request)
    if request.selection_mode == "random":
        # Record every seed used: an omitted --seed draws one, so the selection stays reproducible.
        request = (
            replace(request, seed_source="user")
            if request.seed is not None
            else replace(request, seed=secrets.randbelow(2**31), seed_source="generated")
        )
    args = request
    validate_completed_run_bundle(args.run_dir)
    landscape_source = args.run_dir
    landscape_metadata = read_landscape_metadata(
        landscape_source,
        space=args.space,
        canonical_id=args.canonical_id,
    )
    euler_metadata = _resolve_landscape_euler_metadata(
        args,
        columns=RAW_EULER_ANGLE_COLUMNS,
    )
    policy = _selection_policy_from_args(args, euler_metadata=euler_metadata)
    arrays = None
    landscape_path = Path(str(landscape_metadata["path"]))
    if landscape_path.suffix.casefold() == ".npz":
        arrays = read_landscape_npz_arrays(landscape_path)
        landscape = _selection_landscape_from_arrays(arrays, policy)
    else:
        landscape = read_landscape(
            landscape_source,
            space=args.space,
            canonical_id=args.canonical_id,
        )

    source_metadata = None
    if policy.selection_mode in {"metadata_value", "metadata_group"}:
        source_metadata, metadata_source_file, source_details = _load_run_source_metadata_for_selection(
            args,
            policy=policy,
            landscape=landscape,
        )
        policy = replace(
            policy,
            metadata_source_file=str(metadata_source_file),
            metadata_source_details=source_details,
            metadata_source_row_id_field=f"{policy.metadata_domain}_source_row_id",
        )

    if policy.selection_mode == "metadata_group":
        selection_dirs = _write_metadata_group_selections(
            args=args,
            landscape=landscape,
            landscape_metadata=landscape_metadata,
            euler_metadata=euler_metadata,
            policy=policy,
            source_metadata=source_metadata,
            arrays=arrays,
        )
        return SelectResult(
            output_dir=_metadata_group_output_root(args),
            selection_dirs=tuple(path for path, _count in selection_dirs),
            selected_counts=tuple(count for _path, count in selection_dirs),
        )

    selection = select_particles(landscape, policy=policy, source_metadata=source_metadata)
    selection = replace(selection, parent_run_id=_run_id_from_bundle(args.run_dir))
    _validate_derived_selection_count(args, selection.selected_count)
    output_landscape = _selection_output_landscape(
        landscape,
        arrays=arrays,
        selected_particle_keys=selection.selected_particle_keys,
        materialize_full_rows=bool(args.write_selected_landscape),
    )
    output_dir = _prepare_output_dir(
        _selection_output_dir(args, policy=policy),
        overwrite=args.overwrite,
        label="Selection",
    )
    _write_selection_outputs(
        args=args,
        landscape=output_landscape,
        landscape_metadata=landscape_metadata,
        euler_metadata=euler_metadata,
        policy=policy,
        selection=selection,
        output_dir=output_dir,
    )
    return SelectResult(output_dir=output_dir, selection_dirs=(output_dir,),
                        selected_counts=(len(selection.selected_particle_keys),))


def _write_selection_outputs(
    *,
    args,
    landscape: Landscape,
    landscape_metadata: dict[str, object],
    euler_metadata: dict[str, object],
    policy: SelectionPolicy,
    selection,
    output_dir: Path,
) -> None:
    selected_data = _selected_rows_from_parent_landscape_allow_empty(
        landscape,
        selection.selected_particle_keys,
    )
    selected_landscape_report = None
    if args.write_selected_landscape:
        selected_landscape_report = _write_selected_derived_landscape(
            parent_landscape=landscape,
            selection=selection,
            output_dir=output_dir / "selected_landscape",
            overwrite=args.overwrite,
            recompute_sld=bool(args.recompute_sld),
            density_policy=_density_policy_from_parent_run(args, landscape),
            coordinate_space=args.space,
            parent_landscape_metadata=landscape_metadata,
            euler_metadata=euler_metadata,
        )
    write_selection_artifact(
        SelectionArtifactRequest(
            selection=selection,
            output_dir=output_dir,
            overwrite=args.overwrite,
            selected_rows=SelectedRowProvenance(
                particle_key=selected_data["particle_key"].tolist(),
                ref_source_row_id=(
                    selected_data["ref_source_row_id"].tolist()
                    if "ref_source_row_id" in selected_data.columns
                    else None
                ),
                mov_source_row_id=(
                    selected_data["mov_source_row_id"].tolist()
                    if "mov_source_row_id" in selected_data.columns
                    else None
                ),
            ),
            summary={
            "artifact_type": "selection_summary",
            "schema_version": "1",
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "selection_id": selection.selection_id,
            "overwrite": bool(args.overwrite),
            "mode": _public_selection_mode(selection.selection_mode),
            "run_dir": str(args.run_dir),
            "resolved_by": dict(args.resolved_by or {}),
            "parent_space": args.space,
            "canonical_id": args.canonical_id if args.space == "canonical" else None,
            "parent_landscape_path": landscape_metadata["path"],
            "candidate_count": selection.total_count,
            "landscape_path": landscape_metadata["path"],
            "landscape_artifact_type": landscape_metadata["artifact_type"],
            "landscape_schema_version": landscape_metadata["schema_version"],
            "landscape_row_count": landscape_metadata["row_count"],
            "selection_mode": selection.selection_mode,
            "selection_basis": selection.selection_basis,
            "metric": selection.metric,
            "random_fraction": selection.random_fraction,
            "random_seed": selection.random_seed,
            "random_seed_source": args.seed_source,
            "random_candidate_count": selection.random_candidate_count,
            "metadata_domain": selection.metadata_domain,
            "metadata_source_file": selection.metadata_source_file,
            "metadata_column": selection.metadata_column,
            "metadata_source_details": to_json_safe(policy.metadata_source_details),
            "metadata_values": selection.metadata_values,
            "metadata_source_row_id_field": selection.metadata_source_row_id_field,
            "metadata_candidate_count": selection.metadata_candidate_count,
            "metadata_missing_count": selection.metadata_missing_count,
            "metadata_invalid_count": selection.metadata_invalid_count,
            "density_artifact_policy": selection.density_artifact_policy,
            "density_artifact_flag_field": selection.density_artifact_flag_field,
            "density_artifact_candidate_count_before": (
                selection.density_artifact_candidate_count_before
            ),
            "density_artifact_candidate_count_after": (
                selection.density_artifact_candidate_count_after
            ),
            "density_artifact_excluded_count": selection.density_artifact_excluded_count,
            "radius": selection.radius,
            "radius_unit": selection.radius_unit,
            "sld_min": selection.sld_min,
            "sld_max": selection.sld_max,
            "selected_count": selection.selected_count,
            "total_count": selection.total_count,
            "selection_policy": to_json_safe(policy),
            "policy": to_json_safe(policy),
            "euler_metadata": euler_metadata,
            **euler_metadata,
            "warnings": [],
            "selected_landscape_written": selected_landscape_report is not None,
            "selected_landscape_report": selected_landscape_report,
            },
        )
    )


def _write_metadata_group_selections(
    *,
    args,
    landscape: Landscape,
    landscape_metadata: dict[str, object],
    euler_metadata: dict[str, object],
    policy: SelectionPolicy,
    source_metadata: pd.DataFrame | CsMetadataColumn,
    arrays: LandscapeArrays | None = None,
) -> list[tuple[Path, int]]:
    if not policy.split_by_metadata:
        raise ValueError("metadata split selection requires --split-by-value")
    values = _metadata_group_values_for_landscape(
        landscape,
        source_metadata=source_metadata,
        policy=policy,
    )
    selection_ids = [_metadata_group_selection_id(policy, value) for value in values]
    if len(set(name.casefold() for name in selection_ids)) != len(selection_ids):
        raise ValueError("Metadata split values collide after selection-name sanitization; use --metadata-value with separate names")
    for name in selection_ids:
        path = _metadata_group_output_root(args) / name
        if path.exists() and (not path.is_dir() or not args.overwrite):
            raise FileExistsError(f"Selection output already exists: {path}. Choose another --selection-id or use --overwrite.")
    plans = []
    for value, selection_id in zip(values, selection_ids):
        child_policy = replace(
            policy, selection_mode="metadata_value", metadata_values=(value,),
            selection_id=selection_id, split_by_metadata=False,
        )
        selection = select_particles(landscape, policy=child_policy, source_metadata=source_metadata)
        selection = replace(selection, parent_run_id=_run_id_from_bundle(args.run_dir))
        _validate_derived_selection_count(args, selection.selected_count)
        plans.append((child_policy, selection))
    if not plans:
        raise ValueError("metadata_group selection found no non-missing metadata values")
    output_dirs: list[tuple[Path, int]] = []
    for child_policy, selection in plans:
        output_landscape = _selection_output_landscape(
            landscape, arrays=arrays, selected_particle_keys=selection.selected_particle_keys,
            materialize_full_rows=bool(args.write_selected_landscape),
        )
        output_dir = _prepare_output_dir(
            _metadata_group_output_root(args) / selection.selection_id,
            overwrite=args.overwrite, label="Selection",
        )
        _write_selection_outputs(
            args=args, landscape=output_landscape, landscape_metadata=landscape_metadata,
            euler_metadata=euler_metadata, policy=child_policy, selection=selection, output_dir=output_dir,
        )
        output_dirs.append((output_dir, selection.selected_count))
    return output_dirs


def _validate_derived_selection_count(args: SelectRequest, count: int) -> None:
    if args.write_selected_landscape and count == 0:
        raise ValueError("Selected-derived landscape cannot be empty")
    if args.recompute_sld and count < 2:
        raise ValueError("SLD recomputation on a selected-derived landscape requires at least two particles")



def _selection_landscape_from_arrays(
    arrays: LandscapeArrays,
    policy: SelectionPolicy,
) -> _ArraySelectionLandscape:
    """Expose only fields required by the selected evaluator mode."""

    data: dict[str, object] = {"particle_key": arrays.particle_key}
    if arrays.ref_source_row_id is not None:
        data["ref_source_row_id"] = arrays.ref_source_row_id
    if arrays.mov_source_row_id is not None:
        data["mov_source_row_id"] = arrays.mov_source_row_id
    density_fields = {
        "sld_unfloored": arrays.sld_unfloored,
        "sld_raw": arrays.sld_raw,
        "sld_display": arrays.sld_display,
        "sld_local_k_mean": arrays.sld_local_k_mean,
        "sld_effective_local_k_mean": arrays.sld_effective_local_k_mean,
        "sld_distance_floor": arrays.sld_distance_floor,
    }
    support = policy.density_support_field
    if support in density_fields:
        data[support] = density_fields[support]
    if policy.density_artifact_policy == "exclude_display_outliers":
        data["sld_display_is_outlier"] = arrays.sld_display_is_outlier
    coordinate_needs: set[str] = set()
    if policy.selection_mode == "radius_around_center":
        coordinate_needs.add(policy.evaluation_space)
    if policy.selection_mode == "range_by_coordinates":
        coordinate_needs.add(policy.range_coordinate_source)
    if "analysis" in coordinate_needs:
        data["coordinates_analysis"] = list(arrays.coordinates_analysis)
    if "canonical" in coordinate_needs:
        if arrays.coordinates_canonical is None:
            raise ValueError("canonical coordinate selection requires coordinates_canonical")
        data["coordinates_canonical"] = list(arrays.coordinates_canonical)
    if "canonical_if_available" in coordinate_needs:
        if arrays.coordinates_canonical is not None:
            data["coordinates_canonical"] = list(arrays.coordinates_canonical)
        else:
            data["coordinates_analysis"] = list(arrays.coordinates_analysis)
    return _ArraySelectionLandscape(
        data=pd.DataFrame(data),
        canonical_transform=arrays.canonical_transform,
    )


def _selection_output_landscape(
    landscape,
    *,
    arrays: LandscapeArrays | None,
    selected_particle_keys: Sequence[object],
    materialize_full_rows: bool,
):
    if arrays is None or not materialize_full_rows:
        return landscape
    selected = {str(key) for key in selected_particle_keys}
    indices = np.flatnonzero(
        np.fromiter(
            (key in selected for key in arrays.particle_key),
            dtype=bool,
            count=arrays.n_points,
        )
    )
    if not indices.size:
        raise ValueError("Selected-derived landscape cannot be empty")
    return landscape_from_arrays(arrays, indices)


def _load_run_source_metadata_for_selection(
    args,
    *,
    policy: SelectionPolicy,
    landscape: Landscape,
) -> tuple[pd.DataFrame | CsMetadataColumn, Path, dict[str, object]]:
    if not getattr(args, "run_dir", None):
        raise ValueError("metadata selection requires --run-dir")
    if policy.metadata_domain not in {"ref", "mov"}:
        raise ValueError("metadata selection requires --metadata-domain ref or mov")
    if not policy.metadata_column:
        raise ValueError("metadata selection requires --metadata-column")
    source_info = _resolve_run_source_info_for_metadata(Path(args.run_dir))
    domain_info = source_info.get(str(policy.metadata_domain), {})
    source_path_raw = domain_info.get("path")
    source_type = domain_info.get("source_type")
    if not source_path_raw:
        raise ValueError(f"Cannot locate original {policy.metadata_domain} source metadata path")
    source_path = Path(source_path_raw)
    is_cs = source_type in {"cryosparc", "cryosparc_cs", "cs"} or source_path.suffix.lower() == ".cs"
    summary = _read_json_if_exists(Path(args.run_dir) / "run_summary.json") or {}
    manifest = _read_json_if_exists(Path(args.run_dir) / "run_manifest.json") or {}
    identity = (summary.get("source_identities") or manifest.get("source_identities") or {}).get(policy.metadata_domain)
    details = {"format": "cryosparc_cs" if is_cs else "relion_star", "verification_status": "legacy_unverified"}
    if is_cs and (not identity or not (identity.get("sha256") or identity.get("hash"))):
        raise ValueError("CS metadata selection requires a recorded source SHA-256; rerun the analysis with verified inputs")
    if identity:
        verified = verify_source_identity(identity, operation="metadata selection")
        source_path = Path(verified.resolved_path)
        details.update(verification_status=verified.verification_status, source_identity=identity)
    if is_cs:
        row_field = f"{policy.metadata_domain}_source_row_id"
        if row_field not in landscape.data.columns:
            raise ValueError(f"Landscape missing metadata source-row field: {row_field}")
        values = read_cryosparc_cs_column(source_path, str(policy.metadata_column))
        column = CsMetadataColumn.from_source(values, landscape.data[row_field].to_numpy())
        details.update(dtype=values.dtype.str, comparison="typed_exact", text_encoding="utf-8",
                       missing_value_policy="exclude_empty_strings", invalid_row_policy="error",
                       split_max_groups=CS_METADATA_MAX_GROUPS)
        return column, source_path, details
    if source_type not in {"relion", "relion_star", "star"} and source_path.suffix.lower() != ".star":
        raise ValueError("metadata selection supports run-time CryoSPARC CS and RELION STAR sources")
    star = read_relion_star(source_path)
    _resolve_metadata_column_name(star.particles.columns, str(policy.metadata_column))
    return star.particles, source_path, details


def _resolve_run_source_info_for_metadata(run_dir: Path) -> dict[str, dict[str, str | None]]:
    summary = _read_json_if_exists(run_dir / "run_summary.json") or {}
    manifest = _read_json_if_exists(run_dir / "run_manifest.json") or {}
    input_paths = _domain_string_mapping(summary.get("input_paths"))
    if not input_paths:
        input_paths = _manifest_input_paths_for_metadata(manifest)
    source_types = _domain_string_mapping(summary.get("source_types"))
    if not source_types:
        source_types = _manifest_source_types_for_metadata(manifest, input_paths)
    info: dict[str, dict[str, str | None]] = {}
    for domain in ("ref", "mov"):
        raw_path = input_paths.get(domain)
        if raw_path is None:
            info[domain] = {"path": None, "source_type": None}
            continue
        path = _resolve_run_relative_path(raw_path, run_dir)
        info[domain] = {
            "path": str(path),
            "source_type": source_types.get(domain) or resolve_source_type(str(path)),
        }
    return info


def _domain_string_mapping(value: object) -> dict[str, str]:
    if not isinstance(value, dict):
        return {}
    return {
        domain: str(value[domain])
        for domain in ("ref", "mov")
        if value.get(domain) is not None
    }


def _manifest_input_paths_for_metadata(manifest: dict[str, object]) -> dict[str, str]:
    records = manifest.get("input_provenance", {}).get("inputs", [])
    if not isinstance(records, list) or len(records) < 2:
        return {}
    paths = []
    for record in records[:2]:
        if isinstance(record, dict) and record.get("path") is not None:
            paths.append(str(record["path"]))
    if len(paths) != 2:
        return {}
    return {"ref": paths[0], "mov": paths[1]}


def _manifest_source_types_for_metadata(
    manifest: dict[str, object],
    input_paths: dict[str, str],
) -> dict[str, str]:
    source_types = manifest.get("input_provenance", {}).get("source_types", {})
    if not isinstance(source_types, dict):
        return {}
    by_domain = {}
    for domain, path in input_paths.items():
        if path in source_types:
            by_domain[domain] = str(source_types[path])
    return by_domain


def _resolve_run_relative_path(path: str, run_dir: Path) -> Path:
    candidate = Path(path)
    if candidate.exists():
        return candidate
    relative = run_dir / candidate
    if relative.exists():
        return relative
    return candidate


def _metadata_group_values_for_landscape(
    landscape: Landscape,
    *,
    source_metadata: pd.DataFrame | CsMetadataColumn,
    policy: SelectionPolicy,
) -> tuple[str, ...]:
    if isinstance(source_metadata, CsMetadataColumn):
        return source_metadata.groups()
    row_field = policy.metadata_source_row_id_field or f"{policy.metadata_domain}_source_row_id"
    if row_field not in landscape.data.columns:
        raise ValueError(f"Landscape missing metadata source-row field: {row_field}")
    column = _resolve_metadata_column_name(source_metadata.columns, str(policy.metadata_column))
    values_by_row_id = {
        int(index): _cli_metadata_value_key(value)
        for index, value in source_metadata[column].items()
    }
    values: list[str] = []
    seen: set[str] = set()
    for source_row_id in landscape.data[row_field]:
        row_id = _coerce_cli_source_row_id(source_row_id)
        if row_id is None or row_id not in values_by_row_id:
            continue
        value = values_by_row_id[row_id]
        if value is None or value in seen:
            continue
        seen.add(value)
        values.append(value)
    return tuple(values)


def _metadata_group_selection_id(policy: SelectionPolicy, value: str) -> str:
    return f"{_sanitize_selection_id_component(policy.selection_id)}_{_sanitize_selection_id_component(value)}"


def _sanitize_selection_id_component(value: str) -> str:
    cleaned = "".join(char if char.isalnum() or char in {"-", "_"} else "_" for char in value)
    cleaned = cleaned.strip("_")
    return cleaned or "value"


def _metadata_group_output_root(args) -> Path:
    if args.run_dir:
        return Path(args.run_dir) / "selections"
    raise ValueError("select requires --run-dir")


def _selected_rows_from_parent_landscape_allow_empty(
    landscape: Landscape,
    selected_particle_keys: Sequence[object],
) -> pd.DataFrame:
    if selected_particle_keys:
        return _selected_rows_from_landscape(landscape, selected_particle_keys)
    columns = [
        column
        for column in ("particle_key", "ref_source_row_id", "mov_source_row_id")
        if column in landscape.data.columns
    ]
    selected = landscape.data.iloc[0:0].copy()
    if "sld_display_is_outlier" not in selected.columns:
        selected["sld_display_is_outlier"] = []
    return selected[columns + [column for column in selected.columns if column not in columns]]


def _write_selected_derived_landscape(
    *,
    parent_landscape: Landscape,
    selection,
    output_dir: Path,
    overwrite: bool,
    recompute_sld: bool,
    density_policy: DensityPolicy,
    coordinate_space: str,
    parent_landscape_metadata: dict[str, object],
    euler_metadata: dict[str, object],
) -> dict[str, object]:
    selected_data = _selected_rows_from_landscape(parent_landscape, selection.selected_particle_keys)
    density_source = "recomputed_on_selection" if recompute_sld else "parent_landscape"
    effective_k = None
    display_diagnostics: dict[str, object] = {}
    if recompute_sld:
        selected_data = _with_parent_sld_fields(selected_data)
        coordinates = _density_coordinates_for_selected_space(selected_data, coordinate_space)
        if len(coordinates) < 2:
            raise ValueError(
                "SLD recomputation on a selected-derived landscape requires at least "
                "two selected particles"
            )
        sld = compute_sld_values(
            coordinates,
            k_neighbors=density_policy.k_neighbors,
            distance_floor_fraction=density_policy.distance_floor_fraction,
            sld_metric=density_policy.sld_metric,
        )
        display = compute_sld_display_values(
            sld["sld_raw"],
            mode=density_policy.display_normalization_mode,
            outlier_mode=density_policy.display_outlier_mode,
            tail_search_fraction=density_policy.tail_search_fraction,
            tail_jump_factor=density_policy.tail_jump_factor,
            max_display_outlier_fraction=density_policy.max_display_outlier_fraction,
        )
        selected_data["sld_unfloored"] = sld["sld_unfloored"]
        selected_data["sld_raw"] = sld["sld_raw"]
        selected_data["sld_display"] = display["sld_display"]
        selected_data["sld_display_is_outlier"] = display["sld_display_is_outlier"]
        selected_data["sld_was_floored"] = sld["sld_was_floored"]
        selected_data["sld_local_k_mean"] = sld["sld_local_k_mean"]
        selected_data["sld_effective_local_k_mean"] = sld["sld_effective_local_k_mean"]
        selected_data["sld_distance_floor"] = sld["sld_distance_floor"]
        effective_k = min(density_policy.k_neighbors, len(selected_data) - 1)
        display_diagnostics = {
            key: value
            for key, value in display.items()
            if key != "sld_display" and key != "sld_display_is_outlier"
        }

    selected_active_policies = dict(parent_landscape.active_policies or {})
    if recompute_sld:
        selected_active_policies["density_policy"] = density_policy
    selected_landscape = Landscape(
        data=selected_data,
        canonical_transform=parent_landscape.canonical_transform,
        active_policies=selected_active_policies,
        density_report=None,
        canonicalization_report=parent_landscape.canonicalization_report,
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    npz_path = write_landscape_npz(
        selected_landscape,
        output_dir / "landscape.npz",
        overwrite=overwrite,
        artifact_type="selected_landscape",
    )
    csv_path = _write_selected_landscape_csv(
        selected_landscape,
        output_dir / "landscape.csv",
        overwrite=overwrite,
        euler_sequence=str(euler_metadata["scipy_euler_sequence"]),
    )
    report = {
        "artifact_type": "selected_landscape_report",
        "schema_version": "1",
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "selection_id": selection.selection_id,
        "parent_landscape_path": parent_landscape_metadata.get("path"),
        "parent_landscape_artifact_type": parent_landscape_metadata.get("artifact_type"),
        "parent_landscape_schema_version": parent_landscape_metadata.get("schema_version"),
        "parent_landscape_row_count": parent_landscape_metadata.get("row_count"),
        "selected_count": selection.selected_count,
        "coordinate_space": coordinate_space,
        "density_source": density_source,
        "sld_recomputed": recompute_sld,
        "requested_sld_metric": density_policy.sld_metric if recompute_sld else None,
        "resolved_sld_metric": density_policy.sld_metric if recompute_sld else None,
        "requested_k_neighbors": density_policy.k_neighbors if recompute_sld else None,
        "effective_k_neighbors": effective_k,
        "distance_floor_fraction": (
            density_policy.distance_floor_fraction if recompute_sld else None
        ),
        "display_outlier_policy": (
            {
                "display_outlier_mode": density_policy.display_outlier_mode,
                "tail_search_fraction": density_policy.tail_search_fraction,
                "tail_jump_factor": density_policy.tail_jump_factor,
                "max_display_outlier_fraction": density_policy.max_display_outlier_fraction,
            }
            if recompute_sld
            else None
        ),
        "display_diagnostics": display_diagnostics,
        "parent_sld_fields_preserved": recompute_sld,
        "euler_metadata": euler_metadata,
        "output_paths": {
            "landscape_npz": str(npz_path),
            "landscape_csv": str(csv_path),
            "landscape_report_json": str(output_dir / "landscape_report.json"),
            "selected_landscape_rows_csv": str(output_dir.parent / "selected_landscape_rows.csv"),
        },
    }
    report_path = write_json_artifact(
        report,
        output_dir / "landscape_report.json",
        overwrite=overwrite,
    )
    report["output_paths"]["landscape_report_json"] = str(report_path)
    return report


def _selected_rows_from_landscape(
    landscape: Landscape,
    selected_particle_keys: Sequence[object],
) -> pd.DataFrame:
    if "particle_key" not in landscape.data.columns:
        raise ValueError("Landscape must contain particle_key to write selected landscape")
    data = landscape.data.copy(deep=True)
    key_to_indices: dict[str, list[int]] = {}
    for row_index, particle_key in enumerate(data["particle_key"].map(str)):
        key_to_indices.setdefault(particle_key, []).append(row_index)
    selected_indices = []
    for particle_key in selected_particle_keys:
        key = str(particle_key)
        if key not in key_to_indices:
            raise ValueError(
                "Selection contains particle_key values absent from the parent landscape: "
                f"{key}"
            )
        selected_indices.extend(key_to_indices[key])
    if not selected_indices:
        raise ValueError("Selected-derived landscape cannot be empty")
    selected = data.iloc[selected_indices].reset_index(drop=True)
    if "sld_display_is_outlier" not in selected.columns:
        selected["sld_display_is_outlier"] = False
    return selected


def _with_parent_sld_fields(data: pd.DataFrame) -> pd.DataFrame:
    output = data.copy(deep=True)
    for field in (
        "sld_unfloored",
        "sld_raw",
        "sld_display",
        "sld_display_is_outlier",
        "sld_was_floored",
        "sld_local_k_mean",
        "sld_effective_local_k_mean",
        "sld_distance_floor",
    ):
        if field in output.columns:
            output[f"parent_{field}"] = output[field]
    return output


def _density_coordinates_for_selected_space(data: pd.DataFrame, coordinate_space: str) -> np.ndarray:
    if coordinate_space == "canonical":
        if "coordinates_canonical" not in data.columns:
            raise ValueError("canonical selected SLD recomputation requires coordinates_canonical")
        column = "coordinates_canonical"
    else:
        column = "coordinates_analysis"
    return np.vstack([np.asarray(value, dtype=float) for value in data[column]])


def _write_selected_landscape_csv(
    landscape: Landscape,
    path: Path,
    *,
    overwrite: bool,
    euler_sequence: str,
) -> Path:
    if path.exists() and not overwrite:
        raise FileExistsError(f"Output path already exists: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    data = landscape.data
    raw = np.vstack([np.asarray(value, dtype=float) for value in data["coordinates_analysis"]])
    raw_euler = Rotation.from_rotvec(raw).as_euler(euler_sequence, degrees=True)
    output = pd.DataFrame(
        {
            "particle_key": data["particle_key"].map(str),
            "ref_source_row_id": data.get("ref_source_row_id", -1),
            "mov_source_row_id": data.get("mov_source_row_id", -1),
            "raw_rv_x_rad": raw[:, 0],
            "raw_rv_y_rad": raw[:, 1],
            "raw_rv_z_rad": raw[:, 2],
            "raw_angle_deg": np.degrees(np.linalg.norm(raw, axis=1)),
            "raw_ea_zyx_alpha_deg": raw_euler[:, 0],
            "raw_ea_zyx_beta_deg": raw_euler[:, 1],
            "raw_ea_zyx_gamma_deg": raw_euler[:, 2],
            "sld_unfloored": data["sld_unfloored"],
            "sld_raw": data["sld_raw"],
            "sld_display": data["sld_display"],
            "sld_display_is_outlier": data["sld_display_is_outlier"],
            "sld_was_floored": data["sld_was_floored"],
            "sld_local_k_mean": data["sld_local_k_mean"],
            "sld_effective_local_k_mean": data["sld_effective_local_k_mean"],
            "sld_distance_floor": data["sld_distance_floor"],
        }
    )
    if "coordinates_canonical" in data.columns:
        canonical = np.vstack(
            [np.asarray(value, dtype=float) for value in data["coordinates_canonical"]]
        )
        canonical_euler = Rotation.from_rotvec(canonical).as_euler(
            euler_sequence,
            degrees=True,
        )
        output["canonical_rv_x_rad"] = canonical[:, 0]
        output["canonical_rv_y_rad"] = canonical[:, 1]
        output["canonical_rv_z_rad"] = canonical[:, 2]
        output["canonical_ea_zyx_alpha_deg"] = canonical_euler[:, 0]
        output["canonical_ea_zyx_beta_deg"] = canonical_euler[:, 1]
        output["canonical_ea_zyx_gamma_deg"] = canonical_euler[:, 2]
    for column in (
        "parent_sld_unfloored",
        "parent_sld_raw",
        "parent_sld_display",
        "parent_sld_display_is_outlier",
        "parent_sld_was_floored",
        "parent_sld_local_k_mean",
        "parent_sld_effective_local_k_mean",
        "parent_sld_distance_floor",
    ):
        if column in data.columns:
            output[column] = data[column]
    output.to_csv(path, index=False)
    return path


def _selection_output_dir(args, *, policy: SelectionPolicy | None = None) -> Path:
    if args.run_dir:
        selection_id = policy.selection_id if policy is not None else args.selection_id
        return Path(args.run_dir) / "selections" / str(selection_id)
    raise ValueError("select requires --run-dir")


def _resolve_landscape_euler_metadata(
    args,
    *,
    columns=RAW_EULER_ANGLE_COLUMNS,
    explicit_sequence: str | None = None,
    warn_on_legacy_missing: bool = True,
) -> dict[str, object]:
    if getattr(args, "euler_convention", None):
        resolved = resolve_euler_convention(args.euler_convention, source="cli_override")
        return resolved.metadata(euler_angle_columns=columns)
    if explicit_sequence:
        resolved = resolve_euler_convention(
            scipy_euler_sequence=explicit_sequence,
            source="cli_override",
        )
        return resolved.metadata(euler_angle_columns=columns)
    inherited = _read_parent_euler_metadata(args, columns=columns)
    if inherited is not None:
        return inherited
    resolved = resolve_euler_convention(
        DEFAULT_EULER_CONVENTION,
        source=LEGACY_MISSING_EULER_CONVENTION_SOURCE,
    )
    metadata = resolved.metadata(euler_angle_columns=columns)
    if warn_on_legacy_missing:
        warn_user(
            "parent landscape lacks Euler convention metadata; "
            f"defaulting to {DEFAULT_EULER_CONVENTION} for derived Euler coordinates."
        )
    return metadata


def _read_parent_euler_metadata(args, *, columns) -> dict[str, object] | None:
    run_dir_raw = getattr(args, "run_dir", None)
    if not run_dir_raw:
        return None
    run_dir = Path(run_dir_raw)
    payloads: list[dict[str, object]] = []
    if getattr(args, "space", "raw") == "canonical":
        canonical_id = getattr(args, "canonical_id", "default")
        payload = _read_json_if_exists(
            run_dir / "canonical" / canonical_id / "canonicalize_summary.json"
        )
        if payload is not None:
            payloads.append(payload)
    run_summary = _read_json_if_exists(run_dir / "run_summary.json")
    if run_summary is not None:
        payloads.append(run_summary)
    for payload in payloads:
        convention = payload.get("euler_convention")
        sequence = payload.get("scipy_euler_sequence")
        if convention or sequence:
            resolved = resolve_euler_convention(
                convention if isinstance(convention, str) else None,
                scipy_euler_sequence=sequence if isinstance(sequence, str) else None,
                source="inherited",
            )
            return resolved.metadata(euler_angle_columns=columns)
    return None


def _read_json_if_exists(path: Path) -> dict[str, object] | None:
    if not path.exists():
        return None
    with path.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    if not isinstance(payload, dict):
        return None
    return payload


def _run_id_from_bundle(run_dir: str | Path) -> str | None:
    for filename in ("run_summary.json", "run_manifest.json"):
        payload = _read_json_if_exists(Path(run_dir) / filename)
        if payload is not None and payload.get("run_id") is not None:
            return str(payload["run_id"])
    return None


def _density_policy_from_parent_run(args, landscape: Landscape) -> DensityPolicy:
    active_density_policy = (landscape.active_policies or {}).get("density_policy")
    if isinstance(active_density_policy, DensityPolicy):
        return active_density_policy
    if isinstance(active_density_policy, dict):
        return _density_policy_from_mapping(active_density_policy)

    run_dir_raw = getattr(args, "run_dir", None)
    if run_dir_raw:
        run_dir = Path(run_dir_raw)
        manifest = _read_json_if_exists(run_dir / "run_manifest.json") or {}
        active_policies = manifest.get("active_policies")
        if isinstance(active_policies, dict):
            manifest_density_policy = active_policies.get("density_policy")
            if isinstance(manifest_density_policy, dict):
                return _density_policy_from_mapping(manifest_density_policy)

        summary = _read_json_if_exists(run_dir / "run_summary.json") or {}
        resolved_metric = summary.get("resolved_sld_metric")
        if isinstance(resolved_metric, str):
            values: dict[str, object] = {"sld_metric": resolved_metric}
            k_neighbors = summary.get("k_neighbors")
            if isinstance(k_neighbors, int):
                values["k_neighbors"] = k_neighbors
            return _density_policy_from_mapping(values)

    return DensityPolicy()


def _density_policy_from_mapping(values: dict[str, object]) -> DensityPolicy:
    allowed_fields = {field.name for field in fields(DensityPolicy)}
    return DensityPolicy(
        **{
            key: value
            for key, value in values.items()
            if key in allowed_fields
        }
    )


def _prepare_output_dir(
    output_dir: str | Path,
    *,
    overwrite: bool,
    label: str,
) -> Path:
    path = Path(output_dir)
    if path.exists() and not path.is_dir():
        raise ValueError(f"{label} output path exists and is not a directory: {path}")
    if path.exists() and not overwrite:
        raise FileExistsError(
            f"{label} output directory already exists: {path}. "
            f"Choose another --selection-id to save a separate selection, "
            f"or use --overwrite to replace {path.name!r}."
        )
    path.mkdir(parents=True, exist_ok=True)
    return path


def _selection_policy_from_args(
    args,
    *,
    euler_metadata: dict[str, object] | None = None,
) -> SelectionPolicy:
    radius, radius_unit = _resolve_radius_args(args)
    euler_metadata = euler_metadata or resolve_euler_convention().metadata()
    scipy_sequence = str(euler_metadata["scipy_euler_sequence"])
    euler_convention = str(euler_metadata["euler_convention"])
    internal_mode = _internal_selection_mode(args)
    coordinate_source = _selection_coordinate_source_from_space(args.space)
    range_representation = _range_representation_from_bounds(args.range_bound or ())
    parent_metadata = {
        "euler_metadata": dict(euler_metadata),
        "space": args.space,
        "canonical_id": args.canonical_id if args.space == "canonical" else None,
        "run_dir": str(args.run_dir),
    }
    return SelectionPolicy(
        selection_mode=internal_mode,
        density_support_field="sld_raw",
        density_artifact_policy="include_all",
        top_fraction=0.40,
        random_fraction=getattr(args, "fraction", None),
        random_seed=getattr(args, "seed", None),
        metadata_domain=getattr(args, "metadata_domain", None),
        metadata_column=getattr(args, "metadata_column", None),
        metadata_values=_metadata_values_from_args(args),
        split_by_metadata=bool(getattr(args, "split_by_value", False)),
        threshold=getattr(args, "sld_min", None),
        sld_min=getattr(args, "sld_min", None),
        sld_max=getattr(args, "sld_max", None),
        center_input=tuple(args.center) if args.center is not None else None,
        center_input_representation=args.center_representation,
        center_input_space="evaluation",
        center_euler_sequence=scipy_sequence,
        center_euler_convention=euler_convention,
        center_scipy_euler_sequence=scipy_sequence,
        center_degrees=True,
        evaluation_space=coordinate_source,
        metric=_internal_selection_metric(args.metric),
        radius=radius,
        radius_unit=radius_unit,
        range_coordinate_source=coordinate_source,
        range_representation=range_representation,
        range_euler_sequence=scipy_sequence,
        range_euler_convention=euler_convention,
        range_scipy_euler_sequence=scipy_sequence,
        range_degrees=True,
        range_bounds=dict(args.range_bound or ()),
        parent_landscape_metadata=parent_metadata,
        selection_id=args.selection_id,
    )


def _metadata_values_from_args(args) -> tuple[str, ...]:
    value = getattr(args, "metadata_value", None)
    # Interpret values only after loading the source field's dtype.
    return tuple(str(value).split(",")) if value is not None else ()


def _internal_selection_mode(args) -> str:
    mode = args.selection_mode
    if mode == "radius":
        return "radius_around_center"
    if mode == "threshold":
        return "threshold_by_density"
    if mode == "range":
        return "range_by_coordinates"
    if mode == "random":
        return "random"
    if mode == "metadata":
        if getattr(args, "split_by_value", False):
            if getattr(args, "metadata_value", None) is not None:
                raise ValueError("--split-by-value cannot be used with --metadata-value")
            return "metadata_group"
        return "metadata_value"
    raise ValueError(f"Unsupported selection mode: {mode}")


def _public_selection_mode(internal_mode: str) -> str:
    return {
        "radius_around_center": "radius",
        "threshold_by_density": "threshold",
        "range_by_coordinates": "range",
        "random": "random",
        "metadata_value": "metadata",
        "metadata_group": "metadata",
    }.get(internal_mode, internal_mode)


def _selection_coordinate_source_from_space(space: str) -> str:
    if space == "raw":
        return "analysis"
    if space == "canonical":
        return "canonical"
    raise ValueError(f"Unsupported selection space: {space}")


def _internal_selection_metric(metric: str) -> str:
    if metric == "so3":
        return "so3_geodesic"
    if metric == "rotvec":
        return "rotvec_euclidean"
    raise ValueError(f"Unsupported selection metric: {metric}")


def _range_representation_from_bounds(
    range_bounds: Sequence[tuple[str, tuple[float | None, float | None]]],
) -> str:
    axes = {axis for axis, _bounds in range_bounds}
    euler_axes = {"alpha", "beta", "gamma"}
    rotvec_axes = {"x", "y", "z"}
    if not axes:
        return "euler"
    if axes <= euler_axes:
        return "euler"
    if axes <= rotvec_axes:
        return "rotvec"
    unsupported = sorted(axes - euler_axes - rotvec_axes)
    if unsupported:
        raise ValueError(f"Unsupported range axis names: {unsupported}")
    raise ValueError("Do not mix Euler and rotvec range axes in one selection")


def _resolve_radius_args(args) -> tuple[float | None, str | None]:
    provided = [
        name
        for name, value in (
            ("--radius", args.radius),
            ("--radius-rad", args.radius_rad),
        )
        if value is not None
    ]
    if len(provided) > 1:
        raise ValueError("Use only one radius argument: --radius or --radius-rad")
    if args.radius_rad is not None:
        return args.radius_rad, "radians"
    if args.radius is not None:
        return args.radius, "degrees"
    return None, None
