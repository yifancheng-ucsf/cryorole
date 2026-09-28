"""Implementation of ``cryorole.api``: thin, typed wrappers over the services."""

from __future__ import annotations

from dataclasses import fields, replace
from pathlib import Path
from typing import Any, Mapping, Sequence

from cryorole.errors import CryoroleError
from cryorole.logs import notify_user
from cryorole.workflow.resolve import (
    Resolved,
    canonical_frame_ids,
    resolve_canonical_id,
    resolve_run_dir,
    resolve_selection_id,
)
from cryorole.workflows.progress import CancelToken, ProgressCallback


def _build(cls, values: Mapping[str, Any], *, what: str):
    """Instantiate a request dataclass, rejecting unknown option names clearly."""

    known = {field.name for field in fields(cls)}
    unknown = sorted(set(values) - known)
    if unknown:
        raise CryoroleError(
            "input_invalid",
            f"Unknown {what} option(s): {', '.join(unknown)}.",
            f"Valid options: {', '.join(sorted(known))}.",
        )
    return cls(**values)


def _announce(resolved: Resolved) -> None:
    if resolved.implicit:
        notify_user(f"[cryorole] {resolved.notice()}")


def _run_dir(run_dir: str | Path | None, records: dict[str, dict[str, str]]) -> str:
    resolved = resolve_run_dir(run_dir)
    _announce(resolved)
    records["run_dir"] = resolved.record()
    return resolved.value


def _canonical(run_dir: str, space: str, canonical_id: str | None, records: dict[str, dict[str, str]]) -> str:
    if space != "canonical":
        if canonical_id is not None:
            raise CryoroleError("input_invalid", "canonical_id requires space='canonical'.", "Remove canonical_id or set space='canonical'.")
        return "default"
    resolved = resolve_canonical_id(run_dir, canonical_id)
    _announce(resolved)
    records["canonical_id"] = resolved.record()
    return resolved.value


def _input_policies(ref, mov, *, row_aligned, allow_low_overlap, identity_mode, identity_columns, mapping_file):
    from cryorole.workflows.input_policy import InputPolicyRequest, resolve_input_policies

    return resolve_input_policies(
        InputPolicyRequest(
            ref=ref,
            mov=mov,
            row_aligned=row_aligned,
            allow_low_overlap=allow_low_overlap,
            identity_mode=identity_mode,
            identity_columns=tuple(identity_columns or ()),
            mapping_file=mapping_file,
        )
    )


# --- inputs ------------------------------------------------------------------


def preflight(
    ref: str | Path,
    mov: str | Path,
    *,
    row_aligned: bool = False,
    allow_low_overlap: bool = False,
    identity_mode: str | None = None,
    identity_columns: Sequence[str] = (),
    mapping_file: str | Path | None = None,
    **options: Any,
):
    """Inspect two inputs without writing anything; returns ``PreflightResult``.

    ``result.report["readiness"]`` is ``READY``, ``READY_WITH_WARNINGS`` or
    ``BLOCKED``; ``result.exit_code`` matches the CLI.
    """

    from cryorole.preflight import PreflightRequest, run_preflight

    resolved = _input_policies(
        ref, mov, row_aligned=row_aligned, allow_low_overlap=allow_low_overlap,
        identity_mode=identity_mode, identity_columns=identity_columns, mapping_file=mapping_file,
    )
    request = _build(
        PreflightRequest,
        {
            "ref": ref, "mov": mov, "row_aligned": row_aligned, "allow_low_overlap": allow_low_overlap,
            "identity_mode": identity_mode, "identity_columns": tuple(identity_columns or ()),
            "mapping_file": mapping_file, "resolved_input_policies": resolved, **options,
        },
        what="preflight",
    )
    return run_preflight(request)


def run(
    ref: str | Path,
    mov: str | Path,
    *,
    output_dir: str | Path = "cryorole_outputs",
    row_aligned: bool = False,
    allow_low_overlap: bool = False,
    identity_mode: str | None = None,
    identity_columns: Sequence[str] = (),
    mapping_file: str | Path | None = None,
    overwrite: bool = False,
    progress: ProgressCallback | None = None,
    cancel: CancelToken | None = None,
    stream_progress: bool = False,
    **options: Any,
):
    """Create a run bundle transactionally; returns ``RunResult``.

    Nothing is written to stderr unless ``stream_progress=True``; stage events
    go to ``progress``. If ``cancel`` is triggered the staging directory is
    removed and ``CancelledError`` is raised; no bundle is published.
    """

    from cryorole.workflows.run_service import RunExecutionContext, RunRequest, execute_run

    resolved = _input_policies(
        ref, mov, row_aligned=row_aligned, allow_low_overlap=allow_low_overlap,
        identity_mode=identity_mode, identity_columns=identity_columns, mapping_file=mapping_file,
    )
    request = _build(
        RunRequest,
        {
            "ref": str(ref), "mov": str(mov), "output_dir": str(output_dir), "row_aligned": row_aligned,
            "allow_low_overlap": allow_low_overlap, "identity_mode": identity_mode,
            "identity_column": tuple(identity_columns or ()), "mapping_file": mapping_file,
            "overwrite": overwrite, **options,
        },
        what="run",
    )
    return execute_run(
        replace(request, resolved_input_policies=resolved),
        context=RunExecutionContext(progress=progress, cancel_token=cancel, stream_progress=stream_progress),
    )


# --- downstream --------------------------------------------------------------


def canonicalize(run_dir: str | Path | None = None, *, canonical_id: str = "default", **options: Any):
    """Fit (or with ``use_frame=`` apply) a canonical frame; returns ``CanonicalizeResult``."""

    from cryorole.canonicalize.service import CanonicalizeRequest, canonicalize_bundle

    records: dict[str, dict[str, str]] = {}
    resolved_run = _run_dir(run_dir, records)
    request = _build(
        CanonicalizeRequest,
        {"run_dir": resolved_run, "canonical_id": canonical_id, "resolved_by": records, **options},
        what="canonicalize",
    )
    return canonicalize_bundle(request)


def visualize(
    run_dir: str | Path | None = None,
    *,
    space: str = "raw",
    canonical_id: str | None = None,
    selection_id: str | None = None,
    use_selected_landscape: bool = False,
    **options: Any,
):
    """Render display-only figures; returns ``VisualizationResult``."""

    from cryorole.visualize import VisualizationRequest, visualize as _visualize

    records: dict[str, dict[str, str]] = {}
    resolved_run = _run_dir(run_dir, records)
    resolved_canonical = _canonical(resolved_run, space, canonical_id, records)
    if use_selected_landscape:
        chosen = resolve_selection_id(resolved_run, selection_id)
        _announce(chosen)
        records["selection_id"] = chosen.record()
        selection_id = chosen.value
    request = _build(
        VisualizationRequest,
        {
            "run_dir": resolved_run, "space": space, "canonical_id": resolved_canonical,
            "selection_id": selection_id, "use_selected_landscape": use_selected_landscape,
            "resolved_by": records, **options,
        },
        what="visualize",
    )
    return _visualize(request)


def select(
    run_dir: str | Path | None = None,
    *,
    selection_id: str,
    mode: str = "radius",
    space: str = "raw",
    canonical_id: str | None = None,
    **options: Any,
):
    """Create a named scientific Selection; returns ``SelectResult``.

    ``selection_id`` is always required (a selection is a decision the user
    names). Mode options use the CLI names with underscores, e.g.
    ``center=(0, 0, 0), radius=15`` or ``sld_min=2``.
    """

    from cryorole.select import SelectRequest, create_selection

    records: dict[str, dict[str, str]] = {}
    resolved_run = _run_dir(run_dir, records)
    resolved_canonical = _canonical(resolved_run, space, canonical_id, records)
    if space == "raw" and canonical_frame_ids(resolved_run):
        notify_user(
            "[cryorole] note: selecting in raw space; this bundle also has canonical frame(s) "
            f"{', '.join(canonical_frame_ids(resolved_run))}."
        )
    if "center" in options and options["center"] is not None:
        options["center"] = tuple(float(value) for value in options["center"])
    if "range_bound" in options:
        options["range_bound"] = tuple(options["range_bound"] or ())
    explicit = tuple(options)
    request = _build(
        SelectRequest,
        {
            "run_dir": resolved_run, "selection_id": selection_id, "selection_mode": mode,
            "space": space, "canonical_id": resolved_canonical, "explicit_options": explicit,
            "resolved_by": records, **options,
        },
        what="select",
    )
    return create_selection(request)


def export(
    run_dir: str | Path | None = None,
    *,
    selection_id: str | None = None,
    domain: str = "both",
    format: str = "auto",
    output_dir: str | Path | None = None,
    overwrite: bool = False,
    relocated_ref: str | Path | None = None,
    relocated_mov: str | Path | None = None,
    allow_unverified_source: bool = False,
) -> dict[str, Any]:
    """Write STAR/CS subsets for a selection by recorded source-row IDs; returns the export report."""

    from cryorole.export import export_selection_metadata_subset, read_selection_json
    from cryorole.models.policies import SelectionMetadataExportPolicy

    records: dict[str, dict[str, str]] = {}
    resolved_run = _run_dir(run_dir, records)
    chosen = resolve_selection_id(resolved_run, selection_id)
    _announce(chosen)
    records["selection_id"] = chosen.record()
    selection_path = Path(resolved_run) / "selections" / chosen.value / "selection.json"
    if not selection_path.is_file():
        raise CryoroleError(
            "input_not_found",
            f"Selection JSON file does not exist: {selection_path}.",
            "Check the selection id (cryorole.api.status lists selections).",
        )
    return export_selection_metadata_subset(
        read_selection_json(selection_path),
        policy=SelectionMetadataExportPolicy(
            output_dir=Path(output_dir) if output_dir else None,
            overwrite=overwrite,
            domain=domain,
            format=format,
            run_dir=resolved_run,
            relocated_ref=relocated_ref,
            relocated_mov=relocated_mov,
            allow_unverified_source=allow_unverified_source,
            resolved_by=records,
        ),
        selection_path=selection_path,
    )


# --- RELION alignment --------------------------------------------------------


def align(ref: str | Path, mov: str | Path, **options: Any) -> dict[str, Any]:
    """Align two RELION STAR files exactly (see ``cryorole align``); returns the align report.

    Options use the Python names of ``cryorole.align.align_star_files``, e.g.
    ``key_pairs=["rlnImageName:rlnImageOriginalName"]`` or
    ``coordinate_match="recentered-exact", recenter_shift=(-9, 22, -102)``.
    """

    from cryorole.align import align_star_files

    return align_star_files(ref=ref, mov=mov, **options)


def fix_subtract_coordinates(target: str | Path, **options: Any) -> dict[str, Any]:
    """Write a coordinate-corrected copy of a recentred RELION subtraction; returns the report."""

    from cryorole.align.subtract_fix import fix_subtract_coordinates as _fix

    return _fix(target, **options)


# --- state -------------------------------------------------------------------


def status(run_dir: str | Path | None = None) -> dict[str, Any]:
    """Artifact-derived bundle status (``bundle_status`` is completed/legacy/incomplete/failed/missing)."""

    from cryorole.workflow import inspect_run_status

    records: dict[str, dict[str, str]] = {}
    payload = inspect_run_status(_run_dir(run_dir, records))
    payload["resolved_by"] = records
    return payload


def next_actions(run_dir: str | Path | None = None) -> list[dict[str, Any]]:
    """Conservative next commands derived from the bundle's artifacts."""

    from cryorole.workflow import derive_next_actions

    return derive_next_actions(status(run_dir))
