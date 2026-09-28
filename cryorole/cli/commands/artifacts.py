"""CLI adapters for alignment, manifests, and selection exports."""

from __future__ import annotations

from pathlib import Path

from cryorole.export import (
    export_selection,
    export_selection_metadata_subset,
    read_selection_json,
)
from cryorole.models.policies import (
    RunManifestPolicy,
    SelectionExportPolicy,
    SelectionMetadataExportPolicy,
)
from cryorole.workflows.pipeline_runner import PipelineRunner


def manifest_command(args) -> int:
    runner = PipelineRunner()
    result = runner.write_manifest(
        manifest_policy=RunManifestPolicy(
            output_path=args.output,
            overwrite=args.overwrite,
            workflow_name=args.workflow_name,
            command=args.command_string,
            compute_file_hashes=args.compute_file_hashes,
            hash_algorithm=args.hash_algorithm,
            include_selected_particle_keys=args.include_selected_particle_keys,
        )
    )
    print(result.manifest_report.output_path)
    return 0


def align_command(args) -> int:
    import sys

    if getattr(args, "fix_subtract_coordinates", None):
        from cryorole.align.subtract_fix import fix_subtract_coordinates

        if args.ref or args.mov:
            raise ValueError("--fix-subtract-coordinates works on one Subtract job; do not combine it with --ref/--mov")
        report = fix_subtract_coordinates(
            args.fix_subtract_coordinates,
            subtract_input=args.subtract_input,
            center=args.center,
            model_angpix=args.model_angpix,
            micrograph_angpix=args.micrograph_angpix,
            apply_to=args.apply_to,
            output_dir=args.output_dir,
            overwrite=args.overwrite,
            log=lambda message: print(message, file=sys.stderr),
        )
        print(report["outputs"]["corrected_star"])
        return 0
    if not args.ref or not args.mov:
        raise ValueError("cryorole align needs --ref and --mov (or --fix-subtract-coordinates JOB)")
    for name in ("subtract_input", "center", "model_angpix", "apply_to"):
        if getattr(args, name, None) is not None:
            raise ValueError(f"--{name.replace('_', '-')} is only used with --fix-subtract-coordinates")
    from cryorole.align import align_star_files

    report = align_star_files(
        ref=args.ref,
        mov=args.mov,
        align_id=args.align_id,
        key_columns=args.key,
        key_pairs=args.key_pair,
        float_tolerances=dict(args.float_tol or ()),
        path_mode=args.path_mode,
        duplicate_policy=args.duplicate_policy,
        overwrite=args.overwrite,
        output_dir=args.output_dir,
        coordinate_match=args.coordinate_match,
        recenter_shift=args.recenter_shift,
        ref_angpix=args.ref_angpix,
        micrograph_angpix=args.micrograph_angpix,
        via_extraction=args.via_extraction,
    )
    print(report["output_dir"])
    print(
        f"Aligned {report['matched_count']} particles ({report['strategy']}); "
        f"ref only {report['ref_only_count']}, mov only {report['mov_only_count']}, "
        f"ambiguous {report['ambiguous_ref_row_count']}.",
        file=sys.stderr,
    )
    for warning in report.get("warnings", ()):
        print(f"warning: {warning}", file=sys.stderr)
    print(f"Next: {report['next_command']}", file=sys.stderr)
    return 0


def _resolve_export_ids(args) -> None:
    """Fill in --run-dir / --selection-id unless an explicit selection.json path is used."""

    if getattr(args, "selection", None):
        return
    from cryorole.cli.resolution import resolve_cli_ids

    resolve_cli_ids(args, selection="required")


def export_selection_command(args) -> int:
    _resolve_export_ids(args)
    if getattr(args, "run_dir", None):
        return export_metadata_command(args)
    return _export_selection_artifact_command(args)


def export_metadata_command(args) -> int:
    _resolve_export_ids(args)
    if getattr(args, "selection", None) and not getattr(args, "run_dir", None):
        return _export_selection_artifact_command(args)
    _validate_export_selection_inputs(args)
    selection_path = _resolve_export_selection_path(args)
    selection = read_selection_json(selection_path)
    output_dir = Path(args.output_dir) if getattr(args, "output_dir", None) else None
    report = export_selection_metadata_subset(
        selection,
        policy=SelectionMetadataExportPolicy(
            output_dir=output_dir,
            overwrite=args.overwrite,
            domain=args.domain,
            format=args.format,
            run_dir=args.run_dir,
            relocated_ref=getattr(args, "relocated_ref", None),
            relocated_mov=getattr(args, "relocated_mov", None),
            allow_unverified_source=bool(getattr(args, "allow_unverified_source", False)),
            resolved_by=getattr(args, "resolved_by", None),
        ),
        selection_path=selection_path,
    )
    outputs = report.get("outputs", {})
    selected_keys_path = outputs.get("selected_particle_keys_txt")
    if selected_keys_path:
        print(str(Path(selected_keys_path).parent))
    else:
        print(str(output_dir or Path(args.run_dir) / "exports" / selection.selection_id))
    return 0


def _export_selection_artifact_command(args) -> int:
    _validate_export_selection_inputs(args)
    selection = read_selection_json(args.selection)
    output_dir = _prepare_output_dir(
        args.output_dir,
        overwrite=args.overwrite,
        label="Selection export",
    )
    export_selection(
        selection,
        policy=SelectionExportPolicy(output_dir=output_dir, overwrite=args.overwrite),
    )
    print(str(output_dir))
    return 0


def _resolve_export_selection_path(args) -> Path:
    if getattr(args, "selection", None):
        path = Path(args.selection)
        if not path.exists():
            raise ValueError(
                f"Selection JSON file does not exist: {path}. Provide --run-dir RUN "
                "--selection-id ID, or advanced --selection PATH/to/selection.json."
            )
        return path
    path = Path(args.run_dir) / "selections" / args.selection_id / "selection.json"
    if not path.exists():
        raise ValueError(
            f"Selection JSON file does not exist: {path}. Provide --run-dir RUN "
            "--selection-id ID, or advanced --selection PATH/to/selection.json."
        )
    return path


def _validate_export_selection_inputs(args) -> None:
    has_selection = bool(getattr(args, "selection", None))
    has_selection_id = bool(getattr(args, "selection_id", None))
    has_run_dir = bool(getattr(args, "run_dir", None))
    has_output_dir = bool(getattr(args, "output_dir", None))
    compact_hint = (
        "Provide --run-dir RUN --selection-id ID, or advanced "
        "--selection PATH/to/selection.json."
    )
    if has_selection and has_selection_id:
        raise ValueError("Use either --selection or --selection-id, not both.")
    if has_selection:
        if not has_run_dir and not has_output_dir:
            raise ValueError(
                "Direct --selection export requires --output-dir, or use "
                "--run-dir RUN --selection-id ID."
            )
        return
    if not has_run_dir or not has_selection_id:
        raise ValueError(compact_hint)


def _prepare_output_dir(output_dir: str | Path, *, overwrite: bool, label: str) -> Path:
    path = Path(output_dir)
    if path.exists() and not path.is_dir():
        raise ValueError(f"{label} output path exists and is not a directory: {path}")
    if path.exists() and not overwrite:
        raise FileExistsError(f"{label} output directory already exists: {path}")
    path.mkdir(parents=True, exist_ok=True)
    return path
