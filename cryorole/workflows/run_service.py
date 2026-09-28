"""Typed production run workflow service.

This module owns run orchestration and artifact assembly; the CLI only converts
arguments and publishes the returned result.
"""

from __future__ import annotations

from contextlib import contextmanager
from dataclasses import dataclass, fields, replace
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from cryorole.canonicalize.service import canonicalize_run_artifacts
from cryorole.core.density import compute_landscape_density_arrays
from cryorole.core.input_sanity import assess_ro_angles, same_input_file, sanity_warning_lines
from cryorole.core.display_policy import max_displayed_density
from cryorole.core.euler_conventions import RAW_EULER_ANGLE_COLUMNS, resolve_euler_convention
from cryorole.export import (
    write_json_artifact,
    write_landscape_json,
    write_landscape_npz,
    write_landscape_npz_arrays,
    write_raw_landscape_csv,
    write_report_json,
)
from cryorole.io.writers.landscape_store import write_raw_landscape_csv_from_arrays
from cryorole.models.landscape import Landscape
from cryorole.models.landscape_arrays import LandscapeArrays
from cryorole.models.policies import (
    CanonicalizationPolicy,
    DensityPolicy,
    RepresentationPolicy,
    RunManifestPolicy,
)
from cryorole.preflight import PreflightRequest, run_preflight
from cryorole.provenance import SourceIdentityGuard, build_source_identity
from cryorole.run_bundle import RunBundleWriter
from cryorole.visualize import QuickLookRequest, write_quicklook
from cryorole.workflows.input_policy import (
    InputPolicyRequest,
    ResolvedInputPolicies,
    resolve_input_policies,
)
from cryorole.workflows.pipeline_runner import PipelineRunner
from cryorole.workflows.progress import CancelToken, ProgressCallback, ProgressReporter
from cryorole.workflows.run_pipeline import (
    compute_relative_orientation_arrays,
    preflight_and_normalize_matched_arrays,
)
from cryorole.workflows.run_report import RunReportRequest, write_run_report


RUN_PREVIEW_MAX_POINTS = 50_000
RUN_QUICKLOOK_MAX_POINTS = 500_000


@dataclass(frozen=True)
class RunRequest:
    """Typed, frontend-independent request for one production run."""

    ref: str
    mov: str
    ref_domain: str = "ref"
    mov_domain: str = "mov"
    output_dir: str = "cryorole_outputs"
    row_aligned: bool = False
    allow_low_overlap: bool = False
    sld_metric: str = "rotvec_euclidean"
    run_backend: str = "auto"
    identity_mode: str | None = None
    identity_column: tuple[str, ...] = ()
    mapping_file: str | None = None
    resolved_input_policies: ResolvedInputPolicies | None = None
    k_neighbors: int = 50
    euler_convention: str | None = None
    canonicalize: bool = False
    sign_rule: str = "density_weighted_skewness"
    positive_side: str = "low_density_skew"
    manifest_output: str | None = None
    overwrite: bool = False
    no_visualize: bool = False
    quiet: bool = False
    verbose: bool = False
    profile_time: bool = False
    profile_memory: bool = False
    raw_csv_chunk_size: int = 100_000
    density_query_batch_size: int = 100_000
    no_raw_csv: bool = False
    write_debug_json: bool = False
    sld_display_outlier_mode: str = "tail_jump"
    sld_tail_search_fraction: float = 0.01
    sld_tail_jump_factor: float = 5.0
    sld_max_display_outlier_fraction: float = 0.002

    @classmethod
    def from_namespace(cls, namespace: Any) -> "RunRequest":
        """Copy only declared request fields from an argparse-like object."""

        values = vars(namespace)
        return cls(
            **{
                field.name: values[field.name]
                for field in fields(cls)
                if field.name in values
            }
        )


@dataclass(frozen=True)
class RunExecutionContext:
    """Explicit test/reference dependencies for one run execution."""

    runner: PipelineRunner | None = None
    allow_unverified_test_sources: bool = False
    progress: ProgressCallback | None = None
    cancel_token: CancelToken | None = None
    stream_progress: bool = True


@dataclass(frozen=True)
class RunResult:
    """Typed result returned after staging all run artifacts."""

    output_dir: Path
    output_artifacts: dict[str, str]
    resolved_backend: str
    matched_count: int


@dataclass(frozen=True)
class _RunPhase1:
    ref_row_count: int
    mov_row_count: int
    identity_ref: Any
    identity_mov: Any
    match_table: Any
    match_report: Any
    compatibility_result: Any = None


@dataclass(frozen=True)
class _RunPhase2:
    ro_arrays: Any = None
    compatibility_result: Any = None


@dataclass(frozen=True)
class _RunPhase3:
    landscape: Landscape
    arrays: LandscapeArrays | None = None
    compatibility_result: Any = None


def execute_run(
    request: RunRequest,
    *,
    context: RunExecutionContext | None = None,
) -> RunResult:
    """Create, validate, and atomically publish one complete run bundle."""

    context = context or RunExecutionContext()
    final_output_dir = Path(request.output_dir).expanduser().resolve()
    if request.manifest_output is not None:
        requested_manifest = Path(request.manifest_output).expanduser().resolve()
        required_manifest = final_output_dir / "run_manifest.json"
        if requested_manifest != required_manifest:
            raise ValueError(
                "Transactional run bundles require run_manifest.json inside --output-dir; "
                "external --manifest-output paths are not supported"
            )
    with RunBundleWriter(final_output_dir, overwrite=request.overwrite) as bundle_writer:
        staged_request = replace(
            request,
            output_dir=str(bundle_writer.path),
            manifest_output=None,
            overwrite=True,
        )
        result = _run_command_impl(
            staged_request,
            context=context,
            bundle_writer=bundle_writer,
        )
        if context.cancel_token is not None:
            context.cancel_token.raise_if_cancelled("publishing the run bundle")
        bundle_writer.commit()
    return replace(result, output_dir=final_output_dir)


def _input_policy_request_from_args(args) -> InputPolicyRequest:
    return InputPolicyRequest(
        ref=args.ref,
        mov=args.mov,
        row_aligned=bool(getattr(args, "row_aligned", False)),
        allow_low_overlap=bool(getattr(args, "allow_low_overlap", False)),
        identity_mode=getattr(args, "identity_mode", None),
        identity_columns=tuple(getattr(args, "identity_column", ()) or ()),
        mapping_file=getattr(args, "mapping_file", None),
    )


def _preflight_request_from_args(
    args,
    *,
    resolved_input_policies: ResolvedInputPolicies | None = None,
    same_file_policy: str = "block",
) -> PreflightRequest:
    return PreflightRequest(
        ref=args.ref,
        mov=args.mov,
        ref_domain=getattr(args, "ref_domain", "ref"),
        mov_domain=getattr(args, "mov_domain", "mov"),
        output_dir=getattr(args, "output_dir", "cryorole_outputs"),
        row_aligned=bool(getattr(args, "row_aligned", False)),
        allow_low_overlap=bool(getattr(args, "allow_low_overlap", False)),
        sld_metric=getattr(args, "sld_metric", "rotvec_euclidean"),
        k_neighbors=int(getattr(args, "k_neighbors", 50)),
        run_backend=getattr(args, "run_backend", "auto"),
        density_query_batch_size=int(getattr(args, "density_query_batch_size", 25000)),
        raw_csv=not bool(getattr(args, "no_raw_csv", False)),
        visualize=not bool(getattr(args, "no_visualize", False)),
        identity_mode=getattr(args, "identity_mode", None),
        identity_columns=tuple(getattr(args, "identity_column", ()) or ()),
        mapping_file=getattr(args, "mapping_file", None),
        resolved_input_policies=resolved_input_policies,
        same_file_policy=same_file_policy,
        check_environment=False,
    )


def _run_command_impl(
    args: RunRequest,
    *,
    context: RunExecutionContext,
    bundle_writer: RunBundleWriter,
) -> RunResult:
    """Write all run artifacts into an already-created staging bundle."""

    run_backend_resolved = _resolve_run_backend(args.run_backend)
    if context.runner is not None and run_backend_resolved == "array_native":
        run_backend_resolved = "dataframe_compat"
    reporter = ProgressReporter(
        command="run",
        quiet=args.quiet,
        verbose=args.verbose,
        profile_time=args.profile_time,
        profile_memory=args.profile_memory,
        stream=None if context.stream_progress else False,
        callback=context.progress,
        cancel_token=context.cancel_token,
    )
    reporter.sample_memory("start")
    reporter.info(
        (
            f"backend requested={args.run_backend}, resolved={run_backend_resolved}; "
            f"raw_csv_chunk_size={args.raw_csv_chunk_size}; "
            f"density_query_batch_size={args.density_query_batch_size}"
        ),
        verbose_only=True,
    )
    runner = context.runner or PipelineRunner()
    resolved_inputs = args.resolved_input_policies or resolve_input_policies(
        _input_policy_request_from_args(args)
    )
    source_ref = resolved_inputs.source_ref
    source_mov = resolved_inputs.source_mov
    source_identities = _build_run_source_identities(
        args,
        source_ref=source_ref,
        source_mov=source_mov,
        ref_row_count=-1,
        mov_row_count=-1,
        allow_missing=context.allow_unverified_test_sources,
    )
    identity_policy_ref = resolved_inputs.identity_ref
    identity_policy_mov = resolved_inputs.identity_mov
    convention_policy_ref = resolved_inputs.convention_ref
    convention_policy_mov = resolved_inputs.convention_mov
    match_policy = resolved_inputs.match
    density_policy = DensityPolicy(
        sld_metric=args.sld_metric,
        k_neighbors=args.k_neighbors,
        display_outlier_mode=getattr(args, "sld_display_outlier_mode", "tail_jump"),
        tail_search_fraction=getattr(args, "sld_tail_search_fraction", 0.01),
        tail_jump_factor=getattr(args, "sld_tail_jump_factor", 5.0),
        max_display_outlier_fraction=getattr(
            args,
            "sld_max_display_outlier_fraction",
            0.002,
        ),
    )
    euler_metadata = _resolve_run_euler_metadata(args)
    representation_policy = RepresentationPolicy(
        euler_convention=str(euler_metadata["euler_convention"]),
        scipy_euler_sequence=str(euler_metadata["scipy_euler_sequence"]),
        euler_degrees=True,
    )
    output_dir = _prepare_output_dir(args.output_dir, overwrite=args.overwrite, label="Run")
    data_dir = output_dir / "data"
    reports_dir = output_dir / "reports"
    raw_visualization_dir = output_dir / "visualizations" / "quicklook"
    bundle_writer.set_state("matching")
    with reporter.timed_stage(
        "prepare_output_dirs",
        label="preparing run bundle directories",
        stage_number=1,
        stage_total=8,
        path=str(output_dir),
    ):
        for child in (
            data_dir,
            reports_dir,
            output_dir / "visualizations",
            output_dir / "canonical",
            output_dir / "selections",
            output_dir / "exports",
            output_dir / "debug",
        ):
            child.mkdir(parents=True, exist_ok=True)
    with reporter.timed_stage(
        "read_normalize_match",
        label="reading, normalizing, and matching metadata",
        stage_number=2,
        stage_total=8,
        detail=f"ref={args.ref}; mov={args.mov}",
    ) as stage:
        shared_preflight = None
        if context.runner is None:
            shared_preflight = run_preflight(
                _preflight_request_from_args(
                    args,
                    resolved_input_policies=resolved_inputs,
                    same_file_policy="warn",
                ),
                source_identities=source_identities,
                array_preflight_fn=preflight_and_normalize_matched_arrays,
            )
            if shared_preflight.report["readiness"] == "BLOCKED":
                raise ValueError("; ".join(shared_preflight.report["errors"]))
        if run_backend_resolved == "array_native":
            array_preflight = (
                shared_preflight.array_preflight
                if shared_preflight is not None
                else preflight_and_normalize_matched_arrays(
                    args.ref,
                    args.mov,
                    ref_domain=args.ref_domain,
                    mov_domain=args.mov_domain,
                    identity_policy_ref=identity_policy_ref,
                    identity_policy_mov=identity_policy_mov,
                    convention_policy_ref=convention_policy_ref,
                    convention_policy_mov=convention_policy_mov,
                    match_policy=match_policy,
                )
            )
            phase1 = _RunPhase1(
                ref_row_count=array_preflight.ref_row_count,
                mov_row_count=array_preflight.mov_row_count,
                identity_ref=array_preflight.identity_ref,
                identity_mov=array_preflight.identity_mov,
                match_table=array_preflight.match_table,
                match_report=array_preflight.match_report,
            )
        else:
            compatibility_phase1 = runner.run_phase1(
                args.ref,
                args.mov,
                domain_a_name=args.ref_domain,
                domain_b_name=args.mov_domain,
                identity_policy_a=identity_policy_ref,
                identity_policy_b=identity_policy_mov,
                convention_policy_a=convention_policy_ref,
                convention_policy_b=convention_policy_mov,
                match_policy=match_policy,
            )
            phase1 = _RunPhase1(
                ref_row_count=len(compatibility_phase1.pose_a.data),
                mov_row_count=len(compatibility_phase1.pose_b.data),
                identity_ref=getattr(compatibility_phase1, "identity_a", None),
                identity_mov=getattr(compatibility_phase1, "identity_b", None),
                match_table=getattr(compatibility_phase1, "match_table", None),
                match_report=compatibility_phase1.match_report,
                compatibility_result=compatibility_phase1,
            )
        stage["ref_row_count"] = phase1.ref_row_count
        stage["mov_row_count"] = phase1.mov_row_count
        stage["matched_count"] = phase1.match_report.matched_count
    _assert_run_sources_unchanged(source_identities)
    reporter.sample_memory("source_loaded")
    reporter.info(
        (
            f"matched particles: {phase1.match_report.matched_count} / "
            f"ref={phase1.ref_row_count} / mov={phase1.mov_row_count}"
        )
    )
    for warning in getattr(phase1.match_report, "warnings", ()):
        reporter.warning(str(warning))
    bundle_writer.set_state("computing_ro")
    with reporter.timed_stage(
        "compute_ro",
        label="computing relative orientations",
        stage_number=3,
        stage_total=8,
    ):
        if run_backend_resolved == "array_native":
            ro_arrays = compute_relative_orientation_arrays(
                array_preflight.matched_poses,
                euler_sequence=representation_policy.scipy_euler_sequence,
            )
            phase2 = _RunPhase2(ro_arrays=ro_arrays)
        else:
            compatibility_phase2 = runner.compute_ro_for_phase1_result(
                phase1.compatibility_result,
                euler_sequence=representation_policy.scipy_euler_sequence,
            )
            phase2 = _RunPhase2(compatibility_result=compatibility_phase2)
    reporter.sample_memory("ro_computed")
    input_sanity = assess_ro_angles(
        _ro_angles_for_phase2(phase2),
        same_file=same_input_file(
            source_identities.get("ref"),
            source_identities.get("mov"),
            ref_path=args.ref,
            mov_path=args.mov,
        ),
    )
    reporter.info(str(input_sanity["ro_angle_summary"]["message"]))
    for warning in sanity_warning_lines(input_sanity):
        reporter.warning(warning)
    bundle_writer.set_state("computing_sld")
    raw_arrays: LandscapeArrays | None = None
    with reporter.timed_stage(
        "compute_sld",
        label=(
            f"computing SLD with k={density_policy.k_neighbors}, "
            f"metric={density_policy.sld_metric}"
        ),
        stage_number=4,
        stage_total=8,
        density_query_batch_size=args.density_query_batch_size,
    ):
        if run_backend_resolved == "array_native":
            with _explain_collapsed_sld(input_sanity):
                raw_arrays, density_report = compute_landscape_density_arrays(
                    ro_arrays.particle_key,
                    ro_arrays.rotation_vector,
                    ro_arrays.ref_source_row_id,
                    ro_arrays.mov_source_row_id,
                    policy=density_policy,
                    query_batch_size=args.density_query_batch_size,
                )
            raw_landscape = _landscape_from_run_arrays(
                raw_arrays,
                density_report=density_report,
                density_policy=density_policy,
                max_points=RUN_PREVIEW_MAX_POINTS,
            )
            phase3 = _RunPhase3(landscape=raw_landscape, arrays=raw_arrays)
        else:
            with _explain_collapsed_sld(input_sanity):
                compatibility_phase3 = runner.compute_density_for_ro_result(
                    phase2.compatibility_result,
                    density_policy=density_policy,
                    density_query_batch_size=args.density_query_batch_size,
                )
            phase3 = _RunPhase3(
                landscape=compatibility_phase3.landscape,
                compatibility_result=compatibility_phase3,
            )
            density_report = phase3.landscape.density_report
    if density_report is not None:
        for warning in density_report.warnings:
            reporter.warning(warning)
    reporter.sample_memory("density_computed")
    match_table_data = _match_table_data_for_phase1(phase1, phase3.landscape)
    if run_backend_resolved != "array_native":
        raw_landscape = _landscape_with_match_rows(phase3.landscape, match_table_data)
    raw_active_policies = dict(raw_landscape.active_policies or {})
    raw_active_policies["representation_policy"] = representation_policy
    raw_landscape.active_policies = raw_active_policies
    bundle_writer.set_state("writing")
    with reporter.timed_stage(
        "write_npz",
        label="writing raw landscape NPZ",
        stage_number=5,
        stage_total=8,
        row_count=len(raw_landscape.data),
    ):
        if run_backend_resolved == "array_native":
            raw_landscape_npz_path = write_landscape_npz_arrays(
                raw_arrays,
                data_dir / "raw_landscape.npz",
                overwrite=args.overwrite,
                artifact_type="raw_landscape",
            )
        else:
            raw_landscape_npz_path = write_landscape_npz(
                raw_landscape,
                data_dir / "raw_landscape.npz",
                overwrite=args.overwrite,
                artifact_type="raw_landscape",
            )
    reporter.sample_memory("raw_npz_written")
    raw_landscape_csv_path = None
    if args.no_raw_csv:
        reporter.stage_skipped(
            "write_raw_csv",
            label="raw CSV skipped by --no-raw-csv",
            stage_number=6,
            stage_total=8,
            reason="--no-raw-csv",
        )
    else:
        with reporter.timed_stage(
            "write_raw_csv",
            label="writing raw landscape CSV",
            stage_number=6,
            stage_total=8,
            row_count=len(raw_landscape.data),
            raw_csv_chunk_size=args.raw_csv_chunk_size,
        ):
            if run_backend_resolved == "array_native":
                raw_landscape_csv_path = write_raw_landscape_csv_from_arrays(
                    raw_arrays,
                    data_dir / "raw_landscape.csv",
                    overwrite=args.overwrite,
                    chunk_size=args.raw_csv_chunk_size,
                    euler_sequence=representation_policy.scipy_euler_sequence,
                )
            else:
                raw_landscape_csv_path = write_raw_landscape_csv(
                    raw_landscape,
                    data_dir / "raw_landscape.csv",
                    overwrite=args.overwrite,
                    euler_sequence=representation_policy.scipy_euler_sequence,
                )
        reporter.sample_memory("raw_csv_written")
    with reporter.timed_stage(
        "write_core_reports",
        label="writing core reports",
        row_count=len(raw_landscape.data),
    ):
        match_table_path = _write_csv_artifact(
            match_table_data,
            data_dir / "match_table.csv",
            overwrite=args.overwrite,
        )
        density_report_path = write_report_json(
            density_report,
            reports_dir / "density_report.json",
            overwrite=args.overwrite,
            artifact_type="density_report",
        )
        match_report_path = write_report_json(
            phase1.match_report,
            reports_dir / "match_report.json",
            overwrite=args.overwrite,
            artifact_type="match_report",
        )
        identity_ref_report_path = write_report_json(
            _identity_report_for_phase1(phase1, "ref"),
            reports_dir / "identity_ref_report.json",
            overwrite=args.overwrite,
            artifact_type="identity_report",
        )
        identity_mov_report_path = write_report_json(
            _identity_report_for_phase1(phase1, "mov"),
            reports_dir / "identity_mov_report.json",
            overwrite=args.overwrite,
            artifact_type="identity_report",
        )
        import_ref_report_path = write_json_artifact(
            _simple_import_report(args.ref, source_ref, phase1.ref_row_count),
            reports_dir / "import_ref_report.json",
            overwrite=args.overwrite,
        )
        import_mov_report_path = write_json_artifact(
            _simple_import_report(args.mov, source_mov, phase1.mov_row_count),
            reports_dir / "import_mov_report.json",
            overwrite=args.overwrite,
        )
    source_identities["ref"]["row_count"] = phase1.ref_row_count
    source_identities["mov"]["row_count"] = phase1.mov_row_count
    output_artifacts = {
        "raw_landscape_npz": str(raw_landscape_npz_path),
        "match_table_csv": str(match_table_path),
        "density_report_json": str(density_report_path),
        "match_report_json": str(match_report_path),
        "identity_ref_report_json": str(identity_ref_report_path),
        "identity_mov_report_json": str(identity_mov_report_path),
        "import_ref_report_json": str(import_ref_report_path),
        "import_mov_report_json": str(import_mov_report_path),
    }
    if raw_landscape_csv_path is not None:
        output_artifacts["raw_landscape_csv"] = str(raw_landscape_csv_path)
    if args.write_debug_json:
        debug_landscape = (
            _landscape_from_run_arrays(
                raw_arrays,
                density_report=density_report,
                density_policy=density_policy,
                max_points=None,
            )
            if run_backend_resolved == "array_native"
            else raw_landscape
        )
        landscape_debug_path = write_landscape_json(
            debug_landscape,
            output_dir / "debug" / "landscape_debug.json",
            overwrite=args.overwrite,
            artifact_type="landscape",
        )
        output_artifacts["landscape_debug_json"] = str(landscape_debug_path)
    canonicalization_policy = None
    canonical_landscape = None
    if args.canonicalize:
        with reporter.timed_stage("canonicalize", label="canonicalizing landscape"):
            canonicalization_policy = CanonicalizationPolicy(
                sign_rule=args.sign_rule,
                positive_side=args.positive_side,
            )
            canonical_result = canonicalize_run_artifacts(
                output_dir=output_dir / "canonical" / "default",
                run_backend=run_backend_resolved,
                policy=canonicalization_policy,
                overwrite=args.overwrite,
                csv_chunk_size=args.raw_csv_chunk_size,
                euler_sequence=representation_policy.scipy_euler_sequence,
                raw_arrays=raw_arrays,
                compatibility_result=phase3.compatibility_result,
                runner=runner,
                density_report=density_report,
                density_policy=density_policy,
                preview_max_points=RUN_PREVIEW_MAX_POINTS,
            )
            canonical_landscape = canonical_result.landscape
            canonical_landscape_path = canonical_result.landscape_npz_path
            canonical_landscape_csv_path = canonical_result.landscape_csv_path
            canonicalization_report_path = canonical_result.report_path
        output_artifacts["canonical_landscape_npz"] = str(canonical_landscape_path)
        output_artifacts["canonical_landscape_csv"] = str(canonical_landscape_csv_path)
        output_artifacts["canonicalization_report_json"] = str(canonicalization_report_path)
        if args.write_debug_json:
            canonical_debug_path = write_landscape_json(
                canonical_landscape,
                output_dir / "debug" / "canonical_default_landscape_debug.json",
                overwrite=args.overwrite,
                artifact_type="canonical_landscape",
            )
            output_artifacts["canonical_landscape_debug_json"] = str(canonical_debug_path)

    raw_visualization_performed = False
    quicklook_report: dict[str, Any] | None = None
    if args.no_visualize:
        reporter.stage_skipped(
            "visualize_raw",
            label="raw visualization skipped by --no-visualize",
            stage_number=7,
            stage_total=8,
            reason="--no-visualize",
        )
    else:
        with reporter.timed_stage(
            "visualize_raw",
            label="writing raw quick-look figures",
            stage_number=7,
            stage_total=8,
        ):
            quicklook_report = _write_run_quicklook(
                args=args,
                raw_landscape=(raw_arrays if raw_arrays is not None else raw_landscape),
                raw_visualization_dir=raw_visualization_dir,
                representation_policy=representation_policy,
                euler_metadata=euler_metadata,
                density_report=density_report,
            )
        raw_visualization_performed = True
        for artifact_key, artifact_path in quicklook_report["generated_files"].items():
            output_artifacts[f"raw_quicklook_{artifact_key}"] = artifact_path

    manifest_path = Path(args.manifest_output) if args.manifest_output else output_dir / "run_manifest.json"
    run_summary_path = output_dir / "run_summary.json"
    run_report_path = output_dir / "run_report.md"
    timing_profile_path = reports_dir / "run_timing_profile.json"
    memory_profile_path = reports_dir / "run_memory_profile.json"
    if args.profile_time:
        output_artifacts["run_timing_profile_json"] = str(timing_profile_path)
    if args.profile_memory:
        output_artifacts["run_memory_profile_json"] = str(memory_profile_path)
    output_artifacts["run_summary_json"] = str(run_summary_path)
    output_artifacts["run_manifest_json"] = str(manifest_path)
    output_artifacts["run_report_md"] = str(run_report_path)
    output_artifacts = bundle_writer.publicize(output_artifacts)
    active_policies = {
        "identity_policy_ref": identity_policy_ref,
        "identity_policy_mov": identity_policy_mov,
        "convention_policy_ref": convention_policy_ref,
        "convention_policy_mov": convention_policy_mov,
        "match_policy": match_policy,
        "density_policy": density_policy,
        "representation_policy": representation_policy,
    }
    if canonicalization_policy is not None:
        active_policies["canonicalization_policy"] = canonicalization_policy
    with reporter.timed_stage(
        "manifest",
        label="writing run summary and manifest",
        stage_number=8,
        stage_total=8,
    ):
        run_summary = _run_summary_payload(
                args=args,
                source_ref=source_ref,
                source_mov=source_mov,
                phase1=phase1,
                landscape=phase3.landscape,
                matched_count=phase1.match_report.matched_count,
                k_neighbors=density_policy.k_neighbors,
                canonicalization_performed=args.canonicalize,
                output_artifacts=output_artifacts,
                run_backend_resolved=run_backend_resolved,
                raw_csv_performed=raw_landscape_csv_path is not None,
                raw_csv_backend=(
                    (
                        "array_native_chunked"
                        if run_backend_resolved == "array_native"
                        else "dataframe_compat_full_table"
                    )
                    if raw_landscape_csv_path is not None else "skipped"
                ),
                raw_visualization_performed=raw_visualization_performed,
                quicklook_report=quicklook_report,
                timing_profile_path=timing_profile_path if args.profile_time else None,
                memory_profile_path=memory_profile_path if args.profile_memory else None,
                euler_metadata=euler_metadata,
                run_id=bundle_writer.run_id,
                source_identities=source_identities,
                landscape_row_count=(
                    raw_arrays.n_points
                    if run_backend_resolved == "array_native"
                    else len(raw_landscape.data)
                ),
            )
        run_summary["input_sanity"] = input_sanity
        alignment_provenance = _alignment_provenance_for_run(args, phase1, source_identities)
        run_summary["alignment_provenance"] = alignment_provenance
        if alignment_provenance is not None:
            if not alignment_provenance["attached"]:
                reporter.warning(f"alignment provenance not attached: {alignment_provenance['reason']}")
            else:
                for domain, state in alignment_provenance["original_files"].items():
                    if state["status"] == "mismatch":
                        reporter.warning(
                            f"alignment provenance: the original {domain} file at {state['path']} has changed since "
                            "cryorole align (recorded lineage kept; current file differs)"
                        )
        write_json_artifact(
            run_summary,
            run_summary_path,
            overwrite=args.overwrite,
        )
        write_run_report(
            RunReportRequest(
                output_path=run_report_path,
                summary=run_summary,
                density_report=density_report,
                overwrite=args.overwrite,
            )
        )
        runner.write_manifest(
            manifest_policy=RunManifestPolicy(
                output_path=manifest_path,
                overwrite=args.overwrite,
                workflow_name="run",
                schema_version="3.0",
            ),
            input_paths=(args.ref, args.mov),
            source_types={
                str(args.ref): source_ref,
                str(args.mov): source_mov,
            },
            row_counts={
                str(args.ref): phase1.ref_row_count,
                str(args.mov): phase1.mov_row_count,
            },
            active_policies=active_policies,
            match_report=phase1.match_report,
            density_report=density_report,
            canonicalization_report=(
                canonical_landscape.canonicalization_report
                if canonical_landscape is not None
                else None
            ),
            additional_results={
                "euler_metadata": euler_metadata,
                "alignment_provenance": run_summary.get("alignment_provenance"),
            },
            output_artifacts=output_artifacts,
            run_id=bundle_writer.run_id,
            source_identities=source_identities,
        )
    if args.profile_time:
        write_json_artifact(
            reporter.timing_profile_payload(),
            timing_profile_path,
            overwrite=args.overwrite,
        )
    if args.profile_memory:
        reporter.sample_memory("completed")
        write_json_artifact(
            reporter.memory_profile_payload(),
            memory_profile_path,
            overwrite=args.overwrite,
        )
    reporter.info(f"staging completed successfully: {output_dir}")
    return RunResult(
        output_dir=output_dir,
        output_artifacts=dict(output_artifacts),
        resolved_backend=run_backend_resolved,
        matched_count=int(phase1.match_report.matched_count),
    )


def _alignment_provenance_for_run(args, phase1, source_identities) -> dict[str, Any] | None:
    """Verified ``cryorole align`` lineage for ``--row-aligned`` runs (never blocks the run)."""

    if not getattr(args, "row_aligned", False):
        return None
    from cryorole.align.provenance import discover_alignment_provenance

    try:
        return discover_alignment_provenance(
            args.ref,
            args.mov,
            ref_sha256=(source_identities.get("ref") or {}).get("sha256"),
            mov_sha256=(source_identities.get("mov") or {}).get("sha256"),
            ref_row_count=phase1.ref_row_count,
            mov_row_count=phase1.mov_row_count,
        )
    except (OSError, ValueError) as exc:  # provenance must never stop an explicit --row-aligned run
        return {"attached": False, "reason": f"could not read the align report ({exc})"}


def _ro_angles_for_phase2(phase2) -> np.ndarray | None:
    """Return per-particle RO angles (rad) for either run backend."""

    if phase2.ro_arrays is not None:
        return np.asarray(phase2.ro_arrays.angle_rad, dtype=float)
    ro_result = getattr(phase2.compatibility_result, "ro_result", None)
    data = getattr(ro_result, "data", None)
    if data is None or "rotvec_ro" not in getattr(data, "columns", ()):
        # Injected test runners may not expose an ROResult.
        return None
    if len(data) == 0:
        return np.empty(0, dtype=float)
    rotvecs = np.stack(data["rotvec_ro"].to_numpy())
    return np.linalg.norm(np.asarray(rotvecs, dtype=float), axis=1)


@contextmanager
def _explain_collapsed_sld(input_sanity: dict[str, Any]):
    """Turn the undefined-SLD error into the input-sanity explanation."""

    try:
        yield
    except ValueError as exc:
        lines = sanity_warning_lines(input_sanity)
        if "SLD is undefined" not in str(exc) or not lines:
            raise
        raise ValueError(
            " ".join(lines)
            + " cryoROLE stopped before writing a landscape: every particle has the same relative "
            "orientation, so the local density (SLD) is undefined."
        ) from exc


def _resolve_run_euler_metadata(args) -> dict[str, object]:
    source = "cli_override" if args.euler_convention else "cli_default"
    resolved = resolve_euler_convention(args.euler_convention, source=source)
    return resolved.metadata(euler_angle_columns=RAW_EULER_ANGLE_COLUMNS)


def _run_summary_payload(
    *,
    args,
    source_ref: str,
    source_mov: str,
    phase1,
    landscape,
    matched_count: int,
    k_neighbors: int,
    canonicalization_performed: bool,
    output_artifacts: dict[str, str],
    run_backend_resolved: str,
    raw_csv_performed: bool,
    raw_csv_backend: str,
    raw_visualization_performed: bool,
    quicklook_report: dict[str, Any] | None,
    timing_profile_path: Path | None,
    memory_profile_path: Path | None,
    euler_metadata: dict[str, object],
    run_id: str,
    source_identities: dict[str, dict[str, object]],
    landscape_row_count: int,
) -> dict[str, Any]:
    payload = {
        "artifact_type": "run_summary",
        "schema_version": "1",
        "run_id": run_id,
        "input_paths": {"ref": args.ref, "mov": args.mov},
        "source_identities": source_identities,
        "source_types": {"ref": source_ref, "mov": source_mov},
        "ref_domain": args.ref_domain,
        "mov_domain": args.mov_domain,
        "row_counts": {
            "ref": phase1.ref_row_count,
            "mov": phase1.mov_row_count,
        },
        "matched_count": int(matched_count),
        "match_key": phase1.match_report.match_key,
        "ref_row_count": phase1.match_report.ref_row_count or phase1.ref_row_count,
        "mov_row_count": phase1.match_report.mov_row_count or phase1.mov_row_count,
        "matched_row_count": phase1.match_report.matched_row_count or matched_count,
        "dropped_ref_only_count": phase1.match_report.dropped_ref_only_count,
        "dropped_mov_only_count": phase1.match_report.dropped_mov_only_count,
        "matched_rows_reordered": phase1.match_report.matched_rows_reordered,
        "ref_coverage": phase1.match_report.ref_coverage,
        "mov_coverage": phase1.match_report.mov_coverage,
        "overlap_smaller_input": (
            phase1.match_report.overlap_smaller_input or phase1.match_report.overlap_ratio
        ),
        "low_overlap_allowed": phase1.match_report.low_overlap_allowed,
        "match_warnings": list(phase1.match_report.warnings),
        "ro_coordinate_diagnostics": getattr(getattr(landscape, "density_report", None), "ro_coordinate_diagnostics", None),
        "landscape_row_count": int(landscape_row_count),
        "k_neighbors": k_neighbors,
        "requested_sld_metric": args.sld_metric,
        "resolved_sld_metric": args.sld_metric,
        "canonicalization_performed": canonicalization_performed,
        "selection_performed": False,
        "run_backend_requested": args.run_backend,
        "run_backend_resolved": run_backend_resolved,
        "raw_landscape_npz": output_artifacts.get("raw_landscape_npz"),
        "raw_csv_performed": raw_csv_performed,
        "raw_csv_backend": raw_csv_backend,
        "raw_csv_chunk_size": args.raw_csv_chunk_size,
        "raw_visualization_performed": raw_visualization_performed,
        "raw_visualization_style": "legacy_rainbow",
        "raw_visualization_profile": "quicklook",
        "raw_visualization_files": (
            [
                f"visualizations/quicklook/{filename}"
                for filename in quicklook_report.get("generated_filenames", ())
            ]
            if quicklook_report is not None else []
        ),
        "raw_visualization_representation": "euler_and_rotvec_3view",
        "raw_visualization_filter": "all_sld_ge_1_top_40pct",
        "raw_visualization_display_vmax_cap": 100.0,
        "raw_visualization_max_displayed_sld": (
            quicklook_report.get("max_displayed_sld")
            if quicklook_report is not None else None
        ),
        "raw_visualization_resolved_vmax": (
            quicklook_report.get("resolved_vmax")
            if quicklook_report is not None else None
        ),
        "raw_visualization_rendered_points": (
            (quicklook_report.get("n_points_2d") or {}).get("analysis")
            if quicklook_report is not None else None
        ),
        "raw_visualization_sampling_method": (
            quicklook_report.get("sampling_method")
            if quicklook_report is not None else None
        ),
        "raw_visualization_subset_counts": (
            quicklook_report.get("subset_counts") if quicklook_report is not None else None
        ),
        "raw_visualization_rendered_subset_counts": (
            quicklook_report.get("rendered_subset_counts")
            if quicklook_report is not None else None
        ),
        "raw_visualization_top_40pct_cutoff_sld": (
            quicklook_report.get("top_40pct_cutoff_sld")
            if quicklook_report is not None else None
        ),
        "raw_visualization_top_40pct_tie_policy": (
            quicklook_report.get("top_40pct_tie_policy")
            if quicklook_report is not None else None
        ),
        "raw_visualization_sld_distribution": (
            quicklook_report.get("sld_distribution")
            if quicklook_report is not None else None
        ),
        "density_backend": (
            "array_native_batched"
            if run_backend_resolved == "array_native"
            else "current_dataframe_compat"
        ),
        "density_query_batch_size": args.density_query_batch_size,
        "timing_profile_performed": args.profile_time,
        "timing_profile": str(timing_profile_path) if timing_profile_path is not None else None,
        "memory_profile_performed": args.profile_memory,
        "memory_profile": str(memory_profile_path) if memory_profile_path is not None else None,
        "display_filtering": {
            "display_filter_mode": "fixed_quicklook_ranges",
            "display_top_fraction": 0.40,
            "display_sld_threshold": 1.0,
            "display_density_field": "sld_raw",
            "display_ranges": ["all", "sld_raw >= 1", "top 40% by sld_raw"],
            "top_fraction_tie_policy": "include_all_rows_at_or_above_cutoff",
            "display_max_divisor": None,
            "visual_style": "legacy",
            "sld_display_mode": "identity",
            "sld_display_outlier_mode": getattr(
                args,
                "sld_display_outlier_mode",
                "tail_jump",
            ),
            "sld_tail_search_fraction": getattr(args, "sld_tail_search_fraction", 0.01),
            "sld_tail_jump_factor": getattr(args, "sld_tail_jump_factor", 5.0),
            "sld_max_display_outlier_fraction": getattr(
                args,
                "sld_max_display_outlier_fraction",
                0.002,
            ),
        },
        "output_artifacts": dict(output_artifacts),
    }
    payload.update(euler_metadata)
    return payload


def _landscape_with_match_rows(landscape, match_table: Any):
    data = landscape.data.copy(deep=True)
    if isinstance(match_table, pd.DataFrame):
        match_frame = match_table
    else:
        match_frame = pd.DataFrame(match_table)
    if {"particle_key", "domain_a_row", "domain_b_row"}.issubset(match_frame.columns):
        indexed = match_frame.set_index("particle_key")
        data["ref_source_row_id"] = (
            data["particle_key"].map(indexed["domain_a_row"]).fillna(-1).astype(int)
        )
        data["mov_source_row_id"] = (
            data["particle_key"].map(indexed["domain_b_row"]).fillna(-1).astype(int)
        )
    return Landscape(
        data=data,
        canonical_transform=landscape.canonical_transform,
        active_policies=landscape.active_policies,
        density_report=landscape.density_report,
        canonicalization_report=landscape.canonicalization_report,
    )


def _match_table_data_for_phase1(phase1, landscape) -> pd.DataFrame:
    match_table = getattr(phase1, "match_table", None)
    if match_table is not None and hasattr(match_table, "data"):
        return match_table.data
    particle_keys = list(landscape.data["particle_key"])
    return pd.DataFrame(
        {
            "particle_key": particle_keys,
            "domain_a_row": range(len(particle_keys)),
            "domain_b_row": range(len(particle_keys)),
            "match_status": "matched",
        }
    )


def _landscape_from_run_arrays(
    arrays: LandscapeArrays,
    *,
    density_report,
    density_policy: DensityPolicy,
    max_points: int | None,
    canonicalization_report=None,
) -> Landscape:
    """Inflate only a deterministic bounded preview from production arrays."""

    if max_points is None or arrays.n_points <= max_points:
        indices = np.arange(arrays.n_points, dtype=np.int64)
    else:
        indices = np.linspace(0, arrays.n_points - 1, num=max_points, dtype=np.int64)
        indices = np.unique(indices)
    data = pd.DataFrame(
        {
            "particle_key": arrays.particle_key[indices],
            "coordinates_analysis": list(arrays.coordinates_analysis[indices]),
            "coordinates_display": list(arrays.coordinates_analysis[indices].copy()),
            "sld_unfloored": arrays.sld_unfloored[indices],
            "sld_raw": arrays.sld_raw[indices],
            "sld_display": arrays.sld_display[indices],
            "sld_display_is_outlier": arrays.sld_display_is_outlier[indices],
            "sld_was_floored": arrays.sld_was_floored[indices],
            "sld_local_k_mean": arrays.sld_local_k_mean[indices],
            "sld_effective_local_k_mean": arrays.sld_effective_local_k_mean[indices],
            "sld_distance_floor": arrays.sld_distance_floor[indices],
            "ref_source_row_id": arrays.ref_source_row_id[indices],
            "mov_source_row_id": arrays.mov_source_row_id[indices],
        }
    )
    if arrays.coordinates_canonical is not None:
        data["coordinates_canonical"] = list(arrays.coordinates_canonical[indices])
    preview_report = density_report if len(indices) == arrays.n_points else None
    return Landscape(
        data=data,
        canonical_transform=arrays.canonical_transform,
        active_policies={"density_policy": density_policy},
        density_report=preview_report,
        canonicalization_report=canonicalization_report,
    )


def _identity_report_for_phase1(phase1: _RunPhase1, label: str):
    identity = phase1.identity_ref if label == "ref" else phase1.identity_mov
    if identity is not None and hasattr(identity, "report"):
        return identity.report
    return {
        "identity_mode": "unavailable_in_test_runner",
        "identity_columns": (),
        "column_normalization_rules": {},
        "unique_rate": 1.0,
        "duplicate_count": 0,
        "collision_examples": (),
        "status": f"not_reported_for_{label}",
    }


def _write_csv_artifact(data: pd.DataFrame, path: Path, *, overwrite: bool) -> Path:
    if path.exists() and not overwrite:
        raise FileExistsError(f"Output path already exists: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    data.to_csv(path, index=False)
    return path


def _simple_import_report(path: str, source_type: str, row_count: int) -> dict[str, Any]:
    return {
        "artifact_type": "import_report",
        "schema_version": "1",
        "input_path": path,
        "source_type": source_type,
        "row_count": int(row_count),
        "status": "ok",
    }


def _build_run_source_identities(
    args: RunRequest,
    *,
    source_ref: str,
    source_mov: str,
    ref_row_count: int,
    mov_row_count: int,
    allow_missing: bool,
) -> dict[str, dict[str, object]]:
    records: dict[str, dict[str, object]] = {}
    for domain, path, source_type, row_count in (
        ("ref", args.ref, source_ref, ref_row_count),
        ("mov", args.mov, source_mov, mov_row_count),
    ):
        try:
            records[domain] = build_source_identity(
                path,
                source_type=source_type,
                row_count=row_count,
            ).to_dict()
        except FileNotFoundError:
            if not allow_missing:
                raise
            records[domain] = {
                "schema_version": "1",
                "original_path": str(path),
                "resolved_path": str(Path(path).expanduser().resolve()),
                "source_type": source_type,
                "size_bytes": None,
                "mtime_epoch_sec": None,
                "mtime_ns": None,
                "sha256": None,
                "row_count": int(row_count),
                "unverified_test_source": True,
            }
    return records


def _assert_run_sources_unchanged(source_identities: dict[str, dict[str, object]]) -> None:
    """Fail if any verified source content changed while run was reading it."""

    verified = {
        domain: record
        for domain, record in source_identities.items()
        if not record.get("unverified_test_source")
    }
    if verified:
        SourceIdentityGuard(verified).assert_unchanged()


def _write_run_quicklook(
    *,
    args,
    raw_landscape: Landscape | LandscapeArrays,
    raw_visualization_dir: Path,
    representation_policy: RepresentationPolicy,
    euler_metadata: dict[str, object],
    density_report,
) -> dict[str, Any]:
    preview_values = (
        raw_landscape.sld_raw
        if isinstance(raw_landscape, LandscapeArrays)
        else raw_landscape.data["sld_raw"]
    )
    preview_max_sld = max_displayed_density(preview_values)
    reported_max_sld = (
        float(density_report.max_sld_raw) if density_report is not None else None
    )
    finite_maxima = [
        value
        for value in (preview_max_sld, reported_max_sld)
        if value is not None and np.isfinite(value)
    ]
    max_displayed_sld = max(finite_maxima) if finite_maxima else None
    resolved_vmax = (
        min(max_displayed_sld, 100.0) if max_displayed_sld is not None else None
    )
    result = write_quicklook(
        raw_landscape,
        QuickLookRequest(
            output_dir=raw_visualization_dir,
            overwrite=args.overwrite,
            euler_convention=representation_policy.euler_convention,
            euler_convention_source=str(euler_metadata["euler_convention_source"]),
            color_field="sld_raw",
            display_density_field="sld_raw",
            color_map="rainbow_r",
            color_vmax=resolved_vmax,
            max_points_2d=RUN_QUICKLOOK_MAX_POINTS,
            random_seed=0,
            tail_jump_threshold=(
                density_report.sld_display_outlier_threshold
                if density_report is not None else None
            ),
        ),
    )
    return {
        **result.report,
        "artifact_type": "run_quicklook",
        "schema_version": "1",
        "quicklook_only": True,
        "not_a_selection": True,
        "visual_style": "legacy_rainbow",
        "color_map": "rainbow_r",
        "color_field": "sld_raw",
        "display_filter_mode": "all_sld_ge_1_top_40pct",
        "display_vmax_cap": 100.0,
        "display_vmax_policy": "min(max_displayed_sld, 100)",
        "max_displayed_sld": max_displayed_sld,
        "resolved_vmax": resolved_vmax,
    }


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
        raise FileExistsError(f"{label} output directory already exists: {path}")
    path.mkdir(parents=True, exist_ok=True)
    return path


def _resolve_run_backend(requested: str) -> str:
    if requested in {"auto", "array_native"}:
        return "array_native"
    if requested == "dataframe_compat":
        return requested
    raise ValueError(f"Unsupported run backend: {requested}")
