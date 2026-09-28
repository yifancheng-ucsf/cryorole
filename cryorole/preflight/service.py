"""One shared preflight service used by preflight, run, and guide."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
import shlex
from typing import Any

from cryorole.core.input_sanity import same_input_file
from cryorole.export.serialization import to_json_safe
from cryorole.models.policies import ConventionPolicy
from cryorole.preflight.report import PREFLIGHT_SCHEMA_VERSION, PreflightResult
from cryorole.preflight.resource_estimate import estimate_run_resources
from cryorole.provenance import build_source_identity
from cryorole.workflows.run_pipeline import inspect_source_identity, preflight_and_normalize_matched_arrays
from cryorole.workflows.input_policy import (
    InputPolicyRequest,
    ResolvedInputPolicies,
    resolve_input_policies,
)


@dataclass(frozen=True)
class PreflightRequest:
    ref: str | Path
    mov: str | Path
    ref_domain: str = "ref"
    mov_domain: str = "mov"
    output_dir: str | Path = "cryorole_outputs"
    row_aligned: bool = False
    allow_low_overlap: bool = False
    sld_metric: str = "rotvec_euclidean"
    k_neighbors: int = 50
    run_backend: str = "auto"
    density_query_batch_size: int = 25_000
    raw_csv: bool = True
    visualize: bool = True
    identity_mode: str | None = None
    identity_columns: tuple[str, ...] = ()
    mapping_file: str | Path | None = None
    resolved_input_policies: ResolvedInputPolicies | None = None
    # "block": same ref/mov file is a preflight error (public ``preflight`` and
    # ``run --dry-run``). "warn": recorded as a strong warning (production run).
    same_file_policy: str = "block"


def run_preflight(
    request: PreflightRequest,
    *,
    source_identities: dict[str, dict[str, object]] | None = None,
    array_preflight_fn=preflight_and_normalize_matched_arrays,
) -> PreflightResult:
    """Inspect inputs without creating a run bundle or scientific artifacts."""

    timestamp = datetime.now(timezone.utc).isoformat()
    warnings: list[str] = []
    errors: list[str] = []
    identities: dict[str, dict[str, object]] = source_identities or {}
    source_types: dict[str, str] = {}
    phase = None
    identity_reports: dict[str, object] = {}
    matching: dict[str, object] = {}
    conventions: dict[str, object] = {}
    input_sanity: dict[str, object] = {"same_file": None}
    resolved = request.resolved_input_policies

    try:
        if resolved is None:
            resolved = resolve_input_policies(
                InputPolicyRequest(
                    ref=request.ref,
                    mov=request.mov,
                    row_aligned=request.row_aligned,
                    allow_low_overlap=request.allow_low_overlap,
                    identity_mode=request.identity_mode,
                    identity_columns=request.identity_columns,
                    mapping_file=request.mapping_file,
                )
            )
        source_types = {"ref": resolved.source_ref, "mov": resolved.source_mov}
        for domain, path in (("ref", request.ref), ("mov", request.mov)):
            if domain not in identities:
                identities[domain] = build_source_identity(
                    path,
                    source_type=source_types[domain],
                    row_count=-1,
                ).to_dict()
        same_file = same_input_file(
            identities.get("ref"),
            identities.get("mov"),
            ref_path=request.ref,
            mov_path=request.mov,
        )
        if same_file is not None:
            input_sanity["same_file"] = same_file
            if request.same_file_policy == "block":
                raise ValueError(same_file["message"])
            warnings.append(f"[{same_file['code']}] {same_file['message']}")
        ref_policy = resolved.identity_ref
        mov_policy = resolved.identity_mov
        conventions = {
            "ref": _convention_record(source_types["ref"]),
            "mov": _convention_record(source_types["mov"]),
            "public_ro_euler": {
                "euler_convention": "extrinsic_zyx",
                "scipy_euler_sequence": "zyx",
                "degrees": True,
            },
        }
        phase = array_preflight_fn(
            request.ref,
            request.mov,
            ref_domain=request.ref_domain,
            mov_domain=request.mov_domain,
            identity_policy_ref=ref_policy,
            identity_policy_mov=mov_policy,
            convention_policy_ref=resolved.convention_ref,
            convention_policy_mov=resolved.convention_mov,
            match_policy=resolved.match,
        )
        identities["ref"]["row_count"] = phase.ref_row_count
        identities["mov"]["row_count"] = phase.mov_row_count
        identity_reports = {
            "ref": to_json_safe(phase.identity_ref.report),
            "mov": to_json_safe(phase.identity_mov.report),
        }
        match = phase.match_report
        matching = {
            "identity_key": match.match_key,
            "matched_count": match.matched_count,
            "ref_row_count": match.ref_row_count,
            "mov_row_count": match.mov_row_count,
            "ref_coverage": match.ref_coverage,
            "mov_coverage": match.mov_coverage,
            "overlap_smaller_input": match.overlap_smaller_input,
            "ref_only_count": match.dropped_ref_only_count,
            "mov_only_count": match.dropped_mov_only_count,
            "reordered": match.matched_rows_reordered,
            "low_overlap_allowed": match.low_overlap_allowed,
            "row_aligned": request.row_aligned,
            "status": match.status,
        }
        warnings.extend(str(value) for value in match.warnings)
    except (ValueError, FileNotFoundError, OSError) as exc:
        errors.append(str(exc))
        try:
            if resolved is None:
                resolved = resolve_input_policies(
                    InputPolicyRequest(
                        ref=request.ref,
                        mov=request.mov,
                        row_aligned=request.row_aligned,
                        allow_low_overlap=request.allow_low_overlap,
                        identity_mode=request.identity_mode,
                        identity_columns=request.identity_columns,
                        mapping_file=request.mapping_file,
                    )
                )
            ref_policy = resolved.identity_ref
            mov_policy = resolved.identity_mov
            ref_count, _ref_type, ref_identity = inspect_source_identity(
                request.ref, identity_policy=ref_policy
            )
            mov_count, _mov_type, mov_identity = inspect_source_identity(
                request.mov, identity_policy=mov_policy
            )
            identities["ref"]["row_count"] = ref_count
            identities["mov"]["row_count"] = mov_count
            identity_reports = {
                "ref": to_json_safe(ref_identity.report),
                "mov": to_json_safe(mov_identity.report),
            }
            matching = _diagnostic_matching(ref_identity.keys, mov_identity.keys, request)
        except (ValueError, FileNotFoundError, OSError, KeyError):
            pass

    align_diagnosis = None
    if source_types.get("ref") == "relion" and source_types.get("mov") == "relion":
        try:
            from cryorole.preflight.align_diagnosis import diagnose_star_matching, stale_subtraction_warnings

            for finding in stale_subtraction_warnings((request.ref, request.mov)):
                warnings.append(f"[{finding['code']}] {finding['message']}")
            overlap = float(matching.get("overlap_smaller_input", 0.0) or 0.0)
            if not request.row_aligned and (errors or overlap < 0.5):
                align_diagnosis = diagnose_star_matching(request.ref, request.mov)
        except (ValueError, FileNotFoundError, OSError, KeyError) as exc:
            align_diagnosis = {"error": f"align diagnosis unavailable: {exc}"}

    matched_count = int(matching.get("matched_count", 0))
    resource = estimate_run_resources(
        matched_count=matched_count,
        k_neighbors=request.k_neighbors,
        query_batch_size=request.density_query_batch_size,
        raw_csv=request.raw_csv,
        visualize=request.visualize,
        output_dir=request.output_dir,
    )
    if not resource["disk_sufficient_with_20pct_margin"]:
        errors.append("Estimated output exceeds available disk space with the 20% safety margin")
    if resource["recommend_no_visualize"] and request.visualize:
        warnings.append("large_input_consider_--no-visualize")
    readiness = "BLOCKED" if errors else ("READY_WITH_WARNINGS" if warnings else "READY")
    resolved_command = _resolved_run_command(request)
    report: dict[str, Any] = {
        "artifact_type": "cryorole_preflight_report",
        "schema_version": PREFLIGHT_SCHEMA_VERSION,
        "timestamp": timestamp,
        "readiness": readiness,
        "inputs": {
            "ref": {"original_path": str(request.ref), "source_type": source_types.get("ref")},
            "mov": {"original_path": str(request.mov), "source_type": source_types.get("mov")},
        },
        "source_identities": identities,
        "convention_resolution": conventions,
        "identity_policy": (
            resolved.report_payload()
            if resolved is not None
            else {
                "ref": None,
                "mov": None,
                "row_aligned": request.row_aligned,
                "allow_low_overlap": request.allow_low_overlap,
            }
        ),
        "identity_reports": identity_reports,
        "matching": matching,
        "input_sanity": input_sanity,
        "align_diagnosis": align_diagnosis,
        "run_policy": {
            "resolved_backend": "array_native" if request.run_backend == "auto" else request.run_backend,
            "k_neighbors": request.k_neighbors,
            "sld_metric": request.sld_metric,
            "density_query_batch_size": request.density_query_batch_size,
            "raw_csv": request.raw_csv,
            "visualize": request.visualize,
        },
        "resource_estimate": resource,
        "warnings": warnings,
        "errors": errors,
        "resolved_run_command": resolved_command,
        "recommended_next_command": (
            align_diagnosis["recommended_command"]
            if readiness == "BLOCKED" and align_diagnosis and align_diagnosis.get("recommended_command")
            else _next_command(readiness, errors, request, resolved_command)
        ),
    }
    return PreflightResult(
        report=report,
        array_preflight=phase,
        resolved_input_policies=resolved,
    )


def _convention_record(source_type: str) -> dict[str, object]:
    if source_type == "relion":
        return asdict(ConventionPolicy.relion_default())
    return {
        **asdict(ConventionPolicy.cryosparc_default()),
        "pose_field": "alignments3D/pose",
        "pose_encoding": "rotation_vector_axis_angle",
        "reference_implementation": "pyem csparc2star.py: rot2euler(expmap(pose))",
    }


def _resolved_run_command(request: PreflightRequest) -> str:
    parts = ["cryorole", "run", "--ref", str(Path(request.ref).resolve()), "--mov", str(Path(request.mov).resolve())]
    if str(request.output_dir) != "cryorole_outputs":
        parts.extend(["--output-dir", str(Path(request.output_dir).resolve())])
    if request.row_aligned:
        parts.append("--row-aligned")
    if request.allow_low_overlap:
        parts.append("--allow-low-overlap")
    if request.sld_metric != "rotvec_euclidean":
        parts.extend(["--sld-metric", request.sld_metric])
    if not request.visualize:
        parts.append("--no-visualize")
    if request.identity_mode is not None:
        parts.extend(["--identity-mode", request.identity_mode])
    for column in request.identity_columns:
        parts.extend(["--identity-column", column])
    if request.mapping_file is not None:
        parts.extend(["--mapping-file", str(request.mapping_file)])
    return shlex.join(parts)


def _next_command(readiness: str, errors: list[str], request: PreflightRequest, run_command: str) -> str:
    if readiness != "BLOCKED":
        return run_command
    message = " ".join(errors).casefold()
    if "duplicate" in message:
        return f"cryorole align --ref {shlex.quote(str(request.ref))} --mov {shlex.quote(str(request.mov))}"
    if "overlap" in message or "zero matches" in message:
        return f"cryorole align --ref {shlex.quote(str(request.ref))} --mov {shlex.quote(str(request.mov))}"
    return "Fix the reported input/schema error, then rerun cryorole preflight"


def _diagnostic_matching(ref_keys, mov_keys, request: PreflightRequest) -> dict[str, object]:
    ref_values = ref_keys.astype(str)
    mov_values = mov_keys.astype(str)
    duplicate_ref = int(ref_values.duplicated(keep=False).sum())
    duplicate_mov = int(mov_values.duplicated(keep=False).sum())
    if duplicate_ref or duplicate_mov:
        matched = 0
        reordered = False
    else:
        mov_rows = {key: index for index, key in enumerate(mov_values.tolist())}
        pairs = [(index, mov_rows[key]) for index, key in enumerate(ref_values.tolist()) if key in mov_rows]
        matched = len(pairs)
        reordered = any(ref_row != mov_row for ref_row, mov_row in pairs)
    ref_count, mov_count = len(ref_values), len(mov_values)
    denominator = min(ref_count, mov_count) or 1
    return {
        "identity_key": None,
        "matched_count": matched,
        "ref_row_count": ref_count,
        "mov_row_count": mov_count,
        "ref_coverage": matched / ref_count if ref_count else 0.0,
        "mov_coverage": matched / mov_count if mov_count else 0.0,
        "overlap_smaller_input": matched / denominator,
        "ref_only_count": ref_count - matched,
        "mov_only_count": mov_count - matched,
        "reordered": reordered,
        "duplicate_ref_count": duplicate_ref,
        "duplicate_mov_count": duplicate_mov,
        "low_overlap_allowed": request.allow_low_overlap,
        "row_aligned": request.row_aligned,
        "status": "blocked",
    }
