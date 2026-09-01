"""Resumable guide plans that orchestrate services without becoming a pipeline."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from cryorole.preflight import PreflightRequest, run_preflight
from cryorole.workflow.next_action import derive_next_actions
from cryorole.workflow.status import inspect_run_status
from cryorole.workflows.input_policy import InputPolicyRequest, resolve_input_policies


def build_guide_plan(
    *,
    run_dir: str | Path | None = None,
    ref: str | Path | None = None,
    mov: str | Path | None = None,
    output_dir: str | Path = "cryorole_outputs",
    non_interactive: bool = False,
    row_aligned: bool = False,
    allow_low_overlap: bool = False,
) -> dict[str, Any]:
    """Build a plan from real artifacts or a shared preflight result."""

    if run_dir is not None:
        status = inspect_run_status(run_dir)
        return {
            "artifact_type": "cryorole_guide_plan", "schema_version": "1.0",
            "mode": "resume", "non_interactive": non_interactive,
            "requires_user_input": False, "status": status,
            "actions": derive_next_actions(status),
            "automatic_steps": ["status inspection", "next-command generation"],
            "confirmation_required": ["canonicalize", "create selection", "export", "overwrite"],
        }
    if ref is None or mov is None:
        raise ValueError("guide requires either --run-dir or both --ref and --mov")
    resolved_input_policies = resolve_input_policies(
        InputPolicyRequest(
            ref=ref,
            mov=mov,
            row_aligned=row_aligned,
            allow_low_overlap=allow_low_overlap,
        )
    )
    preflight = run_preflight(
        PreflightRequest(
            ref=ref, mov=mov, output_dir=output_dir,
            row_aligned=row_aligned, allow_low_overlap=allow_low_overlap,
            resolved_input_policies=resolved_input_policies,
        )
    )
    actions = []
    if preflight.report["readiness"] != "BLOCKED":
        actions.append({
            "category": "required", "reason": "Preflight is ready; run creates the first scientific bundle.",
            "command": preflight.report["resolved_run_command"],
        })
    else:
        actions.append({
            "category": "required", "reason": "Preflight is blocked.",
            "command": preflight.report["recommended_next_command"],
        })
    return {
        "artifact_type": "cryorole_guide_plan", "schema_version": "1.0",
        "mode": "new", "non_interactive": non_interactive,
        "requires_user_input": bool(not non_interactive and preflight.exit_code < 2),
        "preflight": preflight.report, "actions": actions,
        "automatic_steps": ["preflight", "resource estimate", "command generation"],
        "confirmation_required": ["run execution", "row-aligned assertion", "low-overlap override", "canonicalize", "create selection", "export", "overwrite"],
    }
