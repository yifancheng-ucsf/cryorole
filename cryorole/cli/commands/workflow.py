"""CLI adapters for run, preflight, and workflow guidance."""

from __future__ import annotations

import json
import sys
import webbrowser
from dataclasses import replace
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable

from cryorole.export import write_json_artifact
from cryorole.export.serialization import to_json_safe
from cryorole.interactive import ExploreSession, create_explore_server
from cryorole.preflight import PreflightRequest, run_preflight
from cryorole.workflow import build_guide_plan, derive_next_actions, inspect_run_status
from cryorole.workflows.input_policy import (
    InputPolicyRequest,
    ResolvedInputPolicies,
    resolve_input_policies,
)
from cryorole.workflows.pipeline_runner import PipelineRunner
from cryorole.workflows.run_service import RunExecutionContext, RunRequest, execute_run


def run_command(args, *, runner: PipelineRunner | None = None) -> int:
    """Create a run bundle transactionally and publish it after validation."""

    _resolve_run_backend(args.run_backend)
    resolved = resolve_input_policies(input_policy_request_from_args(args))
    if getattr(args, "dry_run", False):
        result = run_preflight(
            preflight_request_from_args(args, resolved_input_policies=resolved)
        )
        emit_preflight_result(result, getattr(args, "preflight_json", None))
        return result.exit_code
    if getattr(args, "preflight_json", None) is not None:
        raise ValueError("run --json is valid only together with --dry-run")
    result = execute_run(
        replace(
            RunRequest.from_namespace(args),
            resolved_input_policies=resolved,
        ),
        context=RunExecutionContext(
            runner=runner,
            allow_unverified_test_sources=runner is not None,
        ),
    )
    print(str(result.output_dir))
    return 0


def preflight_command(args) -> int:
    resolved = resolve_input_policies(input_policy_request_from_args(args))
    result = run_preflight(
        preflight_request_from_args(args, resolved_input_policies=resolved)
    )
    emit_preflight_result(result, getattr(args, "preflight_json", None))
    return result.exit_code


def status_command(args) -> int:
    status = inspect_run_status(args.run_dir)
    if emit_optional_json(status, getattr(args, "json_output", None)):
        return 0
    print(f"{status['bundle_status'].upper()} — {status['run_dir']}")
    if status.get("run_id"):
        print(f"Run ID: {status['run_id']}")
    raw = status.get("raw_landscape")
    if raw:
        validity = "valid" if raw.get("valid") else "invalid"
        print(f"Raw landscape: {raw.get('row_count', 'unknown')} particles ({validity})")
    print(
        f"Canonical frames: {len(status.get('canonical_frames', []))}; "
        f"selections: {len(status.get('selections', []))}; "
        f"exports: {len(status.get('exports', []))}."
    )
    for warning in status.get("warnings", []):
        print(f"WARNING: {warning}")
    return 0


def next_command(args) -> int:
    status = inspect_run_status(args.run_dir)
    actions = derive_next_actions(status)
    payload = {
        "artifact_type": "cryorole_next_actions",
        "schema_version": "1.0",
        "status": status,
        "actions": actions,
    }
    if emit_optional_json(payload, getattr(args, "json_output", None)):
        return 0
    print(f"{status['bundle_status'].upper()} — recommended next actions")
    for index, action in enumerate(actions, start=1):
        print(f"{index}. {action['category'].upper()}: {action['reason']}")
        print(f"   {action['command']}")
    return 0


def guide_command(args, *, parser_factory: Callable[[], Any] | None = None) -> int:
    if args.run_dir and (args.ref or args.mov):
        raise ValueError("guide accepts --run-dir or --ref/--mov, not both")
    plan = build_guide_plan(
        run_dir=args.run_dir,
        ref=args.ref,
        mov=args.mov,
        output_dir=args.output_dir,
        non_interactive=bool(args.non_interactive),
        row_aligned=bool(args.row_aligned),
        allow_low_overlap=bool(args.allow_low_overlap),
    )
    if emit_optional_json(plan, getattr(args, "json_output", None)):
        return 0
    print(f"GUIDE {plan['mode'].upper()}")
    if "preflight" in plan:
        print(f"Preflight: {plan['preflight']['readiness']}")
    for index, action in enumerate(plan["actions"], start=1):
        print(f"{index}. {action['category'].upper()}: {action['reason']}")
        print(f"   {action['command']}")

    execute = bool(args.execute_run)
    if plan["mode"] == "new" and not args.non_interactive and not execute:
        if sys.stdin.isatty() and plan.get("requires_user_input"):
            execute = input("Execute the preflight-approved run now? [y/N] ").strip().casefold() in {"y", "yes"}
        else:
            print("Non-TTY detected: no command was executed. Use --execute-run explicitly.")
    if execute:
        if plan["mode"] != "new":
            raise ValueError("--execute-run applies only to guide --ref/--mov")
        if plan["preflight"]["readiness"] == "BLOCKED":
            raise ValueError("guide will not execute a blocked preflight")
        if parser_factory is None:
            from cryorole.cli.main import build_parser

            parser_factory = build_parser
        run_argv = [
            "run", "--ref", str(args.ref), "--mov", str(args.mov),
            "--output-dir", str(args.output_dir),
        ]
        if args.row_aligned:
            run_argv.append("--row-aligned")
        if args.allow_low_overlap:
            run_argv.append("--allow-low-overlap")
        run_command(parser_factory().parse_args(run_argv))
        _record_guide_decision(
            Path(args.output_dir), plan, decision="run_confirmed_and_executed"
        )
    return 0


def explore_command(args) -> int:
    session = ExploreSession(
        args.run_dir,
        space=args.space,
        canonical_id=args.canonical_id,
        selection_id=args.selection_id,
        max_display_points=args.max_display_points,
        display_threshold=args.threshold,
        display_top_fraction=args.top_fraction,
        colormap=args.colormap,
    )
    server = create_explore_server(session, port=args.port)
    payload = session.session_payload()
    print(
        f"EXPLORE READY — {payload['displayed_point_count']} displayed / "
        f"{payload['full_candidate_count']} full candidates"
    )
    print("Draft changes are display-only; Confirm performs an exact full-data SO(3) selection.")
    print(server.url)
    if not args.no_open:
        webbrowser.open(server.url)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("Explore stopped; no background process remains.")
    finally:
        server.shutdown()
    return 0


def emit_optional_json(payload: dict[str, Any], destination: str | None) -> bool:
    if destination is None:
        return False
    if destination == "-":
        print(json.dumps(to_json_safe(payload), indent=2, sort_keys=True))
    else:
        path = write_json_artifact(payload, destination, overwrite=False)
        print(str(path))
    return True


def input_policy_request_from_args(args) -> InputPolicyRequest:
    return InputPolicyRequest(
        ref=args.ref,
        mov=args.mov,
        row_aligned=bool(getattr(args, "row_aligned", False)),
        allow_low_overlap=bool(getattr(args, "allow_low_overlap", False)),
        identity_mode=getattr(args, "identity_mode", None),
        identity_columns=tuple(getattr(args, "identity_column", ()) or ()),
        mapping_file=getattr(args, "mapping_file", None),
    )


def preflight_request_from_args(
    args,
    *,
    resolved_input_policies: ResolvedInputPolicies | None = None,
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
    )


def emit_preflight_result(result, json_destination: str | None) -> None:
    if json_destination == "-":
        print(json.dumps(to_json_safe(result.report), indent=2, sort_keys=True))
        return
    if json_destination:
        path = write_json_artifact(result.report, json_destination, overwrite=False)
        print(f"{result.report['readiness']} — report: {path}")
        return
    report = result.report
    print(report["readiness"])
    matching = report.get("matching", {})
    if matching:
        print(
            "Matched {matched_count} particles; ref coverage {ref_coverage:.1%}; "
            "mov coverage {mov_coverage:.1%}.".format(**matching)
        )
    estimate = report["resource_estimate"]
    print(
        f"Estimated peak memory {estimate['estimated_peak_memory_mib']:.2f} MiB; "
        f"bundle {estimate['estimated_bundle_mib']:.2f} MiB."
    )
    for warning in report["warnings"]:
        print(f"WARNING: {warning}")
    for error in report["errors"]:
        print(f"REQUIRED ACTION: {error}")
    print(f"Next: {report['recommended_next_command']}")


def _record_guide_decision(run_dir: Path, plan: dict[str, Any], *, decision: str) -> None:
    path = run_dir.expanduser().resolve() / "reports" / "guide_history.json"
    history: list[dict[str, object]] = []
    if path.is_file():
        payload = json.loads(path.read_text(encoding="utf-8"))
        history = list(payload.get("history", []))
    history.append({"timestamp": datetime.now(timezone.utc).isoformat(), "decision": decision})
    write_json_artifact(
        {"artifact_type": "cryorole_guide_history", "schema_version": "1.0", "history": history},
        path,
        overwrite=path.exists(),
    )


def _resolve_run_backend(requested: str) -> str:
    if requested in {"auto", "array_native"}:
        return "array_native"
    if requested == "dataframe_compat":
        return requested
    raise ValueError(f"Unsupported run backend: {requested}")
