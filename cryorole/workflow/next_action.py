"""Conservative, command-producing next-action recommendations."""

from __future__ import annotations

import shlex
from typing import Any


def derive_next_actions(status: dict[str, Any]) -> list[dict[str, str]]:
    run_dir = shlex.quote(str(status["run_dir"]))
    bundle = status["bundle_status"]
    if bundle == "missing":
        return [{"category": "required", "reason": "No run bundle exists.", "command": "cryorole preflight --ref REF --mov MOV"}]
    if bundle in {"failed", "incomplete"}:
        return [{"category": "required", "reason": "The run bundle is not complete and cannot be consumed safely.", "command": f"cryorole status --run-dir {run_dir} --json"}]
    if not status.get("raw_landscape") or not status["raw_landscape"].get("valid"):
        return [{"category": "required", "reason": "A valid raw landscape is required.", "command": f"cryorole status --run-dir {run_dir} --json"}]

    actions: list[dict[str, str]] = []
    if not status.get("visualizations"):
        actions.append({
            "category": "recommended", "reason": "Inspect the raw RO landscape before making a scientific selection.",
            "command": f"cryorole visualize --run-dir {run_dir} --space raw",
        })
    if not status.get("canonical_frames"):
        actions.append({
            "category": "optional", "reason": "Canonicalization can align dominant motion axes, but is not scientifically required.",
            "command": f"cryorole canonicalize --run-dir {run_dir}",
        })
    selections = status.get("selections") or []
    if not selections:
        actions.append({
            "category": "recommended", "reason": "Explore a draft neighborhood; display filters remain separate from selection.",
            "command": f"cryorole explore --run-dir {run_dir} --space raw",
        })
    else:
        for selection in selections[:3]:
            selection_id = shlex.quote(str(selection["selection_id"]))
            actions.extend([
                {
                    "category": "recommended", "reason": "Inspect the confirmed scientific selection.",
                    "command": f"cryorole visualize --run-dir {run_dir} --selection-id {selection_id}",
                },
                {
                    "category": "recommended", "reason": "Export uses recorded source-row provenance and verifies source hashes.",
                    "command": f"cryorole export --run-dir {run_dir} --selection-id {selection_id}",
                },
            ])
    if status.get("exports"):
        actions.append({
            "category": "optional", "reason": "Exports exist; review their reports and output integrity.",
            "command": f"cryorole status --run-dir {run_dir} --json",
        })
    return actions
