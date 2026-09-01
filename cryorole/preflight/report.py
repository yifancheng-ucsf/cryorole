"""Preflight result model and identity revalidation."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from cryorole.provenance import SourceIdentityGuard
from cryorole.workflows.input_policy import ResolvedInputPolicies


PREFLIGHT_SCHEMA_VERSION = "1.0"
PREFLIGHT_EXIT_CODES = {"READY": 0, "READY_WITH_WARNINGS": 1, "BLOCKED": 2}


@dataclass(frozen=True)
class PreflightResult:
    """Structured preflight report plus an optional reusable array phase."""

    report: dict[str, Any]
    array_preflight: Any = None
    resolved_input_policies: ResolvedInputPolicies | None = None

    @property
    def exit_code(self) -> int:
        return PREFLIGHT_EXIT_CODES[str(self.report["readiness"])]

    def inputs_unchanged(self) -> bool:
        """Re-hash both inputs before consuming a previously computed result."""

        return SourceIdentityGuard(
            self.report.get("source_identities", {})
        ).inputs_unchanged()
