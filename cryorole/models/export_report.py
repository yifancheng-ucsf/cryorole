"""Export report models."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping


@dataclass(frozen=True)
class ExportReport:
    """Auditable report for derived export artifacts."""

    output_paths: Mapping[str, str]
    row_counts: Mapping[str, int]
    source_selection_id: str
    timestamp: str
    export_policy: Mapping[str, Any]
    export_type: str = "selection_export"
    schema_version: str = "1"
