"""Manifest report model."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping


@dataclass(frozen=True)
class ManifestReport:
    """Report for a written run manifest."""

    output_path: str
    size_bytes: int
    manifest_type: str
    schema_version: str
    timestamp: str
    manifest_policy: Mapping[str, Any]
