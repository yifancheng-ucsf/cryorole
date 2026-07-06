"""Identity resolution report model."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Mapping


@dataclass(frozen=True)
class IdentityReport:
    """Auditable report for particle identity resolution."""

    identity_mode: str
    identity_columns: tuple[str, ...]
    column_normalization_rules: Mapping[str, str]
    unique_rate: float
    duplicate_count: int
    collision_examples: tuple[Mapping[str, Any], ...] = field(default_factory=tuple)
    matched_count: int = 0
    unmatched_a: int = 0
    unmatched_b: int = 0
    overlap_ratio: float = 0.0
    status: str = "ok"

