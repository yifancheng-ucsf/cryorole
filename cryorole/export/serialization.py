"""JSON-safe serialization helpers for export artifacts."""

from __future__ import annotations

from dataclasses import asdict, is_dataclass
from pathlib import Path
from typing import Any, Mapping

import numpy as np


def to_json_safe(value: Any) -> Any:
    """Convert dataclasses and common scientific Python values to JSON-safe data."""

    if is_dataclass(value):
        return to_json_safe(asdict(value))
    if isinstance(value, Mapping):
        return {str(key): to_json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [to_json_safe(item) for item in value]
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    return value
