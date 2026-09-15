"""Typed CS categorical values aligned to the full parent landscape."""

from dataclasses import dataclass
import re
from typing import Sequence

import numpy as np


CS_METADATA_MAX_GROUPS = 100


@dataclass(frozen=True)
class CsMetadataColumn:
    """Compact scalar values; text comparisons never pass through numeric coercion."""

    values: np.ndarray

    def __post_init__(self) -> None:
        if self.values.ndim != 1 or self.values.dtype.kind not in "iubSU":
            raise ValueError("CS metadata requires scalar integer, boolean, or text values")
        if self.values.dtype.kind == "S":
            try:
                decoded = np.char.decode(self.values, "utf-8", errors="strict")
            except UnicodeDecodeError as exc:
                raise ValueError("CS metadata byte strings must be valid UTF-8") from exc
            object.__setattr__(self, "values", decoded)

    @classmethod
    def from_source(cls, values: np.ndarray, row_ids: np.ndarray) -> "CsMetadataColumn":
        rows = np.asarray(row_ids)
        if rows.ndim != 1 or rows.dtype.kind not in "iu" or np.any(rows < 0) or np.any(rows >= len(values)):
            raise ValueError("CS metadata source-row IDs must be valid integer indices into the recorded source")
        return cls(values[rows])

    @property
    def missing(self) -> np.ndarray:
        return self.values == "" if self.values.dtype.kind == "U" else np.zeros(len(self.values), dtype=bool)

    def _requested(self, value: str):
        kind = self.values.dtype.kind
        text = str(value)
        if kind in "iu":
            if not re.fullmatch(r"[+-]?[0-9]+", text.strip()):
                raise ValueError(f"CS metadata value {text!r} must be a decimal integer")
            number = int(text)
            limits = np.iinfo(self.values.dtype)
            if not limits.min <= number <= limits.max:
                raise ValueError(f"CS metadata value {text!r} is outside {self.values.dtype}")
            return number
        if kind == "b":
            normalized = text.strip().lower()
            if normalized not in {"true", "false", "1", "0"}:
                raise ValueError("Boolean CS metadata values must be true, false, 1, or 0")
            return normalized in {"true", "1"}
        if not text:
            raise ValueError("Empty CS metadata strings are missing values and cannot be selected")
        return text

    def select(self, requested: Sequence[str]) -> tuple[np.ndarray, tuple[str, ...]]:
        typed = [self._requested(value) for value in requested]
        # Explicit numeric dtype avoids uint64 values being inferred as float64.
        targets = np.asarray(typed, dtype=self.values.dtype if self.values.dtype.kind != "U" else str)
        resolved = tuple(self._label(value) for value in typed)
        return np.isin(self.values, targets) & ~self.missing, resolved

    def groups(self) -> tuple[str, ...]:
        unique, first = np.unique(self.values[~self.missing], return_index=True)
        if len(unique) > CS_METADATA_MAX_GROUPS:
            raise ValueError(
                f"CS metadata split would create {len(unique)} selections; limit is {CS_METADATA_MAX_GROUPS}. "
                "Use --metadata-value to select explicit values instead."
            )
        return tuple(self._label(value) for value in unique[np.argsort(first)])

    @staticmethod
    def _label(value) -> str:
        return str(value).lower() if isinstance(value, (bool, np.bool_)) else str(value)
