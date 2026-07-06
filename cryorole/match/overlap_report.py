"""Overlap status helpers for match reports."""

from __future__ import annotations


def overlap_status(
    overlap_ratio: float,
    *,
    threshold: float,
    low_overlap_behavior: str,
) -> str:
    """Return an overlap status according to policy."""

    if overlap_ratio >= threshold:
        return "ok"
    if low_overlap_behavior == "warn":
        return "warning_low_overlap"
    if low_overlap_behavior == "fail":
        return "failed_low_overlap"
    raise ValueError(f"Unsupported low overlap behavior: {low_overlap_behavior}")

