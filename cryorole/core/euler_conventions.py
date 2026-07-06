"""Centralized Euler convention policy helpers for derived representations."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable


DEFAULT_EULER_CONVENTION = "extrinsic_zyx"
DEFAULT_EULER_CONVENTION_SOURCE = "cli_default"
LEGACY_MISSING_EULER_CONVENTION_SOURCE = "legacy_missing"
EULER_CONVENTIONS = ("extrinsic_zyx", "intrinsic_zyx")
RAW_EULER_ANGLE_COLUMNS = (
    "raw_ea_zyx_alpha_deg",
    "raw_ea_zyx_beta_deg",
    "raw_ea_zyx_gamma_deg",
)
CANONICAL_EULER_ANGLE_COLUMNS = (
    "canonical_ea_zyx_alpha_deg",
    "canonical_ea_zyx_beta_deg",
    "canonical_ea_zyx_gamma_deg",
)


_SCIPY_SEQUENCE_BY_CONVENTION = {
    "extrinsic_zyx": "zyx",
    "intrinsic_zyx": "ZYX",
}
_CONVENTION_BY_SCIPY_SEQUENCE = {
    value: key for key, value in _SCIPY_SEQUENCE_BY_CONVENTION.items()
}


@dataclass(frozen=True)
class EulerConventionResolution:
    """Resolved user-facing Euler convention and implementation sequence."""

    euler_convention: str
    scipy_euler_sequence: str
    euler_convention_source: str

    def metadata(
        self,
        *,
        euler_angle_columns: Iterable[str] = RAW_EULER_ANGLE_COLUMNS,
    ) -> dict[str, object]:
        """Return the standard report/manifest metadata payload."""

        return {
            "euler_convention": self.euler_convention,
            "scipy_euler_sequence": self.scipy_euler_sequence,
            "euler_angle_columns": list(euler_angle_columns),
            "euler_convention_source": self.euler_convention_source,
        }


def resolve_euler_convention(
    euler_convention: str | None = None,
    *,
    scipy_euler_sequence: str | None = None,
    source: str = DEFAULT_EULER_CONVENTION_SOURCE,
) -> EulerConventionResolution:
    """Resolve public convention names to the SciPy sequence token.

    SciPy uses lowercase sequence strings for extrinsic fixed-axis rotations
    and uppercase sequence strings for intrinsic rotations. Keep that detail in
    this module instead of leaking it into CLI, plotting, or selection code.
    """

    if euler_convention is not None and scipy_euler_sequence is not None:
        resolved_from_sequence = convention_from_scipy_sequence(scipy_euler_sequence)
        if euler_convention != resolved_from_sequence:
            raise ValueError(
                "euler_convention and scipy_euler_sequence disagree: "
                f"{euler_convention!r} vs {scipy_euler_sequence!r}"
            )
    if euler_convention is None:
        if scipy_euler_sequence is None:
            euler_convention = DEFAULT_EULER_CONVENTION
        else:
            euler_convention = convention_from_scipy_sequence(scipy_euler_sequence)
    if euler_convention not in _SCIPY_SEQUENCE_BY_CONVENTION:
        allowed = ", ".join(EULER_CONVENTIONS)
        raise ValueError(f"Unsupported Euler convention {euler_convention!r}; expected {allowed}")
    return EulerConventionResolution(
        euler_convention=euler_convention,
        scipy_euler_sequence=_SCIPY_SEQUENCE_BY_CONVENTION[euler_convention],
        euler_convention_source=source,
    )


def convention_from_scipy_sequence(scipy_euler_sequence: str) -> str:
    """Return the public convention name for a supported SciPy sequence."""

    if scipy_euler_sequence not in _CONVENTION_BY_SCIPY_SEQUENCE:
        raise ValueError(
            "Only ZYX Euler output conventions are supported; expected "
            "'zyx' for extrinsic fixed-axis ZYX or 'ZYX' for intrinsic ZYX"
        )
    return _CONVENTION_BY_SCIPY_SEQUENCE[scipy_euler_sequence]
