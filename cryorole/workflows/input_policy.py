"""Single source of truth for run input, identity, and matching policies."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

from cryorole.match import public_match_policy
from cryorole.models.policies import ConventionPolicy, IdentityPolicy, MatchPolicy


@dataclass(frozen=True)
class InputPolicyRequest:
    """User assertions needed to resolve both run inputs exactly once."""

    ref: str | Path
    mov: str | Path
    row_aligned: bool = False
    allow_low_overlap: bool = False
    identity_mode: str | None = None
    identity_columns: tuple[str, ...] = ()
    mapping_file: str | Path | None = None


@dataclass(frozen=True)
class ResolvedInputPolicies:
    """Fully resolved input policies shared by preflight, run, and guide."""

    source_ref: str
    source_mov: str
    identity_ref: IdentityPolicy
    identity_mov: IdentityPolicy
    convention_ref: ConventionPolicy | None
    convention_mov: ConventionPolicy | None
    match: MatchPolicy
    row_aligned: bool
    allow_low_overlap: bool

    def report_payload(self) -> dict[str, Any]:
        """Return the stable JSON representation recorded by workflow artifacts."""

        return {
            "ref": asdict(self.identity_ref),
            "mov": asdict(self.identity_mov),
            "row_aligned": self.row_aligned,
            "allow_low_overlap": self.allow_low_overlap,
        }


def resolve_input_policies(request: InputPolicyRequest) -> ResolvedInputPolicies:
    """Resolve source, identity, convention, and match policies once."""

    source_ref = resolve_source_type(request.ref)
    source_mov = resolve_source_type(request.mov)
    identity_ref = _resolve_identity_policy(request, source_ref)
    identity_mov = _resolve_identity_policy(request, source_mov)
    return ResolvedInputPolicies(
        source_ref=source_ref,
        source_mov=source_mov,
        identity_ref=identity_ref,
        identity_mov=identity_mov,
        convention_ref=_convention_policy(source_ref),
        convention_mov=_convention_policy(source_mov),
        match=public_match_policy(allow_low_overlap=request.allow_low_overlap),
        row_aligned=request.row_aligned or request.identity_mode == "row_aligned",
        allow_low_overlap=request.allow_low_overlap,
    )


def resolve_source_type(path: str | Path) -> str:
    """Resolve a supported metadata type from its filename suffix."""

    suffix = Path(path).suffix.lower()
    if suffix == ".star":
        return "relion"
    if suffix == ".cs":
        return "cryosparc"
    raise ValueError(f"Unsupported input suffix {suffix!r}; expected .star or .cs")


def _resolve_identity_policy(
    request: InputPolicyRequest,
    source_type: str,
) -> IdentityPolicy:
    if request.row_aligned or request.identity_mode == "row_aligned":
        return IdentityPolicy(identity_mode="row_aligned")
    if request.identity_mode == "explicit_mapping_file":
        if request.mapping_file is None:
            raise ValueError("explicit_mapping_file identity requires --mapping-file")
        if source_type != "relion":
            raise ValueError(
                "Explicit mapping-file identity is currently supported only for RELION STAR"
            )
        return IdentityPolicy(
            identity_mode="explicit_mapping_file",
            mapping_file=str(request.mapping_file),
        )
    if request.identity_columns or request.identity_mode == "relion_user_columns":
        if not request.identity_columns:
            raise ValueError("relion_user_columns identity requires --identity-column")
        if source_type != "relion":
            raise ValueError("Explicit identity columns are currently supported only for RELION STAR")
        return IdentityPolicy.relion_user_columns(request.identity_columns)
    if request.identity_mode not in {None, "cryosparc_uid"}:
        raise ValueError(f"Unsupported identity mode: {request.identity_mode}")
    if request.identity_mode == "cryosparc_uid" and source_type != "cryosparc":
        raise ValueError("cryosparc_uid identity requires CryoSPARC .cs inputs")
    if source_type == "cryosparc":
        return IdentityPolicy.cryosparc_uid()
    return IdentityPolicy.relion_image_name()


def _convention_policy(source_type: str) -> ConventionPolicy | None:
    if source_type == "relion":
        return ConventionPolicy.relion_default()
    if source_type == "cryosparc":
        return ConventionPolicy.cryosparc_default()
    return None
