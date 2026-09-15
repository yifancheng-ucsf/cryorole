"""Stable, content-addressed source metadata identity."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
import hashlib
import hmac
import os
from typing import Mapping


SOURCE_IDENTITY_SCHEMA_VERSION = "1"
HASH_CHUNK_SIZE = 1024 * 1024


@dataclass(frozen=True)
class SourceIdentity:
    """Identity recorded for one immutable run input."""

    original_path: str
    resolved_path: str
    source_type: str
    size_bytes: int
    mtime_epoch_sec: float
    mtime_ns: int
    sha256: str
    row_count: int
    schema_version: str = SOURCE_IDENTITY_SCHEMA_VERSION

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class SourceVerification:
    """Result of export-time source verification."""

    resolved_path: str
    verification_status: str
    relocated: bool
    sha256_verified: bool
    legacy_unverified_allowed: bool = False


@dataclass(frozen=True)
class SourceIdentityGuard:
    """Content-level guard for sources already fingerprinted by a workflow."""

    records: Mapping[str, SourceIdentity | Mapping[str, object]]

    def assert_unchanged(self) -> None:
        """Re-hash every source and fail if its recorded content is no longer present."""

        for domain, record in self.records.items():
            values = record.to_dict() if isinstance(record, SourceIdentity) else dict(record)
            expected_hash = values.get("sha256") or values.get("hash")
            if not expected_hash:
                raise ValueError(f"{domain} source identity has no SHA-256 and cannot be verified")
            path = Path(str(values.get("resolved_path", "")))
            if not path.is_file():
                raise ValueError(f"{domain} source metadata is missing after preflight: {path}")
            observed_hash = stream_sha256(path)
            if not hmac.compare_digest(observed_hash, str(expected_hash)):
                raise ValueError(
                    f"{domain} source metadata content changed after preflight; run aborted"
                )

    def inputs_unchanged(self) -> bool:
        """Return whether all records still resolve to their fingerprinted content."""

        try:
            self.assert_unchanged()
        except (OSError, ValueError):
            return False
        return True


def stream_sha256(path: str | Path, *, chunk_size: int = HASH_CHUNK_SIZE) -> str:
    """Hash a source file without loading it into memory."""

    if chunk_size < 1:
        raise ValueError("hash chunk_size must be positive")
    hasher = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            hasher.update(chunk)
    return hasher.hexdigest()


def build_source_identity(
    original_path: str | Path,
    *,
    source_type: str,
    row_count: int,
) -> SourceIdentity:
    """Resolve and fingerprint one input at run time."""

    original = str(original_path)
    resolved = Path(original_path).expanduser().resolve(strict=True)
    if not resolved.is_file():
        raise FileNotFoundError(f"Source metadata is not a file: {resolved}")
    stat = resolved.stat()
    return SourceIdentity(
        original_path=original,
        resolved_path=str(resolved),
        source_type=source_type,
        size_bytes=int(stat.st_size),
        mtime_epoch_sec=float(stat.st_mtime),
        mtime_ns=int(stat.st_mtime_ns),
        sha256=stream_sha256(resolved),
        row_count=int(row_count),
    )


def verify_source_identity(
    record: SourceIdentity | dict[str, object] | None,
    *,
    relocated_path: str | Path | None = None,
    allow_unverified_source: bool = False,
    operation: str = "export",
) -> SourceVerification:
    """Verify a current or explicitly relocated file against recorded identity."""

    if record is None or not _record_hash(record):
        if not allow_unverified_source:
            raise ValueError(
                f"Legacy run bundle has no verified source SHA-256; {operation} is refused. "
                "Use --allow-unverified-source only after independently validating the source."
            )
        candidate = Path(relocated_path) if relocated_path is not None else _legacy_path(record)
        if candidate is None or not candidate.expanduser().resolve().is_file():
            raise FileNotFoundError("Legacy source file could not be resolved for unverified export")
        return SourceVerification(
            resolved_path=str(candidate.expanduser().resolve()),
            verification_status="legacy_unverified_allowed",
            relocated=relocated_path is not None,
            sha256_verified=False,
            legacy_unverified_allowed=True,
        )

    values = record.to_dict() if isinstance(record, SourceIdentity) else dict(record)
    recorded_path = Path(str(values["resolved_path"]))
    candidate = Path(relocated_path) if relocated_path is not None else recorded_path
    candidate = candidate.expanduser().resolve()
    if not candidate.is_file():
        if relocated_path is None:
            raise FileNotFoundError(
                f"Recorded source is missing: {recorded_path}. Provide an explicit relocated source."
            )
        raise FileNotFoundError(f"Relocated source does not exist: {candidate}")
    expected_type = str(values.get("source_type", ""))
    if expected_type and not _suffix_matches_source_type(candidate, expected_type):
        raise ValueError(
            f"Relocated source type does not match recorded {expected_type!r}: {candidate}"
        )
    observed_hash = stream_sha256(candidate)
    expected_hash = str(values.get("sha256") or values.get("hash"))
    if not hmac.compare_digest(observed_hash, expected_hash):
        raise ValueError(
            f"Source metadata SHA-256 mismatch; {operation} is refused because recorded source-row "
            "provenance may refer to a different file."
        )
    expected_size = values.get("size_bytes")
    if expected_size is not None and int(candidate.stat().st_size) != int(expected_size):
        raise ValueError("Source metadata size does not match the recorded source identity")
    return SourceVerification(
        resolved_path=str(candidate),
        verification_status="verified_sha256",
        relocated=os.path.normcase(str(candidate)) != os.path.normcase(str(recorded_path)),
        sha256_verified=True,
    )


def _record_hash(record: SourceIdentity | dict[str, object]) -> str | None:
    if isinstance(record, SourceIdentity):
        return record.sha256
    value = record.get("sha256") or record.get("hash")
    return str(value) if value else None


def _legacy_path(record: SourceIdentity | dict[str, object] | None) -> Path | None:
    if record is None:
        return None
    values = record.to_dict() if isinstance(record, SourceIdentity) else record
    value = values.get("resolved_path") or values.get("path") or values.get("original_path")
    return Path(str(value)) if value else None


def _suffix_matches_source_type(path: Path, source_type: str) -> bool:
    normalized = source_type.lower()
    if normalized in {"relion", "relion_star", "star"}:
        return path.suffix.lower() == ".star"
    if normalized in {"cryosparc", "cryosparc_cs", "cs"}:
        return path.suffix.lower() == ".cs"
    return True
