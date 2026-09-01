from __future__ import annotations

from pathlib import Path
import os

import pytest

from cryorole.provenance import SourceIdentityGuard, build_source_identity, verify_source_identity


def test_relative_input_identity_resolves_from_original_cwd_and_verifies_elsewhere(
    tmp_path: Path,
    monkeypatch,
) -> None:
    input_dir = tmp_path / "inputs"
    input_dir.mkdir()
    source = input_dir / "particles.star"
    source.write_text("data_particles\n\nloop_\n_rlnImageName\n1@a.mrcs\n", encoding="utf-8")
    monkeypatch.chdir(input_dir)
    identity = build_source_identity(
        "particles.star",
        source_type="relion",
        row_count=1,
    )
    other_dir = tmp_path / "elsewhere"
    other_dir.mkdir()
    monkeypatch.chdir(other_dir)

    verified = verify_source_identity(identity)

    assert identity.original_path == "particles.star"
    assert identity.resolved_path == str(source.resolve())
    assert verified.resolved_path == str(source.resolve())
    assert verified.sha256_verified is True


def test_source_identity_guard_rejects_changed_bytes_with_preserved_stat(tmp_path: Path) -> None:
    source = tmp_path / "particles.star"
    source.write_bytes(b"data_particles\nAAAA\n")
    identity = build_source_identity(source, source_type="relion", row_count=1)
    original_stat = source.stat()
    source.write_bytes(b"data_particles\nBBBB\n")
    os.utime(source, ns=(original_stat.st_atime_ns, original_stat.st_mtime_ns))

    with pytest.raises(ValueError, match="content changed"):
        SourceIdentityGuard({"ref": identity.to_dict()}).assert_unchanged()
