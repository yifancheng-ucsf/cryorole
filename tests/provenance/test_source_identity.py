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
    assert identity.input_path_absolute == str(source)
    assert identity.resolved_path == str(source.resolve())
    assert verified.resolved_path == str(source.resolve())
    assert verified.sha256_verified is True
    SourceIdentityGuard({"ref": identity}).assert_unchanged()


def test_source_identity_guard_rejects_changed_bytes_with_preserved_stat(tmp_path: Path) -> None:
    source = tmp_path / "particles.star"
    source.write_bytes(b"data_particles\nAAAA\n")
    identity = build_source_identity(source, source_type="relion", row_count=1)
    original_stat = source.stat()
    source.write_bytes(b"data_particles\nBBBB\n")
    os.utime(source, ns=(original_stat.st_atime_ns, original_stat.st_mtime_ns))

    with pytest.raises(ValueError, match="content changed"):
        SourceIdentityGuard({"ref": identity.to_dict()}).assert_unchanged()


@pytest.mark.parametrize("legacy", [False, True])
def test_guard_rejects_retargeted_input_link_even_with_identical_target_bytes(tmp_path, legacy):
    first, second, link = (tmp_path / name for name in ("first.star", "second.star", "input.star"))
    first.write_bytes(b"same bytes\n")
    second.write_bytes(first.read_bytes())
    link.symlink_to(first.name)
    record = build_source_identity(link, source_type="relion", row_count=1).to_dict()
    if legacy:
        record.pop("input_path_absolute")
    link.unlink()
    link.symlink_to(second.name)
    with pytest.raises(ValueError, match="source path changed"):
        SourceIdentityGuard({"ref": record}).assert_unchanged()
    assert not SourceIdentityGuard({"ref": record}).inputs_unchanged()
    # Export can still use the recorded, unchanged resolved file.
    assert verify_source_identity(record).sha256_verified


def test_guard_tracks_relative_input_parent_link_from_original_cwd(tmp_path, monkeypatch):
    for name in ("a", "b"):
        folder = tmp_path / name
        folder.mkdir()
        (folder / "particles.star").write_bytes(b"same bytes\n")
    alias = tmp_path / "current"
    alias.symlink_to("a", target_is_directory=True)
    monkeypatch.chdir(tmp_path)
    record = build_source_identity("current/particles.star", source_type="relion", row_count=1)
    monkeypatch.chdir(tmp_path / "b")
    SourceIdentityGuard({"ref": record}).assert_unchanged()
    alias.unlink()
    alias.symlink_to("b", target_is_directory=True)
    with pytest.raises(ValueError, match="source path changed"):
        SourceIdentityGuard({"ref": record}).assert_unchanged()


@pytest.mark.parametrize("broken", ["missing", "loop"])
def test_guard_rejects_broken_input_link_with_original_file_still_present(tmp_path, broken):
    source = tmp_path / "source.star"
    source.write_bytes(b"content\n")
    alias = tmp_path / "input.star"
    alias.symlink_to(source.name)
    record = build_source_identity(alias, source_type="relion", row_count=1)
    alias.unlink()
    alias.symlink_to("missing.star" if broken == "missing" else alias.name)
    assert source.is_file()
    assert not SourceIdentityGuard({"ref": record}).inputs_unchanged()
