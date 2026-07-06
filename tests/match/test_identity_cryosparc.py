import pandas as pd
import numpy as np
import pytest

from cryorole.match import match_resolved_identities, resolve_identity
from cryorole.models.policies import IdentityPolicy
from cryorole.workflows.pipeline_runner import PipelineRunner


def test_cryosparc_identity_defaults_to_uid() -> None:
    table = pd.DataFrame({"uid": [101, 102, 103]})

    resolved = resolve_identity(table, source_type="cryosparc")

    assert list(resolved.keys) == ["101", "102", "103"]
    assert resolved.report.identity_mode == "cryosparc_uid"
    assert resolved.report.identity_columns == ("uid",)


def test_cryosparc_uid_matching_uses_uid_not_row_order() -> None:
    domain_a = pd.DataFrame({"uid": [101, 102, 103]})
    domain_b = pd.DataFrame({"uid": [103, 101, 104]})

    resolved_a = resolve_identity(domain_a, source_type="cryosparc")
    resolved_b = resolve_identity(domain_b, source_type="cryosparc")
    match_table, report = match_resolved_identities(resolved_a, resolved_b)

    assert list(match_table.data["particle_key"]) == ["101", "103"]
    assert list(match_table.data["domain_a_row"]) == [0, 2]
    assert list(match_table.data["domain_b_row"]) == [1, 0]
    assert report.matched_count == 2
    assert report.match_key == "uid"
    assert report.ref_row_count == 3
    assert report.mov_row_count == 3
    assert report.matched_row_count == 2
    assert report.dropped_ref_only_count == 1
    assert report.dropped_mov_only_count == 1
    assert report.matched_rows_reordered is True
    assert "inputs_matched_by_key_and_reordered_before_ro_computation" in report.warnings
    assert "inputs_matched_by_key_with_unmatched_particles_dropped" in report.warnings


def test_cryosparc_row_aligned_preflight_fails_on_row_count_mismatch(tmp_path) -> None:
    ref = tmp_path / "ref.cs"
    mov = tmp_path / "mov.cs"
    _write_cs(ref, [101, 102])
    _write_cs(mov, [101])
    policy = IdentityPolicy(identity_mode="row_aligned")

    with pytest.raises(ValueError, match="row_aligned identity mode requires equal"):
        PipelineRunner().run_phase1(
            str(ref),
            str(mov),
            domain_a_name="ref",
            domain_b_name="mov",
            identity_policy_a=policy,
            identity_policy_b=policy,
        )


def _write_cs(path, uids: list[int]) -> None:
    array = np.zeros(
        len(uids),
        dtype=[
            ("uid", "<u8"),
            ("alignments3D/pose", "<f4", (3,)),
        ],
    )
    array["uid"] = np.asarray(uids, dtype=np.uint64)
    with open(path, "wb") as handle:
        np.save(handle, array)
