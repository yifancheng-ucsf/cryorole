import pandas as pd
import pytest

from cryorole.match import (
    IdentityResolutionError,
    match_resolved_identities,
    resolve_identity,
)
from cryorole.models.policies import IdentityPolicy, MatchPolicy
from cryorole.workflows.pipeline_runner import PipelineRunner


def test_relion_identity_defaults_to_image_name() -> None:
    table = pd.DataFrame({"_rlnImageName": ["1@stack.mrcs", "2@stack.mrcs"]})

    resolved = resolve_identity(table, source_type="relion", policy=IdentityPolicy.relion_image_name())

    assert list(resolved.keys) == ["1@stack.mrcs", "2@stack.mrcs"]
    assert resolved.report.identity_mode == "relion_image_name"
    assert resolved.report.identity_columns == ("_rlnImageName",)


def test_relion_identity_prefers_tomo_particle_name_when_present() -> None:
    table = pd.DataFrame(
        {
            "_rlnTomoParticleName": ["TS_001/1", "TS_001/2"],
            "_rlnImageName": ["1@stack.mrcs", "2@stack.mrcs"],
        }
    )

    resolved = resolve_identity(table, source_type="relion", policy=IdentityPolicy.relion_image_name())

    assert list(resolved.keys) == ["TS_001/1", "TS_001/2"]
    assert resolved.report.identity_mode == "relion_image_name"
    assert resolved.report.identity_columns == ("_rlnTomoParticleName",)


def test_relion_image_name_identity_fails_without_image_name() -> None:
    table = pd.DataFrame({"_rlnParticleName": ["p1", "p2"]})

    with pytest.raises(IdentityResolutionError, match="default identity requires"):
        resolve_identity(table, source_type="relion", policy=IdentityPolicy.relion_image_name())


def test_relion_identity_resolves_user_columns_and_matches_by_key() -> None:
    domain_a = pd.DataFrame({"_rlnParticleName": ["p1", "p2", "p3"]})
    domain_b = pd.DataFrame({"_rlnParticleName": ["p3", "p1", "p4"]})
    policy = IdentityPolicy.relion_user_columns(["_rlnParticleName"])

    resolved_a = resolve_identity(domain_a, source_type="relion", policy=policy)
    resolved_b = resolve_identity(domain_b, source_type="relion", policy=policy)
    match_table, report = match_resolved_identities(
        resolved_a,
        resolved_b,
        policy=MatchPolicy(overlap_threshold=0.5),
    )

    assert list(match_table.data["particle_key"]) == ["p1", "p3"]
    assert report.matched_count == 2
    assert report.unmatched_a == 1
    assert report.unmatched_b == 1
    assert report.status == "ok"


def test_relion_duplicate_keys_fail_by_policy() -> None:
    table = pd.DataFrame({"_rlnParticleName": ["p1", "p1", "p2"]})
    policy = IdentityPolicy.relion_user_columns(["_rlnParticleName"])

    with pytest.raises(IdentityResolutionError, match="failed_duplicates"):
        resolve_identity(table, source_type="relion", policy=policy)


def test_relion_low_overlap_fails_by_match_policy() -> None:
    domain_a = pd.DataFrame({"_rlnParticleName": ["p1", "p2", "p3"]})
    domain_b = pd.DataFrame({"_rlnParticleName": ["p4", "p5", "p1"]})
    policy = IdentityPolicy.relion_user_columns(["_rlnParticleName"])
    resolved_a = resolve_identity(domain_a, source_type="relion", policy=policy)
    resolved_b = resolve_identity(domain_b, source_type="relion", policy=policy)

    with pytest.raises(IdentityResolutionError, match="overlap"):
        match_resolved_identities(
            resolved_a,
            resolved_b,
            policy=MatchPolicy(overlap_threshold=0.5, low_overlap_behavior="fail"),
        )


def test_relion_refinement_columns_are_not_identity_keys() -> None:
    table = pd.DataFrame({"_rlnAngleRot": ["1", "2"]})
    policy = IdentityPolicy.relion_user_columns(["_rlnAngleRot"])

    with pytest.raises(IdentityResolutionError, match="refinement-result"):
        resolve_identity(table, source_type="relion", policy=policy)


def test_relion_coordinate_and_defocus_columns_are_not_identity_keys() -> None:
    table = pd.DataFrame({"_rlnCenteredCoordinateZAngst": ["1", "2"]})
    policy = IdentityPolicy.relion_user_columns(["_rlnCenteredCoordinateZAngst"])

    with pytest.raises(IdentityResolutionError, match="refinement-result"):
        resolve_identity(table, source_type="relion", policy=policy)

    defocus_table = pd.DataFrame({"_rlnDefocusU": ["1", "2"]})
    defocus_policy = IdentityPolicy.relion_user_columns(["_rlnDefocusU"])

    with pytest.raises(IdentityResolutionError, match="refinement-result"):
        resolve_identity(defocus_table, source_type="relion", policy=defocus_policy)


def test_relion_explicit_mapping_file_resolves_keys(tmp_path) -> None:
    table = pd.DataFrame({"_rlnParticleName": ["opaque-a", "opaque-b"]})
    mapping_file = tmp_path / "mapping.csv"
    mapping_file.write_text(
        "source_row_id,particle_key\n0,particle-001\n1,particle-002\n",
        encoding="utf-8",
    )
    policy = IdentityPolicy(
        identity_mode="explicit_mapping_file",
        mapping_file=str(mapping_file),
    )

    resolved = resolve_identity(table, source_type="relion", policy=policy)

    assert list(resolved.keys) == ["particle-001", "particle-002"]
    assert resolved.report.identity_mode == "explicit_mapping_file"


def test_relion_row_aligned_identity_is_explicit_and_records_source_row_id() -> None:
    table = pd.DataFrame({"_rlnParticleName": ["opaque-a", "opaque-b"]})
    policy = IdentityPolicy(identity_mode="row_aligned")

    resolved = resolve_identity(table, source_type="relion", policy=policy)

    assert list(resolved.keys) == ["row:0", "row:1"]
    assert resolved.report.identity_mode == "row_aligned"
    assert resolved.report.identity_columns == ("source_row_id",)


def test_row_aligned_runner_preflight_fails_on_row_count_mismatch(tmp_path) -> None:
    ref = tmp_path / "ref.star"
    mov = tmp_path / "mov.star"
    _write_relion_star(ref, ["p1", "p2"])
    _write_relion_star(mov, ["p1"])
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


def _write_relion_star(path, particle_names: list[str]) -> None:
    rows = [
        f"{index:06d}@stack.mrcs {name} {10 * index} 20 30 0.0 0.0"
        for index, name in enumerate(particle_names, start=1)
    ]
    path.write_text(
        "\n".join(
            [
                "data_particles",
                "",
                "loop_",
                "_rlnImageName #1",
                "_rlnParticleName #2",
                "_rlnAngleRot #3",
                "_rlnAngleTilt #4",
                "_rlnAnglePsi #5",
                "_rlnOriginXAngst #6",
                "_rlnOriginYAngst #7",
                *rows,
            ]
        ),
        encoding="utf-8",
    )
