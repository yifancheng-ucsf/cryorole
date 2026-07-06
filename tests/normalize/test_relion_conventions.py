import numpy as np
import pandas as pd
import pytest
from scipy.spatial.transform import Rotation

from cryorole.models.policies import ConventionPolicy
from cryorole.normalize import ConventionResolver, PoseNormalizer
from cryorole.normalize.conventions import (
    IMPLEMENTED_RELION_CONVERSION_RULE,
    relion_passive_euler_to_active_matrix,
)


def test_relion_uses_uppercase_intrinsic_zyz_and_not_lowercase_zyz() -> None:
    angles = [10.0, 20.0, 30.0]

    passive_upper = Rotation.from_euler("ZYZ", angles, degrees=True).as_matrix()
    passive_lower = Rotation.from_euler("zyz", angles, degrees=True).as_matrix()
    active = relion_passive_euler_to_active_matrix(*angles)

    assert not np.allclose(passive_upper, passive_lower)
    np.testing.assert_allclose(active, passive_upper.T)


def test_convention_resolver_rejects_non_uppercase_relion_sequence() -> None:
    policy = ConventionPolicy(
        source_software="relion",
        source_euler_sequence="zyz",
    )

    with pytest.raises(ValueError, match="uppercase 'ZYZ'"):
        ConventionResolver(policy)


def test_convention_resolver_is_explicitly_relion_only_for_now() -> None:
    policy = ConventionPolicy(
        source_software="cryosparc",
        source_euler_sequence="ZYZ",
    )

    with pytest.raises(ValueError, match="Unsupported source software"):
        ConventionResolver(policy)


def test_convention_policy_rule_must_match_implemented_bridge() -> None:
    policy = ConventionPolicy(
        source_software="relion",
        source_euler_sequence="ZYZ",
        conversion_rule="active_matrix = unexpected_bridge(angles)",
    )

    with pytest.raises(ValueError, match="conversion_rule does not match"):
        ConventionResolver(policy)


def test_default_relion_policy_records_implemented_bridge() -> None:
    resolver = ConventionResolver(ConventionPolicy.relion_default())

    resolution = resolver.resolve()

    assert resolution.conversion_rule_used == IMPLEMENTED_RELION_CONVERSION_RULE


def test_pose_normalizer_centralizes_passive_to_active_conversion() -> None:
    particles = pd.DataFrame(
        [
            {
                "_rlnAngleRot": "10",
                "_rlnAngleTilt": "20",
                "_rlnAnglePsi": "30",
                "_rlnOriginXAngst": "1.5",
                "_rlnOriginYAngst": "-2.0",
                "_rlnParticleName": "p001",
            }
        ]
    )
    resolver = ConventionResolver(ConventionPolicy.relion_default())

    pose_table = PoseNormalizer().normalize_relion(
        particles,
        domain_name="domain-a",
        convention_resolver=resolver,
    )

    expected = relion_passive_euler_to_active_matrix(10, 20, 30)
    np.testing.assert_allclose(
        pose_table.data.loc[0, "rotation_matrix_active"],
        expected,
    )
    assert pose_table.data.loc[0, "raw_metadata"]["_rlnParticleName"] == "p001"
    assert pose_table.data.loc[0, "raw_metadata"]["_cryorole_source_row_id"] == 0
    assert pose_table.data.loc[0, "source_type"] == "relion"
