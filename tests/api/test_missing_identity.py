import pytest

from cryorole import api


@pytest.mark.parametrize("missing", ["?", ".", "''"])
def test_public_run_rejects_missing_relion_identity(tmp_path, missing):
    headers = "\n".join((
        "data_particles", "", "loop_", "_rlnImageName",
        "_rlnAngleRot", "_rlnAngleTilt", "_rlnAnglePsi",
    ))
    ref, mov = tmp_path / "ref.star", tmp_path / "mov.star"
    ref.write_text(f"{headers}\n1@stack.mrcs 0 15 0\n{missing} 12 15 0\n2@stack.mrcs 5 15 0\n")
    mov.write_text(f"{headers}\n1@stack.mrcs 4 15 0\n2@stack.mrcs 15 15 0\n{missing} 80 15 0\n")
    output = tmp_path / "run"

    with pytest.raises(ValueError, match=r"missing identity.*_rlnImageName.*1"):
        api.run(ref, mov, output_dir=output, k_neighbors=1, no_visualize=True)
    assert not output.exists()
    preflight = api.preflight(ref, mov, k_neighbors=1)
    assert preflight.report["readiness"] == "BLOCKED"
