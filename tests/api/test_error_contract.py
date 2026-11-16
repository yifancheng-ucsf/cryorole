"""Actionable service failures share the public API's typed error contract."""
import numpy as np
import pytest

from cryorole import api
import cryorole.io.job_folders as job_folders


@pytest.fixture
def completed_run(tmp_path):
    rng = np.random.default_rng(17)
    paths = [tmp_path / "ref.cs", tmp_path / "mov.cs"]
    for path in paths:
        values = np.zeros(8, dtype=[("uid", "u8"), ("alignments3D/pose", "f8", (3,))])
        values["uid"] = np.arange(100, 108)
        values["alignments3D/pose"] = rng.normal(0., .2, (8, 3))
        with path.open("wb") as handle:
            np.save(handle, values)
    run = tmp_path / "run"
    api.run(*paths, output_dir=run, no_visualize=True, k_neighbors=2)
    return paths[0], paths[1], run


def test_missing_input_inspection_has_stable_api_error(tmp_path):
    with pytest.raises(api.CryoroleError) as exc:
        api.inspect_input(tmp_path / "missing.cs")
    assert exc.value.code == "input_not_found"
    assert exc.value.remedy


def test_missing_run_input_has_stable_api_error(tmp_path):
    with pytest.raises(api.CryoroleError) as exc:
        api.run(tmp_path / "missing_ref.cs", tmp_path / "missing_mov.cs", output_dir=tmp_path / "run")
    assert exc.value.code == "input_not_found"
    assert not (tmp_path / "run").exists()


def test_existing_run_error_is_typed_and_preserves_original_bundle(completed_run):
    ref, mov, run = completed_run
    original = (run / "run_manifest.json").read_bytes()
    with pytest.raises(api.CryoroleError) as exc:
        api.run(ref, mov, output_dir=run, no_visualize=True)
    assert exc.value.code == "output_exists"
    assert (run / "run_manifest.json").read_bytes() == original


def test_invalid_canonical_fit_has_stable_api_error(completed_run):
    _, _, run = completed_run
    with pytest.raises(api.CryoroleError) as exc:
        api.canonicalize(run, fit_top_fraction=0., no_visualize=True)
    assert exc.value.code == "input_invalid"
    assert not (run / "canonical/default").exists()


def test_invalid_selection_radius_has_stable_api_error(completed_run):
    _, _, run = completed_run
    with pytest.raises(api.CryoroleError) as exc:
        api.select(run, selection_id="bad", center=(0., 0., 0.), radius=-1.)
    assert exc.value.code == "input_invalid"
    assert not (run / "selections/bad").exists()


@pytest.mark.parametrize("error", [api.CryoroleError("specific", "known problem", details={"row": 3}),
                                  api.CancelledError("test stage")])
def test_existing_typed_error_and_cancellation_are_preserved(monkeypatch, error):
    def fail(_path):
        raise error
    monkeypatch.setattr(job_folders, "inspect_input", fail)
    with pytest.raises(type(error)) as exc:
        api.inspect_input("unused")
    assert exc.value is error


@pytest.mark.parametrize("error", [TypeError("programming bug"), AssertionError("invariant bug"),
                                  RuntimeError("unexpected backend failure")])
def test_unexpected_programming_errors_are_preserved(monkeypatch, error):
    def fail(_path):
        raise error
    monkeypatch.setattr(job_folders, "inspect_input", fail)
    with pytest.raises(type(error)) as exc:
        api.inspect_input("unused")
    assert exc.value is error
