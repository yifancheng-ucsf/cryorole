from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


def test_architecture_records_workflow_ux_scientific_boundaries() -> None:
    text = (ROOT / "docs" / "architecture.md").read_text(encoding="utf-8")
    for phrase in (
        "cryorole preflight",
        "cryorole explore",
        "READY_WITH_WARNINGS",
        "draft selection",
        "exact full parent landscape",
        "127.0.0.1",
        "cryorole status",
        "cryorole next",
        "cryorole guide",
    ):
        assert phrase in text


def test_workflow_ux_tutorial_explains_required_user_concepts() -> None:
    text = (ROOT / "docs" / "workflow_ux.md").read_text(encoding="utf-8")
    for phrase in (
        "Reference domain",
        "Moving domain",
        "raw",
        "canonical",
        "rotation vector",
        "Euler",
        "SLD",
        "display downsampling",
        "SO(3)",
        "source hash",
    ):
        assert phrase in text


def test_installation_states_explore_has_no_extra_python_dependency() -> None:
    text = (ROOT / "docs" / "installation.md").read_text(encoding="utf-8")
    assert "explore" in text
    assert "no additional Python dependency" in text
