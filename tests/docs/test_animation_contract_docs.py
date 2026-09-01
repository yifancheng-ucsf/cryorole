from pathlib import Path


def test_animation_docs_record_implemented_phase4_contract() -> None:
    root = Path(__file__).resolve().parents[2]
    animation = (root / "docs" / "animation_export.md").read_text(encoding="utf-8")
    architecture = (root / "docs" / "architecture.md").read_text(encoding="utf-8")
    assert "Phase 1–4 and presentation refinements" in animation
    assert "one/two/three structure views, implemented" in animation
    assert "`composite_rendered`" in animation
    assert "`movie_encoded`" in animation
    assert "FFprobe" in animation
    assert "never automatically removed" in animation
    assert "default: stacked" in animation
    assert "coordinate-only" in animation
    assert "default `stacked` layout" in architecture
    assert "`cryorole canonical-views`" in animation
    assert "axes_scene = S @ C" in animation
    assert "camera-only views" in architecture
    assert "Phase 1–4" in architecture
    assert "--secondary-chimerax-session" in architecture
    assert "Multiple structure views" in animation
    assert "Implemented fixed inward multi-view crop" in animation
    assert "--dual-structure-horizontal-crop" in architecture
    assert "concise visible label" in architecture
    assert "`SLD`" in architecture
    assert "underlying field and provenance remain" in architecture
    assert "`sld_display`" in architecture
    assert "--tertiary-chimerax-session" in animation
    assert "one, two, or three" in architecture
    assert "--structure-horizontal-crop" in animation
    assert "--structure-vertical-crop" in animation
    assert "58%" in architecture
