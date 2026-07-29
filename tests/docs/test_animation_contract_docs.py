from pathlib import Path


def test_animation_docs_record_implemented_phase4_contract() -> None:
    root = Path(__file__).resolve().parents[2]
    animation = (root / "docs" / "animation_export.md").read_text(encoding="utf-8")
    architecture = (root / "docs" / "architecture.md").read_text(encoding="utf-8")
    assert "Phase 1–4 and presentation refinements" in animation
    assert "fixed inward dual-view crop, implemented" in animation
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
    assert "Dual structure views" in animation
    assert "Implemented fixed inward dual-view crop" in animation
    assert "--dual-structure-horizontal-crop" in architecture
    assert "concise visible label" in architecture
    assert "`SLD`" in architecture
    assert "underlying field and provenance remain" in architecture
    assert "`sld_display`" in architecture
