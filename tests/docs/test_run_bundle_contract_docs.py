from __future__ import annotations

from pathlib import Path


def test_docs_do_not_restore_landscape_json_as_production_contract() -> None:
    docs_text = "\n".join(
        [
            Path("README.md").read_text(encoding="utf-8"),
            Path("docs/output_files.md").read_text(encoding="utf-8"),
        ]
    ).lower()

    assert "raw_landscape.npz" in docs_text
    assert "landscape.json remains the machine-readable persistence contract" not in docs_text
    assert "landscape.json as the primary production persistence contract" not in docs_text
    assert "full landscape json is debug-only and opt-in" in docs_text
    assert "`landscape.json` is not the\ndefault production landscape" in docs_text
