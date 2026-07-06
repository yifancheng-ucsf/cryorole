from __future__ import annotations

from pathlib import Path


def test_docs_do_not_restore_landscape_json_as_production_contract() -> None:
    docs_text = "\n".join(
        [
            Path("docs/AGENTS.md").read_text(encoding="utf-8"),
            Path("docs/architecture.md").read_text(encoding="utf-8"),
        ]
    ).lower()

    assert "raw_landscape.npz" in docs_text
    assert "landscape.json remains the machine-readable persistence contract" not in docs_text
    assert "landscape.json as the primary production persistence contract" not in docs_text
    assert "do not restore `landscape.json` as the default production persistence contract" in docs_text
