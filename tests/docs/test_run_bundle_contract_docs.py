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
    # The output guide no longer carries a dedicated "full landscape JSON is
    # debug-only" paragraph; it only must not present landscape.json as the
    # production contract (checked above).
