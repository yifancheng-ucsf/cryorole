"""Test configuration for cryoROLE."""

from __future__ import annotations

import logging
from pathlib import Path

import pytest


@pytest.fixture(autouse=True)
def _isolate_cryorole_logger():
    """``cli.main`` installs a stderr handler on the ``cryorole`` logger; undo it after each test."""

    logger = logging.getLogger("cryorole")
    saved = (list(logger.handlers), logger.level, logger.propagate)
    yield
    logger.handlers[:] = saved[0]
    logger.setLevel(saved[1])
    logger.propagate = saved[2]


# --- data-derived fixtures -------------------------------------------------
#
# Some regression tests use small subsets of real CryoSPARC/RELION data under
# tests/fixtures/. These are not published yet (public tutorial datasets will
# replace them), so a public checkout skips those tests instead of failing.
# Mark a module or test with ``pytest.mark.private_fixtures("relative/path", ...)``.

FIXTURES_DIR = Path(__file__).resolve().parent / "fixtures"


def pytest_configure(config):
    config.addinivalue_line(
        "markers",
        "private_fixtures(*paths): skip unless these files exist under tests/fixtures/ "
        "(data-derived fixtures that are not published)",
    )


def pytest_collection_modifyitems(config, items):
    for item in items:
        for marker in item.iter_markers(name="private_fixtures"):
            missing = [name for name in marker.args if not (FIXTURES_DIR / name).exists()]
            if missing:
                item.add_marker(
                    pytest.mark.skip(
                        reason=f"needs unpublished data-derived fixture(s): {', '.join(missing)}"
                    )
                )
