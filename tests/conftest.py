"""Test configuration for cryoROLE."""

from __future__ import annotations

import logging

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
