"""Peak-memory model of ``estimate_run_resources`` (fitted to measured runs, 2026-09-28)."""

from __future__ import annotations

import pytest

from cryorole.preflight.resource_estimate import (
    COMPACT_ARRAY_BYTES_PER_PARTICLE,
    CS_INPUT_RESIDENT_FRACTION,
    FIXED_RUNTIME_BYTES,
    MIB,
    PREVIEW_BYTES_PER_PLOTTED_PARTICLE,
    PREVIEW_FIXED_BYTES,
    PREVIEW_MAX_PLOTTED_PARTICLES,
    STAR_INPUT_RESIDENT_FACTOR,
    estimate_run_resources,
)


def _estimate(tmp_path, **overrides):
    options = dict(
        matched_count=200_000, k_neighbors=50, query_batch_size=100_000, raw_csv=True,
        visualize=False, output_dir=tmp_path / "run", input_paths=(),
    )
    options.update(overrides)
    return estimate_run_resources(**options)


def _sized(path, size):
    with path.open("wb") as handle:
        handle.truncate(size)
    return path


def test_peak_memory_is_the_sum_of_its_documented_terms(tmp_path) -> None:
    ref = _sized(tmp_path / "ref.cs", 40 * MIB)
    mov = _sized(tmp_path / "mov.star", 10 * MIB)
    n = 800_000
    result = _estimate(tmp_path, matched_count=n, visualize=True, input_paths=(ref, mov))
    bounded_batch = min(100_000, 2_000_000 // 51)
    expected = (
        FIXED_RUNTIME_BYTES
        + n * COMPACT_ARRAY_BYTES_PER_PARTICLE
        + bounded_batch * 51 * 16
        + int(40 * MIB * CS_INPUT_RESIDENT_FRACTION + 10 * MIB * STAR_INPUT_RESIDENT_FACTOR)
        + PREVIEW_FIXED_BYTES
        + PREVIEW_MAX_PLOTTED_PARTICLES * PREVIEW_BYTES_PER_PLOTTED_PARTICLE
    )
    assert result["estimated_peak_memory_bytes"] == expected
    assumptions = result["assumptions"]
    assert assumptions["cs_input_bytes"] == 40 * MIB
    assert assumptions["star_input_bytes"] == 10 * MIB
    assert assumptions["preview_memory_bytes"] > 0


def test_figures_and_inputs_raise_the_estimate(tmp_path) -> None:
    base = _estimate(tmp_path)["estimated_peak_memory_bytes"]
    with_figures = _estimate(tmp_path, visualize=True)["estimated_peak_memory_bytes"]
    star = _sized(tmp_path / "in.star", 100 * MIB)
    with_star = _estimate(tmp_path, input_paths=(star,))["estimated_peak_memory_bytes"]
    assert with_figures - base == PREVIEW_FIXED_BYTES + 200_000 * PREVIEW_BYTES_PER_PLOTTED_PARTICLE
    assert with_star - base == int(100 * MIB * STAR_INPUT_RESIDENT_FACTOR)


def test_missing_input_paths_are_ignored(tmp_path) -> None:
    result = _estimate(tmp_path, input_paths=(tmp_path / "absent.cs",))
    assert result["assumptions"]["input_resident_bytes"] == 0


@pytest.mark.parametrize(
    ("n", "visualize", "measured_mib"),
    [(100_000, False, 225.5), (500_000, False, 474.4), (1_000_000, True, 1543.6)],
)
def test_model_stays_within_30_percent_of_measured_synthetic_runs(tmp_path, n, visualize, measured_mib) -> None:
    """Anchors measured with benchmark_full_run inputs (uid + pose, 20 bytes per row per file)."""

    ref = _sized(tmp_path / "ref.cs", 20 * n)
    mov = _sized(tmp_path / "mov.cs", 20 * n)
    estimate = _estimate(tmp_path, matched_count=n, visualize=visualize, input_paths=(ref, mov))
    assert estimate["estimated_peak_memory_mib"] == pytest.approx(measured_mib, rel=0.30)
