"""Preflight must describe the bytes it read and preserve assessed run options."""
from pathlib import Path
import json
import shlex

import numpy as np
import pytest

from cryorole.cli.main import build_parser, run_command
from cryorole.io.writers.landscape_store import read_landscape_npz_arrays
from cryorole.preflight import PreflightRequest, run_preflight
from cryorole.workflows.run_pipeline import preflight_and_normalize_matched_arrays
from cryorole.workflows.run_service import RunExecutionContext, RunRequest, execute_run


def inputs(tmp_path):
    paths = [tmp_path / 'ref.cs', tmp_path / 'mov.cs']
    for seed, path in enumerate(paths):
        values = np.zeros(6, dtype=[('uid', '<u8'), ('alignments3D/pose', '<f8', (3,))])
        values['uid'] = np.arange(1, 7)
        values['alignments3D/pose'] = np.random.default_rng(seed).normal(0, .2, (6, 3))
        with path.open('wb') as stream:
            np.save(stream, values)
    return paths


@pytest.mark.parametrize('change', ['alias_before_read', 'content_before_read', 'content_after_read', 'delete_after_read'])
def test_preflight_blocks_changes_while_loading(tmp_path, change):
    ref, mov = inputs(tmp_path)
    alternate, alias = tmp_path / 'alternate.cs', tmp_path / 'alias.cs'
    values = np.load(ref).copy()[::-1]
    with alternate.open('wb') as stream:
        np.save(stream, values)
    alias.symlink_to(ref)
    chosen = alias if change.startswith('alias') else ref
    called = []

    def read_then_change(*args, **kwargs):
        called.append(True)
        if change == 'alias_before_read':
            alias.unlink()
            alias.symlink_to(alternate)
        elif change == 'content_before_read':
            ref.write_bytes(alternate.read_bytes())
        result = preflight_and_normalize_matched_arrays(*args, **kwargs)
        if change == 'content_after_read':
            ref.write_bytes(alternate.read_bytes())
        elif change == 'delete_after_read':
            ref.unlink()
        return result

    output = tmp_path / 'no-bundle'
    result = run_preflight(PreflightRequest(ref=chosen, mov=mov, output_dir=output, check_environment=False),
                           array_preflight_fn=read_then_change)
    assert called
    assert result.report['readiness'] == 'BLOCKED'
    assert result.exit_code == 2
    assert result.report['errors']
    assert result.array_preflight is None
    assert not output.exists()


@pytest.mark.parametrize('batch', [1, 25_000, 100_000])
def test_suggested_command_preserves_all_assessed_options(tmp_path, batch):
    ref, mov = inputs(tmp_path)
    request = PreflightRequest(ref=ref, mov=mov, output_dir=tmp_path / 'run with spaces',
        ref_domain='fixed domain', mov_domain='mobile domain', k_neighbors=2,
        run_backend='dataframe_compat', density_query_batch_size=batch,
        raw_csv=False, visualize=False, sld_metric='so3_geodesic', row_aligned=True,
        allow_low_overlap=True, check_environment=False)
    result = run_preflight(request)
    assert result.report['readiness'] != 'BLOCKED'
    parsed = build_parser().parse_args(shlex.split(result.report['resolved_run_command'])[1:])
    for name in ['ref_domain', 'mov_domain', 'k_neighbors', 'run_backend', 'density_query_batch_size',
                 'sld_metric', 'row_aligned', 'allow_low_overlap']:
        assert getattr(parsed, name) == getattr(request, name)
    assert parsed.no_raw_csv and parsed.no_visualize
    assert Path(parsed.output_dir) == request.output_dir.resolve()


def test_implicit_row_aligned_policy_is_reported_consistently(tmp_path):
    ref, mov = inputs(tmp_path)
    result = run_preflight(PreflightRequest(ref=ref, mov=mov, identity_mode='row_aligned',
                                           check_environment=False))
    assert result.report['identity_policy']['row_aligned']
    assert result.report['matching']['row_aligned']


def test_executing_recommendation_matches_the_requested_density(tmp_path):
    ref, mov = inputs(tmp_path)
    output, direct = tmp_path / 'suggested', tmp_path / 'direct'
    request = PreflightRequest(ref=ref, mov=mov, output_dir=output, k_neighbors=2,
        run_backend='dataframe_compat', density_query_batch_size=3, raw_csv=False, visualize=False,
        ref_domain='fixed', mov_domain='mobile', check_environment=False)
    result = run_preflight(request)
    assert result.report['readiness'] != 'BLOCKED'
    args = build_parser().parse_args(shlex.split(result.report['recommended_next_command'])[1:])
    assert run_command(args) == 0
    expected = RunRequest(ref=str(ref), mov=str(mov), output_dir=str(direct),
        k_neighbors=2, run_backend='dataframe_compat', density_query_batch_size=3,
        no_raw_csv=True, no_visualize=True, ref_domain='fixed', mov_domain='mobile', quiet=True)
    execute_run(expected, context=RunExecutionContext(stream_progress=False))
    first, second = [read_landscape_npz_arrays(path / 'data/raw_landscape.npz') for path in [output, direct]]
    for name in ['particle_key', 'coordinates_analysis', 'sld_raw', 'ref_source_row_id', 'mov_source_row_id']:
        np.testing.assert_array_equal(getattr(first, name), getattr(second, name))
    summary = json.loads((output / 'run_summary.json').read_text())
    assert summary['k_neighbors'] == 2 and summary['run_backend_resolved'] == 'dataframe_compat'
    assert summary['ref_domain'] == 'fixed' and summary['mov_domain'] == 'mobile'
    assert not (output / 'data/raw_landscape.csv').exists()
