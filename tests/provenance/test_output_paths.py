"""Directory replacement must preserve both source names and resolved targets."""

import pytest

from cryorole.provenance.output_paths import protect_input_directory


@pytest.mark.parametrize('kind', ['inside', 'symlink_outside', 'symlink_inside', 'parent_symlink'])
def test_directory_guard_catches_input_paths_and_targets(tmp_path, kind):
    output = tmp_path / 'run'
    output.mkdir()
    inside, outside = output / 'source.cs', tmp_path / 'outside.cs'
    outside.write_text('source')
    if kind == 'inside':
        inside.write_text('source')
        source = inside
    elif kind == 'symlink_outside':
        inside.write_text('source')
        source = tmp_path / 'alias.cs'
        source.symlink_to(inside)
    elif kind == 'symlink_inside':
        source = inside
        source.symlink_to(outside)
    else:
        alias = tmp_path / 'alias-dir'
        alias.symlink_to(output, target_is_directory=True)
        source = alias / 'source.cs'
        inside.write_text('source')
    with pytest.raises(ValueError, match='contains an input source'):
        protect_input_directory(output, [source])


def test_output_below_source_parent_is_allowed(tmp_path):
    source = tmp_path / 'source.star'
    source.write_text('source')
    protect_input_directory(tmp_path / 'outputs', [source, None])
