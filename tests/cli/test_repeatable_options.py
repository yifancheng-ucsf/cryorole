"""Repeatable CLI flags must preserve immutable defaults across parser reuse."""
import pytest

from cryorole.cli.main import build_parser


def test_repeatable_float_tolerances_parse_and_do_not_change_defaults():
    parser=build_parser()
    request=parser.parse_args(['align','--ref','ref.star','--mov','mov.star',
        '--float-tol','_rlnCoordinateX=0.1','--float-tol','_rlnCoordinateY=0.2'])
    assert request.float_tol==(('_rlnCoordinateX',.1),('_rlnCoordinateY',.2))
    assert parser.parse_args(['align','--ref','ref.star','--mov','mov.star']).float_tol==()


@pytest.mark.parametrize('command',['run','preflight'])
def test_repeatable_identity_columns_parse_and_do_not_change_defaults(command):
    parser=build_parser()
    base=[command,'--ref','ref.star','--mov','mov.star']
    request=parser.parse_args([*base,'--identity-column','_rlnMicrographName',
                              '--identity-column','_rlnCoordinateX'])
    assert request.identity_column==('_rlnMicrographName','_rlnCoordinateX')
    assert parser.parse_args(base).identity_column==()


@pytest.mark.parametrize('value',['nan','inf','-inf'])
def test_nonfinite_float_tolerance_is_a_usage_error(value):
    with pytest.raises(SystemExit)as error:
        build_parser().parse_args(['align','--ref','ref.star','--mov','mov.star','--float-tol',f'_coord={value}'])
    assert error.value.code==2
