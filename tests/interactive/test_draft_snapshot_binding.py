"""Concurrent evaluation and confirmation must use one exact draft snapshot."""
from concurrent.futures import ThreadPoolExecutor
import json
import threading
from urllib.error import HTTPError
from urllib.request import Request, urlopen

import numpy as np
import pytest

from cryorole.errors import CryoroleError
from cryorole.interactive import ExploreSession, create_explore_server
from cryorole.interactive.landscape_view import LandscapeView
from cryorole.select import evaluate_radius_mask
from tests.interactive.test_interactive_selection import _bundle


def test_overlapping_evaluations_return_their_own_complete_snapshots(tmp_path):
    first_written,release_first=threading.Event(),threading.Event()

    class OverlapSession(ExploreSession):
        def __setattr__(self,name,value):
            super().__setattr__(name,value)
            if name=='_draft' and isinstance(value,dict) and value.get('radius')==1. and getattr(self,'armed',False):
                first_written.set()
                assert release_first.wait(5)

    session=OverlapSession(_bundle(tmp_path))
    session.armed=True
    with ThreadPoolExecutor(max_workers=1)as pool:
        future=pool.submit(session.evaluate,center=[0,0,0],representation='rotvec',radius_deg=1.)
        try:
            assert first_written.wait(5)
            second=session.evaluate(center=[0,0,.5],representation='rotvec',radius_deg=2.)
        finally:
            release_first.set()
        first=future.result(timeout=5)
    assert first['radius']==1. and first['center_input']==[0,0,0]
    assert second['radius']==2. and second['center_input']==[0,0,.5]
    expected,_=evaluate_radius_mask(session.coordinates,np.array([0,0,.5]),radius=2.,radius_unit='degrees')
    np.testing.assert_array_equal(session.draft_selected_mask(),expected)
    assert first['draft_id']!=second['draft_id']
    result=session.confirm('latest',expected_draft_id=second['draft_id'])
    assert result['selected_count']==int(expected.sum())


def test_confirm_keeps_the_snapshot_approved_before_a_new_evaluation(tmp_path,monkeypatch):
    session=ExploreSession(_bundle(tmp_path))
    first=session.evaluate(center=[0,0,0],representation='rotvec',radius_deg=1.)
    validate=session._validate_identity
    switched=[]

    def evaluate_during_revalidation(*args,**kwargs):
        validate(*args,**kwargs)
        if kwargs.get('rehash') and not switched:
            switched.append(True)
            session.evaluate(center=[0,0,.5],representation='rotvec',radius_deg=2.)

    monkeypatch.setattr(session,'_validate_identity',evaluate_during_revalidation)
    session.confirm('approved')
    data=json.loads((session.run_dir/'selections/approved/selection.json').read_text())
    assert data['center_input']==first['center_input']
    assert data['radius']==first['radius']
    assert data['selected_count']==first['selected_count']


def test_stale_landscape_view_draft_cannot_confirm_another_session_result(tmp_path):
    view=LandscapeView(_bundle(tmp_path))
    draft=view.evaluate_radius([0,0,0],representation='rotvec',radius_deg=1.)
    view.session.evaluate(center=[0,0,.5],representation='rotvec',radius_deg=2.)
    with pytest.raises(CryoroleError,match='draft|changed'):
        view.confirm('stale',draft)
    assert not (view.run_dir/'selections/stale').exists()


def test_view_refuses_a_mask_replaced_between_evaluation_and_display(tmp_path,monkeypatch):
    view=LandscapeView(_bundle(tmp_path))
    get_mask=view.session.draft_selected_mask

    def replace_before_copy(**kwargs):
        view.session.evaluate(center=[0,0,.5],representation='rotvec',radius_deg=2.)
        return get_mask(**kwargs)

    monkeypatch.setattr(view.session,'draft_selected_mask',replace_before_copy)
    with pytest.raises(CryoroleError,match='draft|changed'):
        view.evaluate_radius([0,0,0],representation='rotvec',radius_deg=1.)
    assert view._draft_token is None


def test_initial_display_top_fraction_retains_cutoff_ties(tmp_path):
    run=_bundle(tmp_path)
    path=run/'data/raw_landscape.npz'
    with np.load(path)as source:
        data={key:source[key]for key in source.files}
    data['sld_display']=np.ones(20)
    np.savez_compressed(path,**data)
    session=ExploreSession(run,display_top_fraction=.4,max_display_points=20)
    assert session.display_indices.tolist()==list(range(20))
    draft=session.evaluate(center=[0,0,.2],representation='rotvec',radius_deg=5.)
    assert draft['full_candidate_count']==20


@pytest.mark.parametrize('supplied',['missing','old'])
def test_http_requires_the_current_evaluated_snapshot(tmp_path,supplied):
    session=ExploreSession(_bundle(tmp_path))
    server=create_explore_server(session)
    server.start_in_thread()
    identity=dict(session_token=session.session_token,run_id=session.run_id,landscape_sha256=session.landscape_sha256)

    def post(path,body):
        request=Request(server.url+'api/'+path,data=json.dumps({**identity,**body}).encode(),
                        headers={'Content-Type':'application/json'},method='POST')
        with urlopen(request,timeout=5)as response:
            return json.load(response)

    try:
        first=post('evaluate',dict(center=[0,0,0],representation='rotvec',radius_deg=1.))
        second=post('evaluate',dict(center=[0,0,.5],representation='rotvec',radius_deg=2.))
        payload={'selection_id':'stale'}
        if supplied=='old':
            payload['draft_id']=first.get('draft_id','old-evaluated-snapshot')
        with pytest.raises(HTTPError)as error:
            post('confirm',payload)
        assert error.value.code==400
        assert not (session.run_dir/'selections/stale').exists()
        result=post('confirm',dict(selection_id='current',draft_id=second['draft_id']))
        assert result['selected_count']==second['selected_count']
    finally:
        server.shutdown()
