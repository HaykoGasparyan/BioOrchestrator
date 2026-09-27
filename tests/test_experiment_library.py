"""Workspace isolation tests use synthetic worker files, never biological results."""
from copy import deepcopy
import json
from pathlib import Path

import pytest

from biocentric import experiment_library as library, discovery_benchmark as benchmark, pipeline_decisions as decisions
from biocentric.discovery_studio import selection


def run_fixture(root, run_id='a'*32, status='running'):
    folder = root / 'evaluation-runs' / run_id
    folder.mkdir(parents=True)
    selected = selection({})
    pair = {'id': 'pair', 'docking': selected['docking'][0], 'scoring': selected['scoring'][0]}
    records = [{'compound_id': f'C{i}', 'p_value': float(i), 'endpoint': 'Ki', 'cohort': 'one'} for i in range(5,8)]
    review = {'records': records, 'selection': selected, 'candidates': [pair], 'reference': {'component': 'LIG'},
              'identity': {'accession': 'P00519', 'gene_names': ['ABL1'], 'organism': 'Homo sapiens'},
              'protocol': {'preparation': 'test preparation', 'pocket': 'test pocket'}, 'created_at': '2026-09-26T00:00:00Z'}
    review['sha256'] = benchmark.digest(review)
    rows = []
    if status == 'completed':
        for r in [*records, {'compound_id': 'crystal:LIG', 'p_value': None, 'endpoint': None}]:
            pose = folder / (r['compound_id'].replace(':','_') + '.sdf')
            pose.write_text('Synthetic integrity fixture')
            rows.append({**r, 'kind': 'reference' if r['p_value'] is None else 'measured', 'combination': 'pair',
                         'docking': 'vina', 'scoring': 'vina_score', 'score': -r['p_value'] if r['p_value'] else -8,
                         'status': 'completed', 'physical_pass': True, 'rmsd_angstrom': 1., 'reference_comparable': True,
                         'pose_path': pose.name, 'pose_sha256': benchmark.science.sha256(pose)})
    state = {'status': status, 'rows': rows, 'completed': len(rows), 'total': 4, 'plan_sha256': review['sha256']}
    (folder/'plan.json').write_text(json.dumps(review)); (folder/'status.json').write_text(json.dumps(state))
    return folder, review, state


def test_draft_roundtrip_and_new_workspace_start_empty(tmp_path):
    one = library.create(tmp_path)
    original = {'studio_draft': {'target': 'ABL1', 'messages': [{'role':'user','text':'ABL1'}], 'evidence_query': ('ABL1','Human')},
                'studio_evidence': {'literature': {'evidence': [{'text':'Saved source'}]}},
                'studio_parameters': {'docking:vina': {'seed': 99}}, 'studio_eval_dataset': {'records': [{'compound_id':'C1'}]},
                'studio_eval_start_requested': True, 'unrelated_secret': 'never-save'}
    assert library.save(tmp_path,one,original) == one
    saved = library.load(tmp_path,one)
    assert 'unrelated_secret' not in saved['state'] and 'studio_eval_start_requested' not in saved['state']
    current = {'studio_identity': {'accession':'OTHER'}}
    library.restore(tmp_path,one,current)
    assert current['studio_evidence'] == original['studio_evidence'] and current['studio_draft']['evidence_query'] == ('ABL1','Human')
    two = library.create(tmp_path)
    library.restore(tmp_path,two,current)
    assert 'studio_draft' not in current and 'studio_parameters' not in current and 'studio_eval_dataset' not in current
    assert library.load(tmp_path,one)['state'] == saved['state']


def test_background_runs_remain_independent_during_navigation(tmp_path,monkeypatch):
    folder, review, state = run_fixture(tmp_path)
    other, _, _ = run_fixture(tmp_path, 'b'*32)
    monkeypatch.setattr(benchmark,'cancel_run',lambda *a: pytest.fail('Navigation must not stop workers'))
    library.adopt_runs(tmp_path)
    assert len(library.summaries(tmp_path)) == 2
    assert all(r['status']=='running' for r in library.summaries(tmp_path))
    before = (folder/'status.json').read_bytes()
    session = {}
    library.restore(tmp_path,folder.name,session)
    library.save(tmp_path,folder.name,session)
    library.restore(tmp_path,library.create(tmp_path),session)
    assert 'studio_eval_run' not in session and (folder/'status.json').read_bytes() == before
    assert (other/'status.json').read_text()
    state.update(status='cancelled')
    (folder/'status.json').write_text(json.dumps(state))
    assert next(r for r in library.summaries(tmp_path) if r['id']==folder.name)['status'] == 'cancelled'


def test_legacy_import_is_idempotent_and_restores_frozen_subset(tmp_path):
    folder,review,_=run_fixture(tmp_path,status='cancelled')
    library.adopt_runs(tmp_path);library.adopt_runs(tmp_path)
    assert len(library.summaries(tmp_path))==1
    session={};library.restore(tmp_path,folder.name,session)
    assert session['studio_eval_review']==review
    assert session['studio_eval_dataset']['loaded']==3
    assert session['studio_eval_dataset']['errors']
    assert session['studio_tools_docking']==[r['tool'] for r in review['selection']['docking']]


def test_locked_chain_and_metrics_survive_reopening(tmp_path):
    folder,review,_=run_fixture(tmp_path,status='completed')
    lock=decisions.lock_chain(folder,'pair')
    library.adopt_runs(tmp_path)
    row=library.summaries(tmp_path)[0]
    assert row['result']=='Final chain locked' and row['rho']==1 and row['chain']=='Crystal pocket → vina → vina_score'
    session={};library.restore(tmp_path,folder.name,session)
    assert decisions.read_lock(session['studio_eval_run'])==lock


def test_rerun_creates_new_row_without_overwriting_old_experiment(tmp_path):
    first,_,_=run_fixture(tmp_path,status='cancelled')
    second,_,_=run_fixture(tmp_path,'b'*32,status='cancelled')
    library.adopt_runs(tmp_path)
    state={};library.restore(tmp_path,first.name,state)
    state['studio_eval_run']=str(second)
    new_id=library.save(tmp_path,first.name,state)
    assert new_id!=first.name
    assert library.load(tmp_path,first.name)['run_id']==first.name
    assert library.load(tmp_path,new_id)['run_id']==second.name


def test_invalid_ids_and_damaged_workspace_do_not_hide_other_runs(tmp_path):
    run_fixture(tmp_path);library.adopt_runs(tmp_path)
    with pytest.raises(ValueError):library.load(tmp_path,'../escape')
    bad=library.create(tmp_path);library.workspace_path(tmp_path,bad).write_text('{}')
    rows=library.summaries(tmp_path)
    assert len(rows)==2 and any(r['status']=='running' for r in rows) and any(r['status']=='unavailable' for r in rows)


def test_delete_hides_running_experiment_without_losing_artifacts(tmp_path):
    folder, review, _ = run_fixture(tmp_path)
    other, _, _ = run_fixture(tmp_path, 'b' * 32)
    library.adopt_runs(tmp_path)
    before = (folder / 'plan.json').read_bytes()
    session = {}; library.restore(tmp_path, folder.name, session)
    library.delete(tmp_path, folder.name)
    assert (folder / 'cancel.request').exists()
    assert not (other / 'cancel.request').exists()
    assert (folder / 'plan.json').read_bytes() == before
    library.adopt_runs(tmp_path)
    assert [r['id'] for r in library.summaries(tmp_path)] == [other.name]
    with pytest.raises(ValueError, match='deleted'):
        library.save(tmp_path, folder.name, session)
    library.undo_delete(tmp_path, folder.name)
    assert len(library.summaries(tmp_path)) == 2
    assert library.load(tmp_path, folder.name)['run_id'] == folder.name
    assert (folder / 'cancel.request').exists()  # Undo never relaunches a worker.


def test_delete_draft_is_recoverable_and_id_is_validated(tmp_path):
    key = library.create(tmp_path)
    library.delete(tmp_path, key)
    assert library.summaries(tmp_path) == []
    with pytest.raises(ValueError):
        library.delete(tmp_path, '../outside')
    library.undo_delete(tmp_path, key)
    assert library.summaries(tmp_path)[0]['id'] == key


def test_bulk_delete_stops_only_selected_runs_and_undo_preserves_evidence(tmp_path):
    first, _, _ = run_fixture(tmp_path)
    second, _, _ = run_fixture(tmp_path, 'b'*32)
    third, _, _ = run_fixture(tmp_path, 'c'*32)
    library.adopt_runs(tmp_path)
    before = (first/'plan.json').read_bytes()
    removed, errors = library.delete_many(tmp_path, [first.name, second.name, first.name])
    assert removed == [first.name, second.name] and not errors
    assert (first/'cancel.request').exists() and (second/'cancel.request').exists()
    assert not (third/'cancel.request').exists()
    assert (first/'plan.json').read_bytes() == before
    assert [r['id'] for r in library.summaries(tmp_path)] == [third.name]
    for key in removed:
        library.undo_delete(tmp_path, key)
    assert len(library.summaries(tmp_path)) == 3
    assert (first/'cancel.request').exists()


def test_bulk_delete_prevalidates_every_id_before_any_mutation(tmp_path):
    one = library.create(tmp_path)
    with pytest.raises(ValueError):
        library.delete_many(tmp_path, [one, '../escape'])
    assert library.load(tmp_path, one)['id'] == one


def test_bulk_selection_confirmation_and_undo_in_ui(tmp_path):
    from streamlit.testing.v1 import AppTest
    ids = [library.create(tmp_path) for _ in range(3)]
    app = AppTest.from_string('''import streamlit as st
from biocentric.experiment_library_ui import render_home
render_home(st, st.session_state['root'])
''')
    app.session_state['root'] = str(tmp_path)
    app.run()
    assert not app.exception
    app.multiselect(key='experiments_bulk_selection').set_value(ids[:2]).run()
    assert not app.exception
    app.button(key='experiments_delete_selected').click().run()
    assert not app.exception and len(library.summaries(tmp_path)) == 3
    confirm = next(b for b in app.button if b.label == 'Delete 2 experiments')
    confirm.click().run()
    assert not app.exception
    assert [r['id'] for r in library.summaries(tmp_path)] == ids[2:]
    app.button(key='experiments_undo_delete').click().run()
    assert not app.exception and len(library.summaries(tmp_path)) == 3
    app.button(key='experiments_select_all').click().run()
    assert set(app.multiselect(key='experiments_bulk_selection').value) == set(ids)
    app.button(key='experiments_delete_selected').click().run()
    next(b for b in app.button if b.label == 'Keep experiments').click().run()
    assert not app.exception and len(library.summaries(tmp_path)) == 3
