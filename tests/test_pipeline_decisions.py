"""Synthetic chain evidence tests; these fixtures make no biological claims."""
from copy import deepcopy
import json

import pytest

from biocentric import discovery_benchmark as benchmark, pipeline_decisions as decisions
from biocentric.pipeline_map import graph_data, graph_html


def fixture():
    pair = {'id': 'one', 'docking': {'tool': 'vina', 'parameters': {'seed': 42}}, 'scoring': {'tool': 'vina_score', 'parameters': {'cpu': 1}}}
    records = [{'compound_id': f'C{i}', 'p_value': float(i), 'endpoint': 'Ki', 'cohort': 'one'} for i in range(5, 13)]
    review = {'records': records, 'reference': {'component': 'LIG'}, 'candidates': [pair],
              'protocol': {'preparation': 'Meeko', 'pocket': 'Crystal box'}, 'identity': {'accession': 'TEST'},
              'selection': {'docking': [pair['docking']], 'scoring': [pair['scoring']]}}
    review['sha256'] = benchmark.digest(review)
    rows = [{'combination': 'one', 'docking': 'vina', 'scoring': 'vina_score', 'kind': 'measured',
             'status': 'completed', 'score': -r['p_value'], 'physical_pass': True, **r} for r in records]
    rows += [{'combination': 'one', 'docking': 'vina', 'scoring': 'vina_score', 'kind': 'reference', 'compound_id': 'crystal:LIG',
              'status': 'completed', 'score': -8., 'physical_pass': True, 'rmsd_angstrom': 1., 'reference_comparable': True}]
    return review, {'status': 'completed', 'rows': rows, 'plan_sha256': review['sha256']}


def test_complete_chain_color_and_orientation():
    review, state = fixture()
    row = decisions.compare(review, state)[0]
    assert row['lockable'] and row['rank'] == 1 and row['color'] == 'green'
    assert row['rho'] == 1 and row['ci_low'] > 0
    state['status'] = 'running'
    assert not decisions.compare(review, state)[0]['lockable']


@pytest.mark.parametrize('case', ['missing', 'duplicate', 'wrong_label', 'wrong_tool', 'failure', 'nonfinite', 'physical', 'rmsd', 'identity', 'constant', 'extra', 'cancelled'])
def test_bad_or_incomplete_evidence_never_qualifies(case):
    review, state = fixture()
    if case == 'missing': state['rows'].pop(0)
    if case == 'duplicate': state['rows'][0] = deepcopy(state['rows'][1])
    if case == 'wrong_label': state['rows'][0]['p_value'] = 1.
    if case == 'wrong_tool': state['rows'][0]['scoring'] = 'gnina_affinity'
    if case == 'failure': state['rows'][0]['status'] = 'failed'
    if case == 'nonfinite': state['rows'][0]['score'] = float('nan')
    if case == 'physical': state['rows'][0]['physical_pass'] = False
    if case == 'rmsd': state['rows'][-1]['rmsd_angstrom'] = 3.
    if case == 'identity': state['rows'][-1]['reference_comparable'] = False
    if case == 'constant':
        for r in state['rows']: r['score'] = -1.
    if case == 'extra': state['rows'].append(deepcopy(state['rows'][0]))
    if case == 'cancelled': state['status'] = 'cancelled'
    row = decisions.compare(review, state)[0]
    assert not row['lockable'] and row['rank'] is None and row['color'] != 'green'


def test_empty_evidence_is_gray_and_small_sample_is_not_green():
    review, state = fixture()
    assert decisions.compare(review, {**state, 'rows': []})[0]['color'] == 'gray'
    review['records'] = review['records'][:3]
    state['rows'] = state['rows'][:3] + [state['rows'][-1]]
    row = decisions.compare(review, state)[0]
    assert row['lockable'] and row['color'] == 'yellow' and row['ci_low'] is None


def persist(folder):
    review, state = fixture()
    for i,row in enumerate(state['rows']):
        path = folder / f'{i}.sdf'
        path.write_text('synthetic test artifact')
        row.update(pose_path=path.name, pose_sha256=benchmark.science.sha256(path))
    (folder/'plan.json').write_text(json.dumps(review))
    (folder/'status.json').write_text(json.dumps(state))
    return review, state


def test_lock_preserves_parameters_and_cannot_be_replaced(tmp_path):
    review, state = persist(tmp_path)
    lock = decisions.lock_chain(tmp_path, 'one')
    assert lock['candidate'] == review['candidates'][0]
    assert decisions.read_lock(tmp_path) == lock
    with pytest.raises(ValueError, match='already'):
        decisions.lock_chain(tmp_path, 'one')
    state['rows'][0]['score'] -= 1
    (tmp_path/'status.json').write_text(json.dumps(state))
    with pytest.raises(ValueError, match='changed'):
        decisions.read_lock(tmp_path)


@pytest.mark.parametrize('case', ['pose', 'outside', 'plan'])
def test_lock_requires_matching_evidence_and_pose_integrity(tmp_path, case):
    _, state = persist(tmp_path)
    if case == 'pose': (tmp_path/'0.sdf').write_text('modified')
    if case == 'outside': state['rows'][0]['pose_path'] = '../outside.sdf'
    if case == 'plan': state['plan_sha256'] = 'wrong'
    (tmp_path/'status.json').write_text(json.dumps(state))
    with pytest.raises(ValueError): decisions.lock_chain(tmp_path, 'one')
    assert not (tmp_path/'decision.json').exists()


def test_all_connections_and_only_executable_path_glows():
    review, state = fixture()
    selected = {**review['selection'], 'structure_preparation': [{'tool': 'pdbfixer'}], 'pocket_identification': [{'tool': 'fpocket'}]}
    selected['docking'] = selected['docking'] + [{'tool': 'custom'}]
    state.update(status='running', active_combination='one', stage='docking')
    data = graph_data(selected, {'custom': '<script>bad()</script>'}, review, state, decisions.compare(review, state))
    assert len(data['edges']) == 1 + 2 + 2 + 1
    assert sum(e['running'] for e in data['edges']) == 4
    assert not any(n['tool'] in {'meeko_fixed', 'pdbfixer', 'fpocket'} for n in data['nodes'])
    assert [group[0]['id'] for group in data['stages'][:2]] == ['protein', 'pocket']
    assert len(data['ligands']) == 9
    assert [n['tool'] for n in data['nodes'] if n['current']] == ['vina']
    rendered = graph_html(data)
    assert '<script>bad()' not in rendered and '&lt;script&gt;' in rendered
    assert 'prefers-reduced-motion' in rendered and 'textContent' in rendered
    assert 'http://' not in rendered and 'https://' not in rendered
    # Finished paths stop pulsing and remain visible after locking.
    state['status'] = 'completed'
    locked = graph_data(selected, {}, review, state, decisions.compare(review, state), {'candidate': review['candidates'][0]})
    assert not any(e['running'] for e in locked['edges'])
    assert sum(e['fixed'] for e in locked['edges']) == 4
    assert len(locked['edges']) == len(data['edges'])


def test_final_ui_locks_chain_and_opens_reference_visualization(tmp_path, monkeypatch):
    from streamlit.testing.v1 import AppTest
    from biocentric import discovery_evaluation_ui
    review, state = persist(tmp_path)
    review['reference']['pdb_id'] = 'SYNTHETIC'
    review.pop('sha256')
    review['sha256'] = benchmark.digest(review)
    state['plan_sha256'] = review['sha256']
    (tmp_path/'plan.json').write_text(json.dumps(review))
    (tmp_path/'status.json').write_text(json.dumps(state))
    monkeypatch.setattr(discovery_evaluation_ui, '_scene', lambda *args: '<div>Synthetic viewer test</div>')
    app = AppTest.from_string('''import streamlit as st
from biocentric.pipeline_results_ui import render_map, render_results, render_final_visualizations
st.session_state['studio_eval_run'] = ''' + repr(str(tmp_path)) + '''
render_map(st)
render_results(st)
render_final_visualizations(st)
''').run()
    assert not app.exception
    app.button(key='studio_lock_final').click().run()
    assert not app.exception
    assert decisions.read_lock(tmp_path)['candidate']['id'] == 'one'
    assert app.selectbox(key='studio_final_complex_' + tmp_path.name).value == 'crystal:LIG'
    assert any('Fixed chain' in s.value for s in app.success)
    assert any('cohort' in frame.value.columns for frame in app.dataframe)


def test_failed_crystal_keeps_full_test_set_metrics_and_best_available(tmp_path):
    review, state = persist(tmp_path)
    state['rows'][-1].update(status='failed', score=None, physical_pass=False, rmsd_angstrom=None)
    (tmp_path / 'status.json').write_text(json.dumps(state))
    result = decisions.compare(review, state)[0]
    assert not result['lockable'] and result['rho'] == 1 and result['completed'] == 8
    lock = decisions.finalize_best_available(tmp_path)
    assert lock['selection_kind'] == 'correlation_first' and not lock['checks_passed']
    assert lock['metrics']['rho'] == 1 and lock['failed_checks']
    assert decisions.read_lock(tmp_path) == lock
    assert decisions.finalize_best_available(tmp_path) == lock


def test_corrupt_labels_never_become_best_available(tmp_path):
    _, state = persist(tmp_path)
    state['rows'][0]['p_value'] = 999
    (tmp_path / 'status.json').write_text(json.dumps(state))
    assert decisions.finalize_best_available(tmp_path) is None


def test_failed_final_ui_displays_fixed_chain_and_unavailable_metrics(tmp_path):
    from streamlit.testing.v1 import AppTest
    _, state = persist(tmp_path)
    for row in state['rows']:
        row.update(status='failed', score=None, rmsd_angstrom=None)
    (tmp_path / 'status.json').write_text(json.dumps(state))
    app = AppTest.from_string('''import streamlit as st
from biocentric.pipeline_results_ui import render_final_chain, render_final_visualizations
st.session_state['studio_eval_run'] = ''' + repr(str(tmp_path)) + '''
render_final_chain(st)
render_final_visualizations(st)
''').run()
    assert not app.exception
    assert any('Best available' in item.value for item in app.warning)
    assert app.metric[0].value == 'Unavailable'


def test_perfect_rho_reports_only_rankable_assay_measurements():
    rows = [{'cohort': 'rankable', 'score': float(i), 'p_value': float(i)} for i in range(4)]
    rows += [{'cohort': 'pair', 'score': 3., 'p_value': 1.}, {'cohort': 'pair', 'score': 2., 'p_value': 2.},
             {'cohort': 'single-a', 'score': 0., 'p_value': 9.}, {'cohort': 'single-b', 'score': 10., 'p_value': 0.}]
    stats = decisions.stratified_statistics(rows, 1)
    assert stats['rho'] == 1.0
    assert stats['ranked_measurements'] == 4 and stats['ranked_cohorts'] == 1 and stats['total_cohorts'] == 4
    assert stats['small_sample'] and stats['ci_low'] is None
    excluded = [row for row in stats['cohort_metrics'] if row['exclusion']]
    assert len(excluded) == 3 and all(row['exclusion'] == 'Fewer than 3 measurements' for row in excluded)


def partial_comparison_fixture():
    review, state = fixture()
    other = deepcopy(review['candidates'][0])
    other.update(id='bio', scoring={'tool': 'biocentric_affinity', 'parameters': {}})
    review['candidates'].append(other)
    copies = []
    for row in state['rows']:
        copied = {**row, 'combination': 'bio', 'scoring': 'biocentric_affinity', 'score': -row['score']}
        copies.append(copied)
    # BioCentric fails one ligand; Vina's first score is deliberately an outlier.
    copies[0].update(status='failed', score=None, error='synthetic ligand failure')
    state['rows'][0]['score'] = -100.
    state['rows'][2]['score'], state['rows'][3]['score'] = state['rows'][3]['score'], state['rows'][2]['score']
    state['rows'].extend(copies)
    return review, state


def test_failed_ligand_is_skipped_and_all_chains_use_identical_subset():
    review, state = partial_comparison_fixture()
    rows = {r['id']: r for r in decisions.compare(review, state)}
    bio, vina = rows['bio'], rows['one']
    assert bio['skipped'] == 1 and bio['completed'] == 7
    assert bio['selection_eligible'] and not bio['lockable']
    assert bio['shared_measurement_ids'] == vina['shared_measurement_ids']
    assert bio['shared_measurements'] == vina['shared_measurements'] == 7
    assert bio['rho'] == pytest.approx(1.) and bio['rho'] > vina['rho']
    assert vina['own_rho'] != vina['rho']
    assert bio['available_rank'] == 1
    from biocentric.bayesian_search import configuration, observe
    review['search'] = configuration(review['candidates'])
    observation = observe(review, state, review['candidates'][1])
    assert observation['feasible'] and observation['objective'] == pytest.approx(1.)
    assert observation['skipped'] == 1


def test_zero_success_chain_does_not_erase_other_chains_metrics():
    review, state = partial_comparison_fixture()
    for row in state['rows']:
        if row['combination'] == 'bio': row.update(status='failed', score=None)
    rows = {r['id']: r for r in decisions.compare(review, state)}
    assert rows['one']['shared_measurements'] == 8 and decisions.finite(rows['one']['rho'])
    assert rows['bio']['rho'] is None and not rows['bio']['selection_eligible']


def test_disjoint_success_sets_do_not_produce_comparable_correlation():
    review, state = partial_comparison_fixture()
    measured = [r for r in state['rows'] if r['kind'] == 'measured']
    for row in measured:
        i = int(row['compound_id'][1:])
        if (row['combination'] == 'one' and i >= 9) or (row['combination'] == 'bio' and i < 9):
            row.update(status='failed', score=None)
    rows = decisions.compare(review, state)
    assert all(r['shared_measurements'] == 0 and r['rho'] is None for r in rows)
    assert all(decisions.finite(r['own_rho']) for r in rows)


@pytest.mark.parametrize('tampered_pose', [False, True])
def test_correlation_wins_despite_pose_flags_and_revision_preserves_evidence(tmp_path, tampered_pose):
    review, state = partial_comparison_fixture()
    next(r for r in state['rows'] if r['combination'] == 'bio' and r['status'] == 'completed')['physical_pass'] = False
    review['sha256'] = benchmark.digest({k: v for k,v in review.items() if k != 'sha256'})
    state['plan_sha256'] = review['sha256']
    for i, row in enumerate(state['rows']):
        file = tmp_path / f'{i}.sdf'
        file.write_text('synthetic pose')
        row.update(pose_path=file.name, pose_sha256=benchmark.science.sha256(file))
    (tmp_path / 'plan.json').write_text(json.dumps(review))
    (tmp_path / 'status.json').write_text(json.dumps(state))
    old = decisions.lock_chain(tmp_path, 'one')
    # A previously saved physical-first decision remains readable.
    old.update(policy=decisions.SHARED_POLICY, selection_kind='shared_subset')
    old['sha256'] = benchmark.digest({k:v for k,v in old.items() if k != 'sha256'})
    (tmp_path / 'decision.json').write_text(json.dumps(old))
    before = {name: (tmp_path / name).read_bytes() for name in ('plan.json', 'status.json', 'decision.json')}
    bio = next(r for r in decisions.compare(review, state) if r['id'] == 'bio')
    assert not bio['selection_eligible'] and bio['available_rank'] == 1
    from biocentric.bayesian_search import configuration, observe
    search_review = {**review, 'search': configuration(review['candidates'])}
    observation = observe(search_review, state, review['candidates'][1])
    assert observation['objective'] == pytest.approx(1) and not observation['physical_and_crystal_checks_passed']
    if tampered_pose:
        row = next(r for r in state['rows'] if r['combination'] == 'bio' and r['status'] == 'completed')
        (tmp_path / row['pose_path']).write_text('changed')
        with pytest.raises(ValueError, match='integrity'):
            decisions.revise_best_selection(tmp_path, 'User requested correlation-first selection', old['sha256'])
        assert (tmp_path / 'decision.json').read_bytes() == before['decision.json']
        return
    new = decisions.revise_best_selection(tmp_path, 'User requested correlation-first selection', old['sha256'])
    assert new['candidate']['id'] == 'bio' and new['rank_on_this_set'] == 1
    assert not new['checks_passed'] and new['failed_checks']
    assert new['metrics']['rho'] == pytest.approx(1)
    assert decisions.read_lock(tmp_path) == new
    assert (tmp_path / 'decision-history' / (old['sha256'] + '.json')).read_bytes() == before['decision.json']
    assert all((tmp_path / name).read_bytes() == before[name] for name in ('plan.json', 'status.json'))
    with pytest.raises(ValueError, match='Previous decision changed'):
        decisions.revise_best_selection(tmp_path, 'stale request', old['sha256'])
