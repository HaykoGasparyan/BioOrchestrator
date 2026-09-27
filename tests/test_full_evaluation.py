"""Full-set and pocket-stage invariants; synthetic data, no predictive claims."""
from copy import deepcopy
import json

import pytest

from biocentric import discovery_benchmark as bench, pipeline_decisions as decisions
from biocentric.pipeline_map import graph_data


def methods():
    return {'pocket_identification': [{'tool': 'reference_box', 'parameters': {'padding': 4}},
                                      {'tool': 'fpocket', 'parameters': {}}],
            'docking': [{'tool': 'vina', 'parameters': {}}],
            'scoring': [{'tool': 'vina_score', 'parameters': {}}]}


def dataset(n=500):
    records = [{'compound_id': f'C{i % 250}', 'source_id': str(i), 'canonical_smiles': 'CCO',
                'cohort': f'assay{i // 250}', 'assay': f'A{i // 250}', 'endpoint': 'Ki',
                'eligible': True, 'p_value': float(i % 10), 'raw': {}, 'organism': 'Homo sapiens'} for i in range(n)]
    return {'records': records, 'identity': {'accession': 'TEST'}, 'source_total': n, 'loaded': n,
            'cap': n, 'complete': True, 'created_at': 'test'}


def setup_runtime(monkeypatch):
    from biocentric import engines
    monkeypatch.setattr(engines, 'engine_catalog', lambda: [{'id': k, 'blockers': [], 'hardware': {'gpu_ready': True}} for k in ['vina', 'gnina', 'unidock', 'gnina_affinity', 'biocentric_affinity']])
    monkeypatch.setattr(bench.pockets, 'executable', lambda: '/synthetic/fpocket')
    monkeypatch.setattr(bench, 'runtime_fingerprint', lambda *args: {})


@pytest.mark.parametrize('count', [500, 1000])
def test_review_keeps_full_eligible_search_and_distinct_assay_measurements(monkeypatch, count):
    setup_runtime(monkeypatch)
    data = dataset(count)
    data['records'][-1]['eligible'] = False
    data['records'][-1]['reasons'] = ['Synthetic source exclusion']
    review = bench.make_review(data, 'assay0', methods(), {'record': {'mutations': []}, 'component': 'LIG'})
    assert len(review['records']) == count - 1
    assert len({bench.measurement_key(r) for r in review['records']}) == count - 1
    assert len(review['candidates']) == 2
    assert len(review['record_exclusions']) == 1
    assert review['protocol']['total_budget_seconds'] is None
    bench.validate_review(review)
    review['records'][1]['measurement_id'] = review['records'][0]['measurement_id']
    review['sha256'] = bench.digest({k: v for k, v in review.items() if k != 'sha256'})
    with pytest.raises(ValueError, match='Duplicate measurement'):
        bench.validate_review(review)


def test_progress_waits_for_entire_set_and_control_validation(monkeypatch):
    setup_runtime(monkeypatch)
    review = bench.make_review(dataset(), None, methods(), {'record': {'mutations': []}, 'component': 'LIG'})
    pair = review['candidates'][0]
    rows = [{'combination': pair['id'], 'status': 'completed', **r} for r in review['records']]
    state = {'status': 'running', 'rows': rows}
    assert bench.combination_progress(review, state)['completed'] == 0
    rows.append({'combination': pair['id'], 'compound_id': 'crystal:LIG', 'status': 'completed'})
    progress = bench.combination_progress(review, state)
    assert (progress['completed'], progress['total']) == (1, 2)
    rows[0]['status'] = 'not_run'
    assert bench.combination_progress(review, state)['completed'] == 0
    rows[0]['status'] = 'failed'
    assert bench.combination_progress(review, state)['statuses'][pair['id']] == 'Failed'


def test_metrics_never_pool_assays_and_do_not_hide_constant_predictions():
    # Both within-assay rankings are perfect, despite contradictory between-assay offsets.
    rows = [{'cohort': 'a', 'p_value': i, 'score': i + 100} for i in range(8)]
    rows += [{'cohort': 'b', 'p_value': i + 20, 'score': i} for i in range(8)]
    stats = decisions.stratified_statistics(rows, 1)
    assert stats['rho'] == 1 and stats['tau'] == pytest.approx(1) and stats['ranked_cohorts'] == 2
    assert stats['ci_low'] == 1
    for r in rows[:8]:
        r['score'] = 0
    assert decisions.stratified_statistics(rows, 1)['rho'] is None


def test_pocket_cartesian_product_and_active_edges(monkeypatch):
    setup_runtime(monkeypatch)
    selected = methods()
    selected['docking'].append({'tool': 'gnina', 'parameters': {}})
    selected['scoring'].append({'tool': 'gnina_affinity', 'parameters': {}})
    candidates, excluded = bench.candidate_protocols(selected)
    assert len(candidates) == 8 and not excluded
    review = {'candidates': candidates, 'records': [], 'reference': {'component': 'LIG'}}
    data = graph_data(selected, {}, review, {'status': 'running', 'active_combination': candidates[-1]['id'], 'stage': 'pocket'})
    assert len(data['stages'][1]) == 2
    assert len(data['edges']) == 2 + 4 + 4 + 2
    assert [n['tool'] for n in data['nodes'] if n['current']] == ['fpocket']
    assert sum(e['running'] for e in data['edges']) == 4
    selected['pocket_identification'] = [{'tool': 'deeppocket', 'parameters': {}}]
    assert not bench.candidate_protocols(selected)[0]
    assert all('not connected' in r['reason'] for r in bench.candidate_protocols(selected)[1])


@pytest.mark.parametrize('pocket_failure', [False, True])
def test_worker_isolates_pocket_caches_and_accounts_repeated_compounds(tmp_path, monkeypatch, pocket_failure):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    setup_runtime(monkeypatch)
    mol = Chem.AddHs(Chem.MolFromSmiles('CCO')); AllChem.EmbedMolecule(mol, randomSeed=42)
    sdf = Chem.MolToMolBlock(mol) + '\n$$$$\n'
    records = dataset(6)['records']
    for i, r in enumerate(records):
        r.update(cohort='a' if i < 3 else 'b', measurement_id=str(i), compound_id=f'C{i % 3}')
    pairs, _ = bench.candidate_protocols(methods())
    review = {'schema': 2, 'records': records, 'candidates': pairs, 'runtime_fingerprint': {},
              'reference': {'component': 'LIG', 'receptor_pdb': 'synthetic', 'reference_sdf': sdf, 'canonical_smiles': 'CCO'},
              'protocol': {'total_budget_seconds': None}}
    review['sha256'] = bench.digest(review)
    (tmp_path / 'plan.json').write_text(json.dumps(review))
    (tmp_path / 'status.json').write_text(json.dumps({'status': 'queued', 'rows': [], 'events': [], 'completed': 0, 'total': 2}))
    prepared = tmp_path / 'prepared.sdf'; prepared.write_text(sdf)
    monkeypatch.setattr(bench.science, 'prepare_receptor', lambda *a, **k: {'output': str(tmp_path / 'receptor.pdbqt')})
    monkeypatch.setattr(bench.science, 'prepare_library', lambda *a, **k: {'compounds': [{}], 'prepared_sdf': str(prepared)})
    monkeypatch.setattr(bench.science, 'validate_pose', lambda *a, **k: {'passed': True, 'status': 'completed'})
    def pocket(method, *args):
        if pocket_failure and method['tool'] == 'fpocket':
            raise ValueError('Synthetic detection failure')
        return {'center': [0 if method['tool'] == 'fpocket' else 10, 0, 0], 'size': [20, 20, 20]}
    monkeypatch.setattr(bench.pockets, 'identify', pocket)
    calls = []
    monkeypatch.setattr(bench, '_dock', lambda pair, *args: calls.append(pair['pocket']['tool']) or [mol])
    monkeypatch.setattr(bench, '_score', lambda *a: [-6.])
    bench.worker(tmp_path)
    state = bench.read_status(tmp_path)
    assert state['completed'] == 2 and state['total'] == 2
    assert len(state['rows']) == 14  # Six measurements plus control for each chain.
    assert len({(r['combination'], bench.measurement_key(r)) for r in state['rows']}) == 14
    assert calls == (['reference_box'] if pocket_failure else ['reference_box', 'fpocket'])
    assert state['status'] == ('partial' if pocket_failure else 'completed')


def test_evaluation_ui_reviews_all_records_without_compound_slider(tmp_path, monkeypatch):
    from streamlit.testing.v1 import AppTest
    setup_runtime(monkeypatch)
    data = dataset(500)
    data['errors'] = []
    data['identity']['organism'] = 'Homo sapiens'
    structure = {'id': 'TEST', 'target_chains': ['A'], 'mutations': [], 'models': 1,
                 'compounds': [{'id': 'LIG', 'formula_weight': 200}]}
    monkeypatch.setattr(bench, 'reference_bundle', lambda *a: {'record': structure, 'component': 'LIG'})
    app = AppTest.from_string('''import streamlit as st
from biocentric.discovery_evaluation_ui import render_data
render_data(st, st.session_state['studio_eval_dataset']['identity'], st.session_state['test_evidence'], '.')
''')
    app.session_state['studio_eval_dataset'] = data
    app.session_state['test_evidence'] = {'structures': {'rows': [structure]}, 'binding': {'rows': [], 'total': 500}}
    app.run()
    assert not app.exception and not app.slider
    assert any('Full eligible set: 500 measurements' in item.value for item in app.caption)
    app.button(key='studio_eval_prepare').click().run()
    assert not app.exception
    review = app.session_state['studio_eval_review']
    assert len(review['records']) == 500 and len(review['candidates']) == 9
    assert not app.session_state.filtered_state.get('studio_eval_run')


def test_short_subset_is_reproducible_score_independent_and_keeps_assay_records(monkeypatch):
    setup_runtime(monkeypatch)
    data = dataset()
    review = bench.make_review(data, None, methods(), {'record': {'mutations': []}, 'component': 'LIG'}, max_compounds=30)
    ids = {r['compound_id'] for r in review['records']}
    assert len(ids) == 30 and len(review['records']) == 60
    assert all(sum(r['compound_id'] == key for r in review['records']) == 2 for key in ids)
    assert len(review['record_exclusions']) == 440
    assert review['ligand_subset']['max_compounds'] == 30
    changed = deepcopy(data)
    changed['records'].reverse()
    for row in changed['records']:
        row['p_value'] = 100 - row['p_value']
    again = bench.make_review(changed, None, methods(), {'record': {'mutations': []}, 'component': 'LIG'}, max_compounds=30)
    assert {r['compound_id'] for r in again['records']} == ids
    bench.validate_review(review)
    assert len(data['records']) == 500


@pytest.mark.parametrize('limit', [True, 0, 2, 3.5, '30'])
def test_ligand_subset_rejects_invalid_limits(limit):
    with pytest.raises(ValueError, match='integer'):
        bench.ligand_subset(dataset()['records'], limit)
