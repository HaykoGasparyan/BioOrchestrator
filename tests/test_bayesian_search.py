"""Synthetic Bayesian scheduling evidence, not scientific validation."""
import json
from itertools import product
from pathlib import Path

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from biocentric import bayesian_search as search, discovery_benchmark as bench
from biocentric import discovery_search_runner as runner, pipeline_decisions as decisions, science


def candidates():
    return [{'id': f'{d}-{s}', 'pocket': bench.pockets.REFERENCE,
             'docking': {'tool': d, 'parameters': {}}, 'scoring': {'tool': s, 'parameters': {}}}
            for d, s in product(['vina', 'gnina', 'unidock'], ['vina_score', 'gnina_cnn', 'biocentric_affinity'])]


def test_gp_posterior_matches_direct_linear_algebra_and_ei_is_finite():
    pairs = candidates()
    config = search.configuration(pairs)
    observations = [{'id': pairs[0]['id'], 'objective': .8}, {'id': pairs[4]['id'], 'objective': -.5}]
    predictions = search.posterior(pairs, observations, config)
    selected = [pairs[0], pairs[4]]
    kernel = lambda a, b: np.exp(-search.distance(a, b) / .7)
    covariance = np.array([[kernel(a, b) for b in selected] for a in selected]) + .0025 * np.eye(2)
    expected_mean = np.array([kernel(pairs[1], x) for x in selected]) @ np.linalg.solve(covariance, [.8, -.5])
    assert predictions[1]['mean'] == pytest.approx(expected_mean)
    assert all(np.isfinite(r['expected_improvement']) and r['expected_improvement'] >= 0 and r['std'] > 0 for r in predictions)
    assert predictions[0]['std'] < predictions[8]['std']
    proposed, trace = search.propose(pairs, observations, config, 2)
    assert len({p['id'] for p in proposed}) == 2
    assert not {p['id'] for p in proposed} & {r['id'] for r in observations}
    assert trace[0]['id'] == max((r for r in predictions if r['id'] not in {o['id'] for o in observations}), key=lambda r:(r['expected_improvement'],r['std'],r['id']))['id']
    assert observations == [{'id': pairs[0]['id'], 'objective': .8}, {'id': pairs[4]['id'], 'objective': -.5}]


def test_initial_design_is_deterministic_diverse_and_budgeted():
    pairs = candidates()
    initial = search.initial_candidates(pairs, 3)
    assert [r['id'] for r in initial] == [r['id'] for r in search.initial_candidates(list(reversed(pairs)), 3)]
    assert len({r['docking']['tool'] for r in initial}) == 3
    assert len({r['scoring']['tool'] for r in initial}) == 3
    assert search.configuration(pairs)['max_trials'] == 6 < len(pairs)
    assert search.configuration(pairs[:3])['max_trials'] == 3


@pytest.mark.parametrize('fail_ligand', [False, True])
@pytest.mark.parametrize('cancel', [False, True])
def test_scheduler_is_ligand_first_adaptive_accounted_and_locks_only_real_evidence(tmp_path, monkeypatch, cancel, fail_ligand):
    mol = Chem.AddHs(Chem.MolFromSmiles('CCO'))
    AllChem.EmbedMolecule(mol, randomSeed=42)
    sdf = Chem.MolToMolBlock(mol) + '\n$$$$\n'
    pairs = candidates()
    records = [{'compound_id': f'C{i}', 'measurement_id': f'M{i}', 'canonical_smiles': 'CCO', 'cohort': 'assay', 'endpoint': 'Ki', 'p_value': i+1.} for i in range(4)]
    review = {'schema': 2, 'records': records, 'candidates': pairs, 'search': search.configuration(pairs),
              'reference': {'component': 'LIG', 'canonical_smiles': 'CCO', 'reference_sdf': sdf, 'receptor_pdb': 'test fixture'},
              'protocol': {'total_budget_seconds': None, 'preparation': 'fixture', 'pocket': 'fixture'},
              'identity': {'accession': 'TEST'}, 'runtime_fingerprint': {}}
    review['sha256'] = bench.digest(review)
    science._write_json(tmp_path / 'plan.json', review)
    science._write_json(tmp_path / 'status.json', {'status': 'queued', 'rows': [], 'events': [], 'plan_sha256': review['sha256']})
    monkeypatch.setattr(bench, 'runtime_fingerprint', lambda *a: {})
    monkeypatch.setattr(science, 'prepare_receptor', lambda *a, **k: {'output': str(tmp_path / 'r.pdbqt')})
    order = []
    def evaluate(pair, record, index, folder, *args):
        order.append((record['compound_id'], pair['id']))
        if cancel and len(order) == 2:
            (folder / 'cancel.request').touch()
            science._cancel(lambda: True)
        file = folder / f"{pair['id']}-{index}.sdf"
        file.write_text(sdf)
        value = record['p_value'] if record['kind'] == 'measured' else 5.
        value *= 1 if bench.SCORES[pair['scoring']['tool']][2] else -1
        failed = fail_ligand and record['compound_id'] == 'C0' and pair['scoring']['tool'] == 'biocentric_affinity'
        return {**record, 'combination': pair['id'], 'docking': pair['docking']['tool'], 'scoring': pair['scoring']['tool'],
                'pocket': 'reference_box', 'status': 'failed' if failed else 'completed', 'score': None if failed else value, 'physical_pass': True,
                'rmsd_angstrom': 0., 'reference_comparable': True, 'pose_path': file.name,
                'pose_sha256': science.sha256(file), 'runtime_seconds': .01}
    monkeypatch.setattr(runner, 'evaluate', evaluate)
    runner.worker(tmp_path, review)
    state = bench.read_status(tmp_path)
    assert len(state['rows']) == 9 * 5
    assert len({(r['combination'],bench.measurement_key(r)) for r in state['rows']}) == 45
    if cancel:
        assert state['status'] == 'cancelled' and state['completed'] == 0
        assert not (tmp_path / 'decision.json').exists()
        assert sum(r['status'] == 'not_run' for r in state['rows']) == 44
        return
    assert state['status'] == ('partial' if fail_ligand else 'completed') and state['completed'] == state['total'] == 6
    if fail_ligand:
        observations = state['search']['observations']
        assert all(o['comparison_measurement_ids'] == ['M1', 'M2', 'M3'] for o in observations)
        assert all(o['feasible'] for o in observations)
        assert all(m['Shared set'] == 3 for m in state['metrics'])
        assert any(o['skipped'] == 1 for o in observations)
    assert len(state['search']['rounds']) == 3 and state['search']['rounds'][1]['acquisition']
    assert len(state['search']['observations']) == 6
    assert sum(r['status'] == 'not_selected' for r in state['rows']) == 15
    initial_ids = state['search']['rounds'][0]['ids']
    assert order[:6] == [('crystal:LIG', key) for key in initial_ids] + [('C0',key) for key in initial_ids]
    assert len({r['id'] for r in state['search']['observations']}) == 6
    lock = decisions.read_lock(tmp_path)
    assert lock['selection_origin'] == 'bayesian_search' and lock['rank_on_this_set'] == 1
    assert state['decision_status'] == 'locked'
    skipped = [r for r in decisions.compare(review,state) if r['id'] not in state['search']['selected_ids']]
    assert all(not r['lockable'] and r['color'] == 'gray' and r['rho'] is None for r in skipped)


def test_final_result_tabs_render_fixed_metrics_and_best_plus_five(tmp_path, monkeypatch):
    # Reuse the synthetic scheduler fixture to exercise the real persisted UI.
    test_scheduler_is_ligand_first_adaptive_accounted_and_locks_only_real_evidence(tmp_path, monkeypatch, False, False)
    from streamlit.testing.v1 import AppTest
    app = AppTest.from_string('''import streamlit as st
from biocentric.pipeline_results_ui import render_results, render_final_chain
render_results(st)
render_final_chain(st)
''')
    app.session_state['studio_eval_run'] = str(tmp_path)
    app.run()
    assert not app.exception
    assert len(app.dataframe[0].value) == 6
    assert app.dataframe[0].value.iloc[0]['Decision'] == 'Locked'
    assert len(app.metric) == 4
    assert any('Fixed' in s.value for s in app.success)


def test_first_round_covers_all_docking_engines_and_old_plans_still_validate():
    # IDs/configurations from the reported case: the old max-min tie left Uni-Dock out.
    pairs = [{'id': key, 'pocket': bench.pockets.REFERENCE, 'docking': {'tool': dock}, 'scoring': {'tool': score}}
             for key, dock, score in [('292c8f823d7c1228', 'vina', 'vina_score'), ('02c653efee1173b6', 'vina', 'biocentric_affinity'),
                 ('8ee0699dec160373', 'gnina', 'vina_score'), ('bcb0c48e8805ed47', 'gnina', 'biocentric_affinity'),
                 ('3d4d28d4f575a375', 'unidock', 'vina_score'), ('bf2c0f4e739227cc', 'unidock', 'biocentric_affinity')]]
    assert {p['docking']['tool'] for p in search.initial_candidates(pairs, 3)} == {'vina', 'gnina', 'unidock'}
    assert len({p['docking']['tool'] for p in search.initial_candidates(pairs, 3, method='categorical_gp_expected_improvement_v1')}) == 2
    review = {'schema': 2, 'records': [], 'candidates': pairs,
              'search': search.configuration(pairs, 'categorical_gp_expected_improvement_v1')}
    review['sha256'] = bench.digest(review)
    bench.validate_review(review)
