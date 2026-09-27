"""Synthetic measurements test invariants, never biological performance."""
from copy import deepcopy
from pathlib import Path
import json

import pytest

from biocentric import discovery_benchmark as bench


def measurement(id=1, **changes):
    return {'compound_id': f'CHEMBL{id}', 'endpoint': 'Ki', 'value': '10', 'units': 'nM', 'relation': '=',
            'organism': 'Homo sapiens', 'assay': 'CHEMBL100', 'document': 'CHEMBL200', 'target_id': 'CHEMBL300',
            'smiles': 'CCO', 'raw': {'activity_id': id}, **changes}


def test_concentration_scale_and_censoring_are_not_pooled():
    values = [measurement(1), measurement(2, value='0.01', units='uM', assay='CHEMBL101'),
              measurement(3, relation='>', assay='CHEMBL102'), measurement(4, endpoint='IC50')]
    rows = bench.normalize(values)
    assert [r['p_value'] for r in rows] == [8, 8, 8, 8]
    assert rows[2]['p_relation'] == '<' and not rows[2]['eligible']
    assert rows[0]['cohort'] != rows[3]['cohort']
    assert values[0]['value'] == '10'


@pytest.mark.parametrize('change', [{'value': 'nan'}, {'value': 'inf'}, {'value': '0'}, {'value': '-1'},
    {'value': True}, {'units': '%'}, {'relation': None}, {'organism': 'Mus musculus'}, {'organism': None},
    {'validity': 'Potential author error'}, {'smiles': 'invalid'}, {'raw': {'activity_id': 1, 'potential_duplicate': 1}}])
def test_bad_measurements_are_visible_exclusions(change):
    row = bench.normalize([measurement(**change)])[0]
    assert not row['eligible'] and row['reasons']


def test_duplicate_identities_and_conflicting_replicates():
    rows = bench.normalize([measurement(1), measurement(2)])
    assert sum(r['eligible'] for r in rows) == 1
    rows = bench.normalize([measurement(1), measurement(2, value='100')])
    assert not any(r['eligible'] for r in rows)


def test_partial_or_failed_protocol_gets_no_correlation():
    pair = {'id': 'pair', 'docking': {'tool': 'vina'}, 'scoring': {'tool': 'vina_score'}}
    review = {'candidates': [pair], 'records': [{}, {}, {}]}
    rows = [{'combination': 'pair', 'kind': 'measured', 'status': 'completed', 'score': s, 'p_value': p} for s, p in [(-6, 6), (-7, 7), (-8, 8)]]
    assert bench.summary_metrics(review, rows)[0]['Spearman ρ'] == 1
    assert bench.summary_metrics(review, rows[:-1])[0]['Spearman ρ'] is None
    rows[1]['status'] = 'failed'
    assert bench.summary_metrics(review, rows)[0]['Spearman ρ'] is None


def test_review_hash_and_target_cohort_guard():
    review = {'records': [{'cohort': 'a', 'p_value': 8}], 'candidates': []}
    review['sha256'] = bench.digest(review)
    bench.validate_review(review)
    review['records'][0]['p_value'] = 7
    with pytest.raises(ValueError, match='changed'):
        bench.validate_review(review)


def test_cancellation_accounts_for_every_unrun_pair(tmp_path, monkeypatch):
    review = {'records': [{'compound_id': 'measured', 'cohort': 'a', 'p_value': 8, 'endpoint': 'Ki'}],
              'candidates': [{'id': 'one', 'docking': {'tool': 'vina'}, 'scoring': {'tool': 'vina_score'}}],
              'reference': {'component': 'LIG', 'receptor_pdb': '', 'reference_sdf': ''},
              'protocol': {'total_budget_seconds': 1800}, 'runtime_fingerprint': {}}
    review['sha256'] = bench.digest(review)
    (tmp_path / 'plan.json').write_text(json.dumps(review))
    (tmp_path / 'status.json').write_text(json.dumps({'status': 'queued', 'rows': [], 'events': [], 'completed': 0, 'total': 2}))
    (tmp_path / 'cancel.request').touch()
    monkeypatch.setattr(bench, 'runtime_fingerprint', lambda *a: {})
    def cancelled(*args, **kwargs):
        assert kwargs['cancelled']()
        raise bench.science.ScienceError('cancelled', 'Cancelled')
    monkeypatch.setattr(bench.science, 'prepare_receptor', cancelled)
    bench.worker(tmp_path)
    state = bench.read_status(tmp_path)
    assert state['status'] == 'cancelled'
    assert len(state['rows']) == 2
    assert all(r['status'] == 'not_run' and r['score'] is None for r in state['rows'])


def test_scene_rejects_changed_or_outside_pose(tmp_path):
    from biocentric.discovery_evaluation_ui import _scene
    pose = tmp_path / 'pose.sdf'
    pose.write_text('changed')
    with pytest.raises(ValueError, match='integrity'):
        _scene(tmp_path, {}, {'pose_path': 'pose.sdf', 'pose_sha256': 'bad'})


def test_worker_accounts_live_results_and_reference_score(tmp_path, monkeypatch):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    # Synthetic test chemistry/labels. This test never dispatches native tools.
    mol = Chem.AddHs(Chem.MolFromSmiles('CCO'))
    AllChem.EmbedMolecule(mol, randomSeed=42)
    sdf = Chem.MolToMolBlock(mol) + '\n$$$$\n'
    pair = {'id': 'one', 'docking': {'tool': 'vina', 'parameters': {}}, 'scoring': {'tool': 'vina_score', 'parameters': {}}}
    review = {'records': [{'compound_id': 'synthetic-1', 'cohort': 'a', 'p_value': 8., 'endpoint': 'Ki', 'canonical_smiles': 'CCO'}],
              'candidates': [pair], 'reference': {'component': 'LIG', 'receptor_pdb': 'fixture', 'reference_sdf': sdf, 'canonical_smiles': 'CCO'},
              'protocol': {'total_budget_seconds': 1800}, 'runtime_fingerprint': {}}
    review['sha256'] = bench.digest(review)
    (tmp_path / 'plan.json').write_text(json.dumps(review))
    (tmp_path / 'status.json').write_text(json.dumps({'status': 'queued', 'rows': [], 'events': [], 'completed': 0, 'total': 2}))
    prepared = tmp_path / 'prepared.sdf'
    prepared.write_text(sdf)
    monkeypatch.setattr(bench, 'runtime_fingerprint', lambda *a: {})
    monkeypatch.setattr(bench.science, 'prepare_receptor', lambda *a, **k: {'output': str(tmp_path / 'receptor.pdbqt')})
    monkeypatch.setattr(bench.science, 'prepare_library', lambda *a, **k: {'compounds': [{}], 'prepared_sdf': str(prepared)})
    monkeypatch.setattr(bench.science, 'validate_pose', lambda *a, **k: {'status': 'completed', 'passed': True})
    calls = []
    monkeypatch.setattr(bench, '_dock', lambda *a, **k: calls.append('dock') or [mol])
    monkeypatch.setattr(bench, '_score', lambda *a, **k: [-6.])
    bench.worker(tmp_path)
    state = bench.read_status(tmp_path)
    assert state['status'] == 'completed' and state['completed'] == 2
    assert calls == ['dock']  # same immutable docking reused across same chemistry
    assert all(r['reference_score'] == -6 and r['rmsd_angstrom'] < .001 for r in state['rows'])
    assert all((tmp_path / r['pose_path']).is_file() for r in state['rows'])


def test_changed_selection_cannot_start_accepted_plan(tmp_path, monkeypatch):
    from types import SimpleNamespace
    from biocentric.discovery_evaluation_ui import _run
    from biocentric.discovery_studio import selection
    state = {'studio_identity': {'accession': 'P00519'}}
    review = {'selection': selection(state), 'identity': {'accession': 'P00519'}}
    state['studio_eval_review'] = review
    state['studio_tools_docking'] = ['gnina']
    calls = []
    monkeypatch.setattr(bench, 'start_run', lambda *a: calls.append(a))
    _run(SimpleNamespace(session_state=state), str(tmp_path))
    assert not calls and 'changed' in state['studio_eval_error']


def test_mutant_assay_cannot_be_compared_to_wild_type_or_different_mutant():
    row = measurement(raw={'activity_id': 1, 'assay_variant_mutation': 'T315I'})
    assert not bench.compatible_reference(row, {'mutations': []})
    assert not bench.compatible_reference(row, {'mutations': ['H396P']})
    assert not bench.compatible_reference(row, {'mutations': ['YES']})
    assert bench.compatible_reference(row, {'mutations': ['T315I']})
    assert bench.compatible_reference(measurement(), {'mutations': []})


def test_assay_species_uses_verified_metadata_not_unchecked_target_defaults(monkeypatch):
    from types import SimpleNamespace
    first = {'rows': [measurement(1, organism=None, raw={'activity_id': 1, 'target_tax_id': 9606, 'target_organism': 'Homo sapiens'})], 'total': 1, 'page': 0}
    identity = {'organism': 'Homo sapiens', 'taxon_id': 9606}
    assay = {'assay_chembl_id': 'CHEMBL100', 'target_chembl_id': 'CHEMBL300', 'confidence_score': 9, 'assay_organism': None}
    monkeypatch.setattr(bench.evidence, 'client', lambda *a: SimpleNamespace(_request=lambda *a, **k: {'data': {'assays': [assay]}}))
    data = bench.fetch_dataset(identity, '.', deepcopy(first))
    assert data['eligible'] == 1
    assert data['records'][0]['organism_basis'].startswith('Direct single-protein')
    assay['confidence_score'] = 8
    data = bench.fetch_dataset(identity, '.', deepcopy(first))
    assert data['eligible'] == 0


@pytest.mark.parametrize('middle', [True, False])
def test_reference_never_drops_internal_or_near_pocket_missing_atoms(monkeypatch, middle):
    import gemmi
    structure = gemmi.Structure()
    model = gemmi.Model('1')
    chain = gemmi.Chain('A')
    names = ['ALA', 'LYS', 'ALA'] if middle else ['ALA', 'ALA', 'LYS']
    for i, name in enumerate(names):
        res = gemmi.Residue(); res.name = name; res.seqid = gemmi.SeqId(i + 1, ' '); res.het_flag = 'A'
        for j, atom_name in enumerate(['N', 'CA', 'C', 'O', 'CB']):
            atom = gemmi.Atom(); atom.name = atom_name; atom.element = gemmi.Element(atom_name[0]); atom.pos = gemmi.Position(i * 3, j, 0)
            res.add_atom(atom)
        chain.add_residue(res)
    ligand = gemmi.Residue(); ligand.name = 'LIG'; ligand.seqid = gemmi.SeqId(20, ' '); ligand.het_flag = 'H'
    atom = gemmi.Atom(); atom.name = 'C1'; atom.element = gemmi.Element('C'); atom.pos = gemmi.Position(6, 1, 2)
    ligand.add_atom(atom); chain.add_residue(ligand); model.add_chain(chain); structure.add_model(model)
    text = structure.make_mmcif_document().as_string()
    monkeypatch.setattr(bench.evidence, 'coordinates', lambda *a: {'text': text})
    with pytest.raises(ValueError, match='Internal missing|within 12'):
        bench.reference_bundle({'id': 'TEST', 'target_chains': ['A'], 'compounds': [{'id': 'LIG'}]}, 'LIG', 'A', '.')


def test_sh2_reference_does_not_accept_kinase_activity_assays():
    reference = {'title': 'CRYSTAL STRUCTURE OF SH2 IN COMPLEX', 'mutations': []}
    assert not bench.compatible_reference({'assay_description': 'Inhibition of human src kinase.'}, reference)
    assert bench.compatible_reference({'assay_description': 'Binding to SRC SH2 domain'}, reference)
