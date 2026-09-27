"""Synthetic orchestration checks; these do not measure predictive accuracy."""
import json
from pathlib import Path

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from biocentric import biocentric_model as model, discovery_affinity as affinity
from biocentric import discovery_benchmark as bench, engines, science


@pytest.fixture
def poses():
    mol = Chem.AddHs(Chem.MolFromSmiles('CCO'))
    AllChem.EmbedMolecule(mol, randomSeed=42)
    mol.SetProp('_Name', 'original')
    other = Chem.Mol(mol)
    other.GetConformer().SetAtomPosition(0, (2, 3, 4))
    return [mol, other]


@pytest.fixture
def native(monkeypatch):
    calls = []
    def run(engine, inputs, parameters, folder, **kwargs):
        assert engine == 'biocentric_affinity'
        assert parameters['input_mode'] == 'pose_records'
        mols = list(Chem.SDMolSupplier(inputs['ligands'], removeHs=False))
        calls.append((parameters, mols))
        provenance = {'weights_sha256': model.SUPPORTED_WEIGHTS_SHA256, 'code_sha256': 'model-code', 'models': ['ensemble']}
        rows = []
        for index, mol in enumerate(mols, 1):
            value = float(index)
            rows.append({'record_index': index, 'compound_id': mol.GetProp('_Name'), 'protein_id': inputs['protein_id'],
                         'smiles': science._canonical(mol), 'status': 'completed', 'score': value,
                         'per_model': {f'{name}@{epoch}': value for name, epoch in model.MODELS},
                         'checkpoint_sha256': provenance['weights_sha256'], 'code_sha256': provenance['code_sha256'],
                         'models': provenance['models'], 'device': parameters['device'], 'ensemble_spread': 0.})
        result = {'status': 'completed', 'results': rows, 'provenance': provenance}
        folder.mkdir(parents=True)
        science._write_json(folder / 'result.json', result)
        return result
    monkeypatch.setattr(engines, 'run_engine', run)
    return calls, run


@pytest.mark.parametrize('failure', [False, True])
def test_batch_binds_all_poses_and_never_selects_only_survivors(tmp_path, monkeypatch, poses, native, failure):
    calls, original = native
    if failure:
        def failed(*args, **kwargs):
            result = original(*args, **kwargs)
            result['status'] = 'partial'
            result['results'][1].update(status='failed', score=None, failure_reason='unsupported chemistry', per_model={})
            return result
        monkeypatch.setattr(engines, 'run_engine', failed)
    receptor = tmp_path / 'r.pdb'
    receptor.write_text('test receptor')
    before = [Chem.MolToMolBlock(mol) for mol in poses]
    result = affinity.score_groups({'measurement-a': poses, 'measurement-b': [poses[0]]}, receptor,
                                  tmp_path / 'batch', {}, lambda: False, None)
    assert len(calls) == 1 and calls[0][0]['max_pairs'] == 3
    assert result['measurement-b']['scores'] == [3.]
    assert result['measurement-a']['scores'] == ([] if failure else [1., 2.])
    assert bool(result['measurement-a']['errors']) == failure
    assert before == [Chem.MolToMolBlock(mol) for mol in poses]
    manifest = json.loads((tmp_path / 'batch/pose-manifest.json').read_text())
    assert [r['pose_index'] for r in manifest['records']] == [0, 1, 0]
    assert len({r['compound_id'] for r in manifest['records']}) == 3


def test_biocentric_catalog_checks_model_and_explicit_device(monkeypatch):
    catalog = [{'id': tool, 'blockers': [], 'hardware': {'gpu_ready': True}} for tool in ('vina', 'biocentric_affinity')]
    monkeypatch.setattr(engines, 'engine_catalog', lambda: catalog)
    selected = {'pocket_identification': [bench.pockets.REFERENCE], 'docking': [{'tool': 'vina', 'parameters': {}}],
                'scoring': [{'tool': 'biocentric_affinity', 'parameters': {}}]}
    assert len(bench.candidate_protocols(selected)[0]) == 1
    catalog[-1]['hardware']['gpu_ready'] = False
    assert 'CUDA' in bench.candidate_protocols(selected)[1][0]['reason']
    selected['scoring'][0]['parameters']['device'] = 'cpu'
    assert len(bench.candidate_protocols(selected)[0]) == 1
    catalog[-1]['blockers'] = ['missing checkpoint']
    assert 'checkpoint' in bench.candidate_protocols(selected)[1][0]['reason']


def test_runtime_binds_model_files_without_requiring_gnina(tmp_path, monkeypatch):
    binary = tmp_path / 'vina'
    binary.write_text('vina runtime')
    checkpoint = tmp_path / 'weights'
    checkpoint.write_text('original')
    def installation():
        return {'blockers': [], 'python': str(binary),
                'files': {'weights': {'path': str(checkpoint), 'sha256': science.sha256(checkpoint)}}}
    monkeypatch.setattr(model, 'installation', installation)
    names = []
    monkeypatch.setattr(engines, '_executable', lambda name: names.append(name) or str(binary))
    pair = {'docking': {'tool': 'vina'}, 'scoring': {'tool': 'biocentric_affinity'}}
    initial = bench.runtime_fingerprint([pair])
    assert names == ['vina']
    checkpoint.write_text('changed')
    assert bench.runtime_fingerprint([pair]) != initial


@pytest.mark.parametrize('cancel', [False, True])
def test_worker_batches_pose_selection_reference_and_full_accounting(tmp_path, monkeypatch, poses, native, cancel):
    pair = {'id': 'bio', 'pocket': bench.pockets.REFERENCE,
            'docking': {'tool': 'vina', 'parameters': {}},
            'scoring': {'tool': 'biocentric_affinity', 'parameters': {'batch_size': 16}}}
    sdf = Chem.MolToMolBlock(poses[0]) + '\n$$$$\n'
    review = {'schema': 2, 'records': [{'compound_id': 'same', 'measurement_id': str(i), 'cohort': 'assay',
               'p_value': float(i), 'endpoint': 'Ki', 'canonical_smiles': 'CCO'} for i in range(3)],
              'candidates': [pair], 'reference': {'component': 'LIG', 'receptor_pdb': 'fixture',
              'reference_sdf': sdf, 'canonical_smiles': 'CCO'},
              'protocol': {'total_budget_seconds': None}, 'runtime_fingerprint': {}}
    review['sha256'] = bench.digest(review)
    science._write_json(tmp_path / 'plan.json', review)
    science._write_json(tmp_path / 'status.json', {'status': 'queued', 'rows': [], 'events': [], 'completed': 0, 'total': 1})
    prepared = tmp_path / 'prepared.sdf'
    prepared.write_text(sdf)
    monkeypatch.setattr(bench, 'runtime_fingerprint', lambda *a: {})
    monkeypatch.setattr(science, 'prepare_receptor', lambda *a, **k: {'output': str(tmp_path / 'r.pdbqt')})
    monkeypatch.setattr(science, 'prepare_library', lambda *a, **k: {'compounds': [{}], 'prepared_sdf': str(prepared)})
    monkeypatch.setattr(science, 'validate_pose', lambda *a, **k: {'status': 'completed', 'passed': True})
    docks = []
    monkeypatch.setattr(bench, '_dock', lambda *a, **k: docks.append(1) or poses)
    if cancel:
        original = native[1]
        def stop(*args, **kwargs):
            result = original(*args, **kwargs)
            (tmp_path / 'cancel.request').touch()
            return result
        monkeypatch.setattr(engines, 'run_engine', stop)
    bench.worker(tmp_path)
    state = bench.read_status(tmp_path)
    assert len(state['rows']) == 4 and len({bench.measurement_key(r) for r in state['rows']}) == 4
    assert len(docks) == len(native[0]) == 1
    if cancel:
        assert state['status'] == 'cancelled' and state['completed'] == 0
        assert all(r['status'] == 'not_run' and r['score'] is None for r in state['rows'])
    else:
        assert state['status'] == 'completed' and state['completed'] == state['total'] == 1
        assert [r['score'] for r in state['rows']] == [2., 4., 6., 8.]
        assert all(r['top_pose_index'] == 2 and r['reference_score'] == 9 for r in state['rows'])
        assert all(r['model_prediction']['experimental_measurement'] is False for r in state['rows'])
        assert all((tmp_path / r['pose_path']).is_file() for r in state['rows'])
