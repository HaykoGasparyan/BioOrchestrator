"""Reuse must preserve exact scientific inputs; synthetic fixtures only."""
from copy import deepcopy
from pathlib import Path
import json

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from biocentric import discovery_benchmark as bench, discovery_studio as studio, science
from biocentric.evaluation_archive import Archive, StageCache, identify_pocket


def molecule():
    mol = Chem.AddHs(Chem.MolFromSmiles('CCO'))
    AllChem.EmbedMolecule(mol, randomSeed=42)
    return mol


def archive(tmp_path, runtime=None):
    return Archive(tmp_path / 'archive', runtime or {}, lambda: False, lambda e: None)


def test_archive_reuses_verifies_and_never_publishes_failure(tmp_path):
    cache = archive(tmp_path)
    calls = []
    def produce(out):
        out.mkdir(parents=True)
        (out / 'result').write_text('real output')
        calls.append(out)
    identity = {'tool': 'fpocket', 'receptor': 'sha', 'rank': 1}
    a = cache.materialize('pocket', identity, tmp_path / 'first', produce)
    b = cache.materialize('pocket', identity, tmp_path / 'second', produce)
    assert not a['hit'] and b['hit'] and len(calls) == 1
    assert (tmp_path / 'second/result').read_text() == 'real output'
    (cache.root / 'pocket' / a['key'] / 'artifacts/result').write_text('corrupt')
    assert not cache.materialize('pocket', identity, tmp_path / 'repair', produce)['hit']
    assert not cache.materialize('pocket', {**identity, 'rank': 2}, tmp_path / 'different', produce)['hit']
    def fail(out):
        out.mkdir(); (out / 'partial').touch()
        raise ValueError('native failure')
    with pytest.raises(ValueError):
        cache.materialize('pocket', {**identity, 'rank': 3}, tmp_path / 'failure', fail)
    assert len(list((cache.root / 'pocket').glob('*/manifest.json'))) == 2


def test_scoring_model_change_does_not_invalidate_pocket_but_pocket_code_does(tmp_path):
    runtime = {'files': {'/source/pocket_identification.py': 'v1', '/model/weights.pt': 'old'}}
    def produce(out):
        out.mkdir(); (out / 'pocket.json').write_text('{}')
    first = archive(tmp_path, runtime).materialize('pocket', {}, tmp_path / 'one', produce)
    runtime['files']['/model/weights.pt'] = 'new'
    assert archive(tmp_path, runtime).materialize('pocket', {}, tmp_path / 'two', produce)['hit']
    runtime['files']['/source/pocket_identification.py'] = 'v2'
    last = archive(tmp_path, runtime).materialize('pocket', {}, tmp_path / 'three', produce)
    assert not last['hit'] and first['key'] != last['key']


def test_pocket_and_docking_reuse_across_scores_runs_and_assay_records(tmp_path, monkeypatch):
    mol = molecule()
    receptor, pdbqt, reference = [tmp_path / name for name in ['protein.pdb', 'protein.pdbqt', 'reference.sdf']]
    receptor.write_text('protein'); pdbqt.write_text('prepared protein')
    bench._write_sdf(reference, mol)
    original = tmp_path / 'original.sdf'; bench._write_sdf(original, mol)
    calls = {'pocket': 0, 'prepare': 0, 'dock': 0}
    box = {'center': [0, 0, 0], 'size': [20, 20, 20]}
    def pocket(method, receptor, reference, out, cancelled):
        calls['pocket'] += 1; out.mkdir(parents=True)
        science._write_json(out / 'pocket.json', {'box': box})
        return box
    def prepare(*a, **k):
        calls['prepare'] += 1
        return {'compounds': [{}], 'prepared_sdf': str(original)}
    def dock(*a, **k):
        calls['dock'] += 1
        return [mol]
    monkeypatch.setattr(bench.pockets, 'identify', pocket)
    monkeypatch.setattr(science, 'prepare_library', prepare)
    monkeypatch.setattr(bench, '_dock', dock)
    pair = {'pocket': bench.pockets.REFERENCE, 'docking': {'tool': 'vina', 'parameters': {'seed': 42}},
            'scoring': {'tool': 'vina_score'}}
    record = {'canonical_smiles': 'CCO'}
    for run in range(2):
        ar = archive(tmp_path); cache = StageCache(ar)
        assert identify_pocket(pair['pocket'], receptor, reference, tmp_path / f'pocket-{run}', lambda: False, ar) == box
        for score in ['vina_score', 'gnina_affinity', 'biocentric_affinity']:
            pair['scoring']['tool'] = score
            poses = bench._prepare_docked_poses(pair, record, tmp_path / f'run-{run}-{score}', receptor, pdbqt, box, cache, lambda: False, lambda e: None)
            assert science._symmetry_rmsd(poses[0], mol)['rmsd_angstrom'] < .001
    assert calls == {'pocket': 1, 'prepare': 1, 'dock': 1}
    changed = deepcopy(pair); changed['docking']['parameters']['seed'] = 43
    bench._prepare_docked_poses(changed, record, tmp_path / 'new-seed', receptor, pdbqt, box, StageCache(archive(tmp_path)), lambda: False, lambda e: None)
    assert calls['dock'] == 2 and calls['prepare'] == 1
    receptor.write_text('different coordinates')
    bench._prepare_docked_poses(pair, record, tmp_path / 'new-receptor', receptor, pdbqt, box, StageCache(archive(tmp_path)), lambda: False, lambda e: None)
    assert calls['dock'] == 3


def test_gpu_defaults_and_explicit_cpu_are_preserved_and_used(tmp_path, monkeypatch):
    selected = studio.selection({'studio_parameters': {'docking:gnina': {'seed': 7}}})
    gnina = next(t for t in selected['docking'] if t['tool'] == 'gnina')
    assert gnina['parameters']['device'] == 'gpu' and gnina['parameters']['seed'] == 7
    assert next(t for t in selected['scoring'] if t['tool'] == 'gnina_cnn')['parameters']['device'] == 'gpu'
    from biocentric import engines, engine_workflow
    seen = []
    monkeypatch.setattr(engines, 'run_engine', lambda engine, inputs, params, *a, **k: seen.append(params) or {'status': 'completed'})
    monkeypatch.setattr(engine_workflow, '_native_molecules', lambda *a: [molecule()])
    for device in ['gpu', 'cpu']:
        gnina['parameters']['device'] = device
        bench._dock({'docking': gnina}, molecule(), tmp_path / device, tmp_path / 'r', tmp_path / 'q',
                    {'center': [0]*3, 'size': [20]*3}, lambda: False, lambda e: None)
    assert [p['device'] for p in seen] == ['gpu', 'cpu']


def test_selected_pool_executable_pool_and_budget_are_distinct(monkeypatch):
    from biocentric import engines, bayesian_search
    monkeypatch.setattr(engines, 'engine_catalog', lambda: [{'id': k, 'blockers': [], 'hardware': {'gpu_ready': True}} for k in ['vina', 'gnina', 'unidock', 'gnina_affinity', 'biocentric_affinity']])
    monkeypatch.setattr(bench.pockets, 'executable', lambda: '/synthetic/fpocket')
    selected = studio.selection({'studio_tools_pocket_identification': ['fpocket', 'deeppocket']})
    pairs, excluded = bench.candidate_protocols(selected)
    review = {'candidates': pairs, 'excluded': excluded, 'records': [], 'reference': {'component': 'LIG'}, 'search': bayesian_search.configuration(pairs)}
    progress = bench.combination_progress(review, {})
    assert (progress['selected_total'], progress['pool_total'], progress['total']) == (18, 9, 6)
    assert bench.combination_caption(review) == '18 combinations · 9 executable · Bayesian budget 6 · 9 unavailable'
    selected['pocket_identification'][1] = bench.pockets.REFERENCE
    pairs, excluded = bench.candidate_protocols(selected)
    assert len(pairs) == 18 and not excluded


def test_runtime_hash_detects_same_size_rewrite_with_restored_mtime(tmp_path, monkeypatch):
    import os
    import time
    runtime = tmp_path / 'weights'
    runtime.write_text('aaa')
    time.sleep(1.05)  # Only stable installed files are hash-cached.
    calls = []
    original = science.sha256
    monkeypatch.setattr(science, 'sha256', lambda p: calls.append(p) or original(p))
    first = science.runtime_sha256(runtime)
    assert science.runtime_sha256(runtime) == first and len(calls) == 1
    old = runtime.stat()
    runtime.write_text('bbb')
    os.utime(runtime, ns=(old.st_atime_ns, old.st_mtime_ns))
    assert science.runtime_sha256(runtime) != first and len(calls) == 2
