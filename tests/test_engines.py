import json
from pathlib import Path

import pytest
from pydantic import ValidationError

from biocentric import engines, science

ROOT = Path(__file__).resolve().parents[1]


def pdbqt(path, ligand=False):
    path.write_text('ATOM      1  C   LIG A   1       1.000   2.000   3.000  1.00  0.00     0.000 C\n' + ('ROOT\nENDROOT\nTORSDOF 0\n' if ligand else ''))
    return str(path)


def dock_input(tmp_path):
    return {'receptor': pdbqt(tmp_path / 'receptor.pdbqt'), 'ligand': pdbqt(tmp_path / 'ligand.pdbqt', True), 'preparation_reviewed': True}


def dock_parameters():
    return {'box': {'center': [1, 2, 3], 'size': [12, 12, 12]}, 'cpu': 1, 'seed': 42, 'exhaustiveness': 1, 'num_modes': 1}


def test_catalog_reports_real_local_blockers_without_running(monkeypatch):
    monkeypatch.setattr(engines, '_executable', lambda name: '/configured/vina' if name == 'vina' else None)
    monkeypatch.setattr(engines.platform, 'system', lambda: 'Darwin')
    monkeypatch.delenv('BIOCENTRIC_PLANET_DIR', raising=False)
    rows = {row['id']: row for row in engines.engine_catalog()}
    assert set(rows) == set(engines.ENGINE_IDS)
    assert rows['vina']['available']
    assert not rows['gnina']['available']
    assert not rows['unidock']['available']
    assert any('Linux NVIDIA' in reason for reason in rows['unidock']['blockers'])
    assert not rows['planet']['available']
    assert rows['gnina_affinity']['output_metrics'][-1]['unit'] == 'pK'
    assert rows['mmpbsa']['output_metrics'][0]['unit'] == 'kcal/mol'


def test_plan_typed_argv_no_execution_or_creation(tmp_path, monkeypatch):
    monkeypatch.setattr(engines, '_executable', lambda name: None)
    inp = dock_input(tmp_path)
    plan = engines.plan_engine_run('vina', inp, dock_parameters(), tmp_path / 'out')
    assert plan['status'] == 'blocked'
    assert plan['commands'][0][:3] == ['vina', '--receptor', inp['receptor']]
    assert not (tmp_path / 'out').exists()
    assert plan['validated_inputs']['ligand']['sha256'] == science.sha256(Path(inp['ligand']))
    for changed in ({'cpu': 99}, {'seed': 0}, {'shell': 'touch BAD'}, {'box': {'center': [float('nan'), 0, 0], 'size': [10, 10, 10]}}):
        with pytest.raises(ValidationError):
            engines.plan_engine_run('vina', inp, {**dock_parameters(), **changed}, tmp_path / 'out')
    with pytest.raises(ValidationError):
        engines.plan_engine_run('vina', {**inp, 'preparation_reviewed': False}, dock_parameters(), tmp_path / 'out')


def test_missing_tools_never_become_success(tmp_path, monkeypatch):
    monkeypatch.setattr(engines, '_executable', lambda name: None)
    monkeypatch.setattr(science, 'run_process', lambda *a, **k: pytest.fail('Must not execute unavailable engine'))
    result = engines.run_engine('vina', dock_input(tmp_path), dock_parameters(), tmp_path / 'out')
    assert result['status'] == 'blocked'
    assert result['executed'] is False
    assert result['results'] == []


def test_vina_parser_names_scores_and_rejects_no_scores(tmp_path):
    path = tmp_path / 'poses.pdbqt'
    path.write_text('REMARK VINA RESULT: -7.2 0.0 0.0\nREMARK VINA RESULT: -6.1 1.0 2.0\n')
    rows = engines.parse_engine_output('vina', path)
    assert [row['metrics'][0]['value'] for row in rows] == [-7.2, -6.1]
    assert rows[0]['metrics'][0]['unit'] == 'kcal/mol'
    assert rows[0]['metrics'][0]['direction'] == 'lower_is_better'
    path.write_text('successful docking without any native score\n')
    with pytest.raises(ValueError, match='No complete'):
        engines.parse_engine_output('vina', path)


def sdf(path, properties=None, is3d=True):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    molecule = Chem.AddHs(Chem.MolFromSmiles('CCO'))
    AllChem.EmbedMolecule(molecule, randomSeed=42)
    molecule.GetConformer().Set3D(is3d)
    molecule.SetProp('_Name', 'ethanol')
    for name, value in (properties or {}).items():
        molecule.SetProp(name, str(value))
    with Chem.SDWriter(str(path)) as writer:
        writer.write(molecule)
    return str(path)


def test_gnina_score_only_keeps_affinity_and_pose_metrics_separate(tmp_path, monkeypatch):
    monkeypatch.setattr(engines, '_executable', lambda name: None)
    inp = {'receptor': pdbqt(tmp_path / 'receptor.pdb'), 'ligand': sdf(tmp_path / 'pose.sdf'), 'preparation_reviewed': True}
    plan = engines.plan_engine_run('gnina_affinity', inp, {}, tmp_path / 'out')
    assert '--score_only' in plan['commands'][0]
    assert '--no_gpu' in plan['commands'][0]
    assert '--exhaustiveness' not in plan['commands'][0]
    result = tmp_path / 'gnina.sdf'
    sdf(result, {'minimizedAffinity': -7, 'CNNscore': .8, 'CNNaffinity': 6.2})
    scores = {m['name']: m for m in engines.parse_engine_output('gnina_affinity', result)[0]['metrics']}
    assert scores['CNNaffinity']['unit'] == 'pK'
    assert scores['CNNscore']['unit'] == 'dimensionless'
    assert scores['minimizedAffinity']['unit'] == 'kcal/mol'
    sdf(result, {'minimizedAffinity': -7, 'CNNscore': .8})
    with pytest.raises(ValueError, match='missing CNNaffinity'):
        engines.parse_engine_output('gnina_affinity', result)


def test_unidock_requires_unique_output_names_and_gpu(tmp_path, monkeypatch):
    monkeypatch.setattr(engines, '_executable', lambda name: None)
    monkeypatch.setattr(engines.platform, 'system', lambda: 'Darwin')
    inp = dock_input(tmp_path)
    plan = engines.plan_engine_run('unidock', {'receptor': inp['receptor'], 'ligands': [inp['ligand']], 'preparation_reviewed': True}, {'box': dock_parameters()['box']}, tmp_path / 'out')
    assert '--gpu_batch' in plan['commands'][0]
    assert plan['status'] == 'blocked'
    with pytest.raises(ValueError, match='basenames'):
        engines.plan_engine_run('unidock', {'receptor': inp['receptor'], 'ligands': [inp['ligand']] * 2, 'preparation_reviewed': True}, {'box': dock_parameters()['box']}, tmp_path / 'out')


def test_autodock4_generated_dpf_binds_maps_and_finite_search(tmp_path, monkeypatch):
    monkeypatch.setattr(engines, '_executable', lambda name: None)
    inp = dock_input(tmp_path)
    header = 'GRID_PARAMETER_FILE receptor.gpf\nGRID_DATA_FILE receptor.maps.fld\nMACROMOLECULE receptor.pdbqt\nSPACING 0.375\nNELEMENTS 40 40 40\nCENTER 1 2 3\n0.0\n'
    for file in ('receptor.C.map', 'receptor.e.map', 'receptor.d.map'):
        (tmp_path / file).write_text(header)
    (tmp_path / 'receptor.maps.fld').write_text('# AVS field\n')
    inp.update(grid_field=str(tmp_path / 'receptor.maps.fld'), atom_maps={'C': str(tmp_path / 'receptor.C.map')}, electrostatic_map=str(tmp_path / 'receptor.e.map'), desolvation_map=str(tmp_path / 'receptor.d.map'), maps_match_receptor_reviewed=True)
    plan = engines.plan_engine_run('autodock4', inp, {'seed': 42, 'ga_runs': 2}, tmp_path / 'out')
    assert 'ga_run 2\n' in plan['generated_files']['protocol.dpf']
    assert 'seed 42 43\n' in plan['generated_files']['protocol.dpf']
    assert '-p' in plan['commands'][0]
    (tmp_path / 'receptor.d.map').write_text(header.replace('CENTER 1 2 3', 'CENTER 9 9 9'))
    with pytest.raises(ValueError, match='grid geometry'):
        engines.plan_engine_run('autodock4', inp, {}, tmp_path / 'out')
    log = tmp_path / 'dock.dlg'
    log.write_text('DOCKED: USER    Estimated Free Energy of Binding    =   -8.12 kcal/mol\n')
    assert engines.parse_engine_output('autodock4', log)[0]['metrics'][0]['value'] == -8.12


def topology(path, atoms):
    names = ''.join('C   ' for _ in range(atoms))
    charges = ''.join(f'{0.:16.8E}' for _ in range(atoms))
    path.write_text(f'%VERSION VERSION_STAMP = V0001.000\n%FLAG POINTERS\n%FORMAT(10I8)\n{atoms:8}\n%FLAG ATOM_NAME\n%FORMAT(20a4)\n{names}\n%FLAG CHARGE\n%FORMAT(5E16.8)\n{charges}\n')
    return str(path)


def mm_input(tmp_path):
    inp = {name: topology(tmp_path / (name + '.prmtop'), n) for name, n in [('solvated_topology', 20), ('complex_topology', 10), ('receptor_topology', 7), ('ligand_topology', 3)]}
    traj = tmp_path / 'trajectory.mdcrd'
    traj.write_text('Reviewed example trajectory\n1 2 3\n')
    inp['trajectory'] = str(traj)
    inp['review'] = {'equilibrated_trajectory': True, 'matching_topologies_and_atom_order': True, 'force_field': 'ff14SB', 'ligand_parameterization': 'GAFF2/AM1-BCC', 'sampling_notes': 'Reviewed production MD with convergence checks.', 'reviewed_input_sha256': {name: science.sha256(Path(path)) for name, path in inp.items()}}
    return inp


def test_mm_requires_hash_bound_review_and_never_uses_docking_pose(tmp_path, monkeypatch):
    monkeypatch.setattr(engines, '_executable', lambda name: None)
    with pytest.raises(ValidationError):
        engines.plan_engine_run('mmpbsa', dock_input(tmp_path), {}, tmp_path / 'out')
    inp = mm_input(tmp_path)
    plan = engines.plan_engine_run('mmpbsa', inp, {'end_frame': 10}, tmp_path / 'out')
    assert '&pb' in plan['generated_files']['mmpbsa.in']
    assert '&gb' not in plan['generated_files']['mmpbsa.in']
    assert '-sp' in plan['commands'][0] and '-y' in plan['commands'][0]
    Path(inp['trajectory']).write_text('Changed after review')
    with pytest.raises(ValueError, match='exact SHA-256'):
        engines.plan_engine_run('mmpbsa', inp, {}, tmp_path / 'out')


def test_mm_parser_uses_correct_method_and_difference_not_complex_total(tmp_path):
    result = tmp_path / 'FINAL_RESULTS.dat'
    result.write_text('All units are reported in kcal/mole.\nGENERALIZED BORN:\nComplex:\nTOTAL -1000 10 1\nDifferences (Complex - Receptor - Ligand):\nDELTA TOTAL -22.1 3.2 0.4\nPOISSON BOLTZMANN:\nComplex:\nTOTAL -2000 10 1\nDifferences (Complex - Receptor - Ligand):\nDELTA G binding = -12.5 +/- 2.0 0.2\n')
    assert engines.parse_engine_output('mmgbsa', result)[0]['metrics'][0]['value'] == -22.1
    row = engines.parse_engine_output('mmpbsa', result)[0]
    assert row['metrics'][0]['value'] == -12.5
    assert row['metrics'][0]['standard_error'] == .2
    assert row['entropy_included'] is False
    result.write_text('GENERALIZED BORN:\nComplex:\nTOTAL -1000 10 1\n')
    with pytest.raises(ValueError, match='complex-minus-components'):
        engines.parse_engine_output('mmgbsa', result)


def test_planet_native_csv_is_prediction_not_pose(tmp_path):
    path = tmp_path / 'result.csv'
    path.write_text(',mol_name,SMILES,PLANET_affinity\n0,ethanol,CCO,4.21\n')
    row = engines.parse_engine_output('planet', path)[0]
    assert row['generated_docking_pose'] is False
    assert row['metrics'][0]['unit'] == 'pK'
    path.write_text('mol_name,SMILES,PLANET_affinity\nethanol,CCO,nan\n')
    with pytest.raises(ValueError, match='Nonfinite'):
        engines.parse_engine_output('planet', path)


def test_runner_records_real_outputs_and_failure(tmp_path, monkeypatch):
    monkeypatch.setattr(engines, '_executable', lambda name: '/configured/vina')
    inp = dock_input(tmp_path)
    def runner(argv, cwd, log, timeout, cancelled):
        assert timeout == 120
        assert callable(cancelled)
        log.write_text('actual captured engine log')
        Path(argv[argv.index('--out') + 1]).write_text('REMARK VINA RESULT: -7.2 0 0\n')
        return {'argv': argv, 'returncode': 0, 'duration_seconds': .1}
    monkeypatch.setattr(science, 'run_process', runner)
    result = engines.run_engine('vina', inp, dock_parameters(), tmp_path / 'out', cancelled=lambda: False)
    assert result['status'] == 'completed'
    assert result['results'][0]['metrics'][0]['value'] == -7.2
    assert json.loads((tmp_path / 'out' / 'result.json').read_text())['executed']
    with pytest.raises(ValueError, match='immutable'):
        engines.run_engine('vina', inp, dock_parameters(), tmp_path / 'out')
    monkeypatch.setattr(science, 'run_process', lambda *a, **kw: (_ for _ in ()).throw(science.ScienceError('cancelled', 'Cancelled')))
    result = engines.run_engine('vina', inp, dock_parameters(), tmp_path / 'cancelled')
    assert result['status'] == 'cancelled'
    assert result['results'] == []


def test_mm_rejects_same_count_but_different_atom_identity(tmp_path, monkeypatch):
    monkeypatch.setattr(engines, '_executable', lambda name: None)
    inp = mm_input(tmp_path)
    ligand = Path(inp['ligand_topology'])
    ligand.write_text(ligand.read_text().replace('C   C   C   ', 'N   C   C   '))
    inp['review']['reviewed_input_sha256']['ligand_topology'] = science.sha256(ligand)
    with pytest.raises(ValueError, match='atom names/charges'):
        engines.plan_engine_run('mmgbsa', inp, {}, tmp_path / 'out')


def test_planet_fixed_native_cli_and_weight_hash_guard(tmp_path, monkeypatch):
    import sys
    install = tmp_path / 'planet'
    install.mkdir()
    for filename in ('PLANET_run.py', 'PLANET_model.py', 'layers.py', 'chemutils.py', 'nnutils.py'):
        (install / filename).write_text('# contract fixture, never executed\n')
    weights = install / 'PLANET.param'
    weights.write_bytes(b'contract fixture, not model weights')
    monkeypatch.setenv('BIOCENTRIC_PLANET_DIR', str(install))
    monkeypatch.setenv('BIOCENTRIC_PLANET_PYTHON', sys.executable)
    monkeypatch.setenv('BIOCENTRIC_PLANET_WEIGHTS_SHA256', '0' * 64)
    inp = {'receptor': pdbqt(tmp_path / 'protein.pdb'), 'pocket_reference': sdf(tmp_path / 'reference.sdf'),
           'ligands': sdf(tmp_path / 'queries.sdf', is3d=False), 'preparation_reviewed': True}
    receptor = Path(inp['receptor'])
    receptor.write_text(receptor.read_text().replace('  C   LIG', '  CA  ALA'))
    plan = engines.plan_engine_run('planet', inp, {}, tmp_path / 'out')
    assert plan['status'] == 'blocked'
    assert any('weights' in message for message in plan['blockers'])
    monkeypatch.setenv('BIOCENTRIC_PLANET_WEIGHTS_SHA256', science.sha256(weights))
    plan = engines.plan_engine_run('planet', inp, {}, tmp_path / 'out')
    assert plan['status'] == 'ready'
    assert plan['commands'][0][1] == str(install / 'PLANET_run.py')
    assert '--prefix' in plan['commands'][0]
    assert plan['validated_inputs']['installation_PLANET.param']['sha256'] == science.sha256(weights)
    assert plan['output_metrics'][0]['unit'] == 'pK'


def test_runner_does_not_accept_input_mutation_during_execution(tmp_path, monkeypatch):
    monkeypatch.setattr(engines, '_executable', lambda name: '/configured/vina')
    inp = dock_input(tmp_path)
    def runner(argv, cwd, log, timeout, cancelled):
        Path(inp['receptor']).write_text('changed during calculation')
        log.write_text('tool output')
        Path(argv[argv.index('--out') + 1]).write_text('REMARK VINA RESULT: -7 0 0\n')
        return {'returncode': 0}
    monkeypatch.setattr(science, 'run_process', runner)
    result = engines.run_engine('vina', inp, dock_parameters(), tmp_path / 'out')
    assert result['status'] == 'failed'
    assert result['results'] == []
    assert 'changed during execution' in result['error']['message']


def native_manifest_fixture(tmp_path):
    from biocentric.service import Service
    service = Service(tmp_path / 'store')
    pid = service.create_project({'name': 'native fixture', 'target': 'P00520', 'organism': 'Mus musculus'})
    source_pdb = Path(pdbqt(tmp_path / 'receptor.pdb'))
    ligand = Path(sdf(tmp_path / 'ligands.sdf', {'minimizedAffinity': -7, 'CNNscore': .8, 'CNNaffinity': 6.1}))
    receptor_artifact = service.store.artifact(pid, 'receptor.pdb', source_pdb.read_bytes(), 'native_input')
    ligand_artifact = service.store.artifact(pid, 'ligands.sdf', ligand.read_bytes(), 'native_input')
    manifest = {'schema_version': 1, 'engine_id': 'gnina', 'phase': 'pilot',
                'input_artifacts': {'receptor': receptor_artifact['id'], 'ligand': ligand_artifact['id']},
                'parameters': dock_parameters(), 'review': {'preparation_reviewed': True}, 'human_reviewed': True,
                'canonical_ligands_artifact': ligand_artifact['id'], 'validation_receptor_artifact': receptor_artifact['id'],
                'reference_artifact': ligand_artifact['id'], 'testset_sha256': 'a' * 64, 'protocol_sha256': 'b' * 64}
    art = service.store.artifact(pid, 'native-pilot.json', json.dumps(manifest).encode(), 'native_engine_pilot')
    return service, pid, art, manifest


def test_native_manifest_project_and_frozen_protocol_guards(tmp_path, monkeypatch):
    from biocentric.engine_workflow import plan_native_manifest
    monkeypatch.setattr(engines, '_executable', lambda name: None)
    service, pid, art, manifest = native_manifest_fixture(tmp_path)
    plan = plan_native_manifest(service.store, pid, art, tmp_path / 'run', expected_engine='gnina', expected_phase='pilot', expected_testset_sha256='a' * 64, expected_protocol_sha256='b' * 64)
    assert plan['status'] == 'blocked'  # dependency absent, valid prepared artifacts
    assert not (tmp_path / 'run').exists()
    with pytest.raises(ValueError, match='frozen test set'):
        plan_native_manifest(service.store, pid, art, tmp_path / 'run', expected_engine='gnina', expected_phase='pilot', expected_testset_sha256='c' * 64)
    with pytest.raises(ValueError, match='engine/phase'):
        plan_native_manifest(service.store, pid, art, tmp_path / 'run', expected_engine='vina', expected_phase='pilot')
    pid2 = service.create_project({'name': 'other', 'target': 'P00520', 'organism': 'Mus musculus'})
    with pytest.raises(ValueError, match='not registered'):
        plan_native_manifest(service.store, pid2, art, tmp_path / 'run', expected_engine='gnina', expected_phase='pilot')


def test_native_workflow_calls_adapter_validates_and_registers_canonical_poses(tmp_path, monkeypatch):
    from biocentric.engine_workflow import run_native_manifest
    monkeypatch.setattr(engines, '_executable', lambda name: '/configured/gnina')
    service, pid, art, manifest = native_manifest_fixture(tmp_path)
    calls = []
    def runner(argv, cwd, log, timeout, cancelled):
        calls.append(argv)
        Path(argv[argv.index('--out') + 1]).write_bytes(Path(argv[argv.index('--ligand') + 1]).read_bytes())
        log.write_text('contract test engine output')
        return {'argv': argv, 'returncode': 0, 'duration_seconds': .1}
    monkeypatch.setattr(science, 'run_process', runner)
    monkeypatch.setattr(science, 'validate_pose', lambda *a, **k: {'status': 'completed', 'passed': True, 'tool': 'PoseBusters'})
    result = run_native_manifest(service.store, pid, art, tmp_path / 'run', expected_engine='gnina', expected_phase='pilot', expected_testset_sha256='a' * 64, expected_protocol_sha256='b' * 64)
    assert len(calls) == 1 and '--cnn_scoring' in calls[0]
    assert result['status'] == 'completed'
    assert len(result['poses']) == 1
    pose = result['poses'][0]
    assert pose['score'] == -7 and pose['validation']['passed'] is True
    assert pose['redocking']['rmsd_angstrom'] == pytest.approx(0, abs=1e-6)
    assert pose['pose_artifact_id'] in {a['id'] for a in service.store.artifacts(pid)}
    assert '/artifacts/' in pose['pose_path']
    assert result['result_artifact_id']


def test_native_screening_requires_frozen_library_and_explicit_exclusions(tmp_path, monkeypatch):
    from biocentric.engine_workflow import plan_native_manifest
    from rdkit import Chem
    from rdkit.Chem import AllChem
    monkeypatch.setattr(engines, '_executable', lambda name: None)
    service, pid, art, manifest = native_manifest_fixture(tmp_path)
    original = tmp_path / 'original.sdf'
    ethanol = Chem.SDMolSupplier(str(tmp_path / 'ligands.sdf'), removeHs=False)[0]
    caffeine = Chem.AddHs(Chem.MolFromSmiles('Cn1c(=O)c2c(ncn2C)n(C)c1=O'))
    AllChem.EmbedMolecule(caffeine, randomSeed=42)
    with Chem.SDWriter(str(original)) as writer:
        writer.write(ethanol)
        writer.write(caffeine)
    library = service.store.artifact(pid, 'original.sdf', original.read_bytes(), 'library')
    from biocentric.testsets import build_testset
    frozen = build_testset(service.snapshot(pid), organism='Mus musculus')
    manifest['testset_sha256'] = frozen['manifest_sha256']
    service.store.update(pid, data={'research': {'testset': frozen}})
    manifest['phase'] = 'screening'
    def register():
        return service.store.artifact(pid, 'native-screening.json', json.dumps(manifest).encode(), 'native_engine_screening')
    with pytest.raises(ValueError, match='explicit exclusion'):
        plan_native_manifest(service.store, pid, register(), tmp_path / 'out', expected_engine='gnina', expected_phase='screening')
    manifest['exclusions'] = [{'record_index': 2, 'canonical_smiles': science._canonical(caffeine), 'reason': 'Excluded for a separately documented chemistry review.'}]
    plan = plan_native_manifest(service.store, pid, register(), tmp_path / 'out', expected_engine='gnina', expected_phase='screening')
    assert plan['frozen_input_binding']['verified'] is True
    assert plan['frozen_input_binding']['selected_chemical_states'] == 1
    assert plan['manifest_artifact']['kind'] == 'native_engine_screening'
    manifest['exclusions'][0]['canonical_smiles'] = 'CC'
    with pytest.raises(ValueError, match='Exclusion chemical identities'):
        plan_native_manifest(service.store, pid, register(), tmp_path / 'out', expected_engine='gnina', expected_phase='screening')


def test_native_protocol_hash_cannot_substitute_receptor_or_pilot_settings(tmp_path, monkeypatch):
    from biocentric.engine_workflow import plan_native_manifest
    monkeypatch.setattr(engines, '_executable', lambda name: None)
    service, pid, art, manifest = native_manifest_fixture(tmp_path)
    receptor = next(a for a in service.store.artifacts(pid) if a['id'] == manifest['input_artifacts']['receptor'])
    service.store.update(pid, data={'protocol': {'input_hashes': {'receptor': 'c' * 64, 'structure': 'd' * 64}}})
    with pytest.raises(ValueError, match='Native receptor differs'):
        plan_native_manifest(service.store, pid, art, tmp_path / 'out', expected_engine='gnina', expected_phase='pilot')
    service.store.update(pid, data={'protocol': {'input_hashes': {'receptor': receptor['sha256'], 'structure': receptor['sha256']}}})
    pilot_plan = plan_native_manifest(service.store, pid, art, tmp_path / 'out', expected_engine='gnina', expected_phase='pilot')
    service.store.update(pid, data={'research': {'pilot': {'native_result': {'engine_id': 'gnina', 'parameters': pilot_plan['engine_plan']['parameters']}}}})
    manifest['phase'] = 'screening'
    manifest['parameters']['exhaustiveness'] = 8
    screen = service.store.artifact(pid, 'screen.json', json.dumps(manifest).encode(), 'native_engine_screening')
    with pytest.raises(ValueError, match='settings differ from the tested pilot'):
        plan_native_manifest(service.store, pid, screen, tmp_path / 'out', expected_engine='gnina', expected_phase='screening')
    manifest['parameters']['exhaustiveness'] = 1
    manifest['parameters']['timeout_seconds'] = 240
    screen = service.store.artifact(pid, 'screen.json', json.dumps(manifest).encode(), 'native_engine_screening')
    plan = plan_native_manifest(service.store, pid, screen, tmp_path / 'out', expected_engine='gnina', expected_phase='screening')
    assert plan['engine_plan']['parameters']['timeout_seconds'] == 240


def test_native_inline_approval_resumes_same_attempt_without_reprompt(tmp_path, monkeypatch):
    from biocentric import engine_workflow
    from biocentric.research_workflow import AwaitingAnswer, ResearchWorkflow
    from biocentric.store import digest, dump
    monkeypatch.setattr(engines, '_executable', lambda name: '/configured/gnina')
    service, pid, art, manifest = native_manifest_fixture(tmp_path)
    reference = next(a for a in service.store.artifacts(pid) if a['id'] == manifest['reference_artifact'])
    protocol = {'testset_sha256': 'a' * 64, 'input_hashes': {'reference': reference['sha256']}}
    manifest['protocol_sha256'] = digest(dump(protocol).encode())
    service.store.artifact(pid, 'pilot-approved.json', json.dumps(manifest).encode(), 'native_engine_pilot')
    workflow = ResearchWorkflow(service.store, pid, 'test-job')
    calls = []
    def execute(store, project_id, manifest_artifact, output_dir, **kwargs):
        calls.append(str(output_dir))
        Path(output_dir).mkdir(parents=True)
        return {'status': 'completed', 'poses': []}
    monkeypatch.setattr(engine_workflow, 'run_native_manifest', execute)
    with pytest.raises(AwaitingAnswer):
        workflow.native_calculation('gnina', 'pilot', protocol)
    assert not calls
    questions = service.store.project(pid)['data']['research']['questions']
    question_id = questions[-1]['id']
    questions[-1].update(status='answered', answer='approve')
    service.store.research_update(pid, questions=questions)
    resumed = ResearchWorkflow(service.store, pid, 'test-job')
    result = resumed.native_calculation('gnina', 'pilot', protocol)
    assert result['status'] == 'completed' and calls[0].endswith('attempt-1')
    assert len(service.store.project(pid)['data']['research']['questions']) == 1
    with pytest.raises(AwaitingAnswer):
        resumed.native_calculation('gnina', 'pilot', protocol)
    pending = service.store.project(pid)['data']['research']['questions'][-1]
    assert pending['id'] != question_id and pending['status'] == 'pending'
    assert len(calls) == 1  # a new execution attempt needs its own concrete review


def test_fast_search_limits_use_engine_specific_native_flags(tmp_path, monkeypatch):
    monkeypatch.setattr(engines, '_executable', lambda name: None)
    vina = engines.plan_engine_run('vina', dock_input(tmp_path), {**dock_parameters(), 'max_evals': 5000}, tmp_path/'v')
    assert vina['commands'][0][-2:] == ['--max_evals', '5000']
    inp = {'receptor': pdbqt(tmp_path/'receptor.pdb'), 'ligand': sdf(tmp_path/'pose.sdf'), 'preparation_reviewed': True}
    gnina = engines.plan_engine_run('gnina', inp, {**dock_parameters(), 'max_mc_steps': 20}, tmp_path/'g')
    argv = gnina['commands'][0]
    assert argv[argv.index('--max_mc_steps')+1] == '20'
    assert '--max_evals' not in argv
    from pydantic import ValidationError
    with pytest.raises(ValidationError):
        engines.plan_engine_run('gnina', inp, {**dock_parameters(), 'max_evals': 5000}, tmp_path/'bad')
