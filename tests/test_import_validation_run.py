"""Import integrity tests use synthetic saved evidence, never run docking or LLMs."""
from copy import deepcopy
import json
from pathlib import Path

import pytest

from biocentric import science
from biocentric.ligand_results import build_ligand_results, export_pose_sdf
from biocentric.service import Service
from biocentric.store import digest
from scripts.import_validation_run import import_run


@pytest.fixture(autouse=True)
def no_engine_execution(monkeypatch):
    monkeypatch.setattr(science, 'run_process', lambda *args, **kwargs: pytest.fail('Imports must never execute scientific tools'))


def write_json(path, data):
    path.write_text(json.dumps(data, indent=2) + '\n')


def artifact(path, role='scientific_output'):
    content = path.read_bytes()
    return {'path': str(path.resolve()), 'sha256': digest(content), 'bytes': len(content), 'role': role}


def saved_bundle(tmp_path):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    base = tmp_path / 'source'
    baseline = base / 'baseline'
    prep_dir = base / 'preparation'
    baseline.mkdir(parents=True)
    prep_dir.mkdir()
    library = base / 'selected_2.smi'
    library.write_text('CCO ethanol\nCCN ethylamine\n')
    library_sha = digest(library.read_bytes())
    prepared_path = prep_dir / 'prepared.sdf'
    compounds, molecules = [], []
    for index, smiles in enumerate(('CCO', 'CCN'), 1):
        molecule = Chem.AddHs(Chem.MolFromSmiles(smiles))
        AllChem.EmbedMolecule(molecule, randomSeed=42)
        canonical = science._canonical(molecule)
        compound_id = 'cmp_' + digest(canonical.encode())[:16]
        name = 'ethanol' if index == 1 else 'ethylamine'
        molecule.SetProp('_Name', name)
        molecule.SetProp('BIOCENTRIC_COMPOUND_ID', compound_id)
        molecules.append(molecule)
        compounds.append({'compound_id': compound_id, 'canonical_smiles': canonical, 'name': name,
                          'source_record': index, 'source_sha256': library_sha, 'metadata': {}})
    with Chem.SDWriter(str(prepared_path)) as writer:
        for molecule in molecules:
            writer.write(molecule)
    receptor = baseline / 'receptor.pdbqt'
    receptor.write_text('Recorded receptor fixture, never sent to docking\n')
    structure = baseline / 'receptor.pdb'
    structure.write_text('Recorded validation receptor fixture\n')
    reference = baseline / 'reference.sdf'
    with Chem.SDWriter(str(reference)) as writer:
        writer.write(molecules[0])
    config = {'seed': 42, 'cpu': 1, 'exhaustiveness': 1, 'n_poses': 1, 'timeout_seconds': 60,
              'receptor_pdb': str(structure), 'require_posebusters': True, 'run_plip': False}
    identity = {'adapter_sha256': 'a' * 64, 'receptor': digest(receptor.read_bytes()),
                'library': digest(prepared_path.read_bytes()), 'reference': digest(reference.read_bytes()),
                'receptor_pdb': digest(structure.read_bytes()), 'config': config, 'vina': 'AutoDock Vina fixture',
                'versions': {'rdkit': 'fixture', 'posebusters': 'fixture'}}
    fingerprint = digest(json.dumps(identity, sort_keys=True).encode())
    request = {**identity, 'fingerprint': fingerprint, 'box': {'center': [1, 2, 3], 'size': [12, 12, 12]}}
    request_path = baseline / 'request.json'
    write_json(request_path, request)
    poses, native_artifacts = [], [artifact(path) for path in (request_path, receptor, structure, reference)]
    for index, (molecule, compound) in enumerate(zip(molecules, compounds), 1):
        folder = baseline / compound['compound_id']
        folder.mkdir()
        score = -8 + index
        molecule.SetDoubleProp('BIOCENTRIC_VINA_SCORE_KCAL_MOL', score)
        pose_path = folder / 'pose.sdf'
        with Chem.SDWriter(str(pose_path)) as writer:
            writer.write(molecule)
        physical = {'status': 'completed', 'tool': 'PoseBusters', 'version': 'synthetic fixture',
                    'checks': {'recorded_check': index == 1}, 'passed': index == 1}
        physical_path = folder / 'posebusters.json'
        write_json(physical_path, physical)
        pose_artifact, physical_artifact = artifact(pose_path), artifact(physical_path, 'physical_validity_report')
        native_artifacts += [pose_artifact, physical_artifact]
        poses.append({'compound_id': compound['compound_id'], 'canonical_smiles': compound['canonical_smiles'],
                      'name': compound['name'], 'pose_id': compound['compound_id'] + '_pose_1', 'pose_number': 1,
                      'score': score, 'score_units': 'kcal/mol', 'score_direction': 'lower_is_better',
                      'engine': 'AutoDock Vina', 'engine_version': identity['vina'], 'receptor_sha256': identity['receptor'],
                      'protocol_fingerprint': fingerprint, 'pose_path': str(pose_path), 'pose_sha256': pose_artifact['sha256'],
                      'validation': {**physical, 'artifact': physical_artifact}, 'descriptors': {},
                      'job_id': 'original-native-job', 'artifact_ids': ['original-native-artifact']})
    preparation = {'status': 'completed', 'compounds': compounds, 'exclusions': [], 'input_records': 2,
                   'prepared_sdf': str(prepared_path), 'artifacts': [artifact(library, 'immutable_original'),
                                                                  artifact(prepared_path, 'canonical_prepared_library')]}
    docking = {'status': 'completed', 'poses': poses, 'failures': [], 'artifacts': native_artifacts,
               'metadata': {'fingerprint': fingerprint, 'config': config, 'box': request['box'],
                            'vina': identity['vina'], 'versions': identity['versions']}}
    report = {'status': 'completed', 'target': 'Curated 1IEP mouse ABL pocket; bound ligand defines box only',
              'adversarial_passed': True, 'checkpoint_reuse': {'status': 'passed'},
              'pinned_inputs': {'count': 2, 'receptor_sha256': identity['receptor'], 'reference_sha256': identity['reference'],
                                'receptor_pdb_sha256': identity['receptor_pdb'], 'science_adapter_sha256': identity['adapter_sha256']},
              'preparation': {'requested': 2, 'accepted': 2, 'excluded': 0}, 'selection': {'requested': 2},
              'protocols': {'baseline': {'status': 'completed', 'config': config, 'protocol_fingerprint': fingerprint,
                                         'requested_ligands': 2, 'ligands_with_poses': 2, 'pose_count': 2, 'valid_poses': 1,
                                         'invalid_or_unverified_poses': 1, 'failed_ligands': 0, 'not_completed_ligands': 0,
                                         'score_units': 'kcal/mol', 'score_direction': 'lower_is_better', 'score_range': [-7, -6]}},
              'resources': {'llm_calls': 0, 'network_requests': 0, 'gpu_used': False, 'cpu_per_vina_process': 1, 'budget_seconds': 1200,
                            'hardware': 'Synthetic local fixture; no calculations executed'}}
    (base / 'events.jsonl').write_text('{"kind":"fixture","message":"Synthetic source event only"}\n')
    for filename in ('selection.json', 'parameter_comparisons.json', 'checkpoint_reuse.json', 'adversarial_verification.json'):
        write_json(base / filename, {'fixture': True})
    bundle = {'base': base, 'report_path': base / 'summary.json', 'preparation_path': prep_dir / 'preparation.json',
              'docking_path': baseline / 'docking.json', 'library_path': library,
              'preparation': preparation, 'docking': docking, 'report': report}
    refresh(bundle)
    return bundle


def refresh(bundle):
    write_json(bundle['preparation_path'], bundle['preparation'])
    write_json(bundle['docking_path'], bundle['docking'])
    bundle['report']['evidence_inventory'] = [artifact(path, 'verification_evidence') for path in sorted(bundle['base'].rglob('*'))
                                             if path.is_file() and path != bundle['report_path']]
    write_json(bundle['report_path'], bundle['report'])


def perform(bundle, root):
    return import_run(bundle['report_path'], bundle['preparation_path'], bundle['docking_path'], bundle['library_path'], root=root)


def test_valid_import_preserves_source_provenance_and_downloads_without_reexecution(tmp_path):
    bundle = saved_bundle(tmp_path)
    original_hashes = {str(path): digest(path.read_bytes()) for path in bundle['base'].rglob('*') if path.is_file()}
    root = tmp_path / 'app'
    pid = perform(bundle, root)
    service = Service(root)
    snap = service.snapshot(pid)
    assert snap['project']['state'] == 'completed'
    assert snap['project']['spec']['allow_model'] is False and snap['project']['spec']['max_model_calls'] == 0
    assert snap['project']['spec']['compute_budget_seconds'] == 1200
    assert snap['data']['result_status'] == 'exploratory_software_validation'
    assert snap['data']['research']['validation_run']['source_bundle_sha256']
    assert snap['data']['research']['validation_run']['read_only'] is True
    assert not [event for event in snap['events'] if event.get('kind') == 'model_call']
    assert len([event for event in snap['events'] if event.get('kind') == 'tool_call']) == 1
    imported = snap['data']['screening']['poses'][0]
    assert imported['job_id'] != 'original-native-job'
    assert imported['source_provenance']['job_id'] == 'original-native-job'
    assert imported['source_provenance']['artifact_ids'] == ['original-native-artifact']
    assert imported['source_provenance']['pose_path'] == bundle['docking']['poses'][0]['pose_path']
    assert '/artifacts/' in imported['pose_path']
    assert imported['pose_artifact_id'] == imported['artifact_ids'][0]
    results = build_ligand_results(snap, root)
    assert len(results['ranked_rows']) == 2 and results['selected_count'] == 1
    selected = [row['selected_pose'] for row in results['ranked_rows'] if row['selected']]
    assert export_pose_sdf(snap, root, pid, selected).count(b'$$$$') == 1
    evidence = snap['data']['source_validation']['evidence_artifact_ids']
    assert set(evidence) == {'events.jsonl', 'selection.json', 'parameter_comparisons.json', 'checkpoint_reuse.json', 'adversarial_verification.json'}
    assert all(evidence.values())
    stage8 = next(step for step in snap['steps'] if step['id'] == '08')
    assert stage8['outputs']['event_log_artifact_id'] == evidence['events.jsonl']
    assert {str(path): digest(path.read_bytes()) for path in bundle['base'].rglob('*') if path.is_file()} == original_hashes
    assert perform(bundle, root) == pid
    assert len(service.list_projects()) == 1


@pytest.mark.parametrize('problem', ['partial_report', 'partial_docking', 'failed_preparation', 'wrong_source_library',
                                   'incomplete_compounds', 'wrong_protocol', 'wrong_counts', 'wrong_units',
                                   'wrong_engine', 'changed_score', 'changed_validity', 'duplicate_pose', 'missing_prepared_artifact'])
def test_invalid_semantic_bundle_rejected_before_app_state_exists(problem, tmp_path):
    bundle = saved_bundle(tmp_path)
    report, docking, prep = bundle['report'], bundle['docking'], bundle['preparation']
    if problem == 'partial_report':
        report['status'] = 'partial'
    elif problem == 'partial_docking':
        docking['status'] = 'budget_stopped'
    elif problem == 'failed_preparation':
        prep['status'] = 'failed'
    elif problem == 'wrong_source_library':
        prep['compounds'][0]['source_sha256'] = 'b' * 64
    elif problem == 'incomplete_compounds':
        docking['poses'].pop()
    elif problem == 'wrong_protocol':
        docking['poses'][0]['protocol_fingerprint'] = 'b' * 64
    elif problem == 'wrong_counts':
        report['protocols']['baseline']['valid_poses'] = 2
    elif problem == 'wrong_units':
        docking['poses'][0]['score_units'] = 'pK'
    elif problem == 'wrong_engine':
        docking['poses'][0]['engine'] = 'gnina'
    elif problem == 'changed_score':
        docking['poses'][0]['score'] = -100
    elif problem == 'changed_validity':
        docking['poses'][0]['validation']['passed'] = False
    elif problem == 'duplicate_pose':
        docking['poses'].append(deepcopy(docking['poses'][0]))
    else:
        prep['artifacts'] = [prep['artifacts'][0]]
    refresh(bundle)
    root = tmp_path / 'app'
    with pytest.raises(ValueError):
        perform(bundle, root)
    assert not root.exists()


@pytest.mark.parametrize('problem', ['tampered_pose', 'tampered_event_log', 'unlisted_preparation', 'missing_inventory', 'foreign_docking'])
def test_unbound_or_tampered_source_files_rejected_before_app_state_exists(problem, tmp_path):
    bundle = saved_bundle(tmp_path)
    if problem == 'tampered_pose':
        Path(bundle['docking']['poses'][0]['pose_path']).write_text('tampered molecule')
    elif problem == 'tampered_event_log':
        (bundle['base'] / 'events.jsonl').write_text('rewritten source history')
    elif problem == 'unlisted_preparation':
        bundle['report']['evidence_inventory'] = [row for row in bundle['report']['evidence_inventory'] if row['path'] != str(bundle['preparation_path'])]
        write_json(bundle['report_path'], bundle['report'])
    elif problem == 'missing_inventory':
        bundle['report']['evidence_inventory'] = []
        write_json(bundle['report_path'], bundle['report'])
    else:
        foreign = bundle['base'] / 'baseline_result.json'
        write_json(foreign, bundle['docking'])
        bundle['docking_path'] = foreign
    root = tmp_path / 'app'
    with pytest.raises(ValueError):
        perform(bundle, root)
    assert not root.exists()


def test_canonical_docking_manifest_survives_stale_wrapper_self_checksum(tmp_path):
    bundle = saved_bundle(tmp_path)
    wrapper = {**deepcopy(bundle['docking']), 'artifacts': [*bundle['docking']['artifacts'],
               {'path': str(bundle['docking_path']), 'sha256': '0' * 64, 'bytes': 1, 'role': 'docking_manifest'}]}
    write_json(bundle['base'] / 'baseline_result.json', wrapper)
    refresh(bundle)
    pid = perform(bundle, tmp_path / 'app')
    assert Service(tmp_path / 'app').snapshot(pid)['project']['state'] == 'completed'


def test_idempotent_import_rechecks_registered_bytes(tmp_path):
    bundle = saved_bundle(tmp_path)
    root = tmp_path / 'app'
    pid = perform(bundle, root)
    service = Service(root)
    pose = service.snapshot(pid)['data']['screening']['poses'][0]
    Path(pose['pose_path']).write_text('tampered registered artifact')
    with pytest.raises(ValueError, match='integrity'):
        perform(bundle, root)
    assert len(service.list_projects()) == 1


def test_explicit_compound_failure_is_complete_accounting_not_a_confirmed_hit(tmp_path):
    bundle = saved_bundle(tmp_path)
    failed = bundle['docking']['poses'].pop()
    bundle['docking']['failures'] = [{'compound_id': failed['compound_id'], 'code': 'fixture_failure', 'reason': 'Recorded synthetic failure'}]
    bundle['report']['protocols']['baseline'].update(ligands_with_poses=1, pose_count=1, invalid_or_unverified_poses=0,
                                                   failed_ligands=1, score_range=[-7, -7])
    refresh(bundle)
    pid = perform(bundle, tmp_path / 'app')
    snap = Service(tmp_path / 'app').snapshot(pid)
    assert snap['project']['state'] == 'completed'
    assert snap['data']['screening']['failures'][0]['compound_id'] == failed['compound_id']
    results = build_ligand_results(snap, tmp_path / 'app')
    assert next(row for row in results['input_rows'] if row['compound_id'] == failed['compound_id'])['status'] == 'Docking failed'
