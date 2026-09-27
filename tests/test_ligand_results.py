from copy import deepcopy
import hashlib
from pathlib import Path

import pytest

from biocentric.ligand_results import (build_ligand_results, export_pose_sdf,
                                      ligand_csv, read_registered, render_ligand_results)

PID = 'c' * 32


def snapshot():
    return {'project': {'id': PID, 'state': 'running'}, 'data': {}, 'candidates': [], 'artifacts': [],
            'events': [], 'jobs': [{'id': 'job-a', 'status': 'running'}]}


def molecule(name='ethanol', smiles='CCO'):
    from rdkit import Chem
    mol = Chem.MolFromSmiles(smiles)
    mol.SetProp('_Name', name)
    return (Chem.MolToMolBlock(mol) + '\n$$$$\n').encode()


def register(snap, root, name, content, kind='scientific_output'):
    sha = hashlib.sha256(content).hexdigest()
    path = root / 'projects' / PID / 'artifacts' / sha / name
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(content)
    row = {'id': hashlib.sha256((sha + name + kind).encode()).hexdigest()[:32], 'project_id': PID,
           'name': name, 'sha256': sha, 'path': str(path), 'kind': kind}
    snap['artifacts'].append(row)
    return row


def pose(score=-7.0, *, compound='cmp_ethanol', passed=True, number=1, engine='AutoDock Vina', direction='lower_is_better', job='job-a'):
    return {'compound_id': compound, 'name': compound.removeprefix('cmp_'), 'pose_id': f'{compound}_{number}',
            'pose_number': number, 'score': score, 'score_units': 'kcal/mol', 'score_direction': direction,
            'engine': engine, 'engine_version': '1.2.7', 'validation': {'passed': passed}, 'job_id': job,
            'receptor_sha256': 'r' * 64, 'pose_sha256': 'a' * 64}


def test_empty_and_prepared_are_visible_without_fabricated_scores():
    snap = snapshot()
    assert build_ligand_results(snap)['ranked_rows'] == []
    snap['data']['preparation'] = {'compounds': [{'source_record': 1, 'compound_id': 'cmp_a', 'name': 'A', 'metadata': {'CONTROL_ROLE': 'untested control'}}],
                                   'exclusions': [{'record_index': 2, 'name': 'bad input', 'reason': 'Unspecified stereochemistry'}]}
    result = build_ligand_results(snap)
    assert result['ranked_rows'] == []
    assert result['input_rows'][0]['selected_for_calculation'] is True
    assert result['input_rows'][0]['status'] == 'Prepared · awaiting docking'
    assert result['input_rows'][1]['status'] == 'Excluded'
    assert 'stereochemistry' in result['input_rows'][1]['reason']


def test_registered_library_is_read_before_preparation_and_draft_stays_distinct(tmp_path):
    snap = snapshot()
    old = register(snap, tmp_path, 'library.sdf', molecule('ethanol'), 'library')
    snap['data']['research'] = {'testset': {'input_artifacts': [{'kind': 'library', 'artifact_id': old['id']}],
                                         'records': [{'compound_id': 'measured', 'partition': 'evaluation', 'endpoint': 'Ki', 'value': 2, 'units': 'nM', 'relation': '='}]}}
    result = build_ligand_results(snap, tmp_path)
    assert result['input_rows'][0]['name'] == 'ethanol'
    assert result['input_rows'][0]['selected_for_calculation'] is None
    assert result['test_rows'][0]['compound_id'] == 'measured'
    draft = register(snap, tmp_path, 'new.sdf', molecule('methane', 'C'), 'library')
    result = build_ligand_results(snap, tmp_path)
    assert result['library_artifact']['id'] == old['id']
    assert result['draft_library']['id'] == draft['id']
    assert result['draft_rows'][0]['Ligand'] == 'methane'
    assert 'not included' in result['draft_rows'][0]['Status']


def test_live_checkpoint_poses_are_provisional_and_pilot_does_not_leak():
    snap = snapshot()
    snap['data']['preparation'] = {'compounds': [{'source_record': 1, 'compound_id': 'cmp_ethanol', 'name': 'ethanol'}]}
    snap['events'] = [{'stage': '06', 'outputs': {'poses': [pose(-100)]}},
                      {'stage': '08', 'kind': 'progress', 'status': 'running', 'details': {'compound_id': 'cmp_ethanol'},
                       'outputs': {'poses': [pose(-7)], 'checkpoint_reused': True}}]
    result = build_ligand_results(snap)
    assert result['provisional']
    assert result['ranked_rows'][0]['best_valid_score'] == -7
    assert result['ranked_rows'][0]['selected'] is None
    assert result['input_rows'][0]['status'] == 'Docked · provisional'


def test_best_valid_pose_beats_failed_low_score_and_selection_is_recorded():
    snap = snapshot()
    invalid, valid = pose(-9, passed=False), pose(-7, number=2)
    snap['data']['screening'] = {'poses': [invalid, valid]}
    snap['candidates'] = [{**valid, 'selected': True, 'rationale': 'Valid pose and diverse scaffold'}]
    result = build_ligand_results(snap)
    row = result['ranked_rows'][0]
    assert row['best_valid_score'] == -7 and row['poses'] == 2
    assert row['selected'] is True and result['selected_count'] == 1
    assert row['selected_pose']['pose_id'] == valid['pose_id']


def test_inconsistent_selected_pose_is_flagged_and_not_exported_as_kept():
    snap = snapshot()
    invalid, valid = pose(-9, passed=False), pose(-7, number=2)
    snap['data']['screening'] = {'poses': [invalid, valid]}
    snap['candidates'] = [{**invalid, 'selected': True}]
    row = build_ligand_results(snap)['ranked_rows'][0]
    assert row['recorded_selection'] is True
    assert row['selected'] is False and row['selected_pose'] is None
    assert 'did not pass' in row['selection_issue']
    snap['candidates'][0] = {**valid, 'pose_sha256': 'wrong', 'selected': True}
    row = build_ligand_results(snap)['ranked_rows'][0]
    assert row['selected'] is False and 'no matching' in row['selection_issue']


def test_unknown_scores_stay_visible_and_failed_execution_has_reason():
    snap = snapshot()
    snap['data']['preparation'] = {'compounds': [{'source_record': 1, 'compound_id': 'cmp_ethanol', 'name': 'ethanol'},
                                               {'source_record': 2, 'compound_id': 'cmp_failed', 'name': 'failed ligand'}]}
    snap['data']['screening'] = {'poses': [pose(-7, direction='unspecified'), pose(float('nan'), compound='cmp_nan')],
                                'failures': [{'compound_id': 'cmp_failed', 'reason': 'Engine timeout'}]}
    result = build_ligand_results(snap)
    assert len(result['ranked_rows']) == 2
    assert all(row['rank'] is None for row in result['ranked_rows'])
    assert all(row['best_valid_score'] is None for row in result['ranked_rows'])
    assert result['input_rows'][0]['status'] == 'Docked · score needs review'
    assert result['input_rows'][1]['reason'] == 'Engine timeout'


def test_runs_and_engines_are_not_pooled_and_cnn_scores_stay_named():
    snap = snapshot()
    vina, gnina = pose(-8), pose(-6, compound='cmp_other', engine='gnina')
    gnina['metrics'] = [{'name': 'minimizedAffinity', 'value': -6, 'unit': 'kcal/mol', 'direction': 'lower_is_better'},
                        {'name': 'CNNaffinity', 'value': 6.5, 'unit': 'pK', 'direction': 'higher_is_better'},
                        {'name': 'CNNscore', 'value': .8, 'unit': 'dimensionless', 'direction': 'higher_is_better'}]
    other_run = pose(-12, compound='cmp_third', job='job-b')
    snap['data']['screening'] = {'poses': [vina, gnina, other_run]}
    result = build_ligand_results(snap)
    assert len(result['groups']) == 3
    assert [row['rank'] for row in result['ranked_rows']] == [1, 1, 1]
    cnn = next(row for row in result['ranked_rows'] if row['engine'] == 'gnina')
    assert cnn['best_valid_score'] == -6 and cnn['cnn_affinity_pK'] == 6.5


def test_exports_only_registered_verified_single_pose_sdf_and_ignores_pose_path(tmp_path):
    snap = snapshot()
    content = molecule()
    artifact = register(snap, tmp_path, 'pose.sdf', content)
    record = {**pose(), 'artifact_ids': [artifact['id']], 'pose_sha256': artifact['sha256'], 'pose_path': '/etc/passwd'}
    assert export_pose_sdf(snap, tmp_path, PID, [record]) == content.rstrip() + b'\n'
    assert export_pose_sdf(snap, tmp_path, PID, [record, record]).count(b'$$$$') == 2
    with pytest.raises(ValueError, match='registered SDF'):
        export_pose_sdf(snap, tmp_path, PID, [{**record, 'artifact_ids': []}])
    multi = register(snap, tmp_path, 'multi.sdf', content * 2)
    with pytest.raises(ValueError, match='more than one'):
        export_pose_sdf(snap, tmp_path, PID, [{**record, 'artifact_ids': [multi['id']], 'pose_sha256': multi['sha256']}])


def test_artifact_cache_detects_changes_and_blocks_cross_project_or_escaped_paths(tmp_path):
    snap = snapshot()
    artifact = register(snap, tmp_path, 'pose.sdf', molecule())
    assert read_registered(snap, tmp_path, PID, artifact['id'])
    path = Path(artifact['path'])
    path.write_bytes(path.read_bytes().replace(b'ethanol', b'methane'))
    with pytest.raises(ValueError, match='hash'):
        read_registered(snap, tmp_path, PID, artifact['id'])
    artifact['project_id'] = 'other'
    with pytest.raises(ValueError, match='not registered'):
        read_registered(snap, tmp_path, PID, artifact['id'])
    artifact['project_id'] = PID
    outside = tmp_path / 'outside.sdf'
    outside.write_bytes(molecule())
    artifact['path'] = str(outside)
    with pytest.raises(ValueError, match='outside'):
        read_registered(snap, tmp_path, PID, artifact['id'])


def test_csv_preserves_numeric_scores_and_neutralizes_text_formula():
    snap = snapshot()
    snap['data']['screening'] = {'poses': [{**pose(), 'name': '=SUM(A1:A2)'}]}
    content = ligand_csv(build_ligand_results(snap)['ranked_rows']).decode()
    assert "'=SUM(A1:A2)" in content
    assert ',-7.0,' in content
    assert 'kcal/mol' in content and 'lower_is_better' in content


class UI:
    def __init__(self):
        self.tables, self.downloads, self.messages = [], [], []
    def __enter__(self): return self
    def __exit__(self, *args): return False
    def dataframe(self, rows, **kwargs): self.tables.append(rows)
    def download_button(self, label, data, **kwargs): self.downloads.append((label, data, kwargs))
    def expander(self, *args, **kwargs): return self
    def columns(self, count): return [self] * count
    def selectbox(self, label, options, **kwargs): return options[0]
    def __getattr__(self, name): return lambda value, **kwargs: self.messages.append((name, value))


def test_renderer_shows_compact_table_and_exact_download_counts_without_mutation(tmp_path):
    snap = snapshot()
    artifact = register(snap, tmp_path, 'pose.sdf', molecule())
    p = {**pose(), 'artifact_ids': [artifact['id']], 'pose_sha256': artifact['sha256']}
    snap['data']['screening'] = {'poses': [p, {**p, 'pose_number': 2, 'pose_id': 'second', 'score': -6}]}
    snap['candidates'] = [{**p, 'selected': True, 'rationale': 'Recorded selection'}]
    before = deepcopy(snap)
    st = UI()
    render_ligand_results(st, snap, tmp_path, PID)
    assert len(st.tables[0][0]) == 5
    downloads = {label: (content, options) for label, content, options in st.downloads}
    assert downloads['Kept ligands · SDF'][0].count(b'$$$$') == 1
    assert downloads['All poses · SDF'][0].count(b'$$$$') == 2
    assert 'Ranked CSV' in downloads
    assert snap == before


def test_native_preflight_and_exclusions_update_full_input_list(tmp_path):
    snap = snapshot()
    register(snap, tmp_path, 'library.sdf', molecule('ethanol') + molecule('methane', 'C'), 'library')
    snap['data']['research'] = {'engine_preflight': {'phase': 'screening', 'validation': {'identities': ['CCO']},
                            'frozen_input_binding': {'exclusions': [{'record_index': 2, 'canonical_smiles': 'C', 'reason': 'Excluded from native protocol after chemistry review.'}]}}}
    result = build_ligand_results(snap, tmp_path)
    assert result['input_rows'][0]['selected_for_calculation'] is True
    assert result['input_rows'][0]['status'].startswith('Native inputs checked')
    assert result['input_rows'][1]['status'] == 'Excluded'
    assert result['input_rows'][1]['selected_for_calculation'] is False


def test_unknown_units_and_unscoped_runs_are_not_ranked():
    snap = snapshot()
    snap['data']['screening'] = {'poses': [{**pose(), 'score_units': 'arbitrary'}, {**pose(compound='cmp_other'), 'job_id': None}]}
    result = build_ligand_results(snap)
    assert len(result['ranked_rows']) == 2
    assert all(row['rank'] is None for row in result['ranked_rows'])
