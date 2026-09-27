"""Focused offline checks for the reproducible diversity-run evidence builder."""
import csv
import hashlib
import json
from pathlib import Path
from unittest.mock import patch

import pytest

from biocentric import science
from scripts.validate_diverse_docking import (
    ADVERSARIAL, EXPECTED_ADVERSARIAL, export_result, immutable_text,
    partial_result, prepare_cached, select_diverse, summarize_result,
)


def test_selection_is_seeded_unique_scaffold_and_accounts_for_rejections(tmp_path):
    source = tmp_path / 'source.smi'
    source.write_text('''O=C(O)c1ccccc1\t1
O=C(O)c1ccncc1\t2
O=C(O)c1ccoc1\t3
O=C(O)c1ccsc1\t4
O=C(O)c1cc[nH]c1\t5
c1ccc2ccccc2c1\t6
OC1CCCCC1\t7
CCC1CCCCC1\t8
O=C(O)c1ccccc1\t9
CCO.[Na+]\t10
CC(O)C(=O)O\t11
[Fe]\t12
not_a_smiles\t13
C\t14
''')
    first = select_diverse(source, count=4, seed=42)
    second = select_diverse(source, count=4, seed=42)
    assert first == second
    assert first['source_sha256'] == science.sha256(source)
    assert first['input_records'] == len(first['source_ledger']) == 14
    assert len({row['scaffold'] for row in first['selected']}) == 4
    assert len({row['canonical_smiles'] for row in first['selected']}) == 4
    assert all(row['status'] == 'selected' for row in first['selected'])
    codes = {row.get('code') for row in first['source_ledger']}
    assert {'exact_duplicate', 'unsupported_salt', 'ambiguous_stereochemistry', 'unsupported_element', 'invalid_chemistry', 'bounded_test_domain'} <= codes
    with pytest.raises(ValueError, match='only .* eligible'):
        select_diverse(source, count=100)


def test_existing_input_evidence_is_never_silently_replaced(tmp_path):
    path = tmp_path / 'fixed.smi'
    immutable_text(path, 'CCO ethanol\n')
    before = path.stat().st_mtime_ns
    immutable_text(path, 'CCO ethanol\n')
    assert path.stat().st_mtime_ns == before
    with pytest.raises(ValueError, match='Existing evidence differs'):
        immutable_text(path, 'CCC propane\n')
    assert path.read_text() == 'CCO ethanol\n'


def test_adversarial_preparation_cache_checks_source_and_artifact_hashes(tmp_path):
    source = tmp_path / 'adversarial.smi'
    source.write_text(ADVERSARIAL)
    result = prepare_cached(source, tmp_path / 'prepared', None, None)
    assert len(result['compounds']) == 2
    assert {row['code'] for row in result['exclusions']} == EXPECTED_ADVERSARIAL
    assert len({row['canonical_smiles'] for row in result['compounds']}) == 2
    with patch.object(science, 'prepare_library', side_effect=AssertionError('cache must not regenerate')):
        cached = prepare_cached(source, tmp_path / 'prepared', None, None)
    assert cached['compounds'] == result['compounds']
    source.write_text(ADVERSARIAL + 'CCO extra\n')
    with pytest.raises(ValueError, match='Preparation input changed'):
        prepare_cached(source, tmp_path / 'prepared', None, None)
    source.write_text(ADVERSARIAL)
    Path(result['prepared_sdf']).write_text('modified evidence')
    with pytest.raises(ValueError, match='Artifact integrity failed'):
        prepare_cached(source, tmp_path / 'prepared', None, None)


def pose(tmp_path, name, smiles, score, valid):
    from rdkit import Chem
    mol = Chem.MolFromSmiles(smiles)
    canonical = science._canonical(mol)
    cid = 'cmp_' + hashlib.sha256(canonical.encode()).hexdigest()[:16]
    mol.SetProp('_Name', name)
    mol.SetProp('BIOCENTRIC_COMPOUND_ID', cid)
    mol.SetDoubleProp('BIOCENTRIC_VINA_SCORE_KCAL_MOL', score)
    path = tmp_path / f'{name}.sdf'
    with Chem.SDWriter(str(path)) as writer:
        writer.write(mol)
    return {'compound_id': cid, 'name': name, 'score': score, 'score_direction': 'lower_is_better',
            'validation': {'passed': valid}, 'pose_path': str(path), 'pose_sha256': science.sha256(path),
            'protocol_fingerprint': 'protocol-a'}, {'compound_id': cid, 'name': name, 'canonical_smiles': canonical}


def test_counts_separate_generated_valid_failed_and_unfinished(tmp_path):
    valid, _ = pose(tmp_path, 'a', 'CCO', -5, True)
    invalid, _ = pose(tmp_path, 'b', 'CCC', -8, False)
    result = {'poses': [valid, invalid], 'failures': [{'compound_id': 'failed', 'code': 'timeout'}]}
    summary = summarize_result(result, [valid['compound_id'], invalid['compound_id'], 'failed', 'unfinished'])
    assert summary['ligands_with_poses'] == 2
    assert summary['valid_poses'] == 1 and summary['invalid_or_unverified_poses'] == 1
    assert summary['failed_ligands'] == 1 and summary['not_completed_ligands'] == 1
    assert summary['score_range'] == [-8, -5]
    assert summary['failure_codes'] == {'timeout': 1}
    with pytest.raises(ValueError, match='Unexpected compound'):
        summarize_result(result, ['different_input'])


def test_exports_are_hash_and_identity_verified_and_preserve_failed_checks(tmp_path):
    from rdkit import Chem
    valid, row_a = pose(tmp_path, 'a', 'CCO', -5, True)
    invalid, row_b = pose(tmp_path, 'b', 'CCC', -8, False)
    row_c = {'compound_id': 'failed', 'name': 'failed ligand', 'canonical_smiles': 'CCCC'}
    result = {'poses': [valid, invalid], 'failures': [{'compound_id': 'failed', 'reason': 'timed out'}]}
    exports = export_result(result, [row_a, row_b, row_c], tmp_path / 'exports')
    assert exports['csv_rows'] == 3 and exports['all_pose_records'] == 2 and exports['valid_pose_records'] == 1
    rows = list(csv.DictReader((tmp_path / 'exports/ligand_results.csv').open()))
    assert rows[-1]['status'] == 'failed' and rows[-1]['score_kcal_mol'] == ''
    recovered = list(Chem.SDMolSupplier(str(tmp_path / 'exports/all_docked_poses.sdf'), removeHs=False))
    assert [mol.GetProp('BIOCENTRIC_PHYSICAL_CHECK_PASSED') for mol in recovered] == ['true', 'false']
    assert [mol.GetProp('BIOCENTRIC_SOURCE_POSE_SHA256') for mol in recovered] == [valid['pose_sha256'], invalid['pose_sha256']]
    wrong = {**row_a, 'canonical_smiles': 'wrong'}
    with pytest.raises(ValueError, match='Pose identity mismatch'):
        export_result(result, [wrong, row_b, row_c], tmp_path / 'wrong_identity')
    Path(valid['pose_path']).write_text('changed')
    with pytest.raises(ValueError, match='Pose hash changed'):
        export_result(result, [row_a, row_b, row_c], tmp_path / 'wrong_hash')


def test_partial_run_recovers_only_matching_hash_verified_checkpoints(tmp_path):
    output = tmp_path / 'run'
    output.mkdir()
    (output / 'request.json').write_text(json.dumps({'fingerprint': 'current', 'config': {'seed': 42}}))
    artifact = output / 'actual_pose.sdf'
    artifact.write_text('canonical saved bytes')
    for name, fingerprint in [('one', 'current'), ('other', 'different')]:
        folder = output / name
        folder.mkdir()
        (folder / 'completed.json').write_text(json.dumps({'fingerprint': fingerprint,
            'poses': [{'compound_id': name}], 'artifacts': [science._artifact(artifact, 'pose')]}))
    recovered = partial_result(output)
    assert recovered['status'] == 'budget_stopped'
    assert recovered['poses'] == [{'compound_id': 'one'}]
    artifact.write_text('tampered')
    with pytest.raises(ValueError, match='Artifact integrity failed'):
        partial_result(output)
