"""Native process interruption preserves authenticated pairs and total accounting."""
import hashlib
import json
from pathlib import Path

import pytest

from biocentric import biocentric_model, engines, science


@pytest.fixture
def planned(tmp_path, monkeypatch):
    pairs = [{'record_index': index, 'compound_id': f'compound_{index}', 'protein_id': 'protein_A', 'smiles': smiles}
             for index, smiles in enumerate(['CCO', 'CCN'], 1)]
    models = [f'{name}@{epoch}' for name, epoch in biocentric_model.MODELS]
    out = tmp_path / 'run'
    plan = {'engine_id': 'biocentric_affinity', 'output_dir': str(out), 'status': 'ready', 'blockers': [],
            'commands': [['native-python', 'bridge.py', '--request', str(out / 'request.json')]],
            'expected_outputs': [str(out / 'predictions.json')],
            'generated_files': {'request.json': json.dumps({'pairs': pairs, 'models': models})},
            'validated_inputs': {}, 'expected_pairs': pairs, 'parameters': {'device': 'cuda', 'timeout_seconds': 10},
            'provenance': {'weights_sha256': biocentric_model.SUPPORTED_WEIGHTS_SHA256, 'code_sha256': 'c' * 64, 'models': models}}
    monkeypatch.setattr(engines, 'plan_engine_run', lambda *args: plan)
    return plan


def row_for(plan, index):
    return {**plan['expected_pairs'][index - 1], 'status': 'completed', 'score': 6.,
            'per_model': {model: 6. for model in plan['provenance']['models']}, 'models': plan['provenance']['models'],
            'checkpoint_sha256': plan['provenance']['weights_sha256'], 'code_sha256': plan['provenance']['code_sha256'],
            'device': 'cuda', 'preprocessing_status': 'passed', 'failure_reason': None}


def checkpoint(plan, index, row=None):
    out = Path(plan['output_dir']); (out / 'pairs').mkdir(exist_ok=True)
    row = row or row_for(plan, index)
    record = {'request_sha256': science.sha256(out / 'request.json'), 'row': row,
              'row_sha256': hashlib.sha256(json.dumps(row, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()}
    path = out / 'pairs' / f'{index:06d}.json'
    path.write_text(json.dumps(record))
    return path


def test_changed_real_progress_emits_once_and_native_resources_are_retained(planned, monkeypatch):
    events = []
    resources = {'peak_cuda_allocated_bytes': 1024, 'peak_rss_kib': 4096, 'duration_seconds': 1.5}
    def runner(argv, cwd, log, timeout, cancelled):
        assert cancelled() is False
        for index in [1, 2]:
            checkpoint(planned, index)
            (cwd / 'progress.json').write_text(json.dumps({'completed_pairs': index, 'total_pairs': 2,
                'latest': row_for(planned, index), 'resource_usage': resources}))
            assert cancelled() is False and cancelled() is False
        (cwd / 'predictions.json').write_text(json.dumps({'schema_version': 1, 'engine_id': 'biocentric_affinity',
            'predictions': [row_for(planned, 1), row_for(planned, 2)], 'resources': resources, 'runtime': {'torch': 'reference'}}))
        return {'returncode': 0}
    monkeypatch.setattr(science, 'run_process', runner)
    result = engines.run_engine('biocentric_affinity', {}, {}, planned['output_dir'], emit=events.append)
    progress = [event for event in events if event['kind'] == 'progress']
    assert [event['data']['completed_pairs'] for event in progress] == [1, 2]
    assert progress[0]['data']['latest']['compound_id'] == 'compound_1'
    assert progress[0]['data']['latest']['metrics'][0]['value'] == 6.
    assert progress[0]['data']['device'] == 'cuda'
    assert result['resource_usage'] == resources and result['native_runtime'] == {'torch': 'reference'}
    assert result['accounting'] == {'requested': 2, 'completed': 2, 'failed': 0}


@pytest.mark.parametrize('code,status', [('timeout', 'partial'), ('cancelled', 'cancelled'), ('tool_failed', 'partial')])
def test_interruption_keeps_verified_pair_and_accounts_for_missing_pair(planned, monkeypatch, code, status):
    original = []
    def runner(argv, cwd, log, timeout, cancelled):
        path = checkpoint(planned, 1); original.append(path.read_bytes())
        raise science.ScienceError(code, 'Native calculation stopped')
    monkeypatch.setattr(science, 'run_process', runner)
    result = engines.run_engine('biocentric_affinity', {}, {}, planned['output_dir'])
    assert result['status'] == status
    assert result['accounting'] == {'requested': 2, 'completed': 1, 'failed': 1}
    assert result['results'][0]['checkpoint_recovered'] and result['results'][0]['score'] == 6.
    assert result['results'][1]['score'] is None and result['results'][1]['metrics'] == []
    assert result['results'][1]['failure_code'] == code
    assert (Path(planned['output_dir']) / 'pairs/000001.json').read_bytes() == original[0]
    assert (Path(planned['output_dir']) / 'partial-predictions.json').exists()


def test_startup_failure_accounts_for_every_requested_pair(planned, monkeypatch):
    def runner(*args):
        raise FileNotFoundError('Missing native dependency')
    monkeypatch.setattr(science, 'run_process', runner)
    result = engines.run_engine('biocentric_affinity', {}, {}, planned['output_dir'])
    assert result['status'] == 'failed'
    assert result['accounting'] == {'requested': 2, 'completed': 0, 'failed': 2}
    assert all(row['score'] is None and 'Missing native dependency' in row['failure_reason'] for row in result['results'])


@pytest.mark.parametrize('change', ['row_hash', 'request_hash', 'identity', 'code', 'device', 'nonfinite', 'outside'])
def test_untrusted_checkpoint_is_retained_but_never_salvaged(planned, monkeypatch, tmp_path, change):
    def runner(argv, cwd, log, timeout, cancelled):
        row = row_for(planned, 1)
        if change == 'identity': row['protein_id'] = 'substitution'
        if change == 'code': row['code_sha256'] = 'other'
        if change == 'device': row['device'] = 'cpu'
        path = checkpoint(planned, 1, row)
        record = json.loads(path.read_text())
        if change == 'row_hash': record['row']['score'] = 9.
        if change == 'request_hash': record['request_sha256'] = 'other'
        if change == 'nonfinite': record['row']['score'] = float('nan')
        path.write_text(json.dumps(record))
        if change == 'outside':
            outside = tmp_path / 'outside.json'; outside.write_text(path.read_text())
            path.unlink(); path.symlink_to(outside)
        raise science.ScienceError('timeout', 'Stopped')
    monkeypatch.setattr(science, 'run_process', runner)
    result = engines.run_engine('biocentric_affinity', {}, {}, planned['output_dir'])
    assert result['accounting'] == {'requested': 2, 'completed': 0, 'failed': 2}
    assert 'checkpoint not accepted' in result['results'][0]['failure_reason']


def test_input_drift_invalidates_even_well_formed_saved_pairs(planned, monkeypatch, tmp_path):
    input_file = tmp_path / 'reviewed.pdb'; input_file.write_text('original')
    planned['validated_inputs']['receptor'] = {'path': str(input_file), 'sha256': science.sha256(input_file)}
    def runner(argv, cwd, log, timeout, cancelled):
        checkpoint(planned, 1)
        input_file.write_text('changed')
        raise science.ScienceError('tool_failed', 'Stopped after input drift')
    monkeypatch.setattr(science, 'run_process', runner)
    result = engines.run_engine('biocentric_affinity', {}, {}, planned['output_dir'])
    assert result['accounting']['completed'] == 0
    assert all('reviewed input' in row['failure_reason'] for row in result['results'])


def test_oversize_or_substituted_progress_never_emits_fake_completion(planned, monkeypatch):
    events = []
    def runner(argv, cwd, log, timeout, cancelled):
        (cwd / 'progress.json').write_bytes(b' ' * (512 * 1024 + 1))
        assert cancelled() is False
        (cwd / 'progress.json').write_text(json.dumps({'completed_pairs': 99, 'total_pairs': 99, 'latest': {}}))
        assert cancelled() is False
        raise science.ScienceError('timeout', 'Stopped')
    monkeypatch.setattr(science, 'run_process', runner)
    result = engines.run_engine('biocentric_affinity', {}, {}, planned['output_dir'], emit=events.append)
    assert not any(event['kind'] == 'progress' for event in events)
    assert len(result['progress_warnings']) == 2 and result['accounting']['failed'] == 2


@pytest.mark.parametrize('change', ['modified', 'missing'])
def test_prestart_input_drift_accounts_for_every_pair_without_starting_native(planned, monkeypatch, tmp_path, change):
    source = tmp_path / 'reviewed.pdb'
    source.write_text('reviewed input')
    planned['validated_inputs']['receptor'] = {'path': str(source), 'sha256': science.sha256(source)}
    if change == 'modified':
        source.write_text('changed input')
    else:
        source.unlink()
    monkeypatch.setattr(science, 'run_process', lambda *args: pytest.fail('Changed input must not reach native execution'))
    events = []
    result = engines.run_engine('biocentric_affinity', {}, {}, planned['output_dir'], emit=events.append)
    assert result['status'] == 'failed' and result['executed'] is False
    assert result['accounting'] == {'requested': 2, 'completed': 0, 'failed': 2}
    assert all(row['score'] is None and row['failure_reason'] for row in result['results'])
    assert not any(event['kind'] == 'tool_started' for event in events)
    saved = json.loads((Path(planned['output_dir']) / 'result.json').read_text())
    assert saved['results'] == result['results']


@pytest.mark.parametrize('field,changed', [('code_sha256', 'other-code'), ('device', 'cpu'), ('models', ['other-model']),
                                          ('code_sha256', None), ('device', None), ('models', None)])
@pytest.mark.parametrize('status', ['completed', 'failed'])
def test_final_rows_require_reviewed_model_code_and_device(planned, monkeypatch, field, changed, status):
    original = []
    def runner(argv, cwd, log, timeout, cancelled):
        rows = [row_for(planned, 1), row_for(planned, 2)]
        if status == 'failed':
            rows[0].update(status='failed', score=None, per_model={}, failure_reason='Actual preparation failure')
        if changed is None:
            rows[0].pop(field)
        else:
            rows[0][field] = changed
        content = json.dumps({'schema_version': 1, 'engine_id': 'biocentric_affinity', 'predictions': rows})
        (cwd / 'predictions.json').write_text(content)
        original.append(content)
        return {'returncode': 0}
    monkeypatch.setattr(science, 'run_process', runner)
    result = engines.run_engine('biocentric_affinity', {}, {}, planned['output_dir'])
    assert result['status'] == 'failed' and result['accounting']['completed'] == 0
    assert all('reviewed plan' in row['failure_reason'] and row['score'] is None for row in result['results'])
    assert (Path(planned['output_dir']) / 'predictions.json').read_text() == original[0]


@pytest.mark.parametrize('change', ['oversize', 'outside', 'malformed_rows', 'malformed_identity'])
def test_final_json_is_bounded_contained_and_typed_before_predictions_are_accepted(planned, monkeypatch, tmp_path, change):
    def runner(argv, cwd, log, timeout, cancelled):
        path = cwd / 'predictions.json'
        if change == 'oversize':
            path.write_bytes(b' ' * (16 * 1024 * 1024 + 1))
        elif change == 'outside':
            outside = tmp_path / 'outside.json'
            outside.write_text('{}')
            path.symlink_to(outside)
        else:
            rows = [row_for(planned, 1), row_for(planned, 2)]
            if change == 'malformed_rows':
                rows[0] = None
            else:
                rows[0]['compound_id'] = ['not-an-identifier']
            path.write_text(json.dumps({'schema_version': 1, 'engine_id': 'biocentric_affinity', 'predictions': rows}))
        return {'returncode': 0}
    monkeypatch.setattr(science, 'run_process', runner)
    monkeypatch.setattr(biocentric_model, 'parse_output', lambda *args: pytest.fail('The engine must parse its bounded in-root bytes'))
    result = engines.run_engine('biocentric_affinity', {}, {}, planned['output_dir'])
    assert result['status'] == 'failed'
    assert result['accounting'] == {'requested': 2, 'completed': 0, 'failed': 2}
    assert all(row['score'] is None for row in result['results'])
