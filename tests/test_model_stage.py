"""The stage result and visible plan must tell the same model-call story."""
import pytest

from biocentric.service import Service
from biocentric.worker import Pipeline


def pipeline(tmp_path, **policy):
    service = Service(tmp_path / 'store')
    pid = service.create_project({'name': 'Model visibility', 'target': 'P00520',
                                  'organism': 'Mus musculus', 'local_only': False,
                                  'allow_model': True, **policy})
    jid = service.start(pid, spawn=False)
    worker = Pipeline(service.store, pid, jid)
    worker.stage = '03'
    return service, pid, worker


@pytest.mark.parametrize('policy,code', [
    ({'local_only': True, 'allow_model': False}, 'local_only'),
    ({'allow_model': False}, 'disabled'),
    ({'max_model_calls': 0}, 'zero_budget'),
])
def test_skipped_model_has_an_explicit_persisted_reason(tmp_path, monkeypatch, policy, code):
    def forbidden(*args, **kwargs):
        pytest.fail('A disabled or unauthorized model must not be called')
    monkeypatch.setattr('biocentric.agent.advise', forbidden)
    service, pid, worker = pipeline(tmp_path, **policy)
    worker.stage_run('03', worker.strategy)
    snapshot = service.snapshot(pid)
    advice = snapshot['data']['model_advice']
    assert advice['status'] == 'not_requested'
    assert advice['skip_reason'] == code and advice['reason']
    assert snapshot['data']['plan']['model_advice'] == advice
    assert snapshot['steps'][2]['outputs']['model_advice'] == advice
    assert not snapshot['data'].get('model_attempted')
    assert not any(event['kind'] == 'model_call' for event in snapshot['events'])
    event = next(event for event in snapshot['events'] if event['kind'] == 'model_status')
    assert event['outputs'] == advice


def test_successful_advice_is_visible_in_both_plan_and_stage_output(tmp_path, monkeypatch):
    expected = {'summary': 'Evidence remains limited.', 'claims': [],
                'proposed_strategy': 'further_experiments', 'unresolved_questions': ['Mechanism?'],
                'limitations': ['Scientist review required'], 'model': 'gpt-4.1-mini', 'calls': 1}
    monkeypatch.setattr('biocentric.agent.advise', lambda *args, **kwargs: expected)
    service, pid, worker = pipeline(tmp_path)
    worker.stage_run('03', worker.strategy)
    snapshot = service.snapshot(pid)
    advice = {**expected, 'status': 'completed'}
    assert snapshot['data']['model_attempted'] is True
    assert snapshot['data']['model_advice'] == advice
    assert snapshot['data']['plan']['model_advice'] == advice
    assert snapshot['steps'][2]['outputs']['model_advice'] == advice


@pytest.mark.parametrize('failure', ['Model response incomplete: incomplete',
                                     'Model cited evidence it did not read'])
def test_model_failure_remains_visible_without_silent_retry(tmp_path, monkeypatch, failure):
    calls = []
    def fail(*args, **kwargs):
        calls.append(True)
        raise ValueError(failure)
    monkeypatch.setattr('biocentric.agent.advise', fail)
    service, pid, worker = pipeline(tmp_path)
    worker.stage_run('03', worker.strategy)
    snapshot = service.snapshot(pid)
    advice = snapshot['data']['model_advice']
    assert advice['status'] == 'failed' and advice['reason'] == failure
    assert snapshot['data']['plan']['model_advice'] == advice
    assert snapshot['steps'][2]['outputs']['model_advice'] == advice
    assert any(event['kind'] == 'model_call' and event['status'] == 'failed'
               for event in snapshot['events'])
    # Even if a new strategy evaluation is needed, it must retain the paid attempt.
    second_plan = worker.strategy()
    assert calls == [True]
    assert second_plan['model_advice'] == advice
    assert 'no automatic paid retry' in second_plan['model_retry_policy']
    snapshot = service.snapshot(pid)
    assert snapshot['data']['plan'] == second_plan
    assert snapshot['events'][-1]['kind'] == 'model_status'
    assert snapshot['events'][-1]['status'] == 'reused'
