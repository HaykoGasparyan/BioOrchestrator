from biocentric.explain import model_activity, stage_explanation
import json


def test_configured_key_cannot_override_project_privacy_or_budget():
    for policy in [{'local_only': True, 'allow_model': False},
                   {'local_only': False, 'allow_model': False},
                   {'local_only': False, 'allow_model': True, 'max_model_calls': 0}]:
        result = model_activity({'project': {'spec': policy}}, key_available=True)
        assert result['state'] == 'disabled'
        assert result['calls'] == 0
        assert result['input_tokens'] is None


def test_public_response_and_usage_survive_missing_key_and_use_actual_requests():
    response = {'response_id': 'r1', 'usage': {'input_tokens': 100, 'output_tokens': 20}}
    snapshot = {'project': {'spec': {'local_only': False, 'allow_model': True}},
                'data': {'model_advice': {'summary': 'Recorded proposal'}},
                'events': [{'kind': 'model_call', 'status': 'running'},
                           {'kind': 'model_call', 'status': 'completed', 'outputs': response},
                           {'kind': 'model_call', 'status': 'completed', 'outputs': response}]}
    result = model_activity(snapshot, key_available=False)
    assert result['state'] == 'completed'
    assert result['calls'] == 1  # Started/completed are two events for one request.
    assert result['input_tokens'] == 100  # Duplicate response events do not double usage.
    assert result['output_tokens'] == 20


def test_failed_advice_is_distinct_from_missing_configuration():
    snapshot = {'project': {'spec': {'local_only': False, 'allow_model': True}},
                'data': {'model_attempted': True, 'model_advice': {'status': 'failed', 'reason': 'Call budget exhausted'}}}
    result = model_activity(snapshot, key_available=False)
    assert result['state'] == 'failed'
    assert result['reason'] == 'Call budget exhausted'


def test_summary_preserves_blocked_reason_and_uses_recorded_strategy_rationale():
    step = {'id': '03', 'status': 'blocked', 'reason': 'Unsupported workflow'}
    assert stage_explanation({}, step)['reason'] == 'Unsupported workflow'
    step.update(status='completed', outputs={'strategy': 'structure_based', 'rationale': 'Recorded basis for choice'})
    explanation = stage_explanation({}, step)
    assert explanation['reason'] == 'Recorded basis for choice'
    assert 'Workflow rules' in explanation['actor']
    assert 'No LLM' in stage_explanation({}, {'id': '08'})['actor']


def test_rejected_public_output_is_visible_but_never_promoted_to_advice():
    public={'summary':'Unaccepted model text','claims':[{'text':'Uncited statement','category':'retrieved_observation','evidence_ids':[]}],
            'proposed_strategy':'further_experiments','unresolved_questions':[],'limitations':[]}
    snapshot={'project':{'spec':{'local_only':False,'allow_model':True}},
              'data':{'model_advice':{'status':'failed','reason':'Missing citations'}},
              'events':[{'kind':'model_call','status':'completed','outputs':{'output':[
                  {'type':'message','content':[{'type':'output_text','text':json.dumps(public)}]}]}}]}
    activity=model_activity(snapshot,key_available=True)
    assert activity['state']=='failed'
    assert 'summary' not in activity['advice']
    assert activity['rejected_advice']['summary']=='Unaccepted model text'
