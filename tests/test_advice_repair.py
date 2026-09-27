import pytest

from biocentric.service import Service


def rejected_project(service, *, used=2, local=False):
    pid=service.create_project({'name':'Repair test','target':'P00520','organism':'Mus musculus',
                                'local_only':local,'allow_model':not local,'max_model_calls':3})
    service.store.update(pid,state='blocked',data={'plan':{'strategy':'structure_based'},
        'model_attempted':True,'model_advice':{'status':'failed','reason':'Unread citation'},
        'evidence':[{'id':'read','text':'Already opened'},{'id':'unread','text':'Never opened'}]})
    for _ in range(used): service.store.event(pid,kind='model_call',status='running')
    service.store.event(pid,kind='tool_call',name='read_evidence',status='completed',inputs={'evidence_id':'read'})
    service.store.step(pid,'03','completed',outputs={'strategy':'structure_based'})
    return pid


def test_repair_uses_only_opened_evidence_and_cannot_repeat(tmp_path,monkeypatch):
    service=Service(tmp_path); pid=rejected_project(service)
    monkeypatch.setenv('OPENAI_API_KEY','synthetic')
    requests=[]
    def repair(spec,evidence,rejected,emit):
        requests.append(evidence)
        emit(kind='model_call',name='test',status='running',purpose='repair')
        return {'summary':'Corrected public advice','proposed_strategy':'further_experiments','claims':[]}
    monkeypatch.setattr('biocentric.agent.repair_advice',repair)
    service.repair_model_advice(pid)
    snapshot=service.snapshot(pid)
    assert requests==[[{'id':'read','text':'Already opened'}]]
    assert snapshot['project']['state']=='blocked'
    assert snapshot['data']['plan']['model_advice']['summary']=='Corrected public advice'
    assert sum(e.get('kind')=='model_call' and e.get('status')=='running' for e in snapshot['events'])==3
    with pytest.raises(ValueError,match='already attempted'): service.repair_model_advice(pid)
    assert len(requests)==1


@pytest.mark.parametrize('used,local,message',[(3,False,'budget is exhausted'),(0,True,'disabled')])
def test_repair_respects_budget_and_privacy(tmp_path,monkeypatch,used,local,message):
    service=Service(tmp_path);pid=rejected_project(service,used=used,local=local)
    monkeypatch.setenv('OPENAI_API_KEY','synthetic')
    with pytest.raises(ValueError,match=message): service.repair_model_advice(pid)
    assert 'model_repair' not in service.snapshot(pid)['data']


def test_failed_repair_keeps_reservation_and_does_not_expose_provider_body(tmp_path,monkeypatch):
    service=Service(tmp_path);pid=rejected_project(service)
    monkeypatch.setenv('OPENAI_API_KEY','synthetic')
    def fail(*args): raise RuntimeError('secret-provider-body')
    monkeypatch.setattr('biocentric.agent.repair_advice',fail)
    with pytest.raises(ValueError,match='RuntimeError'): service.repair_model_advice(pid)
    snapshot=service.snapshot(pid)
    assert snapshot['data']['model_repair']['status']=='failed'
    assert 'secret-provider-body' not in str(snapshot)
    with pytest.raises(ValueError,match='already attempted'): service.repair_model_advice(pid)
