import hashlib
import json
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from unittest.mock import patch
import pytest
from pydantic import ValidationError
from biocentric.models import ProjectSpec, normalize_measurement
from biocentric.service import Service
from biocentric.store import Store
from biocentric.worker import execute, Pipeline, Gate


def project(service,**extra):
    return service.create_project({'name':'Test','target':'P00520','organism':'Mus musculus',**extra})


def identity(service,pid):
    service.decide(pid,'identity',{'accession':'P00520','name':'ABL1','organism':'Mus musculus'},'Reviewed local construct')


def test_schema_limits_and_unknowns():
    with pytest.raises(ValidationError): ProjectSpec(name='x',target='p',organism='m',arbitrary_command='rm')
    with pytest.raises(ValidationError): ProjectSpec(name='x',target='p',organism='m',cpu=999)
    with pytest.raises(ValidationError): ProjectSpec(name='x',target='p',organism='m',seed=0)


def test_endpoint_units_relation_and_original_preserved():
    record=normalize_measurement({'compound_id':'c1','endpoint':'IC50','value':'20','units':'nM','relation':'<','assay_id':'a1','organism':'human'})
    assert record.endpoint=='IC50' and record.value==20 and record.relation=='<'
    assert record.value_molar==pytest.approx(2e-8) and record.original['value']=='20'
    with pytest.raises(ValueError): normalize_measurement({**record.original,'units':'kcal/mol'})
    with pytest.raises(ValueError): normalize_measurement({**record.original,'value':'nan'})


def test_worker_identity_gate_has_no_external_calls(tmp_path):
    service=Service(tmp_path);pid=project(service);jid=service.start(pid,spawn=False)
    with patch('biocentric.retrieval.httpx.Client',side_effect=AssertionError('external call')): execute(tmp_path,pid,jid)
    snap=service.snapshot(pid)
    assert snap['project']['state']=='blocked' and snap['data']['gate']['kind']=='identity'
    assert snap['steps'][0]['status']=='blocked'
    assert all(s['status']=='pending' for s in snap['steps'][1:])


def test_fixture_cannot_run_or_approve(tmp_path):
    service=Service(tmp_path);pid=service.create_demo()
    with pytest.raises(ValueError): service.start(pid,spawn=False)
    with pytest.raises(ValueError): service.decide(pid,'exploratory',{'approved':True},'test')
    exports=service.export(pid)
    assert 'OFFLINE FIXTURE' in Path(exports['markdown']).read_text()


def test_identity_cannot_substitute_species_accession_or_isoform(tmp_path):
    service=Service(tmp_path);pid=project(service,isoform='2')
    for value in [dict(accession='P00519',name='ABL1',organism='Mus musculus',isoform='2'),dict(accession='P00520',name='ABL1',organism='Homo sapiens',isoform='2'),dict(accession='P00520',name='ABL1',organism='Mus musculus',isoform='canonical')]:
        with pytest.raises(ValueError): service.decide(pid,'identity',value,'reason')


def test_atomic_claim_and_idempotent_queue(tmp_path):
    service=Service(tmp_path);pid=project(service)
    with ThreadPoolExecutor(max_workers=8) as pool: ids=list(pool.map(lambda _:service.store.enqueue(pid),range(20)))
    assert len(set(ids))==1
    with ThreadPoolExecutor(max_workers=8) as pool: claimed=list(pool.map(lambda _:service.store.claim(ids[0]),range(20)))
    assert claimed.count(True)==1
    assert service.store.enqueue(pid)==ids[0]
    assert service.store.jobs(pid)[0]['status']=='running'


def test_resume_reuses_completed_stage_and_same_job(tmp_path):
    service=Service(tmp_path);pid=project(service);identity(service,pid)
    jid=service.start(pid,spawn=False);execute(tmp_path,pid,jid)
    snap=service.snapshot(pid)
    assert snap['steps'][0]['status']=='completed' and snap['steps'][3]['status']=='blocked'
    first=[a['id'] for a in snap['artifacts'] if a['name']=='01-result.json']
    assert service.start(pid,spawn=False)==jid
    execute(tmp_path,pid,jid)
    snap=service.snapshot(pid)
    assert [a['id'] for a in snap['artifacts'] if a['name']=='01-result.json']==first
    assert any(e['kind']=='checkpoint' and e['stage']=='01' for e in snap['events'])
    assert snap['jobs'][0]['attempts']==2


def test_stale_lease_recoverable_but_live_lease_not_duplicated(tmp_path):
    service=Service(tmp_path);pid=project(service);jid=service.store.enqueue(pid);service.store.claim(jid)
    assert service.store.enqueue(pid)==jid
    with service.store.connect() as c: c.execute('UPDATE jobs SET lease_until=? WHERE id=?',(time.time()-1,jid))
    assert service.store.enqueue(pid)==jid and service.store.claim(jid)


def test_cancel_queued_job_prevents_claim(tmp_path):
    service=Service(tmp_path);pid=project(service);jid=service.start(pid,spawn=False)
    service.cancel(pid)
    assert not service.store.claim(jid)
    assert service.store.jobs(pid)[0]['status']=='cancelled'


def test_unsupported_workflow_blocks_before_preparation(tmp_path):
    service=Service(tmp_path);pid=project(service,workflow_kind='covalent');identity(service,pid)
    jid=service.start(pid,spawn=False);execute(tmp_path,pid,jid)
    snap=service.snapshot(pid)
    assert snap['data']['gate']['kind']=='unsupported'
    assert snap['steps'][2]['status']=='blocked' and snap['steps'][4]['status']=='pending'


def test_exploratory_is_not_automatic_pass(tmp_path):
    service=Service(tmp_path);pid=project(service);jid=service.start(pid,spawn=False)
    service.store.finish(jid,'blocked')
    pipeline=Pipeline(service.store,pid,jid)
    with pytest.raises(Gate): pipeline.benchmark()
    service.decide(pid,'exploratory',{'approved':True},'Small exploratory software test')
    result=pipeline.benchmark()
    assert result['status']=='exploratory' and not result['metrics']


def test_protocol_requires_concrete_proposal_and_input_review(tmp_path):
    service=Service(tmp_path);pid=project(service)
    with pytest.raises(ValueError): service.decide(pid,'protocol',{'approved':True,'receptor_reviewed':True},'reason')
    service.store.update(pid,data={'protocol':{'engine':'Vina','input_hashes':{'r':'123'}}})
    with pytest.raises(ValueError): service.decide(pid,'protocol',{'approved':True},'reason')
    service.decide(pid,'protocol',{'approved':True,'receptor_reviewed':True,'ph_assumption':7.4},'reason')
    assert service.snapshot(pid)['decisions'][-1]['value']['protocol']['input_hashes']['r']=='123'


def test_artifact_integrity_paths_and_immutable_original(tmp_path):
    service=Service(tmp_path);pid=project(service)
    a=service.upload(pid,'../../example.smi',b'CCO ethanol\n','library')
    assert Path(a['path']).is_relative_to(service.root) and a['name']=='example.smi'
    assert service.store.read_artifact(pid,a)==b'CCO ethanol\n'
    Path(a['path']).write_text('tampered')
    with pytest.raises(ValueError): service.store.read_artifact(pid,a)
    with pytest.raises(ValueError): service.store.directory('../outside')


def test_reports_link_jobs_citations_and_hashes_without_secrets(tmp_path):
    service=Service(tmp_path);pid=project(service)
    service.store.event(pid,kind='tool',inputs={'api_key':'private','authorization':'Bearer private'},outputs={'scores':[]})
    a=service.store.artifact(pid,'pose.sdf',b'original result','scientific_output','completed-job')
    service.store.update(pid,data={'evidence':[{'id':'e1','source':'source','source_id':'id','url':'https://example.org/record','retrieved_at':'2026-09-26','location':'abstract','text':'Untrusted observation','evidence_type':'association'}]})
    result=service.export(pid)
    assert '[e1]' in Path(result['markdown']).read_text()
    manifest=json.loads(Path(result['manifest']).read_text())
    assert any(r['sha256']==hashlib.sha256(b'original result').hexdigest() for r in manifest['artifacts'])
    assert 'private' not in Path(result['events']).read_text()


def test_budget_cannot_be_bypassed_by_resume(tmp_path):
    service=Service(tmp_path);pid=project(service,compute_budget_seconds=10);identity(service,pid)
    service.store.update(pid,data={'compute_seconds':11})
    jid=service.start(pid,spawn=False);execute(tmp_path,pid,jid)
    assert service.snapshot(pid)['project']['state']=='blocked'
    assert all(s['status']=='pending' for s in service.snapshot(pid)['steps'])
