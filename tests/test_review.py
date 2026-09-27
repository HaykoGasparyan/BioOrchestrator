"""Regression checks for identity, resumed work and output provenance."""
import json
import os
import sys
import time
from pathlib import Path
from types import SimpleNamespace

import pytest

from biocentric.service import Service
from biocentric.worker import Pipeline, Cancelled, execute


def project(service, **extra):
    return service.create_project({'name':'Review regression','target':'P00520','organism':'Mus musculus',**extra})


def identify(service,pid,**extra):
    service.decide(pid,'identity',{'accession':'P00520','name':'ABL1','organism':'Mus musculus',**extra},'Reviewed target identity')


def test_imported_identity_cannot_drop_requested_mutation(tmp_path):
    service=Service(tmp_path);pid=project(service,mutations=['T315I'])
    with pytest.raises(ValueError,match='mutations'):
        identify(service,pid)
    identify(service,pid,mutations=['T315I'])
    assert service.snapshot(pid)['decisions'][-1]['value']['mutations']==['T315I']


def test_imported_identity_cannot_substitute_supplied_sequence(tmp_path):
    service=Service(tmp_path);sequence='ACDEFGHIKLMNPQRSTVWY'*3
    pid=project(service,target='>custom\n'+sequence)
    with pytest.raises(ValueError,match='match.*sequence'):
        identify(service,pid,sequence='M'+sequence)
    identify(service,pid,sequence=sequence.lower())
    assert service.snapshot(pid)['decisions'][-1]['value']['sequence']==sequence


def test_invalid_imported_amino_acids_are_rejected(tmp_path):
    service=Service(tmp_path);pid=project(service)
    with pytest.raises(ValueError,match='amino-acid'):
        identify(service,pid,sequence='ACDE12345')


def test_start_of_live_job_preserves_running_state(tmp_path):
    service=Service(tmp_path);pid=project(service);jid=service.start(pid,spawn=False)
    assert service.store.claim(jid)
    service.store.update(pid,state='running')
    assert service.start(pid,spawn=False)==jid
    assert service.snapshot(pid)['project']['state']=='running'


def test_worker_cannot_claim_job_for_another_project(tmp_path):
    service=Service(tmp_path);pid=project(service);other=project(service)
    jid=service.start(pid,spawn=False)
    assert execute(tmp_path,other,jid) is False
    assert service.store.jobs(pid)[0]['status']=='queued'
    assert service.snapshot(other)['project']['state']=='draft'


def test_stale_owner_cannot_heartbeat_or_finish_successor(tmp_path):
    service=Service(tmp_path);pid=project(service);jid=service.store.enqueue(pid)
    assert service.store.claim(jid,owner_token='old')
    with service.store.connect() as c:
        c.execute('UPDATE jobs SET lease_until=? WHERE id=?',(time.time()-1,jid))
    service.store.enqueue(pid)
    assert service.store.claim(jid,owner_token='new')
    assert not service.store.heartbeat(jid,'old')
    assert not service.store.finish(jid,'completed',owner_token='old',project_id=pid,project_data={'compromised':True})
    assert 'compromised' not in service.store.project(pid)['data']
    assert service.store.jobs(pid)[0]['status']=='running'
    old=Pipeline(service.store,pid,jid,owner_token='old')
    with pytest.raises(Cancelled,match='owns'):
        old.check()


def test_expired_lease_does_not_duplicate_a_live_local_worker(tmp_path):
    service=Service(tmp_path);pid=project(service);jid=service.store.enqueue(pid)
    service.store.claim(jid,owner_token='live',worker_pid=os.getpid())
    with service.store.connect() as c:
        c.execute('UPDATE jobs SET lease_until=? WHERE id=?',(time.time()-1,jid))
    service.store.enqueue(pid)
    assert service.store.jobs(pid)[0]['status']=='running'
    assert not service.store.claim(jid,owner_token='duplicate')


def test_budget_block_persists_gate_even_before_stage_entry(tmp_path):
    service=Service(tmp_path);pid=project(service,compute_budget_seconds=10)
    service.store.update(pid,data={'compute_seconds':11})
    jid=service.start(pid,spawn=False);execute(tmp_path,pid,jid)
    snap=service.snapshot(pid)
    assert snap['project']['state']=='blocked'
    assert snap['data']['gate']['kind']=='budget'


def test_resume_rejects_modified_checkpoint_result(tmp_path):
    service=Service(tmp_path);pid=project(service);identify(service,pid)
    jid=service.start(pid,spawn=False);execute(tmp_path,pid,jid)
    snapshot=service.snapshot(pid)
    checkpoint=next(a for a in snapshot['artifacts'] if a['name']=='01-result.json')
    Path(checkpoint['path']).write_text('{}')
    service.start(pid,spawn=False);execute(tmp_path,pid,jid)
    snapshot=service.snapshot(pid)
    assert snapshot['data']['gate']['kind']=='integrity'
    assert snapshot['steps'][0]['status']=='blocked'
    assert 'integrity' in snapshot['steps'][0]['reason'].lower()


def test_next_stage_uses_immutable_prepared_library(tmp_path):
    service=Service(tmp_path);pid=project(service);jid=service.store.enqueue(pid)
    pipeline=Pipeline(service.store,pid,jid)
    working=pipeline.work/'prepared.sdf';working.write_text('original prepared chemistry')
    result={'prepared_sdf':str(working),'artifacts':[{'path':str(working)}]}
    registered=pipeline.register(result);pipeline.bind_prepared(result,registered)
    working.write_text('modified scratch data')
    prepared=pipeline.prepared_input(result)
    assert prepared!=working and prepared.read_text()=='original prepared chemistry'
    prepared.write_text('modified immutable data')
    with pytest.raises(ValueError,match='integrity'):
        pipeline.prepared_input(result)


def test_handoff_exports_selected_poses_with_provenance(tmp_path):
    from rdkit import Chem
    service=Service(tmp_path);pid=project(service);jid=service.store.enqueue(pid)
    molecule=Chem.MolFromSmiles('CCO')
    content=(Chem.MolToMolBlock(molecule)+'\n$$$$\n').encode()
    pose=service.store.artifact(pid,'ethanol.sdf',content,'scientific_output',jid)
    candidate={'compound_id':'ethanol','selected':True,'rank':1,'job_id':jid,'artifact_ids':[pose['id']],'pose_sha256':pose['sha256']}
    service.store.update(pid,data={'candidates':[candidate]})
    pipeline=Pipeline(service.store,pid,jid)
    result=pipeline.handoff()
    exported=service.store.read_artifact(pid,result['shortlist'])
    assert exported==content
    meta=service.store.project(pid)['data']['shortlist']['selection_manifest']
    selection=json.loads(service.store.read_artifact(pid,meta))
    assert selection['selection'][0]['pose_artifact_id']==pose['id']
    assert selection['sdf_sha256']==result['shortlist']['sha256']


def test_existing_immutable_artifact_cannot_be_silently_reregistered(tmp_path):
    service=Service(tmp_path);pid=project(service)
    artifact=service.store.artifact(pid,'example.sdf',b'original','scientific_output')
    Path(artifact['path']).write_bytes(b'tampered')
    with pytest.raises(ValueError,match='integrity'):
        service.store.artifact(pid,'example.sdf',b'original','scientific_output')


@pytest.mark.parametrize('local_only,allow_model,expected_key',[
    (True,False,None), (False,False,None), (False,True,'synthetic-key-for-tests'),
])
def test_worker_inherits_api_key_only_for_explicit_model_access(tmp_path,monkeypatch,local_only,allow_model,expected_key):
    monkeypatch.setenv('OPENAI_API_KEY','synthetic-key-for-tests')
    calls=[]
    monkeypatch.setattr('biocentric.service.subprocess.Popen',lambda argv,**kwargs:calls.append(kwargs))
    service=Service(tmp_path);pid=project(service,local_only=local_only,allow_model=allow_model)
    service.start(pid)
    assert len(calls)==1
    assert calls[0]['env'].get('OPENAI_API_KEY')==expected_key
    assert os.environ['OPENAI_API_KEY']=='synthetic-key-for-tests'


def test_vina_health_probe_receives_no_application_credentials(monkeypatch):
    from biocentric import science
    monkeypatch.setenv('OPENAI_API_KEY','synthetic-key-for-tests')
    monkeypatch.setenv('ANTHROPIC_API_KEY','synthetic-other-key')
    monkeypatch.setenv('AWS_SECRET_ACCESS_KEY','synthetic-cloud-key')
    monkeypatch.setenv('BIOCENTRIC_TEST_ENV','retained')
    monkeypatch.setattr(science,'_vina_executable',lambda:'/fake/vina')
    environments=[]
    def probe(argv,**kwargs):
        assert argv==['/fake/vina','--version']
        environments.append(kwargs['env'])
        return SimpleNamespace(stdout='AutoDock Vina v1.2.7\n')
    monkeypatch.setattr(science.subprocess,'run',probe)
    science.health()
    assert len(environments)==1
    assert not {'OPENAI_API_KEY','ANTHROPIC_API_KEY','AWS_SECRET_ACCESS_KEY'} & environments[0].keys()
    assert environments[0]['BIOCENTRIC_TEST_ENV']=='retained'


def test_docking_version_probe_receives_no_api_key(tmp_path,monkeypatch):
    pytest.importorskip('rdkit')
    pytest.importorskip('meeko')
    from biocentric import science
    monkeypatch.setenv('OPENAI_API_KEY','synthetic-key-for-tests')
    monkeypatch.setattr(science,'_vina_executable',lambda:'/fake/vina')
    class ProbeReached(Exception): pass
    def probe(argv,**kwargs):
        assert argv==['/fake/vina','--version']
        assert 'OPENAI_API_KEY' not in kwargs['env']
        raise ProbeReached
    monkeypatch.setattr(science.subprocess,'run',probe)
    original=Path(__file__).resolve().parents[1]/'examples'/'1iep'/'originals'
    with pytest.raises(ProbeReached):
        science.dock_library(original/'1iep_receptor.pdbqt',original/'1iep_ligand.sdf',original/'1iep_ligand.sdf',tmp_path/'docking',{'require_posebusters':False,'run_plip':False})


def test_scientific_subprocess_does_not_inherit_api_key(tmp_path,monkeypatch):
    from biocentric.science import run_process
    monkeypatch.setenv('OPENAI_API_KEY','synthetic-key-for-tests')
    log=tmp_path/'process.log'
    run_process([sys.executable,'-c',"import os; assert 'OPENAI_API_KEY' not in os.environ; print('credential absent')"],tmp_path,log,5)
    assert log.read_text().strip()=='credential absent'
