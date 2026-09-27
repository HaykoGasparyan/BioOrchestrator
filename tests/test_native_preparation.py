"""Native preparation integrity checks; no docking, model, GPU or network calls.

Synthetic passing-pose metadata is used only to exercise guards. It is not a
physical-validation or scientific-performance result.
"""
from copy import deepcopy
import io
import json
from pathlib import Path

import pytest
from rdkit import Chem

from biocentric import science
from biocentric.engine_workflow import NativeManifest, plan_native_manifest
from biocentric.engines import GninaParameters, UniDockParameters
from biocentric.native_preparation import auto_prepare_native_manifest
from biocentric.research_workflow import ResearchWorkflow
from biocentric.service import Service
from biocentric.store import digest, dump
from biocentric.testsets import build_testset


@pytest.fixture
def context(tmp_path):
    service=Service(tmp_path/'store')
    pid=service.create_research_project('Native preparation integrity test fixture',use_public_example=True)
    jid=service.start_research(pid,spawn=False)
    frozen=build_testset(service.snapshot(pid),organism='Mus musculus')
    service.store.research_update(pid,testset=frozen,tool_plan={'docking_engine':'vina'})
    protocol=ResearchWorkflow(service.store,pid,jid).protocol()
    artifacts={row['kind']:row for row in service.store.artifacts(pid)}
    library=list(Chem.SDMolSupplier(artifacts['library']['path'],removeHs=False))
    mol=library[0]
    smiles=science._canonical(mol)
    buffer=io.StringIO()
    with Chem.SDWriter(buffer) as writer:
        writer.write(mol)
    pose_artifact=service.store.artifact(pid,'synthetic-validated-pose.sdf',buffer.getvalue().encode(),'scientific_output')
    validation_report={'status':'completed','passed':True,'tool':'synthetic_test_fixture_only','mode':'dock',
                       'receptor_checks_run':True,'failed_checks':[], 'scope':'No physical calculation was performed by this fixture.'}
    validation_artifact=service.store.artifact(pid,'synthetic-validation.json',json.dumps(validation_report).encode(),'scientific_output')
    pose={'compound_id':'cmp_'+digest(smiles.encode())[:16],'canonical_smiles':smiles,
          'pose_sha256':pose_artifact['sha256'],'pose_number':1,'receptor_sha256':artifacts['structure']['sha256'],
          'validation':{**validation_report,'artifact':{'sha256':validation_artifact['sha256']}},
          'fixture_scope':'synthetic metadata only; no physical validation'}
    service.store.update(pid,data={'screening':{'poses':[pose]}})
    return service,pid,jid,protocol,pose


def prepare(context,engine='biocentric_affinity',phase='affinity'):
    service,pid,jid,protocol,_=context
    artifact=auto_prepare_native_manifest(service.store,pid,engine,phase,protocol,job_id=jid)
    return artifact,NativeManifest.model_validate_json(service.store.read_artifact(pid,artifact))


def test_affinity_binds_reviewed_identity_receptor_and_exact_compound(context):
    service,pid,_,protocol,pose=context
    artifact,manifest=prepare(context)
    assert manifest.testset_sha256==protocol['testset_sha256']
    assert manifest.protocol_sha256==digest(dump(protocol).encode())
    identity=json.loads(service.store.read_artifact(pid,{'id':manifest.protein_identity_artifact}))
    assert identity['reviewed_identity']['accession']=='P00520'
    assert identity['reviewed_identity']['organism']=='Mus musculus'
    assert identity['receptor_sha256']==protocol['input_hashes']['structure']
    assert manifest.review['protein_id'].startswith('P00520:Mus musculus:chains=')
    assert manifest.review['protein_id'].endswith(protocol['input_hashes']['structure'])
    native=Chem.ForwardSDMolSupplier(io.BytesIO(service.store.read_artifact(pid,{'id':manifest.input_artifacts['ligands']})),removeHs=False)
    molecules=list(native)
    assert len(molecules)==1 and molecules[0].GetProp('_Name')==pose['compound_id']
    assert science._canonical(molecules[0])==pose['canonical_smiles']
    assert len(manifest.exclusions)==2
    assert all('not imputed' in row.reason for row in manifest.exclusions)


@pytest.mark.parametrize('mutation,pattern',[
    ('cid','compound ID'),('receptor','reviewed receptor'),('chemistry','pose identity'),('unregistered','not registered'),('invalid','passing structures')])
def test_affinity_rejects_identity_receptor_or_unverified_pose_substitutions(context,mutation,pattern):
    service,pid,_,_,original=context
    pose=deepcopy(original)
    if mutation=='cid': pose['compound_id']='cmp_wrong'
    elif mutation=='receptor': pose['receptor_sha256']='0'*64
    elif mutation=='chemistry':
        pose['canonical_smiles']='CCO'
        pose['compound_id']='cmp_'+digest(b'CCO')[:16]
    elif mutation=='unregistered': pose['pose_sha256']='f'*64
    else: pose['validation']['passed']=False
    service.store.update(pid,data={'screening':{'poses':[pose]}})
    with pytest.raises(ValueError,match=pattern):
        prepare(context)
    assert not any(row['kind']=='native_engine_affinity' for row in service.store.artifacts(pid))


def test_affinity_structure_artifact_tampering_is_detected_before_model_inputs(context):
    service,pid,_,_,pose=context
    artifact=next(row for row in service.store.artifacts(pid) if row['sha256']==pose['pose_sha256'])
    Path(artifact['path']).write_bytes(b'changed structure')
    with pytest.raises(ValueError,match='integrity'):
        prepare(context)


def test_missing_reviewed_accession_and_species_block_model_preparation(context):
    service,pid,_,_,_=context
    # The nonempty explicit identity must not silently fall back to another
    # historical target when its own scientific identity is incomplete.
    service.store.update(pid,data={'target_identity':{'name':'Unresolved target'}})
    with pytest.raises(ValueError,match='accession and species'):
        prepare(context)


@pytest.mark.parametrize('engine,schema',[('gnina',GninaParameters),('unidock',UniDockParameters)])
def test_native_pilot_box_retains_values_but_strips_nonparameter_metadata(context,engine,schema):
    _,_,_,protocol,_=context
    assert 'reference_sha256' in protocol['box']  # actual derive_box shape
    _,manifest=prepare(context,engine,'pilot')
    parsed=schema.model_validate(manifest.parameters)
    assert parsed.box.model_dump(mode='json')=={key:protocol['box'][key] for key in ('center','size')}
    assert set(manifest.parameters['box'])=={'center','size'}
    assert manifest.parameters['timeout_seconds']==protocol['timeout_seconds']
    invalid=deepcopy(manifest.parameters)
    invalid['box']['size']=[51,20,20]
    with pytest.raises(ValueError,match='box side'):
        schema.model_validate(invalid)


def test_native_preparation_rejects_changed_testset_binding(context):
    service,pid,jid,protocol,_=context
    with pytest.raises(ValueError,match='frozen test set'):
        auto_prepare_native_manifest(service.store,pid,'biocentric_affinity','affinity',
                                     {**protocol,'testset_sha256':'f'*64},job_id=jid)


def test_preflight_links_generated_protein_identity_evidence_without_model_execution(context,tmp_path,monkeypatch):
    service,pid,_,protocol,_=context
    artifact,manifest=prepare(context)
    monkeypatch.setattr('biocentric.engines.plan_engine_run',lambda engine,inputs,parameters,out:
                        {'engine_id':engine,'status':'ready','blockers':[],'parameters':parameters,'commands':[]})
    plan=plan_native_manifest(service.store,pid,artifact,tmp_path/'out',expected_engine='biocentric_affinity',
                              expected_phase='affinity',expected_testset_sha256=protocol['testset_sha256'],
                              expected_protocol_sha256=digest(dump(protocol).encode()))
    assert plan['resolved_artifacts']['protein_identity']['id']==manifest.protein_identity_artifact
    assert plan['frozen_input_binding']['verified']
    assert plan['engine_inputs']['protein_id']==manifest.review['protein_id']


@pytest.mark.parametrize('mutation',['protein_id','missing_evidence','receptor_binding','species','construct','chains',
                                     'isoform','mutations','sequence','domain','extra_field'])
def test_preflight_rejects_forged_protein_identity_evidence(context,tmp_path,monkeypatch,mutation):
    service,pid,_,protocol,_=context
    _,manifest=prepare(context)
    body=manifest.model_dump(mode='json')
    if mutation=='protein_id':
        body['review']['protein_id']='different_protein_claim'
    elif mutation=='missing_evidence':
        body['protein_identity_artifact']=None
    else:
        identity=json.loads(service.store.read_artifact(pid,{'id':body['protein_identity_artifact']}))
        if mutation=='receptor_binding': identity['receptor_sha256']='f'*64
        elif mutation=='species':
            identity['reviewed_identity']['organism']='Homo sapiens'
            body['review']['protein_id']=body['review']['protein_id'].replace('Mus musculus','Homo sapiens')
        elif mutation=='construct': identity['reviewed_identity']['construct']='substituted construct'
        elif mutation=='isoform': identity['reviewed_identity']['isoform']='substituted isoform'
        elif mutation=='mutations': identity['reviewed_identity']['mutations']=['T315I']
        elif mutation=='sequence': identity['reviewed_identity']['sequence']='ACDEFGHIK'
        elif mutation=='domain': identity['reviewed_identity']['domain']='substituted domain'
        elif mutation=='extra_field': identity['reviewed_identity']['unsupported_extra_claim']='not reviewed'
        else: identity['structure_chains']=['Z']
        body['protein_identity_artifact']=service.store.artifact(pid,'forged-test-identity.json',json.dumps(identity).encode(),'native_input')['id']
    altered=service.store.artifact(pid,'forged-test-manifest.json',json.dumps(body).encode(),'native_engine_affinity')
    calls=[]
    monkeypatch.setattr('biocentric.engines.plan_engine_run',lambda *args: calls.append(args))
    with pytest.raises(ValueError,match='protein identity evidence'):
        plan_native_manifest(service.store,pid,altered,tmp_path/'out',expected_engine='biocentric_affinity',
                             expected_phase='affinity',expected_testset_sha256=protocol['testset_sha256'],
                             expected_protocol_sha256=digest(dump(protocol).encode()))
    assert not calls


def test_changed_frozen_manifest_is_rejected_before_preparation_writes(context):
    service,pid,_,_,_=context
    research=service.store.project(pid)['data']['research']
    frozen=deepcopy(research['testset']);frozen['criteria']['pose_recovery']['maximum']=99
    service.store.research_update(pid,testset=frozen)
    before=service.store.artifacts(pid)
    with pytest.raises(ValueError,match='hash mismatch'):
        prepare(context)
    assert service.store.artifacts(pid)==before


def test_missing_optional_identity_fields_are_not_invented_or_required(context,tmp_path,monkeypatch):
    service,pid,_,protocol,_=context
    minimal={'accession':'P00520','name':'ABL test fixture','organism':'Mus musculus'}
    service.store.update(pid,data={'target_identity':minimal})
    artifact,manifest=prepare(context)
    identity=json.loads(service.store.read_artifact(pid,{'id':manifest.protein_identity_artifact}))
    assert identity['reviewed_identity']==minimal
    monkeypatch.setattr('biocentric.engines.plan_engine_run',lambda engine,inputs,parameters,out:
                        {'engine_id':engine,'status':'ready','blockers':[],'parameters':parameters,'commands':[]})
    plan=plan_native_manifest(service.store,pid,artifact,tmp_path/'out',expected_engine='biocentric_affinity',expected_phase='affinity')
    assert plan['engine_inputs']['protein_id'].startswith('P00520:Mus musculus:')


def test_prepared_affinity_sources_retain_registered_pose_and_physical_report(context,tmp_path,monkeypatch):
    service,pid,_,protocol,pose=context
    artifact,manifest=prepare(context)
    provenance=json.loads(service.store.read_artifact(pid,{'id':manifest.input_structure_provenance_artifact}))
    assert provenance['prepared_ligands_artifact_id']==manifest.input_artifacts['ligands']
    assert len(provenance['entries'])==1
    row=provenance['entries'][0]
    assert row['compound_id']==pose['compound_id'] and row['source_pose_sha256']==pose['pose_sha256']
    assert row['source_physical_validation']['passed'] is True
    assert row['source_validation_sha256']==pose['validation']['artifact']['sha256']
    assert row['source_receptor_sha256']==pose['receptor_sha256']
    monkeypatch.setattr('biocentric.engines.plan_engine_run',lambda engine,inputs,parameters,out:
                        {'engine_id':engine,'status':'ready','blockers':[],'parameters':parameters,'commands':[]})
    plan=plan_native_manifest(service.store,pid,artifact,tmp_path/'out',expected_engine='biocentric_affinity',expected_phase='affinity')
    assert plan['input_structure_provenance']==provenance
    assert plan['resolved_artifacts']['input_source_pose_1']['sha256']==pose['pose_sha256']
    assert plan['resolved_artifacts']['input_source_validation_1']['sha256']==row['source_validation_sha256']


def test_unregistered_or_failed_physical_validation_cannot_be_imputed_from_pose_flag(context):
    service,pid,_,_,original=context
    pose=deepcopy(original)
    pose['validation'].pop('artifact')
    service.store.update(pid,data={'screening':{'poses':[pose]}})
    with pytest.raises(ValueError,match='registered physical-validation'):
        prepare(context)
    failed=service.store.artifact(pid,'failed-test-validation.json',json.dumps({'passed':False,'status':'completed'}).encode(),'scientific_output')
    pose['validation']['artifact']={'sha256':failed['sha256']}
    service.store.update(pid,data={'screening':{'poses':[pose]}})
    with pytest.raises(ValueError,match='did not pass'):
        prepare(context)


@pytest.mark.parametrize('mutation',['cid','source_hash','prepared_hash','validation_hash','receptor','coordinates','different_passing_report'])
def test_native_preflight_rejects_substituted_affinity_pose_provenance(context,tmp_path,monkeypatch,mutation):
    service,pid,_,_,_=context
    _,manifest=prepare(context)
    body=manifest.model_dump(mode='json')
    record=json.loads(service.store.read_artifact(pid,{'id':manifest.input_structure_provenance_artifact}))
    if mutation=='cid': record['entries'][0]['compound_id']='cmp_wrong'
    elif mutation=='source_hash': record['entries'][0]['source_pose_sha256']='f'*64
    elif mutation=='prepared_hash': record['prepared_ligands_sha256']='f'*64
    elif mutation=='validation_hash': record['entries'][0]['source_validation_sha256']='f'*64
    elif mutation=='receptor': record['entries'][0]['source_receptor_sha256']='f'*64
    elif mutation=='different_passing_report':
        report=json.loads(service.store.read_artifact(pid,{'id':record['entries'][0]['source_validation_artifact_id']}))
        report['scope']='Another synthetic passing report, not recorded for this source pose.'
        alternate=service.store.artifact(pid,'unrelated-passing-report.json',json.dumps(report).encode(),'scientific_output')
        record['entries'][0]['source_validation_artifact_id']=alternate['id']
        record['entries'][0]['source_validation_sha256']=alternate['sha256']
    else:
        content=service.store.read_artifact(pid,{'id':body['input_artifacts']['ligands']})
        molecules=list(Chem.ForwardSDMolSupplier(io.BytesIO(content),removeHs=False))
        conformer=molecules[0].GetConformer()
        point=conformer.GetAtomPosition(0)
        conformer.SetAtomPosition(0,(point.x+10,point.y,point.z))
        buffer=io.StringIO()
        with Chem.SDWriter(buffer) as writer: writer.write(molecules[0])
        changed=service.store.artifact(pid,'coordinate-substitution.sdf',buffer.getvalue().encode(),'native_input')
        body['input_artifacts']['ligands']=changed['id']
        body['canonical_ligands_artifact']=changed['id']
        record['prepared_ligands_artifact_id']=changed['id']
        record['prepared_ligands_sha256']=changed['sha256']
    altered_provenance=service.store.artifact(pid,'altered-pose-provenance.json',json.dumps(record).encode(),'native_input')
    body['input_structure_provenance_artifact']=altered_provenance['id']
    altered=service.store.artifact(pid,'altered-pose-manifest.json',json.dumps(body).encode(),'native_engine_affinity')
    calls=[]
    monkeypatch.setattr('biocentric.engines.plan_engine_run',lambda *args: calls.append(args))
    with pytest.raises(ValueError,match='provenance|coordinates/chemistry|physical validation'):
        plan_native_manifest(service.store,pid,altered,tmp_path/'out',expected_engine='biocentric_affinity',expected_phase='affinity')
    assert not calls


@pytest.mark.parametrize('status',['completed','partial'])
def test_affinity_result_copies_verified_input_pose_provenance_without_mutating_native(context,tmp_path,monkeypatch,status):
    from biocentric.engine_workflow import run_native_manifest
    service,pid,_,_,pose=context
    artifact,_=prepare(context)
    monkeypatch.setattr('biocentric.engines.plan_engine_run',lambda engine,inputs,parameters,out:
                        {'engine_id':engine,'status':'ready','blockers':[],'parameters':parameters,'commands':[]})
    native={'status':status,'executed':True,'results':[
        {'compound_id':pose['compound_id'],'status':'completed','pK':5.0},
        {'compound_id':'unmapped','status':'failed','error':'synthetic test failure'}]}
    before=deepcopy(native)
    monkeypatch.setattr('biocentric.engines.run_engine',lambda *args,**kwargs:native)
    result=run_native_manifest(service.store,pid,artifact,tmp_path/'out',expected_engine='biocentric_affinity',expected_phase='affinity')
    assert result['status']==status and native==before
    enriched=result['results'][0]
    assert enriched['input_structure_physical_passed'] is True
    assert enriched['input_structure_provenance']['source_pose_sha256']==pose['pose_sha256']
    assert 'shortlist_eligible' not in enriched and 'physical_passed' not in enriched
    assert 'input_structure_physical_passed' not in result['results'][1]
    assert result['results'][1] is not native['results'][1]
    saved=json.loads(service.store.read_artifact(pid,{'id':result['result_artifact_id']}))
    assert saved['results'][0]['input_structure_provenance']==enriched['input_structure_provenance']
    assert saved['native_result']==before
