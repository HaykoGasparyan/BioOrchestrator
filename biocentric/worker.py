"""One separately executing, resumable worker with explicit scientific gates."""
from __future__ import annotations
import argparse
import json
import os
import threading
import time
import uuid
from pathlib import Path
from .models import PipelinePlan, utc
from .store import Store, dump, digest


class Gate(Exception):
    def __init__(self,kind,reason): self.kind=kind; super().__init__(reason)


class Cancelled(Exception): pass


class Pipeline:
    def __init__(self,store:Store,pid:str,jid:str,owner_token=None):
        self.store,self.pid,self.jid=store,pid,jid
        self.owner_token=owner_token
        self.spec=store.project(pid)['spec']; self.stage='01'
        self.started=time.monotonic()
        self.spent=float(store.project(pid)['data'].get('compute_seconds',0))
        self.work=store.work_directory(pid)

    @property
    def data(self): return self.store.project(self.pid)['data']
    def update(self,**data): self.store.update(self.pid,data=data)
    def check(self):
        if self.owner_token and not self.store.owns_lease(self.jid,self.owner_token): raise Cancelled('Worker no longer owns this job; successor state is retained')
        if self.store.cancelled(self.jid): raise Cancelled('Cancellation acknowledged; completed checkpoints retained')
        if self.spent+time.monotonic()-self.started>self.spec['compute_budget_seconds']:
            raise Gate('budget','Authorized compute time exhausted; create a new project with a larger budget. Completed artifacts retained.')
    def cancelled(self):
        return (self.owner_token is not None and not self.store.owns_lease(self.jid,self.owner_token)) or self.store.cancelled(self.jid) or self.spent+time.monotonic()-self.started>self.spec['compute_budget_seconds']
    def emit(self,**event):
        self.check()
        self.store.event(self.pid,stage=self.stage,job_id=self.jid,**event)
    def science_event(self,event):
        kind=event.get('kind','tool_call')
        status={'tool_started':'running','tool_completed':'completed','compound_failed':'failed'}.get(kind,'running')
        self.emit(kind='tool_call' if kind.startswith('tool_') else kind,name=event.get('tool','science'),status=status,inputs=(event.get('data') or {'compound_id':event.get('compound_id'),'protocol':self.data.get('protocol')}) if kind=='tool_started' else event.get('inputs'),outputs=event.get('data') if kind!='tool_started' else None,message=event.get('message'),details={k:v for k,v in event.items() if k not in {'kind','tool','stage','message','data','inputs'}})
    def retrieval_event(self,name,inputs,outputs):
        self.emit(kind='tool_call',name=name,status=outputs.get('status','completed'),inputs=inputs,outputs=outputs)
    def decision(self,kind):
        return next((d['value'] for d in reversed(self.store.decisions(self.pid)) if d['kind']==kind),None)
    def artifact(self,kind):
        a=next((x for x in reversed(self.store.artifacts(self.pid)) if x['kind']==kind),None)
        if a:
            self.store.read_artifact(self.pid,a)
            return Path(a['path'])
        return None
    def register(self,result):
        # Scientific outputs must stay inside this project's worker directory.
        registered=[]
        for artifact in result.get('artifacts',[]):
            self.check()
            path=Path(artifact['path']).resolve()
            if not path.is_relative_to(self.work.resolve()): raise ValueError('Scientific artifact escaped worker directory')
            content=path.read_bytes()
            if artifact.get('sha256') and digest(content)!=artifact['sha256']: raise ValueError('Scientific artifact hash mismatch')
            registered.append(self.store.artifact(self.pid,path.name,content,'scientific_output',self.jid))
        return registered

    def prepared_input(self,prepared):
        artifact_id=prepared.get('prepared_artifact_id')
        if not artifact_id: raise Gate('integrity','Prepared library has no registered immutable input; rerun preparation in a new project.')
        artifact=next((a for a in self.store.artifacts(self.pid) if a['id']==artifact_id),None)
        if not artifact: raise Gate('integrity','Registered prepared library is missing')
        self.store.read_artifact(self.pid,artifact)
        return Path(artifact['path'])

    def bind_prepared(self,result,registered):
        sha=digest(Path(result['prepared_sdf']).read_bytes())
        artifact=next((a for a in registered if a['sha256']==sha and a['name']==Path(result['prepared_sdf']).name),None)
        if not artifact: raise Gate('integrity','Prepared library was not registered successfully')
        result['prepared_artifact_id']=artifact['id']
        result['prepared_sdf']=artifact['path']
    def stage_run(self,sid,fn):
        self.stage=sid;self.check()
        step=next(s for s in self.store.steps(self.pid) if s['id']==sid)
        if step['status']=='completed':
            outputs=step.get('outputs') or {}
            ids=[outputs.get('result_artifact')]+outputs.get('registered_artifacts',[])
            try:
                if not ids[0]: raise ValueError('Completed stage has no immutable result artifact')
                for artifact_id in ids:
                    self.store.read_artifact(self.pid,{'id':artifact_id})
            except (ValueError,OSError) as exc:
                self.store.step(self.pid,sid,'blocked',reason=f'Checkpoint integrity failed: {exc}')
                raise Gate('integrity',f'Stage {sid} checkpoint integrity failed: {exc}') from exc
            self.emit(kind='checkpoint',name=step['label'],status='reused',message='Completed stage checkpoint retained')
            return
        self.store.step(self.pid,sid,'running',reason='Executing validated stage')
        self.emit(kind='stage',name=step['label'],status='running')
        start=time.monotonic()
        try:
            output=fn()
            self.check()
            a=self.store.artifact(self.pid,f'{sid}-result.json',dump(output).encode(),'stage_result',self.jid)
            self.store.step(self.pid,sid,'completed',reason='Outputs parsed and persisted',outputs={**output,'result_artifact':a['id'],'job_id':self.jid})
            self.emit(kind='stage',name=step['label'],status='completed',outputs={'artifact':a,'summary':{k:v for k,v in output.items() if k in {'status','strategy','count','limitations','result_status'}}},duration_ms=round((time.monotonic()-start)*1000))
        except Gate as exc:
            self.store.step(self.pid,sid,'blocked',reason=str(exc))
            self.update(gate={'kind':exc.kind,'reason':str(exc),'stage':sid})
            raise
        except Exception as exc:
            self.store.step(self.pid,sid,'cancelled' if isinstance(exc,Cancelled) else 'failed',reason=str(exc))
            raise

    def target(self):
        from .retrieval import RetrievalClient
        supplied=self.decision('identity')
        if supplied:
            identity={**supplied,'verification':'user-supplied identity; not independently sequence-mapped'}
            self.update(target_identity=identity,target_candidates=[])
        else:
            if self.spec['local_only']:
                raise Gate('identity','Local-only policy blocks remote mapping. Import a reviewed identity JSON or record a local identity in Decisions.')
            client=RetrievalClient(self.store.directory(self.pid)/'cache',allow_external=True,emit=self.retrieval_event)
            result=client.resolve_target(self.spec['target'],self.spec['organism'])
            self.update(target_candidates=result['candidates'],search_log=result.get('search_log',[]),limitations=result.get('limitations',[]))
            selection=self.decision('target_selection')
            candidates=result['candidates']
            if selection:
                matches=[x for x in candidates if x.get('accession')==selection.get('accession')]
                if len(matches)!=1: raise Gate('identity','Recorded selection no longer matches a unique target')
                identity=matches[0]
            elif result['status']=='resolved' and len(candidates)==1: identity=candidates[0]
            else: raise Gate('identity','Target identity unresolved or ambiguous. Review species, accession and isoform before proceeding.')
            if self.spec.get('isoform') and identity.get('isoform')!=self.spec['isoform']:
                raise Gate('identity','Requested isoform/construct needs explicit sequence mapping. Import the reviewed identity; canonical substitution is prohibited.')
            if self.spec.get('mutations'):
                raise Gate('identity','Mutant target requires explicitly reviewed sequence/construct mapping. Import identity including mutations.')
            self.update(target_identity=identity,target_candidates=[])
        identity=self.data['target_identity']
        if self.spec['local_only'] or not identity.get('accession'):
            self.update(limitations=self.data.get('limitations',[])+['External evidence retrieval disabled or target lacks a verified public mapping. User-supplied identity is not independent validation.'])
            return {'identity':identity,'status':'user_supplied','evidence_count':len(self.data.get('evidence',[]))}
        client=RetrievalClient(self.store.directory(self.pid)/'cache',allow_external=True,emit=self.retrieval_event)
        dossier=client.dossier(identity,self.spec.get('context',''))
        dossier['evidence']=self.data.get('evidence',[])+dossier.get('evidence',[])
        self.update(**dossier)
        return {'identity':identity,'evidence_count':len(dossier['evidence']),'limitations':dossier.get('limitations',[])}

    def chemistry(self):
        records=self.data.get('assays',[])
        return {'count':len(records),'assays':records,'control_status':'No labeled held-out active/inactive dataset has been approved','note':'Assay endpoint, units, relations and original source lineage preserved; no pooling across assay contexts.'}

    def strategy(self):
        if self.spec['workflow_kind']!='noncovalent': raise Gate('unsupported','This release executes reversible noncovalent protein–small-molecule workflows only. A separate validated protocol is needed for this target/chemistry.')
        structure_route=all(self.artifact(kind) for kind in ('receptor','library','reference'))
        plan=PipelinePlan(strategy='structure_based' if structure_route else ('ligand_based' if self.data.get('assays') else 'further_experiments'),rationale='Curated receptor and observed-site reference are present; inspect and benchmark before screening.' if structure_route else 'A reviewed receptor PDBQT, site-reference SDF and ligand library are needed for the executable docking path. Research continues without inventing a pocket.',assumptions=['Input protonation and stereochemistry are retained; pH is a recorded assumption, not a calculated assignment.','Only a small smoke/pilot configuration is proposed.'],prerequisites=['Reviewed receptor state and supported binding pocket','Curated receptor PDBQT and matching PDB for pose checks','Bound ligand SDF in the same coordinate frame','Explicit exploratory choice or independent validated benchmark'],stages=[s['label'] for s in self.store.steps(self.pid)]).model_dump()
        self.update(plan=plan)
        skip_reason=None
        if self.spec['local_only']:
            skip_reason=('local_only','This project keeps inputs local, so no remote language model is called.')
        elif not self.spec['allow_model']:
            skip_reason=('disabled','The language model is disabled for this project.')
        elif self.spec.get('max_model_calls',3)==0:
            skip_reason=('zero_budget','The authorized model-call budget is zero, so no language model is called.')
        if skip_reason:
            code,reason=skip_reason
            plan['model_advice']={'status':'not_requested','skip_reason':code,'reason':reason,'model':self.spec['model']}
            self.emit(kind='model_status',name='Language model',status='not_requested',message=reason,outputs=plan['model_advice'])
        else:
            from .agent import advise
            if self.data.get('model_attempted'):
                reason='Prior advisory attempt retained; no automatic paid retry.'
                plan['model_advice']=self.data.get('model_advice',{'status':'prior_attempt','reason':reason,'model':self.spec['model']})
                plan['model_retry_policy']=reason
                self.emit(kind='model_status',name='Language model',status='reused',message=reason,outputs=plan['model_advice'])
            else:
                # Persist before submitting: failed retries cannot silently exceed authorization.
                self.update(model_attempted=True)
                try:
                    advice=advise(self.spec,self.data.get('evidence',[]),self.emit,self.cancelled)
                    plan['model_advice']={**advice,'status':'completed'}
                except Exception as exc:
                    if self.cancelled(): self.check()
                    self.emit(kind='model_call',name='OpenAI Responses',status='failed',message=str(exc))
                    plan['model_advice']={'status':'failed','reason':str(exc),'model':self.spec['model']}
        self.update(plan=plan,model_advice=plan['model_advice'])
        return plan

    def structure(self):
        from .science import inspect_structure, derive_box
        required={k:self.artifact(k) for k in ('receptor','reference','library','structure')}
        if any(v is None for v in required.values()):
            missing=[k for k,v in required.items() if v is None]
            raise Gate('inputs','Upload curated inputs: '+', '.join(missing)+'. Automatic structure prediction is not installed; predicted/uncertain sites require a reviewed modeling plan.')
        if required['structure'].suffix.lower()!='.pdb':
            inspection=inspect_structure(required['structure']);self.update(structure_inspection=inspection)
            raise Gate('inputs','mmCIF original retained. Curated docking/pose validation requires an inspected matching PDB export and receptor PDBQT.')
        inspection=inspect_structure(required['receptor'])
        full=inspect_structure(required['structure'])
        self.update(structure_inspection={'receptor':inspection,'original':full})
        blockers=inspection.get('routine_screening_blockers',[])+full.get('routine_screening_blockers',[])
        if blockers: raise Gate('unsupported',' '.join(blockers))
        # Input hashes bind every approval to an immutable protocol proposal.
        protocol={'engine':'AutoDock Vina','engine_version':'1.2.7','seed':self.spec['seed'],'cpu':self.spec['cpu'],'exhaustiveness':4,'n_poses':3,'box':derive_box(required['reference']),'box_padding_angstrom':4.0,'timeout_seconds':min(120,self.spec['compute_budget_seconds']),'require_posebusters':True,'run_plip':True,'receptor_pdb':str(required['structure']),'input_hashes':{k:digest(v.read_bytes()) for k,v in required.items()},'validation':'unvalidated; pilot and held-out criteria required','expected_artifacts':['prepared SDF','PDBQT poses','canonical SDF poses','PoseBusters flags','available PLIP reports','CSV/JSON shortlist','run manifest'],'limitations':['No automatic missing-loop repair, protonation-state enumeration, or metal protocol.','Receptor assembly, residue mapping and pocket quality need scientist review.']}
        self.update(protocol=protocol)
        decision=self.decision('protocol')
        if not decision or not decision.get('approved') or decision.get('protocol')!=protocol:
            raise Gate('protocol','Review the proposed protocol, receptor components, pocket and supported chemistry; record approval in Plan & Runs.')
        return {'inspection':inspection,'original_inspection':full,'protocol':protocol,'review':decision}

    def preparation(self):
        from .science import prepare_library
        result=prepare_library(self.artifact('library'),self.work/'preparation',seed=self.spec['seed'],emit=self.science_event,cancelled=self.cancelled)
        registered=self.register(result)
        self.bind_prepared(result,registered)
        self.update(preparation=result,exclusions=result['exclusions'])
        if not result['compounds']: raise Gate('chemistry','No supported compounds survived preparation; see explicit exclusions and preserved originals.')
        return {**result,'registered_artifacts':[a['id'] for a in registered]}

    def config(self):
        protocol=self.data['protocol']
        return {k:protocol[k] for k in ('seed','cpu','exhaustiveness','n_poses','box_padding_angstrom','timeout_seconds','require_posebusters','run_plip','receptor_pdb')}

    def pilot(self):
        from .science import prepare_library,dock_library
        prepared=prepare_library(self.artifact('reference'),self.work/'pilot_preparation',self.spec['seed'],emit=self.science_event,cancelled=self.cancelled)
        prepared_artifacts=self.register(prepared)
        self.bind_prepared(prepared,prepared_artifacts)
        if not prepared['compounds']: raise Gate('chemistry','Site-reference ligand is outside the supported preparation domain')
        result=dock_library(self.artifact('receptor'),self.prepared_input(prepared),self.artifact('reference'),self.work/'pilot',self.config(),emit=self.science_event,cancelled=self.cancelled)
        registered=self.register(result);self.update(pilot=result)
        if not result['poses']: raise Gate('pilot','Pilot calculation failed; inspect tool errors. Missing calculations cannot be replaced with model scores.')
        return {'status':'pilot_completed','poses':result['poses'],'failures':result.get('failures',[]),'registered_artifacts':[a['id'] for a in registered],'limitations':['Redocking is not a held-out enrichment benchmark.','RMSD is meaningful only for equivalent ligand identities in the same receptor coordinate frame.']}

    def benchmark(self):
        result={'status':'exploratory','criteria':{'held_out_enrichment':'not established','physical_pose_validity':'required for shortlist'},'metrics':{},'dataset_sha256':None,'limitations':['No independent active/inactive held-out screening dataset in this project.','Pilot pose recovery does not validate ranking; pretrained-model training overlap unknown.']}
        self.update(benchmark=result)
        if not self.decision('exploratory'):
            raise Gate('exploratory','No held-out benchmark establishes predictive performance. Record an explicit exploratory choice to continue; all exports will retain that status.')
        self.update(result_status='exploratory',gate=None)
        return result

    def screening(self):
        from .science import dock_library
        result=dock_library(self.artifact('receptor'),self.prepared_input(self.data['preparation']),self.artifact('reference'),self.work/'screening',self.config(),emit=self.science_event,cancelled=self.cancelled)
        registered=self.register(result)
        for pose in result['poses']:
            pose['job_id']=self.jid
            pose['artifact_ids']=[a['id'] for a in registered if a['sha256']==pose.get('pose_sha256')]
        self.update(screening=result)
        if not result['poses']: raise Gate('screening','Every compound failed. Inspect exclusions and tool failures; no shortlist can be computed.')
        return {'status':'completed','count':len(result['poses']),'failures':result.get('failures',[]),'metadata':result.get('metadata'),'registered_artifacts':[a['id'] for a in registered]}

    def prioritization(self):
        from rdkit import Chem, DataStructs
        from rdkit.Chem import rdFingerprintGenerator
        poses=self.data['screening']['poses']; by_compound={}
        for p in poses:
            old=by_compound.get(p['compound_id'])
            key=(p.get('validation',{}).get('passed') is not True,p['score'])
            if old is None or key<(old.get('validation',{}).get('passed') is not True,old['score']): by_compound[p['compound_id']]=p
        rows=sorted(by_compound.values(),key=lambda p:(p.get('validation',{}).get('passed') is not True,p['score']))
        chosen=[]; candidates=[]; generator=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048)
        for rank,pose in enumerate(rows,1):
            fp=generator.GetFingerprint(Chem.MolFromSmiles(pose['canonical_smiles']))
            similarity=max((DataStructs.TanimotoSimilarity(fp,other) for other in chosen),default=0.0)
            valid=pose.get('validation',{}).get('passed') is True
            selected=valid and len(chosen)<self.spec['shortlist_size'] and similarity<0.65
            if selected: chosen.append(fp)
            rationale=(f"Selected by {pose['engine']} docking rank and fingerprint diversity among physically valid poses." if selected else 'Excluded from shortlist: '+('pose validity failed or incomplete' if not valid else 'similarity to selected representative ≥0.65 or shortlist capacity reached')+'.')
            native_metrics=pose.get('metrics')
            metrics={metric['name']:{**metric,'engine':pose['engine'],'version':pose['engine_version']} for metric in native_metrics} if isinstance(native_metrics,list) else {'vina_score':{'value':pose['score'],'units':'kcal/mol','direction':'lower_is_better','engine':pose['engine'],'version':pose['engine_version']}}
            metrics['descriptors']=pose.get('descriptors',{})
            candidates.append({**pose,'rank':rank,'selected':selected,'rationale':rationale,'max_similarity_to_previously_selected':similarity,'missing_evidence':['measured binding','selectivity','ADMET','compound availability','synthesis feasibility'],'status':'exploratory','metrics':metrics})
        self.update(candidates=candidates)
        return {'count':len(candidates),'selected':len(chosen),'candidates':candidates,'limitations':['Fingerprint diversity is a transparent selection heuristic, not a measured selectivity property.','Top rank cannot override failed physical pose checks.']}

    def handoff(self):
        self.update(gate=None)
        records=[];selection=[]
        known={a['id']:a for a in self.store.artifacts(self.pid)}
        for candidate in self.data.get('candidates',[]):
            if not candidate.get('selected'): continue
            self.check()
            artifact=next((known[aid] for aid in candidate.get('artifact_ids',[]) if aid in known and known[aid]['name'].lower().endswith('.sdf') and known[aid]['sha256']==candidate.get('pose_sha256')),None)
            if not artifact: raise Gate('integrity','Selected compound has no registered canonical SDF pose')
            content=self.store.read_artifact(self.pid,artifact)
            if not content.rstrip().endswith(b'$$$$'): raise Gate('integrity','Selected pose SDF lacks a complete record delimiter')
            records.append(content.rstrip()+b'\n')
            selection.append({'compound_id':candidate['compound_id'],'rank':candidate.get('rank'),'pose_artifact_id':artifact['id'],'pose_sha256':artifact['sha256'],'job_id':candidate['job_id']})
        shortlist=self.store.artifact(self.pid,'shortlist.sdf',b''.join(records),'shortlist',self.jid)
        selection_manifest=self.store.artifact(self.pid,'shortlist-selection.json',dump({'status':'exploratory','selected_count':len(selection),'sdf_artifact_id':shortlist['id'],'sdf_sha256':shortlist['sha256'],'selection':selection}).encode(),'shortlist',self.jid)
        self.update(shortlist={'sdf':shortlist,'selection_manifest':selection_manifest,'selected_count':len(selection)})
        from .reports import export_project
        return {'status':'completed','reports':export_project(self.store,self.pid),'shortlist':shortlist,'registered_artifacts':[shortlist['id'],selection_manifest['id']],'result_status':'exploratory','proposal':'Confirm compound identity/purity, run repeat dose responses, orthogonal binding/function, interference and selectivity counterscreens, then relevant cellular target engagement.'}

    def run(self):
        for sid,method in [('01',self.target),('02',self.chemistry),('03',self.strategy),('04',self.structure),('05',self.preparation),('06',self.pilot),('07',self.benchmark),('08',self.screening),('09',self.prioritization),('10',self.handoff)]:
            self.stage_run(sid,method)


def execute(root,pid,jid):
    store=Store(root)
    owner_token=uuid.uuid4().hex
    if not store.claim(jid,owner_token=owner_token,worker_pid=os.getpid(),project_id=pid): return False
    stop=threading.Event()
    def pulse():
        while not stop.wait(15):
            if not store.heartbeat(jid,owner_token): return
    heartbeat=threading.Thread(target=pulse,daemon=True);heartbeat.start()
    pipeline=Pipeline(store,pid,jid,owner_token)
    store.update(pid,state='running')
    status,error='completed',None
    try: pipeline.run()
    except Gate as exc:
        status,error='blocked',str(exc)
        if store.owns_lease(jid,owner_token): pipeline.update(gate={'kind':exc.kind,'reason':str(exc),'stage':pipeline.stage})
        store.event(pid,kind='gate',stage=pipeline.stage,name=exc.kind,status='blocked',message=str(exc))
    except Cancelled as exc: status,error='cancelled',str(exc)
    except Exception as exc:
        status,error=('cancelled' if getattr(exc,'code',None)=='cancelled' else 'failed'),f'{type(exc).__name__}: {exc}'
        store.event(pid,kind='error',stage=pipeline.stage,name=type(exc).__name__,status=status,message=error)
    finally:
        stop.set();heartbeat.join(timeout=1)
        finished=store.finish(jid,status,error,owner_token=owner_token,project_id=pid,project_data={'compute_seconds':pipeline.spent+time.monotonic()-pipeline.started})
        if finished:
            from .reports import export_project
            export_project(store,pid)
    return True


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',required=True);p.add_argument('--project',required=True);p.add_argument('--job',required=True)
    args=p.parse_args();execute(args.root,args.project,args.job)
if __name__=='__main__': main()
