"""Shared application boundary for CLI and Streamlit; scientific work is out of process."""
from __future__ import annotations
import csv
import io
import json
import os
import platform
import re
import subprocess
import sys
from pathlib import Path
from .models import ProjectSpec, TargetIdentity, normalize_measurement
from .store import Store, dump, digest

STAGE_CONTRACTS = [
 ('01','Target dossier','UniProt / evidence sources','ProjectSpec → species-safe TargetIdentity','Normalize target input; never silently substitute isoforms','Check accession/species; citations retain retrieval dates'),
 ('02','Known chemistry & controls','ChEMBL / PubChem','TargetIdentity → AssayMeasurement[]','Preserve original endpoint, units and relation','Deduplicate source lineage; no IC50 to Kd conversion'),
 ('03','Strategy selection','Deterministic gate + optional OpenAI advisor','EvidenceRecord[] + available inputs → PipelinePlan','Check supported workflow and usable inputs','Separate proposal from authorization; no fabricated runs'),
 ('04','Structure selection','Local structural inspection','PDB / mmCIF / curated PDBQT → inspection','Preserve original chains, waters, metals and coordinate frame','Require reviewed receptor/pocket; unsupported components block'),
 ('05','Protein & ligand preparation','RDKit / Meeko','SMILES / SDF → parent/state SDF + PDBQT','Sanitize, preserve stereochemistry; explicit H and seeded conformers','Record exclusions and properties; no silent protonation assumptions'),
 ('06','Protocol & pilot','AutoDock Vina / PoseBusters','Bound SDF + prepared receptor → pilot poses','Derive search box from observed ligand coordinates','Physical validity and pose recovery are distinct from affinity'),
 ('07','Benchmark gate','Held-out benchmark / recorded exploratory choice','Predeclared criteria + controls → BenchmarkResult','Separate tuning from evaluation; no controls means no pass','Exploratory choice is explicit and survives every export'),
 ('08','Screening & inspection','AutoDock Vina / PoseBusters / PLIP when installed','Prepared ligand states → poses + named raw scores','Reuse approved seed, CPU, receptor and site','Parse scores; preserve per-compound failures and pose checks'),
 ('09','Prioritization','RDKit diversity + deterministic ranking','Valid poses + descriptors → CandidateAssessment[]','Keep score semantics; sort within one engine/protocol','Balance ranking and diversity; missing ADMET remains unknown'),
 ('10','Experimental handoff','Report and manifest exporter','CandidateAssessment[] → Markdown / HTML / CSV / JSON / manifest','Collect artifact hashes, decisions, tool versions and limitations','Scientist-review assay proposals; never claim experiments ran'),
]


def initial_steps():
    return [{"id":i,"label":label,"tool":tool,"inputs":contract,"preprocessing":pre,"postprocessing":post,"status":"pending","reason":"Awaiting prerequisite stages","outputs":None} for i,label,tool,contract,pre,post in STAGE_CONTRACTS]


class Service:
    def __init__(self,root: str | Path | None=None):
        self.root=Path(root or os.getenv('BIOCENTRIC_HOME','.biocentric')).resolve()
        self.store=Store(self.root)

    def list_projects(self): return self.store.list_projects()

    def create_project(self,spec:dict)->str:
        validated=ProjectSpec.model_validate(spec)
        if validated.local_only and validated.allow_model: raise ValueError("Local-only projects cannot submit data to a remote model")
        steps=initial_steps() if validated.workflow_mode=='classic' else [
            {'id':'01','label':'Understand the research question','status':'pending','reason':'Ready to discuss the target and available inputs'}]
        return self.store.create(validated.model_dump(),steps)

    def snapshot(self,pid):
        project=self.store.project(pid); data=project.pop('data')
        # Normalize old pending model questions in this returned copy only. The
        # stored answer contract is upgraded transactionally on submission.
        for question in data.get('research',{}).get('questions',[]):
            if question.get('kind')=='model_clarification' and question.get('status')=='pending':
                original={'prompt':question.get('prompt'),'options':question.get('options',[])}
                question['context']={**(question.get('context') or {}),'original_model_question':original}
                question['kind']='model_clarification_control'
                question['options']=[{'label':'Continue with frozen inputs','value':'continue_frozen'},
                                     {'label':'Revise inputs before continuing','value':'revise_inputs'}]
                question['prompt']=(question.get('prompt') or '')+'\nChoose whether to continue with the current frozen inputs or pause to revise them. Your choice will not silently change model parameters.'
        return {"project":project,"data":data,"steps":self.store.steps(pid),"events":self.store.events(pid),"jobs":self.store.jobs(pid),"decisions":self.store.decisions(pid),"artifacts":self.store.artifacts(pid),"evidence":data.get('evidence',[]),"candidates":data.get('candidates',[]),"relocation":self.store.relocation.describe()}

    def create_research_project(self, question: str, *, use_public_example: bool = False,
                                model: str = 'gpt-4.1-mini', max_model_calls: int = 6) -> str:
        goal=question.strip()
        if not goal or len(goal)>4000: raise ValueError('Describe your research question in 1–4,000 characters')
        if use_public_example:
            pid=self.create_example(allow_model=True,model=model,workflow_mode='research',max_model_calls=max_model_calls)
            target={'target':'P00520','organism':'Mus musculus'}
        else:
            pid=self.create_project({'name':goal[:100],'target':'Unresolved target','organism':'Unspecified',
                                     'context':goal,'workflow_mode':'research','local_only':False,'allow_model':True,
                                     'model':model,'max_model_calls':max_model_calls,'max_model_output_tokens':2200,
                                     'compute_budget_seconds':1800,'cpu':1})
            target=None
        self.store.research_update(pid,goal=goal,phase='ready',questions=[],target=target,
                                   hardware=f'Local {platform.system()} / {platform.machine()}; runtime capability checks determine eligible devices',
                                   public_example=use_public_example)
        self.store.event(pid,kind='research_message',role='user',message=goal,status='recorded')
        self.store.event(pid,kind='research_message',role='assistant',status='recorded',message=
            'I will search the literature, construct a small test set, propose suitable tools, and show the pilot results before a larger run. Any questions will appear here; answering continues the same workflow.')
        return pid

    def start_research(self,pid,*,spawn=True):
        project=self.store.project(pid)
        if project['spec'].get('workflow_mode')!='research': raise ValueError('This is a classic workflow project')
        if project['spec']['mode']=='offline_fixture': raise ValueError('Offline fixtures cannot execute research')
        pending=[q for q in project['data'].get('research',{}).get('questions',[]) if q.get('status')=='pending']
        if pending: raise ValueError('Answer the question in the conversation to continue')
        if project['state']=='completed':
            jobs=self.store.jobs(pid)
            return jobs[0]['id'] if jobs else None
        jid=self.store.enqueue(pid)
        if self.store.jobs(pid)[0]['status']=='running': return jid
        self.store.update(pid,state='queued')
        self.store.research_update(pid,phase='queued')
        self.store.event(pid,kind='control',name='Research queued',status='queued',outputs={'job_id':jid})
        if spawn:
            env=os.environ.copy();env['PYTHONPATH']=str(Path(__file__).parent.parent)
            if project['spec']['local_only'] or not project['spec']['allow_model']: env.pop('OPENAI_API_KEY',None)
            with (self.store.directory(pid)/'worker.log').open('ab') as log:
                subprocess.Popen([sys.executable,'-m','biocentric.research_worker','--root',str(self.root),
                                  '--job',jid,'--project',pid],stdout=log,stderr=subprocess.STDOUT,
                                 cwd=str(Path(__file__).parent.parent),env=env,start_new_session=True)
        return jid

    def answer_research(self,pid,question_id,answer,*,spawn=True):
        from .research_workflow import validate_answer
        resume=True
        with self.store.connect() as connection:
            connection.execute('BEGIN IMMEDIATE')
            row=connection.execute('SELECT spec,data FROM projects WHERE id=?',(pid,)).fetchone()
            if not row: raise ValueError('Unknown project')
            spec,data=json.loads(row['spec']),json.loads(row['data'])
            if spec.get('workflow_mode')!='research': raise ValueError('This is not a research conversation')
            if connection.execute("SELECT 1 FROM jobs WHERE project_id=? AND status IN ('queued','running')",(pid,)).fetchone():
                raise ValueError('Wait for the worker to reach its question before answering')
            questions=data.setdefault('research',{}).setdefault('questions',[])
            question=next((q for q in questions if q['id']==question_id),None)
            if not question or question.get('status')!='pending': raise ValueError('This question is no longer awaiting an answer')
            if question.get('kind')=='model_clarification':
                question['context']={**(question.get('context') or {}),'original_model_question':{
                    'prompt':question.get('prompt'),'options':question.get('options',[])}}
                question['kind']='model_clarification_control'
                question['options']=[{'label':'Continue with frozen inputs','value':'continue_frozen'},
                                     {'label':'Revise inputs before continuing','value':'revise_inputs'}]
            value=validate_answer(question,answer)
            research=data['research']
            completed_pilot=any(research.get(key) for key in ('pilot','comparison','selection_lock')) or connection.execute(
                "SELECT 1 FROM steps WHERE project_id=? AND id>='06' AND status='completed'",(pid,)).fetchone() is not None
            if question.get('kind')=='model_clarification_control' and value=='revise_inputs':
                if completed_pilot:
                    raise ValueError('A completed pilot protects these inputs. This run remains paused; start a new research conversation to revise them, or explicitly continue with the frozen inputs.')
                used=sum(json.loads(row['data']).get('kind')=='model_call' and json.loads(row['data']).get('status')=='running'
                         for row in connection.execute('SELECT data FROM events WHERE project_id=?',(pid,)))
                remaining=max(0,spec['max_model_calls']-used)
                revision=int(research.get('input_revision_version',0))+1
                research['input_revision']={'version':revision,'status':'awaiting_inputs','source_question_id':question_id,
                    'baseline_artifact_ids':[row['id'] for row in connection.execute('SELECT id FROM artifacts WHERE project_id=?',(pid,))],
                    'changed_inputs':[],'remaining_calls_at_request':remaining}
                questions.append({'id':f'input_revision_{revision}','kind':'research_input_revision','status':'pending','stage_id':'04',
                    'prompt':f'Upload revised scientific inputs or matched assay CSV using Add scientific input files. Then rebuild the frozen test set and tool proposal. Rebuilding uses at most 1 GPT call; {remaining} calls remain in this project. Existing sources and completed results are retained.',
                    'options':[{'label':f'Rebuild test set and proposal · 1 GPT call ({remaining} remaining)','value':'rebuild'}],
                    'context':{'source_question_id':question_id,'changed_inputs':[],'model_calls_to_rebuild':1,'remaining_model_calls':remaining}})
                research['phase']='awaiting_revised_inputs'
                resume=False
            elif question.get('kind')=='research_input_revision':
                revision=research.get('input_revision') or {}
                if completed_pilot:
                    raise ValueError('Completed pilot inputs cannot be revised; start a new research conversation')
                if value!='rebuild' or revision.get('status')!='awaiting_inputs' or not revision.get('changed_inputs'):
                    raise ValueError('Upload at least one changed scientific input or assay CSV before rebuilding')
                used=sum(json.loads(row['data']).get('kind')=='model_call' and json.loads(row['data']).get('status')=='running'
                         for row in connection.execute('SELECT data FROM events WHERE project_id=?',(pid,)))
                if used>=spec['max_model_calls']:
                    raise ValueError('No model calls remain. Revised inputs are saved; create a new project with a sufficient budget to rebuild the proposal.')
                history={'version':research.get('input_revision_version',0),'testset':research.get('testset'),
                         'tool_plan':research.get('tool_plan'),'protocol':research.get('protocol'),
                         'model_answers':research.get('model_answers',[]),'stage_results':[]}
                for row in connection.execute("SELECT id,data FROM steps WHERE project_id=? AND id>='04'",(pid,)).fetchall():
                    previous=json.loads(row['data'])
                    history['stage_results'].append(previous)
                    fresh={**previous,'status':'pending','reason':'Revised inputs approved; rebuilding test set and tool proposal','outputs':None}
                    connection.execute('UPDATE steps SET status=?,data=? WHERE project_id=? AND id=?',('pending',dump(fresh),pid,row['id']))
                research.setdefault('revision_history',[]).append(history)
                for key in ('testset','tool_plan','protocol','model_answers','engine_preflight'):
                    research.pop(key,None)
                data['protocol']=None
                for other in questions:
                    if other is not question and other.get('status')=='pending' and str(other.get('stage_id',''))>='04':
                        other['status']='superseded_by_input_revision'
                research['input_revision_version']=revision['version']
                revision['status']='rebuilding'
                revision['approved_model_calls']=1
            question.update(status='answered',answer=value)
            if question.get('kind')=='target_identity': data['research']['target']=value
            connection.execute('UPDATE projects SET data=? WHERE id=?',(dump(data),pid))
        self.store.decision(pid,question_id,value,'User answered an inline research question')
        self.store.event(pid,kind='research_message',role='user',stage=question.get('stage_id'),status='recorded',
                         message=answer if isinstance(answer,str) else json.dumps(answer),question_id=question_id)
        self.store.event(pid,kind='research_message',role='assistant',stage=question.get('stage_id'),status='recorded',
                         message='Answer recorded. I am continuing from the saved step.' if resume else 'Revision requested. Upload the new inputs, then use the explicit rebuild button. No pilot or model request will run while I wait.',question_id=question_id)
        if resume:
            return self.start_research(pid,spawn=spawn)
        jobs=self.store.jobs(pid)
        return jobs[0]['id'] if jobs else None

    def upload(self,pid,name,content:bytes,kind:str):
        if not content or len(content)>50*1024*1024: raise ValueError("Upload must contain 1 byte to 50 MB")
        formats={'receptor':{'.pdbqt'},'library':{'.sdf','.smi','.smiles'},'reference':{'.sdf'},'structure':{'.pdb','.cif','.mmcif'},'identity':{'.json'},'evidence':{'.json'},
                 'native_input':{'.prmtop','.parm7','.nc','.mdcrd','.map','.fld','.pdb','.pdbqt','.sdf','.json'},
                 'native_engine_pilot':{'.json'},'native_engine_screening':{'.json'},'native_engine_affinity':{'.json'}}
        if kind not in formats or Path(name).suffix.lower() not in formats[kind]: raise ValueError(f"Unsupported {kind} file format")
        snapshot=self.snapshot(pid)
        if any(j['status'] in {'running','queued'} for j in snapshot['jobs']): raise ValueError("Cancel the running job before modifying inputs")
        research=snapshot['project']['spec'].get('workflow_mode')=='research'
        if research:
            revision=snapshot['data'].get('research',{}).get('input_revision') or {}
            revising=revision.get('status')=='awaiting_inputs'
            completed_pilot=any(snapshot['data'].get('research',{}).get(key) for key in ('pilot','comparison','selection_lock')) or any(s['id']>='06' and s['status']=='completed' for s in snapshot['steps'])
            if kind.startswith('native_'):
                if any(s['id']>='08' and s['status']=='completed' for s in snapshot['steps']):
                    raise ValueError('Completed calculations are immutable; use a new research project for changed native inputs')
            elif completed_pilot or (any(s['id']>='04' and s['status']=='completed' for s in snapshot['steps']) and not revising):
                raise ValueError('The test set is frozen. Create a new research project to change its scientific inputs or evidence')
        elif any(s['status']=='completed' for s in snapshot['steps'][3:]): raise ValueError("Create a new project for changed scientific inputs after preparation; completed runs are immutable")
        if kind in {'receptor','library','reference','structure'} or kind.startswith('native_engine_'):
            versions=[row for row in snapshot['artifacts'] if row['kind']==kind]
            current=versions[-1] if versions else None
            sha=digest(content)
            if current and current['sha256']==sha:
                self.store.read_artifact(pid,current)
                return {**current,'unchanged':True}
            clean=re.sub(r'[^A-Za-z0-9._-]','_',Path(name).name)[:150] or 'artifact'
            historical=next((row for row in versions if row['sha256']==sha and row['name']==clean),None)
            if historical:
                raise ValueError(f'This exact file is an older {kind} input and would not replace the current input. Use a new filename to register this version as current, or keep the current input. No input was changed.')
        if kind=='identity':
            identity=TargetIdentity.model_validate_json(content).model_dump()
            self.decide(pid,'identity',identity,'User imported target identity; source remains user-supplied')
        if kind=='evidence':
            from .models import EvidenceRecord
            records=json.loads(content)
            if not isinstance(records,list): raise ValueError("Evidence import must be a JSON array")
            records=[EvidenceRecord.model_validate(r).model_dump() for r in records]
            self.store.update(pid,data={'evidence':snapshot['evidence']+records})
        if kind in {'receptor','library','reference','structure'} and not research:
            for sid in ('03','04'):
                self.store.step(pid,sid,'pending',reason='Scientific input changed; proposal must be recomputed',outputs=None)
            self.store.update(pid,data={'protocol':None})
        artifact=self.store.artifact(pid,name,content,kind)
        if research and revising and not kind.startswith('native_') and artifact['id'] not in revision['baseline_artifact_ids']:
            changed={'artifact_id':artifact['id'],'kind':kind,'name':artifact['name']}
            if changed not in revision['changed_inputs']: revision['changed_inputs'].append(changed)
            questions=snapshot['data']['research'].get('questions',[])
            for pending in questions:
                if pending.get('kind')=='research_input_revision' and pending.get('status')=='pending':
                    pending.setdefault('context',{})['changed_inputs']=revision['changed_inputs']
            self.store.research_update(pid,input_revision=revision,questions=questions)
        self.store.event(pid,kind='input',name=f'Import {kind}',status='completed',inputs={'name':name,'bytes':len(content)},outputs=artifact)
        return artifact

    def decide(self,pid,kind,value,reason):
        snap=self.snapshot(pid)
        research=snap['project']['spec'].get('workflow_mode')=='research'
        if research and snap['data'].get('research',{}).get('target'):
            snap['project']['spec']={**snap['project']['spec'],**snap['data']['research']['target']}
        if snap['project']['spec']['mode']=='offline_fixture': raise ValueError("Fixture projects cannot authorize real work")
        if not reason.strip(): raise ValueError("Decision rationale is required")
        if any(j['status'] in {'running','queued'} for j in snap['jobs']): raise ValueError("Wait for a gate or cancel before changing decisions")
        if kind not in {'target_selection','identity','protocol','exploratory'}: raise ValueError("Unknown scientific decision")
        resolved=next((step for step in snap['steps'] if step['id']==('03' if research else '01')), {})
        if kind in {'target_selection','identity'} and resolved.get('status')=='completed': raise ValueError("Resolved identity is immutable; create a new project to change target")
        if kind=='target_selection':
            candidates=snap['data'].get('target_candidates',[])
            if not any(c.get('accession')==value.get('accession') for c in candidates): raise ValueError("Choose an actual resolved candidate")
        if kind=='identity':
            value=TargetIdentity.model_validate(value).model_dump()
            if value['organism'].casefold()!=snap['project']['spec']['organism'].casefold(): raise ValueError("Identity organism differs from project")
            if not value.get('accession') and not value.get('sequence'): raise ValueError("Custom target requires sequence or accession")
            target=snap['project']['spec']['target']
            from .retrieval import ACCESSION, AA
            if ACCESSION.fullmatch(target.upper()) and value.get('accession')!=target.upper(): raise ValueError("Identity accession differs from project target")
            requested_isoform=snap['project']['spec'].get('isoform')
            if requested_isoform and value.get('isoform')!=requested_isoform: raise ValueError("Requested isoform must be preserved")
            requested_mutations=snap['project']['spec'].get('mutations',[])
            if sorted(x.strip().casefold() for x in requested_mutations)!=sorted(x.strip().casefold() for x in value.get('mutations',[])):
                raise ValueError("Imported identity mutations must match the project's requested mutations")
            supplied_sequence=value.get('sequence')
            if supplied_sequence:
                supplied_sequence=''.join(supplied_sequence.split()).upper()
                if not set(supplied_sequence)<=AA: raise ValueError("Identity sequence contains unsupported amino-acid characters")
                value['sequence']=supplied_sequence
            target_lines=target.splitlines()
            target_sequence=''.join(''.join(line.split()) for line in target_lines if not line.lstrip().startswith('>')).upper()
            if target.startswith('>') or (len(target_sequence)>=30 and set(target_sequence)<=AA):
                if sum(line.lstrip().startswith('>') for line in target_lines)>1: raise ValueError("Select one FASTA target; multiple records cannot be concatenated")
                if not target_sequence or not set(target_sequence)<=AA: raise ValueError("Project target contains an invalid amino-acid sequence")
                if supplied_sequence!=target_sequence: raise ValueError("Imported identity sequence must match the project's supplied target sequence")
        if kind=='protocol':
            if not snap['data'].get('protocol'): raise ValueError("Run planning first; approval must reference a concrete proposed protocol")
            if value.get('approved') is not True or value.get('receptor_reviewed') is not True: raise ValueError("Explicit receptor, site and chemistry review required")
            ph=float(value.get('ph_assumption',7.4))
            if not 0<=ph<=14: raise ValueError("pH must be within 0–14")
            value={**value,'ph_assumption':ph,'protocol':snap['data']['protocol']}
        if kind=='exploratory' and value.get('approved') is not True: raise ValueError("Exploratory approval must be explicit")
        self.store.decision(pid,kind,value,reason)
        return True

    def start(self,pid,*,spawn=True):
        p=self.store.project(pid)
        if p['spec'].get('workflow_mode')=='research':
            return self.start_research(pid,spawn=spawn)
        if p['spec']['mode']=='offline_fixture': raise ValueError("Offline fixtures cannot execute scientific jobs")
        jobs=self.store.jobs(pid)
        if jobs and jobs[0]['status']=='completed': return jobs[0]['id']
        jid=self.store.enqueue(pid)
        # Repeated Start must not relabel a live worker as queued or spawn another.
        if self.store.jobs(pid)[0]['status']=='running': return jid
        self.store.update(pid,state='queued')
        self.store.event(pid,kind='control',name='Pipeline queued',status='queued',outputs={'job_id':jid})
        if spawn:
            log=self.store.directory(pid)/'worker.log'
            env=os.environ.copy(); env['PYTHONPATH']=str(Path(__file__).parent.parent)
            if p['spec']['local_only'] or not p['spec']['allow_model']:
                env.pop('OPENAI_API_KEY',None)
            with log.open('ab') as output:
                subprocess.Popen([sys.executable,'-m','biocentric.worker','--root',str(self.root),'--job',jid,'--project',pid],stdout=output,stderr=subprocess.STDOUT,cwd=str(Path(__file__).parent.parent),env=env,start_new_session=True)
        return jid

    def resume(self,pid): return self.start(pid)
    def cancel(self,pid): self.store.cancel(pid); return True

    def repair_model_advice(self, pid):
        """Explicit one-shot correction using one remaining authorized model call."""
        if not os.getenv('OPENAI_API_KEY'):
            raise ValueError('This app process has no API key configured')
        # Reserve once under a write transaction before any submission. Concurrent
        # clicks cannot create a second repair or exceed the project call budget.
        with self.store.connect() as connection:
            connection.execute('BEGIN IMMEDIATE')
            row=connection.execute('SELECT spec,data FROM projects WHERE id=?',(pid,)).fetchone()
            if not row: raise ValueError('Unknown project')
            spec,data=json.loads(row['spec']),json.loads(row['data'])
            if spec['local_only'] or not spec['allow_model'] or spec['mode']=='offline_fixture':
                raise ValueError('Remote model use is disabled for this project')
            if connection.execute("SELECT 1 FROM jobs WHERE project_id=? AND status IN ('running','queued','cancelling')",(pid,)).fetchone():
                raise ValueError('Wait for the workflow to stop before correcting advice')
            if data.get('model_repair'):
                raise ValueError('The one-shot correction was already attempted; no automatic paid retry')
            if (data.get('model_advice') or {}).get('status')!='failed':
                raise ValueError('Only rejected or failed model advice can be corrected')
            events=[json.loads(row['data']) for row in connection.execute('SELECT data FROM events WHERE project_id=? ORDER BY id',(pid,))]
            used=sum(event.get('kind')=='model_call' and event.get('status')=='running' for event in events)
            if used>=spec['max_model_calls']:
                raise ValueError('The authorized model-call budget is exhausted')
            read_ids={event['inputs']['evidence_id'] for event in events if event.get('name')=='read_evidence' and event.get('status')=='completed' and isinstance(event.get('inputs'),dict) and 'evidence_id' in event['inputs']}
            evidence=[record for record in data.get('evidence',[]) if record['id'] in read_ids]
            rejected=[]
            for event in events:
                if event.get('kind')=='model_call' and event.get('status')=='completed':
                    texts=[part['text'] for item in (event.get('outputs') or {}).get('output',[]) if item.get('type')=='message' for part in item.get('content',[]) if part.get('type')=='output_text']
                    if texts: rejected=texts
            data['model_repair']={'status':'running','reserved_calls':1,'previous_calls':used}
            connection.execute('UPDATE projects SET data=? WHERE id=?',(dump(data),pid))
        def emit(**event):
            event.setdefault('purpose','advisory_repair')
            self.store.event(pid,stage='03',**event)
        try:
            from .agent import repair_advice
            advice={**repair_advice(spec,evidence,rejected,emit),'status':'completed'}
        except Exception as exc:
            # SDK response bodies can contain request content; preserve only a safe error class.
            reason=str(exc) if isinstance(exc,ValueError) else f'Correction failed ({type(exc).__name__})'
            emit(kind='model_status',name='Advisory correction',status='failed',message=reason)
            failed={'status':'failed','reason':reason,'model':spec['model']}
            self.store.update(pid,data={'model_repair':{'status':'failed','reason':reason},
                                        'model_advice':failed,'plan':{**data.get('plan',{}),'model_advice':failed}})
            raise ValueError(reason) from None
        plan={**data.get('plan',{}),'model_advice':advice}
        artifact=self.store.artifact(pid,'03-corrected-advice.json',dump(plan).encode(),'stage_result')
        step=next(step for step in self.store.steps(pid) if step['id']=='03')
        self.store.step(pid,'03',step['status'],outputs={**(step.get('outputs') or {}),'model_advice':advice,'result_artifact':artifact['id']})
        self.store.update(pid,data={'plan':plan,'model_advice':advice,'model_repair':{'status':'completed','artifact_id':artifact['id']}})
        emit(kind='model_status',name='Advisory correction',status='completed',outputs={'artifact_id':artifact['id'],'validation':advice.get('validation')})
        return advice

    def adapters(self):
        from .registry import adapters
        return adapters()

    def export(self,pid):
        from .reports import export_project
        return export_project(self.store,pid)

    def import_measurements(self,pid,content:bytes):
        project=self.store.project(pid)
        revision=project['data'].get('research',{}).get('input_revision') or {}
        revising=revision.get('status')=='awaiting_inputs'
        if project['spec'].get('workflow_mode')=='research':
            if any(job['status'] in {'queued','running'} for job in self.store.jobs(pid)):
                raise ValueError('Wait for an input question before importing test-set measurements')
            steps=self.store.steps(pid)
            completed_pilot=any(project['data'].get('research',{}).get(key) for key in ('pilot','comparison','selection_lock')) or any(step['id']>='06' and step['status']=='completed' for step in steps)
            if completed_pilot or (any(step['id']=='04' and step['status']=='completed' for step in steps) and not revising):
                raise ValueError('The research test set is frozen; import new measurements into a new research project')
        if len(content)>10*1024*1024: raise ValueError("Assay CSV exceeds 10 MB")
        rows=list(csv.DictReader(io.StringIO(content.decode('utf-8-sig'))))
        if not rows: raise ValueError("No assay rows")
        parsed=[normalize_measurement(row).model_dump() for row in rows]
        compounds={str(x.get('compound_id',x.get('id'))) for x in self.store.project(pid)['data'].get('candidates',[])}
        if compounds and any(row['compound_id'] not in compounds for row in parsed): raise ValueError("Measured compound IDs must map to this run's candidates")
        versions=[row for row in self.store.artifacts(pid) if row['kind']=='measured_results']
        sha=digest(content)
        current=versions[-1] if versions else None
        if current and current['sha256']==sha:
            self.store.read_artifact(pid,current)
            if project['data'].get('measurements')!=parsed:
                raise ValueError('Stored measurements differ from the current immutable assay CSV; create a new project to resolve the input mismatch')
            return {'count':len(parsed),'artifact':current,'measurements':parsed,'unchanged':True}
        if any(row['sha256']==sha and row['name']=='measured-results.csv' for row in versions):
            raise ValueError('This is an older assay CSV and cannot replace the current measurements under the fixed import filename. Keep the current measurements, or create a new research project to restore this historical CSV. No measurements were changed.')
        a=self.store.artifact(pid,'measured-results.csv',content,'measured_results')
        self.store.update(pid,data={'measurements':parsed})
        if project['spec'].get('workflow_mode')=='research' and revising and a['id'] not in revision['baseline_artifact_ids']:
            changed={'artifact_id':a['id'],'kind':'measured_results','name':a['name']}
            if changed not in revision['changed_inputs']: revision['changed_inputs'].append(changed)
            questions=project['data']['research'].get('questions',[])
            for pending in questions:
                if pending.get('kind')=='research_input_revision' and pending.get('status')=='pending':
                    pending.setdefault('context',{})['changed_inputs']=revision['changed_inputs']
            self.store.research_update(pid,input_revision=revision,questions=questions)
        self.store.event(pid,kind='measurement_import',name='User experimental results',status='completed',outputs={'count':len(parsed),'artifact_id':a['id']})
        return {'count':len(parsed),'artifact':a,'measurements':parsed}

    def create_demo(self):
        pid=self.create_project({'name':'Offline workbench tour — fixture','target':'Fixture target','organism':'Example organism','mode':'offline_fixture'})
        for sid in ('01','02','03'):
            self.store.step(pid,sid,'completed',reason='Offline UI fixture — no external request or calculation',outputs={'fixture':True,'calculation_performed':False})
        self.store.step(pid,'04','blocked',reason='Fixture illustrates a receptor-review gate; no prepared structure supplied')
        self.store.update(pid,state='blocked',data={'fixture':True,'gate':{'kind':'fixture','reason':'Offline UI demonstration only'},'plan':{'strategy':'further_experiments','status':'fixture','limitations':['No scientific calculations or external model calls were made.']}})
        self.store.event(pid,kind='fixture',stage='03',name='Illustrative model-call location',status='fixture',inputs={'model':'not called'},outputs={'fixture':True,'message':'Real model requests and responses appear here after explicit model opt-in.'})
        return pid

    def create_example(self, *, allow_model: bool = False, model: str = "gpt-4.1-mini", workflow_mode: str = 'classic', max_model_calls: int = 3):
        base=Path(__file__).parent.parent/'examples'/'1iep'
        if not (base/'originals/1iep_ligand.sdf').exists(): raise ValueError('Example files unavailable; see examples/1iep/README.md')
        pid=self.create_project({'name':'1IEP · GPT-enabled public research' if allow_model else '1IEP · real local smoke calculation','target':'P00520','organism':'Mus musculus','context':'Verify software operation on curated ABL–imatinib example; not predictive validation.','intended_effect':'binding characterization','shortlist_size':3,'cpu':1,'compute_budget_seconds':1200,'local_only':not allow_model,'allow_model':allow_model,'model':model,'workflow_mode':workflow_mode,'max_model_calls':max_model_calls,'max_model_output_tokens':2200 if workflow_mode=='research' else 1500})
        for name,kind in [('1iep_ligand.sdf','reference'),('1iep_receptor.pdbqt','receptor'),('1iep_receptorH.pdb','structure')]:
            self.upload(pid,name,(base/'originals'/name).read_bytes(),kind)
        library=base/'library.sdf'
        self.upload(pid,library.name if library.exists() else '1iep_ligand.sdf',library.read_bytes() if library.exists() else (base/'originals/1iep_ligand.sdf').read_bytes(),'library')
        identity_document=json.loads((base/'identity.json').read_text())
        self.decide(pid,'identity',identity_document['identity'],'Curated AutoDock Vina 1IEP example. Bundled UniProt P00520 reference identity was retrieved and species-checked; its full-length sequence is separate from the experimental kinase-domain construct. Construct mapping and pocket suitability still require review.')
        self.store.artifact(pid,'identity.json',(base/'identity.json').read_bytes(),'identity')
        source=identity_document['source']
        self.store.update(pid,data={'example':True,'evidence':[{'id':'example-1iep','source':'RCSB PDB / AutoDock Vina tutorial','source_id':'1IEP','url':'https://www.rcsb.org/structure/1IEP','retrieved_at':'2026-09-26','evidence_type':'structure','text':'Curated ABL kinase–imatinib tutorial input; see bundled provenance and license.','location':'PDB 1IEP and official Vina basic_docking example','claim_kind':'retrieved_observation'},{'id':'example-uniprot-P00520','source':'UniProt','source_id':'P00520','url':source['url'],'retrieved_at':source['retrieved_at'],'evidence_type':'reference_target_identity','text':f"Mus musculus ABL1 reference sequence; {source['sequence_length']} residues; sequence SHA256 {source['sequence_sha256']}. Full-length reference is distinct from the 1IEP receptor construct.",'location':source['location'],'claim_kind':'retrieved_observation'}]})
        return pid
