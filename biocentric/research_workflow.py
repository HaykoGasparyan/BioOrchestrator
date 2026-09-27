"""Durable literature → test set → tool proposal → pilot → evaluation workflow."""
from __future__ import annotations

import json
import math
import re
import time
from pathlib import Path

from .store import Store, digest, dump
from .worker import Cancelled, Pipeline


class AwaitingAnswer(Exception):
    pass


def validate_answer(question: dict, answer):
    value=answer.get('choice',answer) if isinstance(answer,dict) else answer
    if question.get('kind')=='target_identity':
        if isinstance(value,dict):
            target,organism=value.get('target',''),value.get('organism','')
        else:
            pieces=re.split(r'[;,\n]',str(value),maxsplit=1)
            if len(pieces)!=2: raise ValueError('Enter a target accession and organism, for example P00520; Mus musculus')
            target,organism=pieces
        if not target.strip() or not organism.strip(): raise ValueError('Both target and organism are required')
        if len(target)>300 or len(organism)>160: raise ValueError('Target or organism is too long')
        return {'target':target.strip(),'organism':organism.strip()}
    options=question.get('options',[])
    allowed=[option.get('value',option.get('label')) if isinstance(option,dict) else option for option in options]
    if options and value not in allowed: raise ValueError('Choose one of the displayed answers')
    if not isinstance(value,str) or not value.strip() or len(value)>4000:
        raise ValueError('Provide a short answer')
    return value.strip()


def pilot_checks(pilot: dict, criteria: dict, expected_references: int = 1) -> dict:
    """No LLM may modify these measurements or promote a redock to a benchmark."""
    poses=pilot.get('poses',[])
    tops=[pose for pose in poses if pose.get('pose_number')==1]
    rmsds=[pose.get('redocking',{}).get('rmsd_angstrom') for pose in tops]
    finite=bool(tops) and all(isinstance(value,(int,float)) and not isinstance(value,bool) and math.isfinite(value) and value>=0 for value in rmsds)
    maximum=float(criteria.get('maximum',2.0))
    checks={'poses_produced':bool(poses),'all_reference_compounds_scored':len(tops)==expected_references and expected_references>0 and not pilot.get('failures'),
            'top_pose_physical_validity':bool(tops) and all(pose.get('validation',{}).get('passed') is True for pose in tops),
            'top_pose_recovery':finite and all(value<=maximum for value in rmsds)}
    return {'status':'software_pilot_passed' if all(checks.values()) else 'pilot_failed','checks':checks,
            'metrics':{'pose_count':len(poses),'reference_compounds':len(tops),'top_pose_rmsd_angstrom':rmsds},
            'criteria':criteria,'predictive_validation':'not_established','affinity_validation':'not_established',
            'comparison_scope':'One installed engine; this is not evidence that it outperforms other engines.',
            'limitations':['Pose recovery and physical validity do not establish affinity accuracy or screening enrichment.',
                           'An independently held-out assay dataset is required for predictive validation.']}


class ResearchWorkflow:
    def __init__(self,store:Store,pid:str,jid:str,owner_token=None):
        self.store,self.pid,self.jid,self.owner_token=store,pid,jid,owner_token
        self.stage='01';self.started=time.monotonic()
        self.spent=float(self.data.get('compute_seconds',0))
        self.work=store.work_directory(pid)

    @property
    def data(self): return self.store.project(self.pid)['data']
    @property
    def research(self): return self.data.get('research',{})
    @property
    def spec(self):
        original=self.store.project(self.pid)['spec']
        return {**original,**(self.research.get('target') or {}),'context':self.research.get('goal',original['context'])}
    def check(self):
        if self.owner_token and not self.store.owns_lease(self.jid,self.owner_token): raise Cancelled('Research worker lease changed')
        if self.store.cancelled(self.jid): raise Cancelled('Research cancelled; recorded results retained')
        if self.spent+time.monotonic()-self.started>self.spec['compute_budget_seconds']:
            raise ValueError('The authorized compute-time budget is exhausted; saved results remain available')
    def cancelled(self):
        try: self.check()
        except (Cancelled,ValueError): return True
        return False
    def emit(self,**event):
        self.check()
        if event.get('kind')=='model_call' and event.get('status')=='running':
            used=sum(e.get('kind')=='model_call' and e.get('status')=='running' for e in self.store.events(self.pid))
            if used>=self.spec['max_model_calls']:
                raise ValueError('The authorized model-call budget is exhausted; no further request was sent')
        event.setdefault('stage',self.stage)
        self.store.event(self.pid,job_id=self.jid,**event)
    def note(self,message,**extra):
        self.emit(kind='research_message',role='assistant',status='recorded',message=message,**extra)
    def artifact(self,kind):
        artifact=next((a for a in reversed(self.store.artifacts(self.pid)) if a['kind']==kind),None)
        if artifact:
            self.store.read_artifact(self.pid,artifact)
            return Path(artifact['path'])
        return None
    def ask(self,identifier,prompt,options=None,*,kind='choice',context=None):
        qid=identifier if context is None else identifier+'_'+digest(dump(context).encode())[:10]
        questions=self.research.get('questions',[])
        question=next((q for q in questions if q['id']==qid),None)
        if question and question['status']=='answered': return question['answer']
        if not question:
            question={'id':qid,'prompt':prompt,'options':options or [],'status':'pending','kind':kind,
                      'stage_id':self.stage,'context':context}
            questions.append(question)
            self.store.research_update(self.pid,questions=questions,phase='awaiting_answer')
            self.emit(kind='question',name=identifier,status='awaiting_answer',message=prompt,outputs=question)
        raise AwaitingAnswer(prompt)
    def stage_run(self,sid,label,fn,*,model=False):
        self.stage=sid;self.check()
        self.store.ensure_step(self.pid,sid,label,reason='Waiting for preceding research steps')
        step=next(step for step in self.store.steps(self.pid) if step['id']==sid)
        if step['status']=='completed':
            artifact_id=(step.get('outputs') or {}).get('result_artifact')
            if not artifact_id: raise ValueError(f'Stage {sid} has no immutable result')
            self.store.read_artifact(self.pid,{'id':artifact_id})
            return step.get('outputs')
        if model:
            previous=sum(e.get('kind')=='model_call' and e.get('status')=='running' and e.get('stage')==sid for e in self.store.events(self.pid))
            if previous and step['status'] in {'failed','running','waiting'}:
                action=self.ask('retry_'+sid,'The previous model step did not finish successfully. Retry within the remaining project budget?',
                                [{'label':'Retry this step','value':'retry'},{'label':'Keep results and stop','value':'stop'}],context={'attempts':previous})
                if action=='stop': raise Cancelled('Stopped at the user’s request; no model retry was made')
        self.store.research_update(self.pid,phase=label,current_stage=sid)
        self.store.step(self.pid,sid,'running',reason=label)
        self.note(label+'.')
        start=time.monotonic()
        try:
            result=fn()
            self.check()
            artifact=self.store.artifact(self.pid,f'research-{sid}.json',dump(result).encode(),'research_stage',self.jid)
            output={**result,'result_artifact':artifact['id']}
            self.store.step(self.pid,sid,'completed',reason=result.get('summary',result.get('reason',label+' completed')),outputs=output)
            self.emit(kind='stage',name=label,status='completed',duration_ms=round((time.monotonic()-start)*1000),outputs=output)
            return output
        except AwaitingAnswer as exc:
            self.store.step(self.pid,sid,'waiting',reason=str(exc));raise
        except Exception as exc:
            self.store.step(self.pid,sid,'failed',reason=str(exc));raise

    def brief(self):
        if self.spec.get('workflow_kind')!='noncovalent':
            raise ValueError('This local workflow supports noncovalent protein–small-molecule studies. Other chemistry needs a separate validated protocol.')
        if not self.research.get('target'):
            self.ask('target_identity','Which target and species should I study? Enter an accession and organism, for example P00520; Mus musculus.',kind='target_identity')
        from .engines import engine_catalog
        catalog=engine_catalog()
        self.store.research_update(self.pid,capabilities=catalog)
        return {'summary':'Research question and target recorded. Literature search comes before selecting or running scientific tools.',
                'target':self.research['target'],'goal':self.research['goal'],'capabilities':catalog}

    def literature(self):
        from .research_agent import research_literature
        result=research_literature(self.spec,self.research['capabilities'],self.emit,cancelled=self.cancelled)
        records=result.get('records',[])
        if not records: raise ValueError('No sourced literature was returned; tool selection is blocked')
        self.store.research_update(self.pid,literature=result,sources=records)
        self.note(result.get('summary','Literature review recorded with citations.'))
        return {**result,'sources':records}

    def evidence(self):
        from .retrieval import RetrievalClient
        client=RetrievalClient(self.store.directory(self.pid)/'cache',allow_external=not self.spec['local_only'],
                               emit=lambda name,inputs,outputs:self.emit(kind='tool_call',name=name,status=outputs.get('status','completed'),inputs=inputs,outputs=outputs))
        # Keep an interactive public-source outage bounded. Missing coverage is
        # recorded in the dossier; it must not silently become assay evidence.
        client.timeout=8.0
        client.max_attempts=1
        supplied=next((decision['value'] for decision in reversed(self.store.decisions(self.pid)) if decision['kind']=='identity'),None)
        if supplied:
            identity={**supplied,'verification':'user-supplied curated identity; construct mapping requires review'}
        else:
            resolution=client.resolve_target(self.spec['target'],self.spec['organism'])
            candidates=resolution.get('candidates',[])
            if not candidates:
                self.store.update(self.pid,data={'target_candidates':[]})
                raise ValueError('No species-matched target resolved; create a new research question with a precise accession')
            if len(candidates)>1:
                chosen=self.ask('select_identity','Several target records match. Which accession should this run use?',
                                [{'label':f"{row['accession']} · {row.get('name','')}",'value':row['accession']} for row in candidates],context=candidates)
                identity=next(row for row in candidates if row['accession']==chosen)
            else: identity=candidates[0]
        self.store.update(self.pid,data={'target_identity':identity})
        dossier=client.dossier(identity,'')
        dossier['evidence']=self.data.get('evidence',[])+dossier.get('evidence',[])
        self.store.update(self.pid,data=dossier)
        return {'summary':f"Collected {len(dossier['evidence'])} source records and {len(dossier.get('assays',[]))} assay measurements; original endpoints are preserved.",
                'identity':identity,'evidence_count':len(dossier['evidence']),'assay_count':len(dossier.get('assays',[])),
                'limitations':dossier.get('limitations',[])}

    def testset(self):
        from .service import Service
        from .testsets import build_testset
        missing=[kind for kind in ('receptor','reference','library','structure') if not self.artifact(kind)]
        if missing:
            self.ask('scientific_inputs','Please upload the missing scientific inputs ('+', '.join(missing)+') using Files, then continue.',
                     [{'label':'Inputs uploaded; check again','value':'check'}],context={'missing':missing,'artifacts':[a['id'] for a in self.store.artifacts(self.pid)]})
            missing=[kind for kind in missing if not self.artifact(kind)]
            if missing: raise ValueError('Required inputs are still missing: '+', '.join(missing))
        manifest=build_testset(Service(self.store.root).snapshot(self.pid),seed=self.spec['seed'],organism=self.spec['organism'])
        self.store.research_update(self.pid,testset=manifest)
        self.note('The test-set manifest and pilot criteria are frozen before any docking score is computed. Untested compounds are not labeled inactive.')
        return manifest

    def choose_tools(self):
        from .research_agent import recommend_tools
        result=recommend_tools(self.spec,self.research['literature'],self.research['testset'],
                               self.research['capabilities'],self.emit,cancelled=self.cancelled)
        self.store.research_update(self.pid,tool_plan=result)
        self.note(result.get('reason',result.get('summary','Tool proposal recorded.')))
        return result

    def protocol(self):
        from .science import derive_box,inspect_structure
        from .testsets import validate_testset
        if self.spec.get('workflow_kind')!='noncovalent':
            raise ValueError('Noncovalent preparation cannot be used for this workflow kind')
        validate_testset(self.research['testset'])
        inputs={kind:self.artifact(kind) for kind in ('receptor','reference','library','structure')}
        if any(path is None for path in inputs.values()): raise ValueError('The pilot requires all four curated scientific inputs')
        if inputs['structure'].suffix.lower()!='.pdb': raise ValueError('Provide an inspected matching PDB export for the curated pilot')
        inspections={kind:inspect_structure(inputs[kind]) for kind in ('receptor','structure')}
        blockers=[item for inspection in inspections.values() for item in inspection.get('routine_screening_blockers',[])]
        if blockers: raise ValueError('Unsupported receptor components: '+'; '.join(blockers))
        result={'engine':self.research['tool_plan'].get('docking_engine'),'seed':self.spec['seed'],'cpu':self.spec['cpu'],
                'exhaustiveness':4,'n_poses':3,'box':derive_box(inputs['reference']),'box_padding_angstrom':4.0,
                'timeout_seconds':min(120,self.spec['compute_budget_seconds']),'require_posebusters':True,'run_plip':True,
                'receptor_pdb':str(inputs['structure']),'input_hashes':{kind:digest(path.read_bytes()) for kind,path in inputs.items()},
                'testset_sha256':self.research['testset']['manifest_sha256'],'validation':'pilot only; no established predictive validity',
                'ph_assumption':7.4,'inspection':inspections,'criteria':self.research['testset']['criteria'],
                'answers_to_model_questions':self.research.get('model_answers',[])}
        self.store.update(self.pid,data={'protocol':result})
        self.store.research_update(self.pid,protocol=result)
        return result

    def model_questions(self, questions, prefix):
        answers=[]
        for index,question in enumerate(questions):
            context={'question':question}
            # Generation zero deliberately preserves legacy question IDs. After
            # an explicit input revision, old answers cannot approve new inputs.
            revision=self.research.get('input_revision_version',0)
            if revision:
                context['input_revision_version']=revision
                context['testset_sha256']=(self.research.get('testset') or {}).get('manifest_sha256')
            controls=[{'label':'Continue with frozen inputs','value':'continue_frozen'},
                      {'label':'Revise inputs before continuing','value':'revise_inputs'}]
            prompt=question['question']+'\nChoose whether to continue with the current frozen inputs or pause to revise them. A model question cannot silently change the protocol.'
            answer=self.ask(f'{prefix}_{index}',prompt,controls,
                            kind='model_clarification_control',context=context)
            if answer not in {'continue_frozen','revise_inputs'}:
                # Old free-form answers are retained as evidence, never treated
                # as implicit permission or executable parameter changes.
                answer=self.ask(f'{prefix}_{index}_confirm',prompt,controls,
                                kind='model_clarification_control',context={**context,'legacy_answer':answer})
            if answer!='continue_frozen':
                raise AwaitingAnswer('Waiting for revised inputs and explicit permission to rebuild the proposal; no pilot will run')
            answers.append({'question':question['question'],'answer':answer,
                            'input_revision_version':revision,'original_model_options':question.get('options',[])})
        self.store.research_update(self.pid,model_answers=self.research.get('model_answers',[])+[
            answer for answer in answers if answer not in self.research.get('model_answers',[])])
        return answers

    def legacy(self):
        pipeline=Pipeline(self.store,self.pid,self.jid,self.owner_token);pipeline.stage=self.stage
        pipeline.spec=self.spec
        # Every child shares this worker attempt's clock and persisted prior
        # usage; creating the next pilot must not reset its cancellation budget.
        pipeline.started=self.started
        pipeline.spent=self.spent
        return pipeline

    def comparison_candidates(self, protocol):
        """Bounded operational candidates; unavailable/native-unprepared choices stay visible."""
        from .science import _version
        capabilities={row['id']:row for row in self.research.get('capabilities',[])}
        engine=protocol['engine']
        proposed=[engine]+[row['engine'] for row in self.research['tool_plan'].get('candidates',[]) if row.get('role')=='docking']
        eligible=list(dict.fromkeys(item for item in proposed if item==engine or (
            item in {'vina','gnina','unidock'} and capabilities.get(item,{}).get('available') is True)))
        parameters=[{**protocol,'engine':item,'exhaustiveness':amount} for item in eligible for amount in ((4,8) if item=='vina' else (4,))]
        candidates=[]
        for settings in parameters:
            from .engines import GninaParameters,UniDockParameters,AutoDockParameters
            if settings['engine']=='gnina':
                settings['device']='gpu'
                settings['native_parameters']=GninaParameters.model_validate({key:settings[key] for key in
                    ('cpu','seed','exhaustiveness','timeout_seconds','device')} | {'num_modes':settings['n_poses'],'box':{key:settings['box'][key] for key in ('center','size')}}).model_dump(mode='json')
            elif settings['engine']=='unidock':
                settings['search_mode']='balance'
                settings['native_parameters']=UniDockParameters.model_validate({key:settings[key] for key in
                    ('seed','timeout_seconds','search_mode')} | {'num_modes':settings['n_poses'],'box':{key:settings['box'][key] for key in ('center','size')}}).model_dump(mode='json')
            elif settings['engine']=='autodock4':
                settings['native_parameters']=AutoDockParameters.model_validate({'seed':settings['seed'],'timeout_seconds':settings['timeout_seconds']}).model_dump(mode='json')
            capability=capabilities.get(settings['engine'],{})
            executable=capability.get('executable')
            provenance={'executable':executable,'executable_sha256':digest(Path(executable).read_bytes()) if executable and Path(executable).is_file() else None,
                        'dependencies':{name:_version(name) for name in ('rdkit','meeko','posebusters','vina')},
                        'source_ids':self.research['tool_plan'].get('citations',[]),
                        'validation_status':capability.get('validation_status','runtime execution still pending')}
            adapter_files=('science.py','worker.py') if settings['engine']=='vina' else ('science.py','engines.py','engine_workflow.py','native_preparation.py')
            provenance['adapter_code_sha256']={str(Path(__file__).with_name(name)):digest(Path(__file__).with_name(name).read_bytes()) for name in adapter_files}
            runtime_manifest=Path(__file__).resolve().parents[1]/'.tools'/'linux'/'gpu-runtime-manifest.json'
            if settings['engine'] in {'gnina','unidock'} and runtime_manifest.is_file():
                provenance['runtime_manifest']={'path':str(runtime_manifest),'sha256':digest(runtime_manifest.read_bytes()),
                                                'record':json.loads(runtime_manifest.read_text())}
            ident=settings['engine']+'-exhaustiveness-'+str(settings['exhaustiveness'])
            if settings['engine']=='unidock': ident='unidock-'+settings['search_mode']
            candidates.append({'protocol_id':ident,'engine_id':settings['engine'],
                               'metric':'top_1_symmetry_aware_rmsd_angstrom','units':'angstrom','direction':'lower_is_better',
                               'parameters':settings,'input_hashes':settings['input_hashes'],'provenance':provenance,
                               'training_overlap':'not_applicable' if settings['engine'] in {'vina','autodock4','unidock'} else 'unknown'})
        excluded=[]
        for row in self.research['tool_plan'].get('candidates',[]):
            if row.get('role')=='docking' and row['engine'] not in eligible:
                capability=capabilities.get(row['engine'],{})
                excluded.append({'engine_id':row['engine'],'reason':'; '.join(capability.get('blockers',[])) or
                                 'Not selected for this bounded pilot; native prepared inputs and a reviewed protocol are required before inclusion.'})
        return candidates,excluded

    def check_candidate_runtime(self, candidate):
        from .science import sha256,_version
        provenance=candidate['provenance']
        executable=provenance.get('executable')
        if executable and provenance.get('executable_sha256')!=sha256(Path(executable)):
            raise ValueError('The scientific executable changed after candidate freeze')
        if any(sha256(Path(path))!=expected for path,expected in provenance.get('adapter_code_sha256',{}).items()):
            raise ValueError('A scientific adapter changed after candidate freeze')
        runtime=provenance.get('runtime_manifest')
        if runtime and sha256(Path(runtime['path']))!=runtime['sha256']:
            raise ValueError('The GPU runtime manifest changed after candidate freeze')
        if any(_version(name)!=version for name,version in provenance.get('dependencies',{}).items()):
            raise ValueError('A scientific preprocessing dependency changed after candidate freeze')

    def _comparison_file(self, name):
        frozen=self.research['testset']['manifest_sha256']
        return self.work/'comparison'/frozen/name

    def _comparison_artifact(self, name, value):
        artifact=self.store.artifact(self.pid,name,dump(value).encode(),'research_comparison',self.jid)
        return artifact['id']

    def pilot(self):
        from .comparative import freeze_comparison,validate_comparison
        from .science import ScienceError
        from .worker import Gate
        self.model_questions(self.research['tool_plan'].get('questions',[]),'tool_question')
        if not self.research['tool_plan'].get('docking_engine'):
            raise ValueError('The proposal found no executable docking engine. Configure a compatible local engine before starting a new research run.')
        protocol=self.protocol()
        candidates,excluded=self.comparison_candidates(protocol)
        criteria={'mode':'pose_software','maximum_pose_rmsd_angstrom':self.research['testset']['criteria']['pose_recovery']['maximum'],
                  'scientist_reviewed':True,'seed':self.spec['seed']}
        review_context={'protocol':protocol,'candidate_protocols':candidates,'comparison_criteria':criteria,'excluded_methods':excluded}
        self.store.step(self.pid,self.stage,'running',outputs=review_context)
        action=self.ask('protocol_review','Review the candidate protocols, receptor, observed binding site, chemistry and pH assumption. Run these bounded comparative pose pilots?',
                        [{'label':'Reviewed — run the pilots','value':'approve'},{'label':'Keep the plan; stop for input review','value':'stop'}],context=review_context)
        if action=='stop': raise Cancelled('Stopped for receptor and chemistry review at the user’s request')
        path=self._comparison_file('frozen-comparison.json')
        if path.exists():
            frozen=json.loads(path.read_text())
            validate_comparison(self.research['testset'],frozen)
            # A rerun may not silently substitute a binary or change a candidate.
            from .comparative import CandidateProtocol
            current=sorted([CandidateProtocol.model_validate(row).model_dump() for row in candidates],key=lambda row:row['protocol_id'])
            if current!=frozen['candidates']:
                raise ValueError('Candidate protocol/runtime provenance changed after freeze; start a new comparison')
        else:
            frozen=freeze_comparison(self.research['testset'],candidates,criteria,path)
        frozen_artifact=self._comparison_artifact('frozen-comparison.json',frozen)
        self.store.research_update(self.pid,comparison=frozen,comparison_artifact=frozen_artifact,comparison_exclusions=excluded)
        completed=self.research.get('comparison_results',{})
        if any(row['engine_id']=='vina' for row in frozen['candidates']) and not self.data.get('preparation'):
            self.legacy().preparation()
        for candidate in frozen['candidates']:
            self.check()
            ident=candidate['protocol_id']
            if ident in completed:
                # Durable per-protocol checkpoints are read from immutable artifacts.
                completed[ident]=json.loads(self.store.read_artifact(self.pid,{'id':completed[ident]['artifact_id']})) | {'artifact_id':completed[ident]['artifact_id']}
                self.emit(kind='checkpoint',name=ident,status='reused',message='Completed comparative pilot retained; no new engine call')
                continue
            self.check_candidate_runtime(candidate)
            settings=candidate['parameters']
            self.store.update(self.pid,data={'protocol':settings})
            self.note(f"Running {ident} on the frozen reference; its scores cannot establish affinity accuracy.")
            started=time.monotonic()
            try:
                if candidate['engine_id']=='vina':
                    pipeline=self.legacy()
                    pipeline.work=self.work/'comparative-pilots'/frozen['comparison_sha256']/ident
                    pipeline.work.mkdir(parents=True,exist_ok=True)
                    raw=pipeline.pilot()
                else:
                    raw=self.native_calculation(candidate['engine_id'],'pilot',settings)
            except (ScienceError,Gate,ValueError) as exc:
                if getattr(exc,'code',None)=='cancelled' or getattr(exc,'kind',None)=='budget': raise
                raw={'poses':[],'failures':[{'reason':str(exc),'code':getattr(exc,'code',getattr(exc,'kind','calculation_failed'))}],
                     'status':'failed'}
            raw={**raw,'engine':candidate['engine_id'],'protocol_id':ident,
                 'elapsed_seconds':time.monotonic()-started,'comparison_sha256':frozen['comparison_sha256'],
                 'testset_sha256':protocol['testset_sha256'],'protocol_sha256':digest(dump(settings).encode())}
            raw['artifact_id']=self._comparison_artifact('pilot-'+ident+'.json',raw)
            completed[ident]=raw
            self.store.research_update(self.pid,comparison_results=completed)
        # Preserve legacy public pilot fields until deterministic selection in step 07.
        first=completed[frozen['candidates'][0]['protocol_id']]
        result={**first,'comparison':frozen,'comparison_results':completed,'excluded_methods':excluded,
                'summary':f"Completed {len(completed)} frozen pose protocols. Selection is pending; affinity and screening remain unvalidated."}
        self.store.research_update(self.pid,pilot=result)
        return result

    def lock_comparison(self):
        from .comparative import select_and_lock,evaluate_locked,validate_comparison,_verify,content_hash
        frozen=self.research.get('comparison')
        if not frozen:
            return None  # Readable historical single-pilot records retain their old path.
        validate_comparison(self.research['testset'],frozen)
        results=[]
        references=self.research['testset']['pose_records']
        for candidate in frozen['candidates']:
            stored=self.research['comparison_results'][candidate['protocol_id']]
            raw=json.loads(self.store.read_artifact(self.pid,{'id':stored['artifact_id']}))
            tops=[row for row in raw.get('poses',[]) if row.get('pose_number')==1]
            # Current curated input contract permits exactly one observed reference.
            if len(references)!=1:
                raise ValueError('Comparative pose pilots require one exact observed reference per protocol')
            good=len(tops)==1 and not raw.get('failures')
            value=tops[0].get('redocking',{}).get('rmsd_angstrom') if good else None
            good=good and isinstance(value,(int,float)) and not isinstance(value,bool) and math.isfinite(value) and value>=0
            results.append({'protocol_id':candidate['protocol_id'],'compound_id':references[0]['id'],'partition':'tuning',
                            'status':'completed' if good else 'failed','value':float(value) if good else None,
                            'metric':candidate['metric'],'units':candidate['units'],'direction':candidate['direction'],
                            'physical_validity':tops[0].get('validation',{}).get('passed') is True if good else None,
                            'reason':None if good else 'Reference did not produce one finite top-pose RMSD without execution failures; inspect immutable raw result.',
                            'job_id':self.jid,'artifact_sha256':digest(dump(raw).encode()),'comparison_sha256':frozen['comparison_sha256']})
        path=self._comparison_file('selection-lock.json')
        if path.exists():
            lock=json.loads(path.read_text());_verify(lock,'selection_sha256')
            from .comparative import _results
            normalized=_results(results,frozen,validate_comparison(self.research['testset'],frozen),'tuning',frozen['candidates'])
            if lock['comparison_sha256']!=frozen['comparison_sha256'] or lock['tuning_results_sha256']!=content_hash(normalized):
                raise ValueError('Frozen selection cannot be changed after tuning or final evaluation')
        else:
            lock=select_and_lock(self.research['testset'],results,frozen,path)
        self.store.research_update(self.pid,selection_lock=lock,selection_lock_artifact=self._comparison_artifact('selection-lock.json',lock))
        if lock['status']=='locked':
            selected=lock['selected_protocol']
            self.store.update(self.pid,data={'protocol':selected['parameters']})
            plan={**self.research['tool_plan'],'docking_engine':selected['engine_id']}
            pilot=self.research['comparison_results'][selected['protocol_id']]
            self.store.research_update(self.pid,tool_plan=plan,protocol=selected['parameters'],pilot=pilot)
            evaluation_path=self._comparison_file('heldout-evaluation.json')
            if evaluation_path.exists():
                evaluation=json.loads(evaluation_path.read_text());_verify(evaluation,'evaluation_sha256')
                if evaluation['selection_sha256']!=lock['selection_sha256']:
                    raise ValueError('Held-out evaluation belongs to another selection lock')
            else:
                evaluation=evaluate_locked(self.research['testset'],[],frozen,lock,evaluation_path)
        else:
            evaluation={'status':'not_evaluated','predictive_validation':'not_established','reason':'No candidate passed tuning checks; no held-out data were consumed.'}
        self.store.research_update(self.pid,heldout_evaluation=evaluation,
                                   heldout_evaluation_artifact=self._comparison_artifact('heldout-evaluation.json',evaluation))
        return {'selection_lock':lock,'heldout_evaluation':evaluation,'comparison':frozen,
                'comparison_results':lock['candidate_results']}

    def external_pilot(self,engine,protocol):
        result=self.native_calculation(engine,'pilot',protocol)
        if not result.get('poses'): raise ValueError('The native pilot returned no chemically verified poses')
        self.store.research_update(self.pid,pilot=result)
        return result

    def native_calculation(self,engine,phase,protocol):
        """Only registered, hash-bound native inputs can reach an engine runner."""
        from .engine_workflow import plan_native_manifest,run_native_manifest,inspect_native_resume
        kind='native_engine_'+phase
        manifest=None
        for artifact in reversed(self.store.artifacts(self.pid)):
            if artifact['kind']!=kind: continue
            body=json.loads(self.store.read_artifact(self.pid,artifact))
            if body.get('engine_id')==engine and body.get('phase')==phase and body.get('protocol_sha256')==digest(dump(protocol).encode()) and body.get('testset_sha256')==protocol['testset_sha256']:
                manifest=artifact;break
        if not manifest and engine in {'gnina','unidock','gnina_affinity','biocentric_affinity'}:
            from .engine_workflow import auto_prepare_native_manifest
            manifest=auto_prepare_native_manifest(self.store,self.pid,engine,phase,protocol,job_id=self.jid,
                                                  emit=self.legacy().science_event,cancelled=self.cancelled)
        if not manifest:
            self.ask('native_inputs_'+phase,
                     f'{engine} needs its native prepared inputs. Upload the files and a {kind} JSON manifest using Advanced native engine inputs, then continue.',
                     [{'label':'Native inputs uploaded; inspect them','value':'check'}],context={
                         'missing':kind,'engine':engine,'phase':phase,
                         'testset_sha256':protocol['testset_sha256'],'protocol_sha256':digest(dump(protocol).encode()),
                         'artifacts':[a['id'] for a in self.store.artifacts(self.pid)]})
            raise ValueError('The native input manifest is still missing')
        parent=self.work/'native'/phase/manifest['sha256'][:16]
        attempt_key=engine+':'+phase+':'+manifest['sha256']
        saved=self.research.get('native_attempts',{}).get(attempt_key) if engine=='biocentric_affinity' else None
        resume_plan=None
        if saved:
            reviewed=next((item for item in self.store.artifacts(self.pid) if item['id']==saved['reviewed_plan_artifact_id']),None)
            if not reviewed or reviewed['kind']!='native_engine_review':
                raise ValueError('The saved native attempt has no registered reviewed plan; preserve its outputs.')
            resume_plan=json.loads(self.store.read_artifact(self.pid,reviewed))
            output=Path(resume_plan['output_dir'])
            if output.parent.resolve()!=parent.resolve() or not re.fullmatch(r'attempt-[1-9][0-9]*',output.name):
                raise ValueError('Saved native resume output is outside the original reviewed experiment.')
        else:
            attempts=[int(path.name.split('-')[1]) for path in parent.glob('attempt-*')
                      if re.fullmatch(r'attempt-[1-9][0-9]*',path.name)]
            output=parent/f'attempt-{max(attempts,default=0)+1}'
        kwargs={'expected_engine':engine,'expected_phase':phase,
                'expected_testset_sha256':protocol['testset_sha256'],
                'expected_protocol_sha256':digest(dump(protocol).encode())}
        resume_kwargs={'resume_plan':resume_plan} if resume_plan is not None else {}
        plan=plan_native_manifest(self.store,self.pid,manifest,output,**kwargs,**resume_kwargs)
        self.check()
        parameters=plan['engine_plan']['parameters']
        if phase in {'pilot','screening'} and protocol.get('native_parameters') and parameters!=protocol['native_parameters']:
            raise ValueError('Native parameters differ from the frozen comparative candidate')
        if parameters.get('cpu',1)>self.spec['cpu']:
            raise ValueError('The native manifest requests more CPUs than the project budget')
        if parameters.get('timeout_seconds',120)>self.spec['compute_budget_seconds']:
            raise ValueError('The native timeout exceeds the project compute budget')
        if phase=='pilot':
            reference=plan['resolved_artifacts'].get('reference')
            if not reference or reference['sha256']!=protocol['input_hashes']['reference']:
                raise ValueError('The native pilot must use the frozen observed reference ligand')
        self.store.research_update(self.pid,engine_preflight=plan)
        self.store.step(self.pid,self.stage,'running',outputs={'protocol':protocol,'engine_preflight':plan})
        if plan.get('blockers'):
            raise ValueError('Native engine preflight blocked: '+'; '.join(plan['blockers']))
        action=self.ask('native_review_'+phase,
                        f'Review the {engine} input hashes, preprocessing attestations, exact command and resource limits. Run this {phase} calculation?',
                        [{'label':'Reviewed — run this calculation','value':'approve'},{'label':'Keep the proposal and stop','value':'stop'}],context=plan)
        if action=='stop': raise Cancelled('Stopped before the native calculation at the user’s request')
        retry_failed=False
        resume_state=None
        if resume_plan is not None:
            resume_state=inspect_native_resume(plan)
            if resume_state['failed']:
                action=self.ask('native_resume_'+phase,
                                f"This reviewed BioCentric run retains {resume_state['completed']} completed, {resume_state['failed']} failed and {resume_state['missing']} unfinished pairs. Resume unfinished pairs, or explicitly retry the recorded failures too? Original evidence will be archived.",
                                [{'label':'Resume unfinished pairs; retain recorded failures','value':'resume'},
                                 {'label':'Retry failed and unfinished pairs','value':'retry_failed'},
                                 {'label':'Keep all results and stop','value':'stop'}],
                                context={'reviewed_plan_sha256':reviewed['sha256'],'checkpoints':resume_state})
                if action=='stop': raise Cancelled('Stopped before retrying the native calculation')
                retry_failed=action=='retry_failed'
            self.note(f"Resuming the reviewed BioCentric request: {resume_state['completed']} completed pairs retained, {resume_state['missing']} unfinished pairs, {resume_state['failed']} recorded failures.")
        elif engine=='biocentric_affinity':
            reviewed=self.store.artifact(self.pid,'reviewed-biocentric-workflow-plan.json',dump(plan).encode(),'native_engine_review',self.jid)
            attempts={**self.research.get('native_attempts',{}),attempt_key:{'reviewed_plan_artifact_id':reviewed['id'],
                       'manifest_sha256':manifest['sha256'],'output_dir':str(output),'job_id':self.jid}}
            self.store.research_update(self.pid,native_attempts=attempts)
        result=run_native_manifest(self.store,self.pid,manifest,output,**kwargs,job_id=self.jid,
                                   emit=self.legacy().science_event,cancelled=self.cancelled,
                                   **({'resume_plan':resume_plan,'reviewed_resume_state':resume_state,
                                       'retry_failed_pairs':retry_failed} if resume_plan is not None else {}))
        if result.get('status')!='completed':
            self.store.research_update(self.pid,native_failure=result)
            raise ValueError('Native calculation did not complete: '+str(result.get('failures') or result.get('blockers') or result['status']))
        for pose in result.get('poses',[]):
            pose['job_id']=self.jid
            pose['artifact_ids']=[pose['pose_artifact_id']]
        return result

    def verify(self):
        from .research_agent import review_pilot
        comparative=self.lock_comparison()
        checks=pilot_checks(self.research['pilot'],self.research['testset']['criteria']['pose_recovery'],
                            expected_references=len(self.research['testset']['pose_records']))
        if comparative:
            checks={**checks,**comparative,'comparison_scope':comparative['selection_lock']['comparison_scope']}
            if comparative['selection_lock']['status']!='locked': checks['status']='pilot_failed'
        self.store.research_update(self.pid,verification=checks)
        self.note(f"Measured pilot result: {checks['status'].replace('_',' ')}. These checks establish neither binding-affinity accuracy nor screening enrichment.")
        # Send measured summaries, not duplicated atom inspections and raw pose
        # reports for every candidate. Full evidence remains in immutable artifacts.
        measured=self.research['pilot']
        review_checks={key:value for key,value in checks.items() if key not in {'comparison','comparison_results','selection_lock'}}
        if comparative:
            lock=comparative['selection_lock']
            review_checks['selection_lock']={key:lock[key] for key in ('selected_protocol_id','selection_sha256','selection_scope','status','reason','candidate_results')}
            review_checks['comparison_sha256']=comparative['comparison']['comparison_sha256']
        review_evidence={key:measured.get(key) for key in ('engine','protocol_id','protocol_sha256','testset_sha256','failures','artifact_id')}
        review_evidence['poses']=[{key:pose.get(key) for key in ('compound_id','pose_number','score','score_units','engine','engine_version','redocking')} |
                                 {'physical_validity_passed':pose.get('validation',{}).get('passed')} for pose in measured.get('poses',[])]
        review_evidence['verification']=review_checks
        review=review_pilot(self.spec,self.research['literature'],self.research['testset'],self.research['tool_plan'],
                           review_evidence,self.emit,cancelled=self.cancelled)
        result={**checks,'model_review':review}
        self.store.research_update(self.pid,verification=result)
        self.note(review.get('summary','Pilot interpretation recorded.'))
        return result

    def production(self):
        verification=self.research['verification']
        self.model_questions((verification.get('model_review') or {}).get('questions',[]),'pilot_question')
        if verification['status']!='software_pilot_passed':
            self.ask('failed_pilot','The pilot failed its frozen checks. The larger run is blocked. Keep these results for revision?',
                     [{'label':'Keep results and stop','value':'stop'}],context=verification)
            raise Cancelled('Pilot checks failed; larger computation was not run')
        lock=self.research.get('selection_lock')
        if lock:
            from .comparative import validate_comparison,_verify
            validate_comparison(self.research['testset'],self.research['comparison'])
            _verify(lock,'selection_sha256')
            if lock['status']!='locked' or lock['comparison_sha256']!=self.research['comparison']['comparison_sha256'] or self.data['protocol']!=lock['selected_protocol']['parameters'] or self.research['tool_plan']['docking_engine']!=lock['selected_protocol']['engine_id']:
                raise ValueError('Production must use the locked protocol and engine without retuning')
            self.check_candidate_runtime(lock['selected_protocol'])
        recommendation=(verification.get('model_review') or {}).get('recommendation')
        if recommendation in {'revise','need_data'}:
            self.ask('review_pilot','The advisor recommends revision or more data. Preserve these results before a new protocol is designed?',
                     [{'label':'Keep results and stop','value':'stop'}],context=verification)
            raise Cancelled('Advisor requested revision; no larger computation was run')
        action=self.ask('exploratory_execution','The pose-recovery pilot passed, but screening and affinity remain unvalidated. Run the supplied library as an explicitly exploratory calculation?',
                        [{'label':'Run the exploratory library','value':'run_exploratory'},{'label':'Stop and obtain measured validation data','value':'stop_for_data'}],
                        context={'protocol':self.data['protocol'],'verification':verification})
        if action=='stop_for_data': raise Cancelled('Waiting for measured validation data; no larger computation was run')
        self.store.update(self.pid,data={'result_status':'exploratory','benchmark':{'status':'exploratory','limitations':verification['limitations']}})
        pipeline=self.legacy()
        # An affinity-input question must not rerun completed docking on resume.
        docking=self.research.get('production_docking')
        engine=self.research['tool_plan']['docking_engine']
        if docking:
            docking=json.loads(self.store.read_artifact(self.pid,{'id':docking['result_artifact']}))
        if not docking:
            if engine=='vina':
                screening=pipeline.screening()
            else:
                native=self.native_calculation(engine,'screening',self.data['protocol'])
                self.store.update(self.pid,data={'screening':native})
                screening={'count':len(native['poses']),'native_result':native}
            ranking=pipeline.prioritization()
            docking={'screening':screening,'ranking':ranking}
            artifact=self.store.artifact(self.pid,'research-production-docking.json',dump(docking).encode(),'research_stage',self.jid)
            self.store.research_update(self.pid,production_docking={**docking,'result_artifact':artifact['id']})
        screening,ranking=docking['screening'],docking['ranking']
        affinity=self.research['tool_plan'].get('affinity_engine')
        affinity_result=self.native_calculation(affinity,'affinity',self.data['protocol']) if affinity else {
                'requested_engine':None,'status':'not_run',
                'reason':'No configured affinity method with the required inputs was selected. Docking scores are not measured affinities.'}
        result={'summary':f"Completed exploratory {engine} screening and selected {ranking['selected']} candidates for review.",
                'screening':screening,'ranking':ranking,'affinity':affinity_result,
                'metrics':{'poses':screening['count'],'selected_candidates':ranking['selected']},'limitations':verification['limitations']}
        self.store.research_update(self.pid,production=result)
        return result

    def report(self):
        result=self.legacy().handoff()
        self.store.research_update(self.pid,phase='completed')
        self.note('The exploratory run is complete. The cards preserve the literature, tool proposal, frozen pilot, measured checks, decisions, and downloadable results. Affinity and screening predictive validity remain unestablished.')
        return result

    def run(self):
        for sid,label,method,model in [
            ('01','Understand the research question',self.brief,False),
            ('02','Search literature and tool documentation',self.literature,True),
            ('03','Resolve target evidence and measured controls',self.evidence,False),
            ('04','Construct and freeze the small test set',self.testset,False),
            ('05','Propose tools from the literature and test set',self.choose_tools,True),
            ('06','Freeze protocols and run comparative pilots',self.pilot,False),
            ('07','Lock the protocol and review independent evaluation',self.verify,True),
            ('08','Run the approved exploratory calculation',self.production,False),
            ('09','Collect results and experimental follow-up',self.report,False)]:
            self.stage_run(sid,label,method,model=model)
