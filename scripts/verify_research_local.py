"""Local scientific integration smoke test; all advisors and retrieval are mocked.

Run: .venv-linux/bin/python -m scripts.verify_research_local --output runs/linux-research-new
This spends no API credit. It runs actual molecular preparation, docking and
PoseBusters, preserves optional-tool availability, and
answers only the bundled software-test approvals inside an isolated test store.
"""
import argparse
import json
import time
from pathlib import Path
from unittest.mock import patch

from biocentric.research_worker import execute
from biocentric.research_workflow import ResearchWorkflow
from biocentric.service import Service


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('runs/linux-research-conversation-verification'))
    parser.add_argument('--engines', nargs='+', choices=['vina','gnina','unidock'], default=['vina'])
    parser.add_argument('--affinity', choices=['gnina_affinity','biocentric_affinity'])
    parser.add_argument('--pause-before-screening', action='store_true', help='Leave the actual exploratory review pending for browser continuation')
    args = parser.parse_args()
    root = (args.output / 'store').resolve()
    report_path = args.output.resolve() / 'local-science-verification.json'
    if report_path.exists():
        raise SystemExit('Choose a new --output directory; historical verification records are immutable.')
    started = time.monotonic()
    answered = []
    service = Service(root)
    pid = service.create_research_project('LOCAL INTEGRATION TEST: mocked advisors; actual 1IEP Vina calculations.', use_public_example=True)
    jid = service.start_research(pid, spawn=False)

    def literature(*args, **kwargs):
        return {'status':'test_double','records':[{'id':'TEST-ONLY','title':'Mock advisor; no literature claim','url':None}],
                'summary':'OFFLINE TEST DOUBLE: no literature search and no model request.', 'search_calls':0}

    def proposal(spec, sources, testset, capabilities, emit, **kwargs):
        return {'status':'test_double','docking_engine':args.engines[0],'affinity_engine':args.affinity,'questions':[],
                'candidates':[{'engine':engine,'role':'docking','reason':'Synthetic advisor request for actual local integration testing','citations':[]} for engine in args.engines],
                'reason':'TEST DOUBLE: exercise the installed local Vina route; no comparison claim.',
                'criteria_sha256':testset['criteria_sha256'],'dataset_sha256':testset['dataset_sha256']}

    def evidence(self):
        identity=next((row['value'] for row in reversed(self.store.decisions(self.pid)) if row['kind']=='identity'),None)
        if identity is None:
            raise ValueError('The local verification requires the bundled reviewed identity record')
        self.store.update(self.pid,data={'target_identity':identity})
        return {'summary':'TEST DOUBLE: retain the actual bundled curated identity and reference evidence; no remote retrieval.',
                'identity':identity}


    def review(*args, **kwargs):
        return {'status':'test_double','recommendation':'proceed_exploratory','questions':[],
                'summary':'TEST DOUBLE: exercise the deterministic pilot gate; no model interpretation.'}

    with patch('biocentric.research_agent.research_literature', literature), \
         patch('biocentric.research_agent.recommend_tools', proposal), \
         patch('biocentric.research_agent.review_pilot', review), \
         patch.object(ResearchWorkflow, 'evidence', evidence):
        for _ in range(16):
            execute(root, pid, jid)
            snapshot = service.snapshot(pid)
            if snapshot['project']['state'] == 'completed':
                break
            pending = [q for q in snapshot['data']['research']['questions'] if q['status']=='pending']
            if len(pending)!=1:
                raise RuntimeError(json.dumps(snapshot['jobs'], indent=2))
            question = pending[0]
            if question['id'].startswith(('protocol_review_','native_review_')):
                answer='approve'
            elif question['id'].startswith('exploratory_execution_'):
                if args.pause_before_screening:
                    handoff={'status':'awaiting_browser_review','project_id':pid,'job_id':jid,'root':str(root),
                             'question_id':question['id'],'completed_stages':[s['id'] for s in snapshot['steps'] if s['status']=='completed'],
                             'advisor':'MOCKED — completed advisor stages only; remaining real worker stages need no model calls',
                             'selection_sha256':snapshot['data']['research']['selection_lock']['selection_sha256'],
                             'elapsed_seconds':time.monotonic()-started}
                    (root.parent/'browser-review-ready.json').write_text(json.dumps(handoff,indent=2))
                    print(json.dumps(handoff,indent=2))
                    return
                answer='run_exploratory' 
            else:
                raise RuntimeError('Unexpected integration question: '+question['id'])
            answered.append({'question_id':question['id'],'answer':answer})
            assert service.answer_research(pid, question['id'], answer, spawn=False)==jid

    result=service.snapshot(pid)
    assert result['project']['state']=='completed', result['jobs']
    assert len(result['steps'])==9 and all(step['status']=='completed' for step in result['steps'])
    assert result['data']['research']['verification']['status']=='software_pilot_passed'
    assert result['data']['research']['verification']['affinity_validation']=='not_established'
    assert len(result['jobs'])==1 and result['jobs'][0]['attempts']==len(answered)+1
    research=result['data']['research']
    assert research['selection_lock']['status']=='locked'
    assert research['heldout_evaluation']['status']=='not_evaluated'
    assert result['data']['protocol']==research['selection_lock']['selected_protocol']['parameters']
    recorded_calls=len(result['events'])
    assert not execute(root,pid,jid)  # Completed same-job resume executes no new work.
    assert len(service.snapshot(pid)['events'])==recorded_calls
    for artifact in result['artifacts']:
        service.store.read_artifact(pid, artifact)
    report={'status':'passed','project_id':pid,'job_id':jid,'root':str(root),'stages_completed':9,
            'advisor':'MOCKED — no network and no paid model calls',
            'scientific_execution':'Actual registered molecular preparation, comparative docking and physical validation; optional tool status remains in raw outputs',
            'requested_engines':args.engines,'elapsed_seconds':time.monotonic()-started,
            'comparison':research['comparison'],'selection_lock':research['selection_lock'],
            'heldout_evaluation':research['heldout_evaluation'],'comparison_exclusions':research.get('comparison_exclusions',[]),
            'same_job_completed_resume_new_calls':0,'review_answers':answered,
            'verification':result['data']['research']['verification'],
            'screening_pose_count':len(result['data']['screening']['poses']),
            'affinity':result['data']['research']['production']['affinity'],
            'resume_attempts':result['jobs'][0]['attempts'],'artifact_count':len(result['artifacts'])}
    path=report_path
    path.write_text(json.dumps(report, indent=2))
    print(json.dumps({key:report[key] for key in ('status','project_id','job_id','root','elapsed_seconds','screening_pose_count','resume_attempts','artifact_count')}, indent=2))
    print('Full immutable evidence:',path)


if __name__=='__main__':
    main()
