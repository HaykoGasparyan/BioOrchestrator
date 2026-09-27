"""Deterministic report exports: observations, computations and proposals remain distinct."""
from __future__ import annotations
import csv
import html
import importlib.metadata
import io
import json
import platform
import sys
from .models import RunManifest, utc
from .store import dump


def export_project(store,pid):
    p=store.project(pid); data=p['data']; candidates=data.get('candidates',[])
    artifacts=store.artifacts(pid); jobs=store.jobs(pid); decisions=store.decisions(pid)
    fixture=p['spec']['mode']=='offline_fixture'
    status='OFFLINE FIXTURE — NO CALCULATION' if fixture else data.get('result_status','exploratory / incomplete')
    env={'python':sys.version,'platform':platform.platform(),'machine':platform.machine(),'cpu':p['spec']['cpu'],'gpu_requested':p['spec']['has_gpu'],'seed':p['spec']['seed'],'prompt_version':'research-v2' if p['spec'].get('workflow_mode')=='research' else 'runtime-v1'}
    for name in ['rdkit','meeko','posebusters','openai','streamlit']:
        try: env[name]=importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError: env[name]='not installed'
    env['scientific_tools']=data.get('screening',{}).get('metadata',{})
    limitations=list(dict.fromkeys(data.get('limitations',[])+['Computational ranking is a hypothesis, not confirmed binding or a drug.','No held-out enrichment validation is established unless a separately audited benchmark is attached.','ADMET, selectivity, availability and synthesis feasibility are unknown unless evidenced.','Reproducible configuration and hashes do not imply bit-identical results across hardware.']))
    manifest=RunManifest(run_id=jobs[0]['id'] if jobs else 'not_run',project_id=pid,created_at=utc(),status=status,spec=p['spec'],input_hashes={f"{a['kind']}/{a['name']}/{a['id']}":a['sha256'] for a in artifacts if a['kind'] in {'receptor','library','structure','reference','identity','immutable_input'}},artifacts=[a for a in artifacts if a['kind']!='report'],jobs=jobs,decisions=decisions,environment=env,limitations=limitations).model_dump()
    lines=[f"# {__import__('os').getenv('BIOCENTRIC_PRODUCT_NAME','BioCentricGPT')} — {p['spec']['name']}",f"\n**Execution: {p['state']} · Scientific status: {status}**",f"\nProject `{pid}` · Generated {utc()}","\n## Objective and target",p['spec'].get('context') or 'Biological context unresolved.',"```json",dump(data.get('target_identity',{'status':'unresolved'})),"```","\n## Evidence and coverage"]
    for e in data.get('evidence',[]):
        lines.append(f"- [{e.get('id')}] [{e.get('source')} / {e.get('source_id')}]({e.get('url')}) — {e.get('evidence_type')}; retrieved {e.get('retrieved_at')}; location: {e.get('location')}. {e.get('text','')}")
    if not data.get('evidence'): lines.append('No retrieved external evidence. No literature claim is inferred.')
    if data.get('research'):
        lines+=['\n## Research conversation and source review','```json',dump(data['research']),'```']
    lines+=['\n## Pipeline and scientific gates','```json',json.dumps({'plan':data.get('plan'),'protocol':data.get('protocol'),'benchmark':data.get('benchmark'),'gate':data.get('gate'),'steps':store.steps(pid)},indent=2,default=str),'```','\n## Candidate hypotheses']
    if candidates:
        lines += ['| Compound | Docking score (see engine and units in candidates.json) | Selected | Rationale |','|---|---:|---|---|']
        for c in candidates:
            lines.append('| '+ ' | '.join(str(x).replace('|',' / ').replace('\n',' ') for x in [c['compound_id'],c.get('score','unknown'),c.get('selected'),c.get('rationale')])+' |')
    else: lines.append('No candidate results. Failed or blocked computations are not replaced with estimates.')
    lines += ['\nEvery computed candidate includes a job ID and artifact IDs in candidates.json. Original SDF, pose SDF/PDBQT, preparation records and logs are retained.','\n## Experimental handoff — proposed for scientist review','Confirm compound identity and purity, then repeat concentration–response testing in a mechanism-appropriate assay. Use an orthogonal binding or functional measurement, interference/aggregation counterscreens, selectivity assays, and relevant cellular target-engagement/function assays. No physical experiments were performed by this software.']
    if not p['spec'].get('intended_effect'): lines.append('Intended mechanism is unresolved: select the functional assay only after clarifying the intended biological effect.')
    if data.get('measurements'):
        lines += ['\n## Experimental results imported by the user','Original endpoints, units, relations and assay context are preserved in measured-results.csv. They are not inferred from docking.','```json',json.dumps(data['measurements'],indent=2),'```']
    lines += ['\n## Limitations']+[f'- {x}' for x in limitations]
    lines += ['\n## Reproducibility','See manifest.json for input/output hashes, versions, CPU, seeds, job IDs and decisions. See events.json for recorded tool and model requests and outputs.']
    markdown='\n'.join(lines)+'\n'
    document='<!doctype html><html><head><meta charset="utf-8"><title>BioCentricGPT report</title><style>body{max-width:1100px;margin:48px auto;font:15px/1.6 system-ui;color:#173b3b}pre{white-space:pre-wrap;overflow-wrap:anywhere;background:#f5f8f8;padding:24px}h1{color:#087e83}</style></head><body><h1>BioCentricGPT report</h1><pre>'+html.escape(markdown)+'</pre></body></html>'
    output=io.StringIO(); fields=['compound_id','rank','score','selected','rationale','status','job_id','artifact_ids','missing_evidence']
    writer=csv.DictWriter(output,fieldnames=fields);writer.writeheader()
    for c in candidates:
        row={k:json.dumps(c[k]) if isinstance(c.get(k),(dict,list)) else c.get(k,'') for k in fields}
        row={k:("'"+v if isinstance(v,str) and v.startswith(('=','+','-','@')) else v) for k,v in row.items()}
        writer.writerow(row)
    documents={'markdown':('report.md',markdown),'html':('report.html',document),'csv':('candidates.csv',output.getvalue()),'json':('candidates.json',json.dumps(candidates,indent=2,default=str)),'manifest':('manifest.json',json.dumps(manifest,indent=2,default=str)),'events':('events.json',json.dumps(store.events(pid),indent=2,default=str))}
    from .affinity_results import affinity_rows, affinity_csv
    affinity = affinity_rows({'data': data})
    if affinity:
        documents['affinity_csv'] = ('affinity_predictions.csv', affinity_csv(affinity).decode())
        documents['affinity_json'] = ('affinity_predictions.json', json.dumps(affinity, indent=2, allow_nan=False))
    return {key:store.artifact(pid,name,content.encode(),'report',jobs[0]['id'] if jobs else None)['path'] for key,(name,content) in documents.items()}
