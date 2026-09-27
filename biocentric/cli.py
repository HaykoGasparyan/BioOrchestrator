"""CLI over exactly the same service and durable jobs as the interactive interface."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
from .service import Service


def main():
    parser=argparse.ArgumentParser(prog='biocentric',description='BioCentricGPT inspectable discovery workflow')
    parser.add_argument('--root',default=None,help='Project store (default BIOCENTRIC_HOME or .biocentric)')
    sub=parser.add_subparsers(dest='command',required=True)
    sub.add_parser('adapters');sub.add_parser('list');sub.add_parser('demo');sub.add_parser('example')
    create=sub.add_parser('create');create.add_argument('spec',type=Path)
    for name in ['status','run','resume','cancel','export']:
        cmd=sub.add_parser(name);cmd.add_argument('project')
        if name in ['run','resume']: cmd.add_argument('--foreground',action='store_true',help='Run the separate-worker entry point synchronously for automation')
    upload=sub.add_parser('upload');upload.add_argument('project');upload.add_argument('kind',choices=['receptor','structure','library','reference','identity','evidence']);upload.add_argument('file',type=Path)
    decide=sub.add_parser('decide');decide.add_argument('project');decide.add_argument('kind',choices=['identity','target_selection','protocol','exploratory']);decide.add_argument('--value',required=True,help='JSON object');decide.add_argument('--reason',required=True)
    measured=sub.add_parser('measurements');measured.add_argument('project');measured.add_argument('file',type=Path)
    bench=sub.add_parser('benchmark');bench.add_argument('dataset',type=Path);bench.add_argument('config',type=Path);bench.add_argument('scores',type=Path)
    freeze=sub.add_parser('freeze-benchmark');freeze.add_argument('dataset',type=Path);freeze.add_argument('parameters',type=Path);freeze.add_argument('output',type=Path)
    streamed=sub.add_parser('export-streamed');streamed.add_argument('project');streamed.add_argument('output',type=Path)
    streamed.add_argument('--poses-per-part',type=int,default=100)
    args=parser.parse_args();service=Service(args.root)
    try:
        if args.command=='adapters': result=service.adapters()
        elif args.command=='list': result=service.list_projects()
        elif args.command=='demo': result={'project_id':service.create_demo()}
        elif args.command=='example': result={'project_id':service.create_example()}
        elif args.command=='create': result={'project_id':service.create_project(json.loads(args.spec.read_text()))}
        elif args.command=='status': result=service.snapshot(args.project)
        elif args.command in ['run','resume']:
            jid=service.start(args.project,spawn=not args.foreground)
            if args.foreground:
                if service.store.project(args.project)['spec'].get('workflow_mode')=='research':
                    from .research_worker import execute
                else:
                    from .worker import execute
                execute(service.root,args.project,jid)
            result={'job_id':jid,'status':service.store.project(args.project)['state']}
        elif args.command=='cancel': result={'cancel_requested':service.cancel(args.project)}
        elif args.command=='export': result=service.export(args.project)
        elif args.command=='export-streamed':
            from .stream_exports import export_streamed
            result=export_streamed(service.snapshot(args.project),service.root,args.project,args.output,poses_per_part=args.poses_per_part)
        elif args.command=='upload': result=service.upload(args.project,args.file.name,args.file.read_bytes(),args.kind)
        elif args.command=='decide': result={'recorded':service.decide(args.project,args.kind,json.loads(args.value),args.reason)}
        elif args.command=='measurements': result=service.import_measurements(args.project,args.file.read_bytes())
        elif args.command=='benchmark':
            from .benchmark import evaluate
            result=evaluate(json.loads(args.dataset.read_text()),json.loads(args.scores.read_text()),json.loads(args.config.read_text()))
        elif args.command=='freeze-benchmark':
            from .benchmark import freeze_protocol
            result=freeze_protocol(json.loads(args.dataset.read_text()),json.loads(args.parameters.read_text()),args.output)
        print(json.dumps(result,indent=2,default=str))
    except Exception as exc:
        parser.exit(1,f'{type(exc).__name__}: {exc}\n')

if __name__=='__main__': main()
