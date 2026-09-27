"""Bounded real GPU interruption/resume through the registered research path.

Reuses actual saved Linux Vina poses and their PoseBusters reports. No docking,
training, model responses, network calls or scientific fixture outputs are used.
"""
from __future__ import annotations

import argparse
import io
import json
import os
from pathlib import Path
import time

from rdkit import Chem

from biocentric import science
from biocentric.research_workflow import AwaitingAnswer, ResearchWorkflow
from biocentric.service import Service
from biocentric.testsets import build_testset


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--baseline', type=Path, default=Path('runs/ubuntu-diverse-docking-verification/baseline/docking.json'))
    parser.add_argument('--count', type=int, default=16)
    parser.add_argument('--model-dir', type=Path, default=Path('/home/sargis/dti'))
    parser.add_argument('--model-python', type=Path, default=Path('/home/sargis/anaconda3/bin/python'))
    args = parser.parse_args()
    if not 4 <= args.count <= 32:
        parser.error('count must be between 4 and 32')
    out = args.output.resolve()
    if out.exists():
        raise ValueError('Use a fresh verification output directory')
    out.mkdir(parents=True)
    os.environ['BIOCENTRIC_MODEL_DIR'] = str(args.model_dir.resolve())
    os.environ['BIOCENTRIC_MODEL_PYTHON'] = str(args.model_python.resolve())
    baseline = json.loads(args.baseline.read_text())
    selected = [row for row in baseline['poses'] if row['validation']['passed'] is True][:args.count]
    assert len(selected) == args.count
    service = Service(out / 'store')
    pid = service.create_research_project('Bounded actual BioCentric GPU checkpoint resume verification', use_public_example=True)
    buffer = io.StringIO()
    with Chem.SDWriter(buffer) as writer:
        for row in selected:
            path = Path(row['pose_path'])
            assert science.sha256(path) == row['pose_sha256']
            mol = next(iter(Chem.SDMolSupplier(str(path), removeHs=False)))
            assert mol is not None and science._canonical(mol) == row['canonical_smiles']
            writer.write(mol)
    service.upload(pid, 'actual-linux-vina-resume-inputs.sdf', buffer.getvalue().encode(), 'library')
    for row in selected:
        path = Path(row['pose_path'])
        service.store.artifact(pid, path.name, path.read_bytes(), 'scientific_output')
        evidence = row['validation']['artifact']
        path = Path(evidence['path'])
        assert science.sha256(path) == evidence['sha256']
        service.store.artifact(pid, path.name, path.read_bytes(), 'scientific_output')
    jid = service.start_research(pid, spawn=False)
    frozen = build_testset(service.snapshot(pid), organism='Mus musculus')
    service.store.research_update(pid, testset=frozen, tool_plan={'docking_engine': 'vina'})
    service.store.update(pid, data={'screening': {'poses': selected}})
    workflow = ResearchWorkflow(service.store, pid, jid)
    workflow.stage = '08'
    service.store.ensure_step(pid, '08', 'Real GPU interrupted request resume')
    protocol = workflow.protocol()
    try:
        workflow.native_calculation('biocentric_affinity', 'affinity', protocol)
    except AwaitingAnswer:
        pass
    else:
        raise AssertionError('Native execution bypassed concrete review')
    questions = workflow.research['questions']
    pending = [row for row in questions if row['status'] == 'pending']
    assert len(pending) == 1 and pending[0]['id'].startswith('native_review_affinity')
    pending[0].update(status='answered', answer='approve')
    service.store.research_update(pid, questions=questions)
    engine_out = Path(workflow.research['engine_preflight']['engine_plan']['output_dir'])
    ordinary_cancelled = workflow.cancelled
    def interrupt_after_checkpoint():
        return any((engine_out / 'pairs').glob('*.json')) or ordinary_cancelled()
    workflow.cancelled = interrupt_after_checkpoint
    started = time.monotonic()
    try:
        workflow.native_calculation('biocentric_affinity', 'affinity', protocol)
    except ValueError as error:
        assert 'did not complete' in str(error)
    else:
        raise AssertionError('The bounded native process did not interrupt')
    interrupted = workflow.research['native_failure']
    assert interrupted['status'] == 'cancelled'
    before = {path.name: {'sha256': science.sha256(path), 'mtime_ns': path.stat().st_mtime_ns,
                          'row': json.loads(path.read_text())['row']}
              for path in (engine_out / 'pairs').glob('*.json')}
    assert 0 < len(before) < args.count
    assert all(item['row']['status'] == 'completed' for item in before.values())
    request_sha = science.sha256(engine_out / 'request.json')
    later = ResearchWorkflow(service.store, pid, jid)
    later.stage = '08'
    result = later.native_calculation('biocentric_affinity', 'affinity', protocol)
    assert result['status'] == 'completed' and len(result['results']) == args.count
    assert result['native_result']['accounting'] == {'requested': args.count, 'completed': args.count, 'failed': 0}
    assert not result['native_result']['progress_warnings']
    assert science.sha256(engine_out / 'request.json') == request_sha
    for name, item in before.items():
        path = engine_out / 'pairs' / name
        assert science.sha256(path) == item['sha256'] and path.stat().st_mtime_ns == item['mtime_ns']
    reused = sum(row.get('reused') is True for row in result['results'])
    assert reused == len(before)
    assert {row['compound_id'] for row in result['results']} == {row['compound_id'] for row in selected}
    assert len({row['compound_id'] for row in result['results']}) == args.count
    archive = json.loads(service.store.read_artifact(pid, {'id': result['metadata']['checkpoint_resume']['archive_artifact_id']}))
    for item in archive['artifacts']:
        service.store.read_artifact(pid, {'id': item['artifact_id']})
    summary = {'status': 'passed', 'scope': 'Actual GPU application interruption/resume; no predictive accuracy established',
               'project_id': pid, 'job_id': jid, 'store': str(out / 'store'), 'requested': args.count,
               'source_baseline': str(args.baseline.resolve()), 'source_baseline_sha256': science.sha256(args.baseline),
               'interrupted_status': interrupted['status'], 'completed_before_interrupt': len(before),
               'reused_after_resume': reused, 'new_predictions_after_resume': args.count - reused,
               'completed_final': len(result['results']), 'failed_final': 0,
               'completed_checkpoint_hashes_and_mtimes_unchanged': True,
               'same_request_sha256': request_sha, 'same_output_directory': str(engine_out),
               'preserved_before_resume': before, 'archive_artifact_id': result['metadata']['checkpoint_resume']['archive_artifact_id'],
               'archived_artifacts_verified': len(archive['artifacts']), 'final_result_artifact_id': result['result_artifact_id'],
               'runtime': result['native_result']['native_runtime'], 'resources': result['native_result']['resource_usage'],
               'elapsed_seconds': time.monotonic() - started, 'paid_model_calls': 0, 'network_calls': 0}
    science._write_json(out / 'verification.json', summary)
    print(json.dumps({key: value for key, value in summary.items() if key != 'preserved_before_resume'}, indent=2))


if __name__ == '__main__':
    main()
