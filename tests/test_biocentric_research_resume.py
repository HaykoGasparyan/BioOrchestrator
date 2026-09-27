"""Application resume with a synthetic subprocess; no learned predictions here."""
import io
import json
from pathlib import Path

import pytest
from rdkit import Chem

from biocentric import biocentric_model, science
from biocentric.biocentric_predictor import read_checkpoint, row_sha256
from biocentric.research_workflow import AwaitingAnswer, ResearchWorkflow
from biocentric.service import Service
from biocentric.store import digest
from biocentric.testsets import build_testset


@pytest.fixture
def research_run(tmp_path, monkeypatch):
    service = Service(tmp_path / 'store')
    pid = service.create_research_project('Synthetic interruption control', use_public_example=True)
    jid = service.start_research(pid, spawn=False)
    frozen = build_testset(service.snapshot(pid), organism='Mus musculus')
    service.store.research_update(pid, testset=frozen, tool_plan={'docking_engine': 'vina'})
    workflow = ResearchWorkflow(service.store, pid, jid)
    protocol = workflow.protocol()
    sources = {row['kind']: row for row in service.store.artifacts(pid)}
    poses = []
    for i, mol in enumerate(Chem.SDMolSupplier(sources['library']['path'], removeHs=False), 1):
        from rdkit.Chem import AllChem
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        buffer = io.StringIO()
        with Chem.SDWriter(buffer) as writer:
            writer.write(mol)
        artifact = service.store.artifact(pid, f'fixture-{i}.sdf', buffer.getvalue().encode(), 'scientific_output')
        validation = service.store.artifact(pid, f'fixture-{i}-validation.json',
                            json.dumps({'passed': True, 'fixture_scope': 'Synthetic metadata only'}).encode(), 'scientific_output')
        smiles = science._canonical(mol)
        cid = 'cmp_' + digest(smiles.encode())[:16]
        poses.append({'compound_id': cid, 'pose_id': cid + '_pose_1', 'canonical_smiles': smiles,
                      'pose_sha256': artifact['sha256'], 'pose_number': 1,
                      'receptor_sha256': sources['structure']['sha256'], 'validation': {'passed': True, 'artifact': validation},
                      'fixture_scope': 'Synthetic passing metadata for workflow guards only'})
    service.store.update(pid, data={'screening': {'poses': poses}})
    workflow.stage = '08'
    service.store.ensure_step(pid, '08', 'Synthetic resume control')
    model_file = tmp_path / 'fixture-model.py'
    model_file.write_text('# synthetic runtime binding\n')
    models = [f'{name}@{epoch}' for name, epoch in biocentric_model.MODELS]
    def installed():
        return {'root': str(tmp_path), 'python': None, 'blockers': [], 'models': models,
                'weights_sha256': biocentric_model.SUPPORTED_WEIGHTS_SHA256,
                'code_sha256': science.sha256(model_file),
                'files': {str(model_file): {'path': str(model_file), 'sha256': science.sha256(model_file), 'bytes': model_file.stat().st_size}}}
    monkeypatch.setattr(biocentric_model, 'installation', installed)
    return service, pid, jid, protocol, model_file


def invoke(context):
    service, pid, jid, protocol, _ = context
    workflow = ResearchWorkflow(service.store, pid, jid)
    workflow.stage = '08'
    return workflow.native_calculation('biocentric_affinity', 'affinity', protocol)


def answer_pending(context, answer):
    service, pid, _, _, _ = context
    questions = service.store.project(pid)['data']['research']['questions']
    question = next(row for row in questions if row['status'] == 'pending')
    assert answer in {row['value'] for row in question['options']}
    question.update(status='answered', answer=answer)
    service.store.research_update(pid, questions=questions)


def fake_runner(monkeypatch, *, first_failure=False):
    calls, fresh, requests = [], [], []
    def run(argv, cwd, log, timeout, cancelled):
        calls.append(str(cwd))
        request_path = Path(argv[-1])
        request = json.loads(request_path.read_text())
        request_sha = science.sha256(request_path)
        requests.append(request_sha)
        (cwd / 'pairs').mkdir(exist_ok=True)
        rows = []
        log.write_text(f'Synthetic process attempt {len(calls)}\n')
        for pair in request['pairs']:
            path = cwd / 'pairs' / f"{pair['record_index']:06d}.json"
            if path.exists():
                row = {**read_checkpoint(path, request_sha, pair), 'reused': True}
            else:
                fresh.append(pair['compound_id'])
                row = {**pair, 'status': 'completed', 'score': 6., 'failure_reason': None,
                       'per_model': {model: 6. for model in request['models']}, 'models': request['models'],
                       'checkpoint_sha256': request['weights_sha256'], 'code_sha256': request['code_sha256'],
                       'device': request['parameters']['device'], 'preprocessing_status': 'passed', 'reused': False}
                if first_failure and len(calls) == 1 and pair['record_index'] == 2:
                    row.update(status='failed', score=None, per_model={}, failure_reason='Synthetic transient pair failure')
                path.write_text(json.dumps({'request_sha256': request_sha, 'row_sha256': row_sha256(row), 'row': row}))
            rows.append(row)
            if len(calls) == 1 and pair['record_index'] == (2 if first_failure else 1):
                raise science.ScienceError('timeout', 'Synthetic interruption after checkpoint')
        (cwd / 'predictions.json').write_text(json.dumps({'schema_version': 1, 'engine_id': 'biocentric_affinity', 'predictions': rows}))
        return {'returncode': 0}
    monkeypatch.setattr(science, 'run_process', run)
    return calls, fresh, requests


def start_interrupted(context):
    with pytest.raises(AwaitingAnswer):
        invoke(context)
    answer_pending(context, 'approve')
    with pytest.raises(ValueError, match='did not complete'):
        invoke(context)


def test_application_resume_reuses_completed_pair_and_exact_request(research_run, monkeypatch):
    calls, fresh, requests = fake_runner(monkeypatch)
    start_interrupted(research_run)
    checkpoint = Path(calls[0]) / 'pairs/000001.json'
    original = checkpoint.read_bytes()
    result = invoke(research_run)
    assert result['status'] == 'completed'
    assert len(calls) == 2 and calls[0] == calls[1]
    assert requests[0] == requests[1]
    assert len(fresh) == 3 and len(set(fresh)) == 3
    assert result['results'][0]['reused'] is True
    assert checkpoint.read_bytes() == original
    assert result['metadata']['checkpoint_resume']['completed'] == 1
    service, pid, *_ = research_run
    archive = service.store.read_artifact(pid, {'id': result['metadata']['checkpoint_resume']['archive_artifact_id']})
    receipt = json.loads(archive)
    prior_log = next(item for item in receipt['artifacts'] if item['relative_path'] == 'engine/engine.log')
    assert service.store.read_artifact(pid, {'id': prior_log['artifact_id']}) == b'Synthetic process attempt 1\n'


@pytest.mark.parametrize('retry_failed', [False, True])
def test_recorded_failures_need_explicit_review_and_are_archived(research_run, monkeypatch, retry_failed):
    calls, fresh, requests = fake_runner(monkeypatch, first_failure=True)
    start_interrupted(research_run)
    failed_path = Path(calls[0]) / 'pairs/000002.json'
    failed_bytes = failed_path.read_bytes()
    with pytest.raises(AwaitingAnswer):
        invoke(research_run)
    assert len(calls) == 1
    answer_pending(research_run, 'retry_failed' if retry_failed else 'resume')
    if retry_failed:
        result = invoke(research_run)
        assert result['status'] == 'completed'
        assert len(fresh) == 4 and len(set(fresh)) == 3
        assert json.loads(failed_path.read_text())['row']['status'] == 'completed'
        assert len({row['compound_id'] for row in result['results']}) == 3
    else:
        with pytest.raises(ValueError, match='did not complete'):
            invoke(research_run)
        assert len(fresh) == 3 and failed_path.read_bytes() == failed_bytes
    assert calls[0] == calls[1] and requests[0] == requests[1]
    service, pid, *_ = research_run
    prior = [item for item in service.store.artifacts(pid) if item['kind'] == 'native_engine_prior_attempt' and item['name'] == '000002.json']
    assert len(prior) == 1 and service.store.read_artifact(pid, prior[0]) == failed_bytes


@pytest.mark.parametrize('mutation', ['model', 'checkpoint', 'reviewed_plan'])
def test_application_resume_rejects_changed_provenance_before_execution(research_run, monkeypatch, mutation):
    calls, _, _ = fake_runner(monkeypatch)
    start_interrupted(research_run)
    service, pid, _, _, model_file = research_run
    if mutation == 'model':
        model_file.write_text('# changed runtime\n')
    elif mutation == 'checkpoint':
        path = Path(calls[0]) / 'pairs/000001.json'
        record = json.loads(path.read_text()); record['row']['score'] = 42
        path.write_text(json.dumps(record))
    else:
        reviewed = next(item for item in service.store.artifacts(pid) if item['kind'] == 'native_engine_review')
        Path(reviewed['path']).write_text('{}')
    with pytest.raises(ValueError):
        invoke(research_run)
    assert len(calls) == 1
