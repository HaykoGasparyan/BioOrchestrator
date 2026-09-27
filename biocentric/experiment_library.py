"""Durable experiment workspaces; navigation never owns or cancels workers."""
from copy import deepcopy
import json
from pathlib import Path
import re
import uuid

from . import discovery_benchmark as benchmark, science

KEYS = {'studio_docking_preset', 'studio_draft', 'studio_target', 'studio_identity', 'studio_resolution', 'studio_evidence',
        'studio_eval_dataset', 'studio_eval_review', 'studio_eval_config', 'studio_parameters',
        'studio_custom_plugins', 'studio_parameter_errors', 'studio_eval_limit', 'studio_eval_cohort',
        'studio_eval_reference', 'studio_structure_page', 'studio_binding_page', 'studio_gallery_batch',
        'studio_eval_inspect', 'studio_eval_requested', 'studio_eval_max_compounds'}
PREFIXES = ('studio_tools_', 'studio_eval_count_', 'studio_eval_ligand_', 'studio_eval_chain_', 'studio_final_complex_')


def valid_id(value):
    return isinstance(value, str) and bool(re.fullmatch(r'[a-f0-9]{32}', value))


def workspace_path(root, experiment_id):
    if not valid_id(experiment_id):
        raise ValueError('Invalid experiment ID')
    return Path(root) / 'experiments' / experiment_id / 'workspace.json'


def load(root, experiment_id):
    value = json.loads(workspace_path(root, experiment_id).read_text())
    if value.get('id') != experiment_id or not isinstance(value.get('state'), dict):
        raise ValueError('Invalid experiment workspace')
    if value.get('deleted_at'):
        raise ValueError('This experiment was deleted')
    return value


def is_deleted(root, experiment_id):
    return bool(json.loads(workspace_path(root, experiment_id).read_text()).get('deleted_at'))


def delete(root, experiment_id):
    """Remove from the library, stop its worker, and retain recoverable local files."""
    saved = load(root, experiment_id)
    run_id = saved.get('run_id')
    if run_id:
        if not valid_id(run_id):
            raise ValueError('Invalid evaluation ID')
        folder = Path(root) / 'evaluation-runs' / run_id
        if (folder / 'status.json').is_file() and benchmark.read_status(folder)['status'] in {'queued', 'running'}:
            benchmark.cancel_run(folder)
    # Also hide any recovery alias for the same worker, so it cannot reappear.
    targets = [(workspace_path(root, experiment_id), saved)]
    if run_id:
        for path in (Path(root) / 'experiments').glob('*/workspace.json'):
            if path == targets[0][0]:
                continue
            try:
                value = json.loads(path.read_text())
                if value.get('run_id') == run_id:
                    targets.append((path, value))
            except (OSError, ValueError):
                continue
    for path, value in targets:
        value['deleted_at'] = science._now()
        science._write_json(path, value)


def undo_delete(root, experiment_id):
    path = workspace_path(root, experiment_id)
    saved = json.loads(path.read_text())
    if saved.get('deleted_at'):
        saved.pop('deleted_at')
        science._write_json(path, saved)
    # A stopped worker stays stopped; undo never restarts scientific work.


def delete_many(root, experiment_ids):
    """Validate the complete selection before stopping or hiding any experiment."""
    selected = list(dict.fromkeys(experiment_ids))
    for experiment_id in selected:
        saved = load(root, experiment_id)
        if saved.get('run_id') and not valid_id(saved['run_id']):
            raise ValueError('Invalid evaluation ID')
    deleted, errors = [], {}
    for experiment_id in selected:
        try:
            # Deleting another selected alias may already have hidden this row.
            if not is_deleted(root, experiment_id):
                delete(root, experiment_id)
            deleted.append(experiment_id)
        except (OSError, ValueError) as exc:
            errors[experiment_id] = str(exc)
    return deleted, errors


def create(root):
    experiment_id = uuid.uuid4().hex
    value = {'schema': 1, 'id': experiment_id, 'created_at': science._now(), 'updated_at': science._now(), 'run_id': None, 'state': {}}
    science._write_json(workspace_path(root, experiment_id), value)
    return experiment_id


def snapshot(state):
    return deepcopy({k: v for k, v in state.items() if k in KEYS or k.startswith(PREFIXES)})


def save(root, experiment_id, state):
    previous = load(root, experiment_id)
    run_id = Path(state['studio_eval_run']).name if state.get('studio_eval_run') else None
    if run_id and not valid_id(run_id):
        raise ValueError('Invalid evaluation ID')
    # A rerun or a different target becomes its own row; never overwrite a past run.
    if previous.get('run_id') and previous['run_id'] != run_id:
        experiment_id = create(root)
        previous = load(root, experiment_id)
    current = snapshot(state)
    if current != previous['state'] or run_id != previous.get('run_id'):
        previous.update(state=current, run_id=run_id, updated_at=science._now())
        science._write_json(workspace_path(root, experiment_id), previous)
    return experiment_id


def clear_workspace(state):
    for key in list(state):
        if key.startswith('studio_'):
            del state[key]
    state['studio_version'] = 3


def restore(root, experiment_id, state):
    saved = load(root, experiment_id)
    clear_workspace(state)
    state.update(snapshot(saved['state']))
    if state.get('studio_draft', {}).get('evidence_query'):
        state['studio_draft']['evidence_query'] = tuple(state['studio_draft']['evidence_query'])
    run_id = saved.get('run_id')
    if run_id:
        if not valid_id(run_id):
            raise ValueError('Invalid evaluation ID in workspace')
        folder = Path(root) / 'evaluation-runs' / run_id
        review = json.loads((folder / 'plan.json').read_text())
        benchmark.validate_review(review)
        state['studio_eval_run'] = str(folder.resolve())
        state['studio_eval_review'] = review
        state['studio_eval_max_compounds'] = (review.get('ligand_subset') or {}).get('max_compounds') or 0
        state['studio_identity'] = review['identity']
        state['studio_resolution'] = {'status': 'resolved', 'candidates': [review['identity']]}
        target = (review['identity'].get('gene_names') or [review['identity']['accession']])[0]
        draft = state.setdefault('studio_draft', {})
        draft.update(target=target, request=target, organism='Human', status='review', evidence_query=(target, 'Human'))
        draft.setdefault('messages', [{'role': 'user', 'text': target}])
        state['studio_target'] = target
        state['studio_eval_requested'] = True
        for phase, rows in review['selection'].items():
            state[f'studio_tools_{phase}'] = [r['tool'] for r in rows]
            for row in rows:
                if row['tool'].startswith('custom_'):
                    state.setdefault('studio_custom_plugins', {}).setdefault(row['tool'],
                        {'name': row['tool'], 'phase': phase, 'parameters': row['parameters']})
                state.setdefault('studio_parameters', {})[f"{phase}:{row['tool']}"] = row['parameters']
        if 'studio_eval_dataset' not in state:
            # Older runs retained only their frozen subset, not the full source search.
            records = [{**r, 'eligible': True, 'reasons': [], 'p_relation': '='} for r in review['records']]
            coverage = review.get('source_coverage', {})
            state['studio_eval_dataset'] = {'identity': review['identity'], 'records': records, 'source_total': coverage.get('source_total', len(records)),
                'loaded': len(records), 'cap': len(records), 'complete': False, 'created_at': review.get('created_at'),
                'errors': ['Restored frozen evaluation subset. The full search was not saved in this older experiment.']}
    state['experiments_current'] = experiment_id
    return saved


def adopt_runs(root):
    """Make pre-library and orphaned worker runs discoverable without changing artifacts."""
    known = set()
    for path in (Path(root) / 'experiments').glob('*/workspace.json'):
        try:
            known.add(json.loads(path.read_text()).get('run_id'))
        except (OSError, ValueError):
            continue
    for path in (Path(root) / 'evaluation-runs').glob('*/plan.json'):
        run_id = path.parent.name
        if run_id in known or not valid_id(run_id):
            continue
        try:
            review = json.loads(path.read_text())
            benchmark.validate_review(review)
            target = (review['identity'].get('gene_names') or [review['identity']['accession']])[0]
            value = {'schema': 1, 'id': run_id, 'created_at': review.get('created_at', science._now()), 'updated_at': science._now(),
                     'run_id': run_id, 'state': {'studio_draft': {'target': target}, 'studio_identity': review['identity']}}
            destination = workspace_path(root, run_id)
            if not destination.exists():
                science._write_json(destination, value)
        except (OSError, ValueError, KeyError):
            continue


def summaries(root):
    from . import pipeline_decisions as decisions
    rows = []
    for path in (Path(root) / 'experiments').glob('*/workspace.json'):
        if not valid_id(path.parent.name):
            continue
        try:
            if json.loads(path.read_text()).get('deleted_at'):
                continue
        except (OSError, ValueError):
            pass
        row = {'id': path.parent.name, 'target': 'New experiment', 'status': 'draft', 'created': '', 'progress': '—',
               'measurements': None, 'chain': '—', 'rho': None, 'result': 'Not started', 'run_id': None, 'saved_fields': 0}
        try:
            saved = load(root, path.parent.name)
            row.update(target=saved['state'].get('studio_draft', {}).get('target') or 'New experiment', created=saved['created_at'],
                       run_id=saved.get('run_id'), saved_fields=len(saved['state']))
            if saved.get('run_id'):
                folder = Path(root) / 'evaluation-runs' / saved['run_id']
                review, status = decisions.context(folder)
                lock = decisions.read_lock(folder, review, status)
                progress = benchmark.combination_progress(review, status)
                completed = progress['completed']
                row.update(status=status['status'], progress=f"{completed}/{progress['total']} evaluated · {progress['selected_total']} total", measurements=len(review['records']))
                if lock:
                    selected = next(r for r in decisions.compare(review, status) if r['id'] == lock['candidate']['id'])
                    row.update(chain=benchmark.chain_label(lock['candidate']),
                               rho=lock['metrics']['rho'], result='Final chain locked',
                               ranked_measurements=selected['ranked_measurements'], small_sample=selected['small_sample'])
                elif status['status'] in {'queued', 'running'}:
                    row['result'] = 'Evaluating selected tools'
                elif status['status'] in {'cancelled', 'interrupted'}:
                    row['result'] = 'Saved partial results'
                else:
                    metrics = decisions.compare(review, status)
                    best = next((r for r in metrics if r['rank'] == 1), None)
                    if best:
                        row.update(chain=benchmark.chain_label(best['pair']), rho=best['rho'], result='Ready to choose final chain')
                    else:
                        row['result'] = 'No eligible final chain' if completed else 'No evaluated combinations'
        except (OSError, ValueError, KeyError) as exc:
            row.update(status='unavailable', result=str(exc))
        rows.append(row)
    # A different browser can discover a new worker just before its workspace is saved.
    # Prefer the full workspace over that recovery entry, keeping one row per run.
    unique = {}
    for row in sorted(rows, key=lambda r: r['saved_fields']):
        unique[row['run_id'] or row['id']] = row
    return sorted(unique.values(), key=lambda r: (r['created'], r['id']), reverse=True)
