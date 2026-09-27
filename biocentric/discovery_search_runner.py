"""Ligand-first search rounds, using real docking/scoring and immutable evidence."""
import json
import math
import time
from pathlib import Path

from . import discovery_benchmark as bench, science, bayesian_search as search
from .evaluation_archive import Archive, StageCache, identify_pocket


def summary_metrics(review, rows):
    if not review['search']['method'].endswith(('_v3', '_v4')):
        return bench.summary_metrics(review, rows)
    from .pipeline_decisions import compare
    return [{'Combination': bench.chain_label(r['pair']), 'Measurements': r['completed'],
             'Expected': r['expected'], 'Skipped': r['skipped'], 'Shared set': r['shared_measurements'],
             'Correlation measurements': r['ranked_measurements'], 'Assays ranked': r['ranked_cohorts'],
             'Spearman ρ': r['rho'], 'Kendall τ': r['tau']}
            for r in compare(review, {'status': 'running', 'rows': rows})]


def evaluate(pair, record, index, folder, receptor, pdbqt, reference, pocket_cache,
             dock_cache, reference_scores, state, save, cancelled, emit):
    from rdkit import Chem
    archive = getattr(dock_cache, 'archive', None)
    receipt_start = len(archive.receipts) if archive else 0
    method = pair.get('pocket', bench.pockets.REFERENCE)
    row = {'combination': pair['id'], 'docking': pair['docking']['tool'], 'scoring': pair['scoring']['tool'],
           'compound_id': record['compound_id'], 'measurement_id': bench.measurement_key(record),
           'cohort': record.get('cohort'), 'pocket': method['tool'], 'kind': record['kind'],
           'p_value': record.get('p_value'), 'endpoint': record.get('endpoint'), 'score': None,
           'score_unit': bench.SCORES[pair['scoring']['tool']][1], 'rmsd_angstrom': None,
           'physical_pass': None, 'status': 'failed'}
    began = time.monotonic()
    state.update(active_combination=pair['id'], active_compound=record['compound_id'], stage='pocket')
    save()
    try:
        key = bench.digest(method)
        if key not in pocket_cache:
            try:
                pocket_cache[key] = identify_pocket(method, receptor, folder / 'reference.sdf', folder / 'pockets' / key[:16], cancelled, archive)
            except Exception as exc:
                science._cancel(cancelled)
                pocket_cache[key] = str(exc)
        box = pocket_cache[key]
        if isinstance(box, str):
            raise ValueError('Pocket identification failed: ' + box)
        row['pocket_box'] = box
        work = folder / 'calculations' / pair['id'] / str(index)
        work.mkdir(parents=True, exist_ok=True)
        state['stage'] = 'docking'
        save()
        poses = bench._prepare_docked_poses(pair, record, work, receptor, pdbqt, box, dock_cache, cancelled, emit)
        state['stage'] = 'scoring'
        save()
        same = science._canonical(poses[0]) == science._canonical(reference)
        score_key = bench.digest([method, box, pair['scoring']])
        prediction = None
        if pair['scoring']['tool'] == 'biocentric_affinity':
            from .discovery_affinity import score_groups
            groups = {'poses': poses}
            if same and score_key not in reference_scores:
                groups['reference'] = [Chem.AddHs(reference, addCoords=True)]
            output = score_groups(groups, receptor, work / 'score', pair['scoring']['parameters'], cancelled, emit, crop_reference=reference)
            prediction = output['poses']
            if prediction['errors']:
                raise ValueError('; '.join(prediction['errors']))
            scores = prediction['scores']
            if 'reference' in output:
                if output['reference']['errors']:
                    raise ValueError('Crystal scoring failed: ' + '; '.join(output['reference']['errors']))
                reference_scores[score_key] = output['reference']['scores'][0]
        else:
            scores = bench._score(pair, poses, work / 'score', receptor, pdbqt, cancelled, emit, box, archive=archive)
            if same and score_key not in reference_scores:
                reference_scores[score_key] = bench._score(pair, [Chem.AddHs(reference, addCoords=True)],
                    work / 'reference-score', receptor, pdbqt, cancelled, emit, box, archive=archive)[0]
        if len(scores) != len(poses) or any(not math.isfinite(value) for value in scores):
            raise ValueError('Missing or nonfinite pose scores')
        best = (max if bench.SCORES[pair['scoring']['tool']][2] else min)(range(len(scores)), key=lambda i: scores[i])
        pose = poses[best]
        pose_path = Path(bench._write_sdf(work / 'best-pose.sdf', pose))
        state['stage'] = 'validation'
        save()
        validity = science.validate_pose(pose_path, receptor, work / 'validation', 120, cancelled)
        if prediction:
            row['model_prediction'] = {'device': prediction['device'], 'provenance': prediction['provenance'],
                'best_pose': prediction['predictions'][best],
                'artifact': str(Path(prediction['artifact']).relative_to(folder)),
                'artifact_sha256': science.sha256(Path(prediction['artifact'])),
                'ensemble_spread_calibrated': False, 'experimental_measurement': False}
        row.update(status='completed', score=scores[best], pose_path=str(pose_path.relative_to(folder)),
                   pose_sha256=science.sha256(pose_path), physical_pass=validity.get('passed'),
                   physical_status=validity.get('status'),
                   rmsd_angstrom=science._symmetry_rmsd(pose, reference).get('rmsd_angstrom') if same else None,
                   reference_comparable=same, reference_score=reference_scores.get(score_key) if same else None,
                   top_pose_index=best+1, all_scores=scores)
    except Exception as exc:
        science._cancel(cancelled)
        row['error'] = str(exc)
    if archive:
        row['archive'] = archive.receipts[receipt_start:]
    row['runtime_seconds'] = time.monotonic() - began
    return row


def worker(folder, review):
    from rdkit import Chem
    from . import pipeline_decisions as decisions
    folder = Path(folder).resolve()
    state = bench.read_status(folder)
    began = time.monotonic()
    config = review['search']
    tasks = [{'compound_id': 'crystal:' + review['reference']['component'], 'canonical_smiles': review['reference']['canonical_smiles'],
              'kind': 'reference', 'p_value': None, 'endpoint': None},
             *sorted([{**r, 'kind': 'measured'} for r in review['records']], key=lambda r: (r['compound_id'], bench.measurement_key(r)))]
    state.update(status='running', started_at=science._now(), stage='preparation',
                 search={'round': 0, 'observations': [], 'rounds': [], 'selected_ids': [], 'budget': config['max_trials']})
    def save():
        state.update(updated_at=science._now(), elapsed_seconds=time.monotonic()-began)
        progress = bench.combination_progress(review, state)
        state.update(completed=progress['completed'], total=progress['total'], progress_unit='combinations')
        science._write_json(folder / 'status.json', state)
    def cancelled():
        budget = review['protocol'].get('total_budget_seconds')
        return (folder / 'cancel.request').exists() or (budget is not None and time.monotonic()-began > budget)
    def emit(event):
        if event.get('stage') == 'archive':
            counts = state.setdefault('reuse', {'hits': 0, 'stored': 0})
            counts['hits' if event.get('message', '').startswith('Reused') else 'stored'] += 1
        state['events'] = (state.get('events', []) + [{k: event.get(k) for k in ('tool', 'stage', 'message')}])[-40:]
        save()
    try:
        if bench.runtime_fingerprint(review['candidates']) != review['runtime_fingerprint']:
            raise ValueError('Runtime changed after approval')
        save()
        receptor, reference_path = folder / 'receptor.pdb', folder / 'reference.sdf'
        receptor.write_text(review['reference']['receptor_pdb'])
        reference_path.write_text(review['reference']['reference_sdf'])
        pdbqt = Path(science.prepare_receptor(receptor, folder / 'receptor-preparation', 'meeko', reviewed=True, cancelled=cancelled)['output'])
        reference = next(iter(Chem.SDMolSupplier(str(reference_path), removeHs=False)))
        archive = (Archive(folder.parent.parent / 'evaluation-archive', review['runtime_fingerprint'], cancelled, emit)
                   if review['runtime_fingerprint'].get('files') else None)
        pockets, dock_cache, reference_scores = {}, StageCache(archive), {}
        chosen = search.initial_candidates(review['candidates'], min(config['initial_trials'], config['max_trials']), config['seed'], config['method'])
        acquisition = []
        while chosen:
            science._cancel(cancelled)
            state['search']['round'] += 1
            state['search']['selected_ids'].extend(p['id'] for p in chosen)
            state['search']['rounds'].append({'round': state['search']['round'], 'ids': [p['id'] for p in chosen],
                                             'acquisition': acquisition, 'selection': 'diverse initialization' if not acquisition else 'expected improvement'})
            save()
            for index, record in enumerate(tasks):
                for pair in chosen:
                    science._cancel(cancelled)
                    state.update(active_record_index=index, records_per_chain=len(tasks))
                    row = evaluate(pair, record, index, folder, receptor, pdbqt, reference,
                                   pockets, dock_cache, reference_scores, state, save, cancelled, emit)
                    state['rows'].append(row)
                    save()
            if config['method'].endswith(('_v3', '_v4')):
                # Every observation must refer to the current common successful set.
                selected = {p['id']: p for p in review['candidates']}
                state['search']['observations'] = [search.observe(review, state, selected[key])
                                                   for key in state['search']['selected_ids']]
            else:
                for pair in chosen:
                    state['search']['observations'].append(search.observe(review, state, pair))
            state['metrics'] = summary_metrics(review, state['rows'])
            remaining = config['max_trials'] - len(state['search']['observations'])
            if remaining <= 0:
                break
            state['stage'] = 'optimizing'
            save()
            chosen, acquisition = search.propose(review['candidates'], state['search']['observations'], config,
                                                  min(config['round_size'], remaining))
        failed = sum(r['status'] != 'completed' for r in state['rows'])
        state.update(status='completed' if not failed else 'partial' if failed < len(state['rows']) else 'failed', stage='finished')
    except Exception as exc:
        state.update(status='cancelled' if cancelled() else 'failed', stage='stopped', error=str(exc))
    finally:
        accounted = {(r['combination'], bench.measurement_key(r)) for r in state['rows']}
        search_finished = state['status'] in {'completed', 'partial', 'failed'} and len(state['search']['observations']) == config['max_trials']
        for pair in review['candidates']:
            skipped = search_finished and pair['id'] not in state['search']['selected_ids']
            for record in tasks:
                if (pair['id'], bench.measurement_key(record)) in accounted:
                    continue
                state['rows'].append({'combination': pair['id'], 'docking': pair['docking']['tool'], 'scoring': pair['scoring']['tool'],
                    'compound_id': record['compound_id'], 'measurement_id': bench.measurement_key(record), 'cohort': record.get('cohort'),
                    'pocket': pair.get('pocket', bench.pockets.REFERENCE)['tool'], 'kind': record['kind'], 'p_value': record.get('p_value'),
                    'endpoint': record.get('endpoint'), 'status': 'not_selected' if skipped else 'not_run', 'score': None,
                    'error': 'Not selected within Bayesian search budget' if skipped else state.get('error', 'Not run')})
        state['metrics'] = summary_metrics(review, state['rows'])
        state['search']['finished'] = search_finished
        if search_finished and config['auto_lock']:
            state['decision_status'] = 'pending'
        save()
    if state['search']['finished'] and config['auto_lock']:
        try:
            decision = decisions.finalize_best_available(folder, selection_origin='bayesian_search')
            state['decision_status'] = 'locked' if decision else 'no_results'
        except (OSError, ValueError) as exc:
            state['decision_error'] = str(exc)
            state['decision_status'] = 'failed'
        save()
