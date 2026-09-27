"""Evidence-derived pipeline comparison and immutable, explicitly chosen decisions."""
from __future__ import annotations

from copy import deepcopy
from functools import lru_cache
import json
import math
from pathlib import Path

from . import discovery_benchmark as benchmark, science

POLICY = {
    'version': 1,
    'ranking': 'Eligible complete chains first; descending Spearman rho, then Kendall tau; ties by candidate ID. No runtime tie-break because poses may be reused.',
    'lock_gate': 'Terminal run, every expected compound and crystal control scored exactly once, all physical checks passed, cognate crystal RMSD <= 2 Å, finite rank correlations.',
    'green': 'Lock gate plus >= 5 measurements and positive lower bound of the 95% stratified-bootstrap interval for mean within-assay Spearman (2000 resamples, seed 42).',
    'yellow': 'Some results exist, but the green evidence gate is not met. Inspect failures, pose checks and sample size.',
    'gray': 'No result, excluded combination, or UI-only tool.',
    'interpretation': 'Evidence strength on this evaluation set (within-assay metrics averaged equally across rankable assays), not calibrated probability of binding or independent predictive validation. Bootstrap does not remove analog dependence or pretrained training overlap.',
}


AVAILABLE_POLICY = {
    'version': 1,
    'ranking': 'Passing chains first; then complete measured-set correlations, measurement coverage, Spearman rho, Kendall tau, candidate ID.',
    'scope': 'Best available tested chain. Failed checks remain failures; selection does not establish validation.',
}


SHARED_POLICY = {
    'version': 2,
    'ranking': 'Chains with passing scored-pose and crystal checks first; descending Spearman and Kendall on the identical shared successful measurement set, then coverage and candidate ID.',
    'missing': 'Failed ligands are skipped, never imputed. Zero-success chains do not remove measurements from other comparisons. Fewer than three varying labels per assay cannot produce a correlation.',
    'scope': 'Exploratory comparison on a success-conditioned common subset; missingness can bias it. Report every chain coverage and skipped count; independent validation unavailable.',
}


CORRELATION_POLICY = {
    **SHARED_POLICY,
    'version': 3,
    'ranking': 'Descending Spearman, then Kendall on the identical shared successful measurement set; then coverage and candidate ID. Physical and crystal checks are reported separately and do not override correlation ranking.',
    'scope': 'Best observed correlation among evaluated chains, not a validated pipeline. Pose failures and missingness remain explicit; independent validation unavailable.',
}


def _trusted_measured(review, pair, measured):
    expected = {benchmark.measurement_key(r): r for r in review['records']}
    ids = [benchmark.measurement_key(r) for r in measured]
    if len(ids) != len(expected) or set(ids) != set(expected):
        return False
    for row in measured:
        label = expected[benchmark.measurement_key(row)]
        if (row.get('docking') != pair['docking']['tool'] or row.get('scoring') != pair['scoring']['tool']
                or ('pocket' in pair and row.get('pocket') != pair['pocket']['tool'])
                or row.get('p_value') != label['p_value'] or row.get('endpoint') != label['endpoint']):
            return False
        if review.get('schema', 1) >= 2 and (row.get('cohort') != label.get('cohort') or row.get('compound_id') != label['compound_id']):
            return False
    return True


def finite(value):
    return isinstance(value, (float, int)) and not isinstance(value, bool) and math.isfinite(value)


@lru_cache(maxsize=512)
def rank_statistics(scores, measurements):
    import numpy as np
    from scipy.stats import rankdata, spearmanr, kendalltau
    if len(scores) < 3 or len(set(scores)) < 2 or len(set(measurements)) < 2:
        return {'rho': None, 'tau': None, 'ci_low': None, 'ci_high': None}
    result = {'rho': float(spearmanr(scores, measurements).statistic),
              'tau': float(kendalltau(scores, measurements).statistic), 'ci_low': None, 'ci_high': None}
    if len(scores) >= 5:
        indices = np.random.default_rng(42).integers(0, len(scores), (2000, len(scores)))
        x = rankdata(np.asarray(scores)[indices], axis=1)
        y = rankdata(np.asarray(measurements)[indices], axis=1)
        x -= x.mean(axis=1, keepdims=True)
        y -= y.mean(axis=1, keepdims=True)
        denominator = np.sqrt((x*x).sum(axis=1) * (y*y).sum(axis=1))
        valid = denominator > 0
        values = (x[valid] * y[valid]).sum(axis=1) / denominator[valid]
        if len(values) >= 1900:
            result['ci_low'], result['ci_high'] = map(float, np.quantile(values, [.025, .975]))
    return result


@lru_cache(maxsize=256)
def _stratified_cached(groups, intervals):
    import numpy as np
    from scipy.stats import spearmanr, kendalltau, rankdata
    details, bootstraps = [], []
    missing_prediction = False
    rng = np.random.default_rng(42)
    for cohort, x, y in groups:
        valid_labels = len(y) >= 3 and len(set(y)) >= 2
        rho = tau = None
        if valid_labels and len(set(x)) >= 2:
            rho, tau = float(spearmanr(x, y).statistic), float(kendalltau(x, y).statistic)
            if intervals:
                idx = rng.integers(0, len(x), (2000, len(x)))
                xx, yy = rankdata(np.asarray(x)[idx], axis=1), rankdata(np.asarray(y)[idx], axis=1)
                xx -= xx.mean(axis=1, keepdims=True); yy -= yy.mean(axis=1, keepdims=True)
                denom = np.sqrt((xx*xx).sum(axis=1) * (yy*yy).sum(axis=1))
                values = np.full(2000, np.nan)
                np.divide((xx*yy).sum(axis=1), denom, out=values, where=denom > 0)
                bootstraps.append(values)
        elif valid_labels:
            missing_prediction = True
        reason = ('Fewer than 3 measurements' if len(y) < 3 else 'Constant experimental values' if len(set(y)) < 2
                  else 'Constant model scores' if len(set(x)) < 2 else '')
        details.append({'cohort': cohort, 'measurements': len(y), 'rho': rho, 'tau': tau, 'exclusion': reason})
    ranked = [d for d in details if finite(d['rho']) and finite(d['tau'])]
    result = {'rho': None, 'tau': None, 'ci_low': None, 'ci_high': None,
              'ranked_cohorts': len(ranked), 'total_cohorts': len(groups),
              'ranked_measurements': sum(d['measurements'] for d in ranked),
              'small_sample': any(d['measurements'] < 5 for d in ranked), 'cohort_metrics': details}
    if ranked and not missing_prediction:
        result.update(rho=float(np.mean([d['rho'] for d in ranked])), tau=float(np.mean([d['tau'] for d in ranked])))
        if intervals and bootstraps and sum(d['measurements'] for d in ranked) >= 5:
            boot = np.mean(np.array(bootstraps), axis=0)
            boot = boot[np.isfinite(boot)]
            if len(boot) >= 1900:
                result['ci_low'], result['ci_high'] = map(float, np.quantile(boot, [.025, .975]))
    return result


def stratified_statistics(rows, sign, intervals=True):
    from collections import defaultdict
    groups = defaultdict(list)
    for row in rows:
        groups[row.get('cohort') or 'legacy'].append(row)
    frozen = tuple((key, tuple(sign * r['score'] for r in values), tuple(r['p_value'] for r in values))
                   for key, values in sorted(groups.items()))
    return _stratified_cached(frozen, intervals)


def evidence_digest(review, state):
    return benchmark.digest({'plan': review['sha256'], 'status': state['status'], 'rows': state.get('rows', [])})


def compare(review, state):
    """Keep strict validation eligibility separate from best-available selection."""
    expected = {benchmark.measurement_key(r): r for r in review['records']}
    reference_id = 'crystal:' + review['reference']['component']
    terminal = state.get('status') in {'completed', 'partial', 'failed'}
    result = []
    successful_sets = []
    for pair in review['candidates']:
        measured = [r for r in state.get('rows', []) if r.get('combination') == pair['id'] and r.get('kind') == 'measured']
        if _trusted_measured(review, pair, measured):
            successful = {benchmark.measurement_key(r) for r in measured if r.get('status') == 'completed' and finite(r.get('score'))}
            if successful:
                successful_sets.append(successful)
    shared = set.intersection(*successful_sets) if successful_sets else set()
    for pair in review['candidates']:
        rows = [r for r in state.get('rows', []) if r.get('combination') == pair['id']]
        measured = [r for r in rows if r.get('kind') == 'measured']
        control = [r for r in rows if r.get('kind') == 'reference']
        ids = [benchmark.measurement_key(r) for r in measured]
        reasons = []
        if set(ids) != set(expected) or len(ids) != len(expected):
            reasons.append('Incomplete or duplicate measured-compound accounting')
        if len(control) != 1 or control[0].get('compound_id') != reference_id:
            reasons.append('Missing or duplicate crystal control')
        if len(rows) != len(expected) + 1:
            reasons.append('Incomplete or unexpected result rows')
        finished = [r for r in rows if r.get('status') == 'completed' and finite(r.get('score'))]
        if len(finished) != len(expected) + 1:
            reasons.append('Failed, unrun or nonfinite scores')
        if any(r.get('docking') != pair['docking']['tool'] or r.get('scoring') != pair['scoring']['tool'] or ('pocket' in pair and r.get('pocket') != pair['pocket']['tool']) for r in rows):
            reasons.append('Result tool identity differs from frozen protocol')
        if any(r.get('p_value') != expected[benchmark.measurement_key(r)]['p_value']
               or r.get('endpoint') != expected[benchmark.measurement_key(r)]['endpoint']
               or (review.get('schema', 1) >= 2 and (r.get('cohort') != expected[benchmark.measurement_key(r)]['cohort']
                   or r.get('compound_id') != expected[benchmark.measurement_key(r)]['compound_id']))
               for r in measured if benchmark.measurement_key(r) in expected):
            reasons.append('Experimental labels differ from frozen test set')
        integrity_ok = not any(reason in reasons for reason in ('Incomplete or duplicate measured-compound accounting', 'Result tool identity differs from frozen protocol', 'Experimental labels differ from frozen test set'))
        stats = {'rho': None, 'tau': None, 'ci_low': None, 'ci_high': None, 'ranked_cohorts': 0,
                 'total_cohorts': len({r.get('cohort') or 'legacy' for r in expected.values()}),
                 'ranked_measurements': 0, 'small_sample': False, 'cohort_metrics': []}
        own_stats = dict(stats)
        if integrity_ok:
            sign = 1 if benchmark.SCORES[pair['scoring']['tool']][2] else -1
            ordered = sorted((r for r in measured if r.get('status') == 'completed' and finite(r.get('score'))), key=benchmark.measurement_key)
            frozen = [{**r, 'p_value': expected[benchmark.measurement_key(r)]['p_value'],
                       'cohort': expected[benchmark.measurement_key(r)].get('cohort')} for r in ordered]
            own_stats = stratified_statistics(frozen, sign)
            stats = stratified_statistics([r for r in frozen if benchmark.measurement_key(r) in shared], sign)
            stats = {**stats, 'total_cohorts': len({r.get('cohort') or 'legacy' for r in expected.values()})}
        reference_record = review['reference'].get('record')
        if reference_record and any(not benchmark.compatible_reference(record, reference_record) for record in review['records']):
            reasons.append('Test set includes incompatible or unresolved crystal/assay domains or mutations')
        physical = sum(r.get('physical_pass') is True for r in finished)
        if physical != len(expected) + 1:
            reasons.append('Not all physical checks passed')
        rmsd = control[0].get('rmsd_angstrom') if len(control) == 1 else None
        if not finite(rmsd) or rmsd < 0 or rmsd > 2 or len(control) != 1 or control[0].get('reference_comparable') is not True:
            reasons.append('Crystal pose recovery missing or above 2 Å')
        if not finite(stats['rho']) or not finite(stats['tau']):
            reasons.append('Ranking correlation unavailable')
        selection_eligible = (terminal and integrity_ok and len(rows) == len(expected) + 1 and finite(stats['rho']) and finite(stats['tau'])
                              and physical == len(finished) and bool(finished) and len(control) == 1
                              and control[0].get('compound_id') == reference_id and finite(control[0].get('score'))
                              and control[0].get('status') == 'completed' and control[0].get('reference_comparable') is True
                              and finite(rmsd) and 0 <= rmsd <= 2
                              and not (reference_record and any(not benchmark.compatible_reference(r, reference_record) for r in review['records'])))
        lockable = terminal and not reasons
        positive_evidence = terminal and not reasons and stats['ranked_measurements'] >= 5 and finite(stats['ci_low']) and stats['ci_low'] > 0
        result.append({'id': pair['id'], 'pocket': pair.get('pocket', benchmark.pockets.REFERENCE)['tool'], 'docking': pair['docking']['tool'], 'scoring': pair['scoring']['tool'],
                       'expected': len(expected), 'completed': sum(r.get('status') == 'completed' and finite(r.get('score')) for r in measured),
                       'failed': sum(r.get('status') in {'failed', 'not_run'} for r in rows), 'physical_passed': physical,
                       'physical_expected': len(expected) + 1, 'rmsd': rmsd if finite(rmsd) else None,
                       'seconds': sum(r['runtime_seconds'] for r in rows if finite(r.get('runtime_seconds'))),
                       **stats, 'own_rho': own_stats['rho'], 'own_tau': own_stats['tau'],
                       'shared_measurements': len(shared), 'shared_measurement_ids': sorted(shared),
                       'skipped': len(expected) - sum(r.get('status') == 'completed' and finite(r.get('score')) for r in measured),
                       'selection_eligible': selection_eligible, 'lockable': lockable, 'available': terminal and integrity_ok and any(r.get('status') in {'completed', 'failed'} for r in rows), 'color': 'green' if positive_evidence else 'yellow' if any(r.get('status') not in {'not_selected', 'not_run'} for r in rows) else 'gray',
                       'reasons': list(dict.fromkeys(reasons)), 'pair': pair})
    result.sort(key=lambda r: (not r['lockable'], -(r['rho'] if finite(r['rho']) else -2), -(r['tau'] if finite(r['tau']) else -2), r['id']))
    rank = 0
    for row in result:
        if row['lockable']:
            rank += 1
        row['rank'] = rank if row['lockable'] else None
    available = sorted((r for r in result if r['available']), key=available_key)
    for row in result:
        row['available_rank'] = next((i for i, r in enumerate(available, 1) if r['id'] == row['id']), None)
    return result


def available_key(row):
    return (not finite(row['rho']),
            -(row['rho'] if finite(row['rho']) else -2), -(row['tau'] if finite(row['tau']) else -2), -row['completed'], row['id'])


def context(folder):
    folder = Path(folder)
    review = json.loads((folder / 'plan.json').read_text())
    benchmark.validate_review(review)
    state = benchmark.read_status(folder)
    if state.get('plan_sha256') != review['sha256']:
        raise ValueError('Results do not belong to this frozen evaluation plan')
    return review, state


def read_lock(folder, review=None, state=None):
    path = Path(folder) / 'decision.json'
    if not path.exists():
        return None
    if review is None or state is None:
        review, state = context(folder)
    lock = json.loads(path.read_text())
    if benchmark.digest({k: v for k, v in lock.items() if k != 'sha256'}) != lock.get('sha256'):
        raise ValueError('Final decision integrity check failed')
    if lock.get('plan_sha256') != review['sha256'] or lock.get('evidence_sha256') != evidence_digest(review, state):
        raise ValueError('Evaluation changed after the chain was locked; the old decision cannot be reused')
    if lock.get('policy') != (CORRELATION_POLICY if lock.get('selection_kind') == 'correlation_first' else SHARED_POLICY if lock.get('selection_kind') == 'shared_subset' else AVAILABLE_POLICY if lock.get('selection_kind') == 'best_available' else POLICY):
        raise ValueError('Decision policy changed; inspect the archived decision')
    candidate = next((p for p in review['candidates'] if p['id'] == lock['candidate']['id']), None)
    if candidate != lock['candidate']:
        raise ValueError('Locked tools differ from the evaluated candidate')
    return lock


def _make_lock(folder, candidate_id, selection_origin='user', allow_best_available=False):
    folder = Path(folder)
    review, state = context(folder)
    rows = compare(review, state)
    selected = next((r for r in rows if r['id'] == candidate_id), None)
    fallback = bool(selected and allow_best_available and selected['available'] and not selected['lockable'])
    if not selected or not (selected['lockable'] or fallback):
        raise ValueError('Only a completed chain passing the ranking, physical and reference-pose gates can be locked')
    for row in state['rows']:
        if row.get('combination') != candidate_id or row.get('status') != 'completed':
            continue
        path = (folder / row.get('pose_path', '')).resolve()
        if not path.is_relative_to(folder.resolve()) or not path.is_file() or science.sha256(path) != row.get('pose_sha256'):
            raise ValueError('Selected chain pose files failed integrity verification')
    lock = {'schema': 1, 'plan_sha256': review['sha256'], 'evidence_sha256': evidence_digest(review, state),
            'candidate': deepcopy(selected['pair']), 'rank_on_this_set': selected['available_rank'], 'policy': CORRELATION_POLICY,
            'selection_kind': 'correlation_first',
            'checks_passed': selected['lockable'], 'failed_checks': selected['reasons'],
            'metrics': {k: selected[k] for k in ('rho', 'tau', 'ci_low', 'ci_high', 'rmsd', 'physical_passed', 'expected', 'completed', 'color', 'ranked_measurements', 'ranked_cohorts', 'total_cohorts', 'small_sample', 'shared_measurements', 'skipped', 'own_rho')},
            'comparison_measurement_ids': selected['shared_measurement_ids'],
            'preparation': review['protocol']['preparation'], 'pocket': review['protocol']['pocket'],
            'reference_sha256': benchmark.digest(review['reference']), 'identity': review['identity'],
            'selection_origin': selection_origin,
            'scope': ('Bayesian-search-selected' if selection_origin == 'bayesian_search' else 'User-selected') + ' exploratory chain on this evaluation set; independent predictive validation not established',
            'created_at': science._now()}
    lock['sha256'] = benchmark.digest(lock)
    return lock


def lock_chain(folder, candidate_id, selection_origin='user', allow_best_available=False):
    folder = Path(folder)
    if read_lock(folder):
        raise ValueError('This evaluation already has a fixed final chain. A new evaluation is required to change it.')
    lock = _make_lock(folder, candidate_id, selection_origin, allow_best_available)
    # Exclusive create: a second click/process cannot replace an existing decision.
    with (folder / 'decision.json').open('x') as handle:
        json.dump(lock, handle, indent=2, allow_nan=False)
    return lock


def finalize_best_available(folder, selection_origin='automatic'):
    """Fix an immutable result without hiding failed gates or inventing metrics."""
    review, state = context(folder)
    existing = read_lock(folder, review, state)
    if existing:
        return existing
    choices = sorted((r for r in compare(review, state) if r['available']), key=available_key)
    if not choices:
        return None
    try:
        return lock_chain(folder, choices[0]['id'], selection_origin, allow_best_available=True)
    except FileExistsError:
        return read_lock(folder)


def revise_best_selection(folder, reason, expected_previous_sha256):
    """Explicit correction only: archive the previous decision; preserve all evidence."""
    import fcntl
    folder = Path(folder)
    if not reason or not reason.strip():
        raise ValueError('A decision revision requires an explicit reason')
    with (folder / '.decision-revision.lock').open('a') as guard:
        fcntl.flock(guard, fcntl.LOCK_EX)
        review, state = context(folder)
        old = read_lock(folder, review, state)
        if not old or old['sha256'] != expected_previous_sha256:
            raise ValueError('Previous decision changed; inspect it before revising')
        choices = sorted((r for r in compare(review, state) if r['available']), key=available_key)
        if not choices or not finite(choices[0]['rho']):
            raise ValueError('No trustworthy comparable correlation is available')
        new = _make_lock(folder, choices[0]['id'], 'user_requested_correction', True)
        new['revision'] = {'previous_sha256': old['sha256'], 'previous_candidate_id': old['candidate']['id'],
                           'reason': reason.strip(), 'evidence_unchanged': True}
        new['sha256'] = benchmark.digest({k: v for k, v in new.items() if k != 'sha256'})
        # Recheck evidence and decision immediately before publishing the revision.
        latest_review, latest_state = context(folder)
        if evidence_digest(latest_review, latest_state) != old['evidence_sha256'] or read_lock(folder)['sha256'] != old['sha256']:
            raise ValueError('Evidence or decision changed during revision')
        history = folder / 'decision-history'
        history.mkdir(exist_ok=True)
        archived = history / (old['sha256'] + '.json')
        original = (folder / 'decision.json').read_bytes()
        if archived.exists():
            if archived.read_bytes() != original:
                raise ValueError('Archived decision differs from the original')
        else:
            with archived.open('xb') as handle:
                handle.write(original)
        science._write_json(folder / 'decision.json', new)
        return new
