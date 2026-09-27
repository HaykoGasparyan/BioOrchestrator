"""Small, deterministic GP/EI search over frozen categorical tool configurations.

The surrogate objective is within-assay rank correlation after feasibility
checks, never raw scores from different affinity models. No accuracy claim is
made for this surrogate or for selection on the same evaluation set.
"""
import math

from . import discovery_benchmark as benchmark


def configuration(candidates, method="categorical_gp_expected_improvement_v4"):
    if method not in {"categorical_gp_expected_improvement_v1", "categorical_gp_expected_improvement_v2", "categorical_gp_expected_improvement_v3", "categorical_gp_expected_improvement_v4"}:
        raise ValueError("Unsupported Bayesian search method")
    count = len(candidates)
    return {'method': method, 'seed': 42,
            'max_trials': min(count, max(6, min(12, math.ceil(2 * math.sqrt(count))))),
            'initial_trials': min(3, count), 'round_size': 2,
            'length_scale': .7, 'noise_variance': .0025, 'improvement_margin': .01,
            'objective': ('Mean within-assay Spearman on the identical shared successful measurements; physical and crystal checks are reported separately, not feasibility penalties. Missing correlations receive -1.25; failed ligands are skipped. Recompute earlier observations after each round.' if method.endswith('_v4') else 'Mean within-assay Spearman on the identical shared successful measurements, with passing scored-pose and crystal checks; skip failed ligands and report coverage. Otherwise feasibility penalty -1.25. Recompute earlier observations after each round.' if method.endswith('_v3') else 'Mean within-assay Spearman rho for a complete chain passing all physical and crystal RMSD gates; otherwise feasibility penalty -1.25, never a fabricated affinity.'),
            'execution_order': 'ligand_first_within_search_round', 'auto_lock': True,
            'scope': 'Selected tool/configuration combinations only; no continuous hyperparameter tuning or held-out predictive validation.'}


def tokens(pair):
    return tuple(benchmark.digest(pair.get(key, benchmark.pockets.REFERENCE))
                 for key in ('pocket', 'docking', 'scoring'))


def distance(left, right):
    return sum(a != b for a, b in zip(tokens(left), tokens(right))) / 3


def initial_candidates(candidates, count, seed=42, method="categorical_gp_expected_improvement_v4"):
    ordered = sorted(candidates, key=lambda p: benchmark.digest([seed, p['id']]))
    chosen = [ordered.pop(0)] if ordered and count else []
    while ordered and len(chosen) < count:
        def diversity(pair):
            distance_score = min(distance(pair, c) for c in chosen)
            if method == 'categorical_gp_expected_improvement_v1':
                return (distance_score,)
            # Show every docking engine early, then cover pockets and scorers.
            unseen = tuple(int(pair.get(k, benchmark.pockets.REFERENCE)['tool'] not in
                               {c.get(k, benchmark.pockets.REFERENCE)['tool'] for c in chosen})
                           for k in ('docking', 'pocket', 'scoring'))
            return (*unseen, distance_score)
        next_pair = max(ordered, key=diversity)
        ordered.remove(next_pair)
        chosen.append(next_pair)
    return chosen


def posterior(candidates, observations, config):
    """Exact finite-domain GP posterior and expected improvement (maximization)."""
    import numpy as np
    from scipy.special import ndtr
    mapping = {p['id']: p for p in candidates}
    observed = [mapping[row['id']] for row in observations]
    y = np.asarray([row['objective'] for row in observations], dtype=float)
    if not len(y) or not np.isfinite(y).all():
        raise ValueError('The Bayesian surrogate needs finite recorded objective observations')
    kernel = lambda left, right: np.exp(-np.asarray([[distance(a, b) for b in right] for a in left]) / config['length_scale'])
    covariance = kernel(observed, observed) + np.eye(len(observed)) * config['noise_variance']
    cholesky = np.linalg.cholesky(covariance)
    cross = kernel(candidates, observed)
    mean = cross @ np.linalg.solve(cholesky.T, np.linalg.solve(cholesky, y))
    solved = np.linalg.solve(cholesky, cross.T)
    std = np.sqrt(np.maximum(1 - np.sum(solved * solved, axis=0), 1e-12))
    improvement = mean - np.max(y) - config['improvement_margin']
    z = improvement / std
    ei = improvement * ndtr(z) + std * np.exp(-.5 * z*z) / math.sqrt(2*math.pi)
    return [{'id': pair['id'], 'mean': float(mu), 'std': float(sigma), 'expected_improvement': float(max(0, acq))}
            for pair, mu, sigma, acq in zip(candidates, mean, std, ei)]


def propose(candidates, observations, config, count):
    """Greedy batch EI with posterior-mean fantasies; fantasies never become evidence."""
    pending, selected, trace = list(observations), [], []
    seen = {row['id'] for row in observations}
    for _ in range(min(count, len(candidates) - len(seen))):
        predictions = [r for r in posterior(candidates, pending, config) if r['id'] not in seen]
        choice = max(predictions, key=lambda r: (r['expected_improvement'], r['std'], r['id']))
        selected.append(next(p for p in candidates if p['id'] == choice['id']))
        trace.append(choice)
        seen.add(choice['id'])
        pending.append({'id': choice['id'], 'objective': choice['mean']})
    return selected, trace


def observe(review, state, pair):
    from .pipeline_decisions import compare
    # A finished search round can be evaluated before the overall run terminates.
    evidence = next(r for r in compare(review, {**state, 'status': 'completed'}) if r['id'] == pair['id'])
    eligible = evidence['selection_eligible'] if review.get('search', {}).get('method', '').endswith('_v3') else evidence['lockable']
    if review.get('search', {}).get('method', '').endswith('_v4'):
        eligible = evidence['available'] and evidence['rho'] is not None and evidence['tau'] is not None
    return {'id': pair['id'], 'objective': evidence['rho'] if eligible else -1.25,
            'feasible': eligible, 'physical_and_crystal_checks_passed': evidence['selection_eligible'], 'completed': evidence['completed'], 'skipped': evidence['skipped'],
            'comparison_measurement_ids': evidence['shared_measurement_ids'], 'rho': evidence['rho'], 'tau': evidence['tau'],
            'reasons': evidence['reasons']}
