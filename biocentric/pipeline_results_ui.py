"""Live chain map and explicit final decisions for the approved evaluation."""
from pathlib import Path
import json

from . import discovery_benchmark as benchmark, pipeline_decisions as decisions
from .pipeline_map import graph_data, graph_html, graph_height


def _labels(state):
    from .discovery_studio import PHASES, catalog
    return {key: item[0] for phase in PHASES for key, item in catalog(state, phase).items()}


def _issues(lock, rows):
    current = next((r['reasons'] for r in rows if r['id'] == lock['candidate']['id']), [])
    return '; '.join(dict.fromkeys([*lock.get('failed_checks', []), *current]))


def _correlation_coverage(st, result):
    st.caption(f"Shared successful set: {result['shared_measurements']}/{result['expected']} measurements · this chain skipped {result['skipped']}. All chains use the same comparison set.")
    if not decisions.finite(result['rho']):
        return
    message = (f"Correlation uses {result['ranked_measurements']}/{result['expected']} measurements "
               f"in {result['ranked_cohorts']}/{result['total_cohorts']} assays.")
    if result['small_sample']:
        message += ' Small sample: fewer than 5 measurements in a ranked assay; even ρ = 1 does not establish predictive accuracy.'
    if result['small_sample'] or result['ranked_measurements'] < result['expected']:
        st.warning(message)
    else:
        st.caption(message)


def _load(st):
    folder = st.session_state.get('studio_eval_run')
    if not folder:
        return None
    try:
        review, state = decisions.context(folder)
        lock = decisions.read_lock(folder, review, state)
        rows = decisions.compare(review, state)
        if not lock and not any(r['lockable'] for r in rows) and state.get('decision_status') != 'pending':
            lock = decisions.finalize_best_available(folder)
        return Path(folder), review, state, rows, lock
    except (OSError, ValueError, KeyError) as exc:
        st.error(str(exc))
        return None


def render_map(st):
    from .discovery_studio import selection
    import streamlit.components.v1 as components
    folder = st.session_state.get('studio_eval_run')
    initial = benchmark.read_status(folder) if folder else {}
    active = bool(folder and (initial['status'] in {'queued', 'running'} or initial.get('decision_status') == 'pending'))
    @st.fragment(run_every=2 if active else None)
    def panel():
        loaded = _load(st) if folder else None
        if folder and not loaded:
            return
        if loaded:
            _, review, state, rows, lock = loaded
            if active and state['status'] not in {'queued', 'running'} and state.get('decision_status') != 'pending':
                st.rerun()
            selected = review['selection']
            progress = benchmark.combination_progress(review, state)
            st.markdown('**' + benchmark.combination_caption(review) + '**')
            st.caption(('Final chain locked' if lock else state['status'].capitalize()) + f" · {progress['completed']} / {progress['total']} planned evaluations finished")
            if state['status'] == 'running':
                st.caption(f"{state.get('active_compound', '')} · {state.get('stage', '')}")
            if review.get('search'):
                st.caption(f"Search round {state.get('search', {}).get('round', 0)}")
                if review.get('excluded'):
                    st.caption('Unavailable: ' + '; '.join(sorted({r['reason'] for r in review['excluded']})))
                if state.get('records_per_chain') and state['status'] == 'running':
                    st.caption(f"Record {state.get('active_record_index', 0) + 1} / {state['records_per_chain']} · all active chains before the next ligand")
        else:
            review = st.session_state.get('studio_eval_review')
            selected = review['selection'] if review else selection(st.session_state)
            state, rows, lock = {}, [], None
            st.caption('Ready to start' if review else 'Candidate connections')
        data = graph_data(selected, _labels(st.session_state), review, state, rows, lock, st.session_state.get('studio_identity'))
        with st.container(key='studio_pipeline_map'):
            graph = components.declare_component('pipeline_graph', path=Path(__file__).with_name('components') / 'pipeline_graph')
            graph(html=graph_html(data), key='studio_pipeline_graph_' + (Path(folder).name if folder else 'draft'), default=None)
        with st.expander('How connection colors work'):
            st.caption('Colors describe complete paths on this evaluation set. Shared upstream links use the weakest covered path; they do not estimate individual tool quality. Cyan identifies activity, not confidence.')
            for key in ('green', 'yellow', 'gray', 'interpretation'):
                st.write(f"{key.capitalize()}: {decisions.POLICY[key]}")
            st.caption('Protein preparation is fixed. Each pocket method connects to the selected docking and scoring tools. Only executed paths gain evidence colors; unavailable adapters stay gray. Older experiments retain their original crystal-defined pocket.')
    panel()


def render_results(st):
    loaded = _load(st)
    if not loaded:
        if not st.session_state.get('studio_eval_run'):
            st.caption('Results appear after you click Start and the evaluation runs.')
        return
    folder, review, state, rows, lock = loaded
    st.caption('Ranked by shared-set Spearman, then Kendall. Pose checks remain separate · independent predictive validation unavailable')
    best_id = lock['candidate']['id'] if lock else next((r['id'] for r in rows if r['available_rank'] == 1), None)
    if lock:
        st.success(f"Fixed chain: Protein → {benchmark.chain_label(lock['candidate'])}")
        if not lock.get('checks_passed', True):
            st.warning('Selected by correlation · checks incomplete or failed. ' + _issues(lock, rows))
        _correlation_coverage(st, next(r for r in rows if r['id'] == best_id))
    rows = sorted(rows, key=lambda r: (r['available_rank'] is None, r['available_rank'] or 0))
    table = [{'Rank': r['available_rank'], 'Chain': benchmark.chain_label(r['pair']),
              'Decision': 'Locked' if lock and lock['candidate']['id'] == r['id'] else 'Best on shared set' if r['available_rank'] == 1 else 'Comparable' if r['available'] and decisions.finite(r['rho']) else 'Not selected' if any(x.get('combination') == r['id'] and x.get('status') == 'not_selected' for x in state['rows']) else 'Not eligible',
              'Measurements': f"{r['completed']}/{r['expected']}", 'Skipped': r['skipped'], 'Shared set': r['shared_measurements'],
              'Correlation measurements': f"{r['ranked_measurements']}/{r['expected']}",
              'Assays ranked': f"{r['ranked_cohorts']}/{r['total_cohorts']}", 'Spearman ρ': r['rho'],
              'ρ 95% interval': f"{r['ci_low']:.2f} … {r['ci_high']:.2f}" if r['ci_low'] is not None else 'Unavailable',
              'Kendall τ': r['tau'], 'Crystal RMSD (Å)': r['rmsd'], 'Physical passes': f"{r['physical_passed']}/{r['physical_expected']}",
              'Seconds': round(r['seconds'], 1), 'Evidence': r['color'], 'Issue': '; '.join(r['reasons'])} for r in rows]
    for pair in review.get('excluded', []):
        table.append({'Chain': benchmark.chain_label(pair), 'Decision': 'Excluded', 'Evidence': 'gray', 'Issue': pair['reason']})
    import pandas as pd
    evaluated = [i for i, row in enumerate(rows) if any(r.get('combination') == row['id'] and r.get('status') in {'completed', 'failed'} for r in state['rows'])]
    finalists = ([next(i for i, row in enumerate(rows) if row['id'] == best_id)] if best_id else [])
    finalists += [i for i in evaluated if i not in finalists][:5]
    if finalists:
        best_label = benchmark.chain_label(next(r['pair'] for r in rows if r['id'] == best_id)) if best_id else None
        frame = pd.DataFrame([table[i] for i in finalists])
        first = ['Rank', 'Chain', 'Spearman ρ', 'Kendall τ', 'Measurements', 'Skipped', 'Correlation measurements', 'Decision']
        frame = frame[first + [column for column in frame.columns if column not in first]]
        styled = frame.style.apply(lambda row: ['background-color: #e7f8ee; color: #176a47; font-weight: 600' if row['Chain'] == best_label else '' for _ in row], axis=1)
        st.caption('Fixed chain and up to five evaluated alternatives' if best_id else 'Evaluated candidates · no eligible final chain yet')
        st.dataframe(styled, hide_index=True, width='stretch')
    with st.expander('All candidates'):
        st.dataframe(table, hide_index=True, width='stretch')
    if review.get('search'):
        with st.expander('Bayesian search'):
            st.caption(f"{len(state.get('search', {}).get('observations', []))} evaluated / {len(review['candidates'])} available · same frozen evaluation set for every trial")
            st.latex(r'k(x,x\prime)=\exp[-d_{\mathrm{Hamming}}(x,x\prime)/0.7],\qquad EI(x)=(\mu-f^*-\xi)\Phi(z)+\sigma\phi(z)')
            st.caption('z = (μ − f* − ξ) / σ; ξ = 0.01. Surrogate uncertainty guides exploration; it is not biological confidence. Only the selected tool configurations are searched.')
            st.write(review['search']['objective'])
            st.dataframe(state.get('search', {}).get('observations', []), hide_index=True, width='stretch')
            st.caption('Method: [Gaussian-process Bayesian optimization and expected improvement](https://arxiv.org/abs/1012.2599). Selection on this set is exploratory.')
    if lock:
        with st.expander('Locked parameters & decision'):
            st.json(lock, expanded=False)
        st.download_button('Download final decision', json.dumps(lock, indent=2), 'bio-coordinator-final-chain.json', 'application/json', key='studio_final_download')
    elif state['status'] in {'queued', 'running'}:
        st.caption('The full evaluation must finish before a final chain can be locked.')
    else:
        eligible = {r['id']: r for r in rows if r['available'] and decisions.finite(r['rho'])}
        if eligible:
            selected = st.selectbox('Final chain', list(eligible), format_func=lambda key: f"#{eligible[key]['available_rank']} · {benchmark.chain_label(eligible[key]['pair'])}", key='studio_final_choice_' + folder.name)
            st.caption('Locking saves these tools and parameters for this evaluation. Changing them requires a new evaluation.')
            if st.button('Lock final chain', type='primary', key='studio_lock_final'):
                try:
                    decisions.lock_chain(folder, selected, allow_best_available=True)
                    st.rerun()
                except (OSError, ValueError) as exc:
                    st.error(str(exc))
        else:
            st.info('No chain meets the completion, pose and ranking checks needed to lock a final decision.')
    with st.expander('Scores on each chain’s own successful subset'):
        st.caption('Diagnostic only: subsets may differ, so these numbers do not select the winner.')
        st.dataframe([{'Chain': benchmark.chain_label(r['pair']), 'Scored': r['completed'], 'Skipped': r['skipped'], 'Own-subset Spearman': r['own_rho']} for r in rows], hide_index=True)
    with st.expander('Per-assay metrics'):
        assay_labels = {r.get('cohort') or 'legacy': {'Assay': r.get('assay'), 'Endpoint': r.get('endpoint')} for r in review['records']}
        st.dataframe([{'Chain': benchmark.chain_label(r['pair']), **assay_labels.get(m['cohort'], {}), **m}
                      for r in rows for m in r.get('cohort_metrics', [])], hide_index=True)
    with st.expander('Comparison policy'):
        st.json(decisions.CORRELATION_POLICY, expanded=False)
        st.caption('This decision policy ranks observed results. It is not a prespecified held-out validation. Runtime includes reused docking poses and is not a fair speed comparison.')


def render_final_chain(st):
    loaded = _load(st)
    if not loaded:
        return
    _, review, state, rows, lock = loaded
    if not lock:
        st.info('The best eligible chain will be fixed here after search and validation finish.' if state['status'] in {'queued', 'running'} or state.get('decision_status') == 'pending' else 'No usable result is available yet. Inspect execution failures in Final results.')
        if state.get('decision_error'):
            st.error(state['decision_error'])
        return
    import streamlit.components.v1 as components
    pair = lock['candidate']
    selected = {key: [pair[key]] for key in ('pocket', 'docking', 'scoring') if key in pair}
    selected['pocket_identification'] = selected.pop('pocket', [benchmark.pockets.REFERENCE])
    data = graph_data(selected, _labels(st.session_state), {**review, 'candidates': [pair]}, state, rows, lock)
    st.success('Fixed · ' + benchmark.chain_label(pair))
    if not lock.get('checks_passed', True):
        st.warning('Best available by correlation · ' + _issues(lock, rows))
    components.html(graph_html(data), height=graph_height(data) + 110, scrolling=True)
    result = next(r for r in rows if r['id'] == pair['id'])
    columns = st.columns(4)
    columns[0].metric('Spearman ρ', f"{result['rho']:.3f}" if decisions.finite(result['rho']) else "Unavailable")
    columns[1].metric('Kendall τ', f"{result['tau']:.3f}" if decisions.finite(result['tau']) else "Unavailable")
    columns[2].metric('Crystal RMSD', f"{result['rmsd']:.2f} Å" if decisions.finite(result['rmsd']) else "Unavailable")
    columns[3].metric('Physical passes', f"{result['physical_passed']}/{result['physical_expected']}")
    st.caption(f"{result['completed']}/{result['expected']} measurements · {result['ranked_cohorts']} ranked assays · independent predictive validation unavailable")
    _correlation_coverage(st, result)
    with st.expander('Fixed tools & parameters'):
        st.json(pair, expanded=False)


def render_final_visualizations(st):
    loaded = _load(st)
    if not loaded:
        if not st.session_state.get('studio_eval_run'):
            st.caption('Lock a completed chain to compare its poses with the reference.')
        return
    folder, review, state, _, lock = loaded
    if not lock:
        st.caption('Choose and lock a chain in Final results to open its visualizations.')
        return
    pair = lock['candidate']
    selected_rows = [r for r in state['rows'] if r['combination'] == pair['id']]
    selected_rows = [r for r in selected_rows if r.get('status') == 'completed' and r.get('pose_path') and decisions.finite(r.get('score'))]
    if not selected_rows:
        st.info('This chain has no successfully scored pose to visualize.')
        return
    selected_rows.sort(key=lambda r: (r['kind'] != 'reference', r['compound_id']))
    by_id = {benchmark.measurement_key(r): r for r in selected_rows}
    compound = st.selectbox('Complex', list(by_id), format_func=lambda key: by_id[key]['compound_id'] + (' · ' + str(by_id[key].get('endpoint') or '') + ' · ' + key if key != by_id[key]['compound_id'] else ''), key='studio_final_complex_' + folder.name)
    row = by_id[compound]
    st.caption(f"{benchmark.chain_label(pair)} · {review['reference']['pdb_id']} reference")
    import streamlit.components.v1 as components
    from .discovery_evaluation_ui import _scene
    try:
        components.html(_scene(folder, review, row), height=525, scrolling=False)
    except (OSError, ValueError) as exc:
        st.error(str(exc))
    st.caption('Orange: docked ligand · Cyan: experimental crystal ligand')
    if row['kind'] != 'reference':
        st.caption('Different ligands: binding-site comparison only. No experimental pose RMSD is inferred.')
    other = [r for r in state['rows'] if benchmark.measurement_key(r) == compound]
    st.dataframe([{'Chain': r.get('pocket', 'reference_box') + ' → ' + r['docking'] + ' → ' + r['scoring'], 'Final': r['combination'] == pair['id'],
                   'Docked score': r.get('score'), 'Score units': r.get('score_unit'), 'Crystal-pose score': r.get('reference_score'),
                   'Experimental pX': r.get('p_value'), 'Endpoint': r.get('endpoint'), 'RMSD (Å)': r.get('rmsd_angstrom'), 'Status': r['status']} for r in other], hide_index=True, width='stretch')
    st.caption('Raw scores are comparable only within the same scoring method. Experimental pX is a measured endpoint, not a docking energy.')
    measured = [r for r in selected_rows if r['kind'] == 'measured']
    if measured:
        assays = list(dict.fromkeys(r.get('cohort') or 'legacy' for r in measured))
        assay = st.selectbox('Plot assay', assays, key='studio_final_plot_assay_' + folder.name)
        measured = [r for r in measured if (r.get('cohort') or 'legacy') == assay]
        sign = 1 if benchmark.SCORES[pair['scoring']['tool']][2] else -1
        st.scatter_chart([{'Experimental pX': r['p_value'], 'Oriented model score': sign*r['score'], 'Compound': r['compound_id']} for r in measured], x='Experimental pX', y='Oriented model score', size=75, color='#8374dd')
        st.caption('Model scores are oriented so higher is preferred by the scoring method; their units remain distinct from experimental pX.')
