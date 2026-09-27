"""Evaluation review, explicit chat acceptance, and live benchmark presentation."""
from __future__ import annotations

import html
import json
from pathlib import Path
import re

from . import discovery_benchmark as benchmark
from . import discovery_evidence as evidence
from . import science


def request_evaluation(state):
    state['studio_eval_requested'] = True


def restore_run(st, root):
    run_id = st.query_params.get('evaluation', '')
    if run_id and st.session_state.get('studio_eval_detached') == run_id:
        st.query_params.pop('evaluation', None)
        return
    if st.session_state.get('studio_eval_run') or not re.fullmatch(r'[a-f0-9]{32}', run_id):
        return
    folder = Path(root) / 'evaluation-runs' / run_id
    if not (folder / 'plan.json').is_file():
        return
    review = json.loads((folder / 'plan.json').read_text())
    benchmark.validate_review(review)
    target = review['identity'].get('requested_pdb') or (review['identity'].get('gene_names') or [review['identity']['accession']])[0]
    st.session_state.update(studio_eval_run=str(folder), studio_eval_review=review,
        studio_eval_max_compounds=(review.get('ligand_subset') or {}).get('max_compounds') or 0,
        studio_eval_requested=True, studio_identity=review['identity'], studio_resolution={'status': 'resolved', 'candidates': [review['identity']]},
        studio_draft={'messages': [{'role': 'user', 'text': target}], 'target': target, 'request': target,
                      'organism': 'Human', 'status': 'review', 'evidence_query': (target, 'Human')})
    for phase, rows in review['selection'].items():
        st.session_state[f'studio_tools_{phase}'] = [r['tool'] for r in rows if not r['tool'].startswith('custom_')]
        for row in rows:
            st.session_state.setdefault('studio_parameters', {})[f"{phase}:{row['tool']}"] = row['parameters']


def _run(st, root):
    from .discovery_studio import selection
    review = st.session_state.get('studio_eval_review')
    if not review or review['selection'] != selection(st.session_state) or review['identity']['accession'] != st.session_state['studio_identity']['accession'] or ((review.get('ligand_subset') or {}).get('max_compounds') or 0) != st.session_state.get('studio_eval_max_compounds', 0):
        st.session_state['studio_eval_error'] = 'The target or tools changed. Review the evaluation plan again.'
        return
    if st.session_state.get('studio_eval_run'):
        status = benchmark.read_status(st.session_state['studio_eval_run'])
        if status['status'] in {'queued', 'running'}:
            return
    try:
        if st.session_state.get('experiments_current'):
            from .experiment_library import load
            load(root, st.session_state['experiments_current'])
        st.session_state['studio_eval_run'] = benchmark.start_run(root, review)
        st.query_params['evaluation'] = Path(st.session_state['studio_eval_run']).name
        st.session_state['studio_eval_requested'] = True
        st.session_state.pop('studio_eval_error', None)
        from .experiment_library_ui import save_current
        save_current(st, root)
    except (ValueError, OSError) as exc:
        st.session_state['studio_eval_error'] = str(exc)


def prepare_start(state, root, progress=lambda message: None):
    """Resolve the default benchmark in one action, retaining advanced choices."""
    from .discovery_studio import selection, selection_issues
    issues = selection_issues(state)
    if issues:
        raise ValueError('; '.join(issues))
    identity = state.get('studio_identity')
    if not identity:
        raise ValueError('Select an exact protein match in Experiments, then click Start.')
    selected = selection(state)
    requested_pdb = identity.get('requested_pdb')
    target = (state.get('studio_draft', {}).get('target') or '').strip()
    if re.fullmatch(r'[0-9][A-Za-z0-9]{3}', target):
        requested_pdb = target.upper()
    progress('Finding experimental measurements…')
    data = state.get('studio_evidence') or evidence.overview(identity, root)
    state['studio_evidence'] = data
    reference_rows = data.get('structures', {}).get('rows', [])
    if requested_pdb:
        reference_rows = evidence.structure_page(identity, root, entry_ids=[requested_pdb])['rows']
        if not reference_rows:
            raise ValueError(f'{requested_pdb} does not contain the resolved target. No other structure was substituted.')
        if not any(any((c.get('formula_weight') or 0) > 150 for c in r.get('compounds', [])) for r in reference_rows):
            raise ValueError(f'{requested_pdb} has no supported complete small-molecule crystal ligand. Peptide/polymer binders require a curated complex; a residue fragment or another structure will not be substituted.')
    dataset = state.get('studio_eval_dataset')
    if not dataset or dataset.get('cap') != state.get('studio_eval_limit', 500):
        dataset = benchmark.fetch_dataset(identity, root, data.get('binding', {}), state.get('studio_eval_limit', 500),
            progress=lambda n, total: progress(f'Finding measurements · {n} / {total} records'))
        state['studio_eval_dataset'] = dataset
    groups = {key: rows for key, rows in benchmark.cohorts(dataset).items() if len(rows) >= 3}
    order = sorted(groups, key=lambda key: (key != state.get('studio_eval_cohort'), benchmark.assay_mutations(groups[key][0]) != set(), -len(groups[key]), key))
    if not order:
        raise ValueError('No assay with at least three usable measurements was found. Open Evaluation data to increase the search limit.')
    progress('Choosing a compatible crystal structure…')
    errors = []
    attempts = 0
    for cohort in order:
        group = groups[cohort]
        structures = [r for r in reference_rows if r.get('target_chains') and r.get('models', 1) == 1
                      and any((c.get('formula_weight') or 0) > 150 for c in r.get('compounds', [])) and benchmark.compatible_reference(group[0], r)]
        structures.sort(key=lambda r: r['id'] != state.get('studio_eval_reference'))
        for record in structures:
            components = [c['id'] for c in record['compounds'] if (c.get('formula_weight') or 0) > 150]
            component = state.get('studio_eval_ligand_' + record['id'])
            component = component if component in components else components[0]
            chain = state.get('studio_eval_chain_' + record['id'])
            chain = chain if chain in record['target_chains'] else record['target_chains'][0]
            attempts += 1
            try:
                progress(f"Preparing {record['id']} · {component}…")
                reference = benchmark.reference_bundle(record, component, chain, root)
            except (ValueError, OSError, evidence.RetrievalError, science.ScienceError) as exc:
                errors.append(f"{record['id']}: {exc}")
                if attempts >= 8:
                    break
                continue
            limit = state.get('studio_eval_max_compounds', 0)
            review = benchmark.make_review(dataset, cohort, selected, reference, max_compounds=limit or None)
            state['studio_eval_review'] = review
            state['studio_eval_config'] = benchmark.digest([dataset['cap'], record['id'], component, chain, selected, limit])
            return review
        if attempts >= 8:
            break
    if requested_pdb:
        raise ValueError(f'No usable reference and compatible assay set for {requested_pdb}. No different PDB was substituted. ' + (errors[-1] if errors else 'Check construct/domain and ligand compatibility.'))
    raise ValueError('No usable ligand-bound reference was found. Choose another structure in Evaluation data.' + (' ' + errors[-1] if errors else ''))


def _request_start(state):
    state['studio_eval_start_requested'] = True
    state.pop('studio_eval_error', None)


def chat_review(st, root):
    from .discovery_studio import selection_issues
    run = st.session_state.get('studio_eval_run')
    if run:
        state = benchmark.read_status(run)
        st.caption('Evaluation · ' + state['status'])
        if state['status'] in {'queued', 'running'}:
            st.session_state.pop('studio_eval_start_requested', None)
            st.button('Stop', key='studio_eval_stop', width='stretch', on_click=benchmark.cancel_run, args=(run,))
            return
    if st.session_state.pop('studio_eval_start_requested', False):
        with st.status('Starting evaluation…', expanded=True) as progress:
            try:
                prepare_start(st.session_state, root, lambda message: progress.update(label=message))
                progress.update(label='Starting tools…')
                _run(st, root)
                if st.session_state.get('studio_eval_error'):
                    raise ValueError(st.session_state['studio_eval_error'])
                progress.update(label='Evaluation started', state='complete', expanded=False)
                st.rerun()
            except (ValueError, OSError, evidence.RetrievalError, science.ScienceError) as exc:
                st.session_state['studio_eval_error'] = str(exc)
                progress.update(label='Could not start', state='error')
    st.button('Start', key='studio_eval_start', type='primary', width='stretch', disabled=bool(selection_issues(st.session_state)),
              on_click=_request_start, args=(st.session_state,))
    if st.session_state.get('studio_eval_error'):
        st.error(st.session_state['studio_eval_error'])


def _table(records):
    return [{'Compound': r.get('compound_id'), 'Endpoint': r.get('endpoint'), 'Measured': f"{r.get('relation') or '?'} {r.get('value')} {r.get('units') or ''}",
             'p scale': r.get('p_value'), 'p relation': r.get('p_relation'), 'Assay': r.get('assay'),
             'Eligible': r.get('eligible'), 'Reason': '; '.join(r.get('reasons', [])), 'Paper': r.get('publication_url'), 'Source': r.get('source_url')} for r in records]


def render_data(st, identity, data, root):
    from .discovery_studio import selection
    run = st.session_state.get('studio_eval_run')
    if run and benchmark.read_status(run)['status'] in {'queued', 'running'}:
        st.info('The accepted test set is frozen while evaluation runs. Follow progress in Evaluation run.')
        return
    dataset = st.session_state.get('studio_eval_dataset')
    initial = data.get('binding', {})
    records = dataset['records'] if dataset else benchmark.normalize(initial.get('rows', []), identity['organism'])
    total = dataset['source_total'] if dataset else initial.get('total')
    eligible = sum(r['eligible'] for r in records)
    with st.container(key='studio_evaluation_data'):
        a, b, c = st.columns(3)
        a.metric('Source records', f'{total:,}' if total is not None else 'Unavailable')
        b.metric('Inspected', len(records))
        c.metric('Eligible exact values', eligible if dataset else '—')
        counts = {endpoint: sum(r.get('endpoint') == endpoint for r in records) for endpoint in ('Kd', 'Ki', 'IC50')}
        st.caption(' · '.join(f'{endpoint}: {count}' for endpoint, count in counts.items()) + ' · counts refer to inspected records')
        st.caption('pX = −log₁₀(X in mol/L). Higher means stronger measured activity. Assays and endpoints stay separate.')
        if initial.get('error') and not dataset:
            st.warning('Binding data source unavailable: ' + initial['error'])
        st.session_state.setdefault('studio_eval_limit', 500)
        amount = st.number_input('Evaluation data search · maximum records', min_value=1, step=1,
                                 key='studio_eval_limit')
        ligand_limit = st.selectbox('Evaluation ligands', [0, 30, 50, 100], key='studio_eval_max_compounds',
                                    format_func=lambda n: 'All eligible ligands' if not n else str(n))
        if st.button('Find evaluation data', key='studio_eval_find', type='primary'):
            with st.status('Reading literature-linked measurements…', expanded=True) as progress:
                try:
                    dataset = benchmark.fetch_dataset(identity, root, initial, amount,
                        progress=lambda n, total: progress.update(label=f'Inspecting {n} / {total} source records…'))
                    st.session_state['studio_eval_dataset'] = dataset
                    st.session_state.pop('studio_eval_review', None)
                    progress.update(label='Evaluation data ready', state='complete', expanded=False)
                    st.rerun()
                except (evidence.RetrievalError, OSError, ValueError) as exc:
                    progress.update(label='Source lookup failed', state='error')
                    st.error(str(exc))
        if dataset and not dataset['complete']:
            st.caption(f"Partial search: {dataset['loaded']:,} inspected of {total:,} source records; not an exhaustive literature test set.")
        if dataset and dataset['errors']:
            with st.expander('Source failures'):
                st.write(dataset['errors'])
        st.dataframe(_table(records), hide_index=True, width='stretch', height=220,
                     column_config={'Paper': st.column_config.LinkColumn('Paper'), 'Source': st.column_config.LinkColumn('Source')})
        if not dataset:
            return
        groups = {key: rows for key, rows in benchmark.cohorts(dataset).items() if len(rows) >= 3}
        if not groups:
            st.info('No single assay has at least three eligible compounds in this search. Increase the search limit or choose another target.')
            return
        structures = [r for r in data.get('structures', {}).get('rows', []) if r.get('target_chains')
                      and any((c.get('formula_weight') or 0) > 150 for c in r.get('compounds', []))
                      and r.get('models', 1) == 1
                      and any(benchmark.compatible_reference(rows[0], r) for rows in groups.values())]
        structures.sort(key=lambda r: bool(r.get('mutations')))
        if not structures:
            st.info('No ligand-bound reference with compatible mutation annotations on this structure page. Browse additional structures.')
            return
        mapped = {r['id']: r for r in structures}
        pdb_id = st.selectbox('Crystal reference', list(mapped), key='studio_eval_reference', format_func=lambda key: f"{key} · {mapped[key].get('resolution') or '—'} Å")
        record = mapped[pdb_id]
        left, right = st.columns(2)
        component = left.selectbox('Reference ligand', [r['id'] for r in record['compounds'] if (r.get('formula_weight') or 0) > 150], key=f'studio_eval_ligand_{pdb_id}')
        chain = right.selectbox('Protein chain', record['target_chains'], key=f'studio_eval_chain_{pdb_id}')
        included, excluded = benchmark.evaluation_records(dataset, record)
        try:
            included, subset_exclusions, subset = benchmark.ligand_subset(included, ligand_limit or None)
        except ValueError as exc:
            st.warning(str(exc))
            return
        excluded.extend(subset_exclusions)
        cohort_counts = {}
        for row in included:
            key = (row['assay'], row['endpoint'])
            cohort_counts[key] = cohort_counts.get(key, 0) + 1
        label = 'Evaluation subset' if ligand_limit else 'Full eligible set'
        st.caption(f"{label}: {len(included):,} measurements · {len({r['compound_id'] for r in included}):,} compounds · {len(cohort_counts)} assay/endpoint groups · {len(excluded):,} excluded")
        if ligand_limit:
            st.caption('Fixed subset shared by all combinations; crystal control added separately.')
        with st.expander('Assays & exclusions'):
            st.dataframe([{'Assay': key[0], 'Endpoint': key[1], 'Measurements': n} for key, n in cohort_counts.items()], hide_index=True)
            st.dataframe(excluded, hide_index=True)
        current_config = benchmark.digest([dataset['cap'], pdb_id, component, chain, selection(st.session_state), ligand_limit])
        if st.session_state.get('studio_eval_review') and st.session_state.get('studio_eval_config') != current_config and not run:
            st.session_state.pop('studio_eval_review', None)
            st.rerun()
        st.caption('Exploratory ranking: separate assays, one receptor conformation. Crystal/assay construct equivalence and learned-model training overlap remain unverified.')
        with st.expander('Evaluation protocol', expanded=False):
            st.write('All assay records for the selected ligands are evaluated. A ligand limit uses a fixed, score-independent subset shared by every combination. Spearman ρ and Kendall τ are calculated within each assay/endpoint group, then averaged equally across rankable groups. No cross-assay affinity pooling. All selected measurements must finish before a chain receives ranking metrics.')
            st.write('Reference redocking: symmetry-aware heavy-atom RMSD in the receptor frame, plus PoseBusters checks. Other ligands have no pose RMSD unless their identity matches the crystal ligand.')
            st.write('Fixed preprocessing: Meeko on the selected protein chain; waters, other chains and hetero residues omitted. Incomplete terminal residues are omitted only beyond 12 Å from the ligand, with exact changes listed before approval. No missing-atom repair. RDKit seed 42 retains input charge and tautomer. The selected pocket method defines each chain’s docking box.')
            st.write('Preparation remains fixed. Pocket methods run once per configured method and are reused across its docking/scoring chains. fpocket uses the selected predicted rank without the crystal ligand. Crystal ligand pocket is a known-site baseline. Unconnected adapters are explicitly excluded. Docking/validation commands have a 120 s timeout; BioCentric uses its configured timeout. Ligand-first search rounds continue until the search budget is complete or stopped.')
            st.write('All combinations are exploratory comparisons on this evaluation set. Selecting or tuning with these results requires a new independent holdout before claiming predictive validation.')
        if st.button('Use these settings', key='studio_eval_prepare'):
            with st.spinner('Checking experimental coordinates and executable combinations…'):
                try:
                    reference = benchmark.reference_bundle(record, component, chain, root)
                    review = benchmark.make_review(dataset, None, selection(st.session_state), reference, max_compounds=ligand_limit or None)
                    st.session_state['studio_eval_review'] = review
                    st.session_state['studio_eval_config'] = current_config
                    st.session_state.pop('studio_eval_run', None)
                    st.session_state.pop('studio_eval_error', None)
                    st.rerun()
                except Exception as exc:
                    st.error(str(exc))
        review = st.session_state.get('studio_eval_review')
        if review:
            st.success("Ready to start · " + benchmark.combination_caption(review))
            if review.get('search'):
                st.caption('Bayesian search · ligand-first rounds · best eligible evaluated chain fixed automatically')
            if review['reference'].get('omitted_terminal_residues'):
                omitted = review['reference']['omitted_terminal_residues']
                st.warning('Proposed receptor omits incomplete terminal residues farther than 12 Å from the ligand: ' + ', '.join(r['chain'] + ':' + r['residue'] + ' ' + r['name'] for r in omitted) + '. These changes are part of the approval.')
            st.dataframe([{'Pocket': r.get('pocket', benchmark.pockets.REFERENCE)['tool'], 'Docking': r['docking']['tool'], 'Scoring': r['scoring']['tool'], 'Status': 'Ready'} for r in review['candidates']] +
                         [{'Pocket': r.get('pocket', benchmark.pockets.REFERENCE)['tool'], 'Docking': r['docking']['tool'], 'Scoring': r['scoring']['tool'], 'Status': r['reason']} for r in review['excluded']], hide_index=True, width='stretch')
            with st.expander('Frozen plan'):
                st.json({k: v for k, v in review.items() if k not in {'reference', 'records'}}, expanded=False)
                st.download_button('Download evaluation plan', json.dumps(review, indent=2), 'evaluation-plan.json', 'application/json')


def _scene(folder, review, row):
    from .viewer import viewer_html
    path = (Path(folder) / row['pose_path']).resolve()
    if not path.is_relative_to(Path(folder).resolve()) or science.sha256(path) != row['pose_sha256']:
        raise ValueError('Pose file failed integrity verification')
    reference = review['reference']
    models = [{'role': 'receptor', 'format': 'pdb', 'text': reference['receptor_pdb'], 'name': reference['pdb_id']},
              {'role': 'reference', 'format': 'sdf', 'text': reference['reference_sdf'], 'name': 'Crystal ligand'},
              {'role': 'ligand', 'format': 'sdf', 'text': path.read_text(), 'name': row['compound_id']}]
    return viewer_html({'models': models, 'coordinate_frame': reference['pdb_id'] + ' deposited receptor frame',
                        'frame_source': 'Frozen experimental coordinates', 'interpretation': 'Orange: docked ligand. Cyan: crystal ligand. Different compounds are a binding-site comparison, not pose recovery.'})


def render_live(st, folder):
    initial = benchmark.read_status(folder)
    was_active = initial['status'] in {'queued', 'running'} or initial.get('decision_status') == 'pending'
    @st.fragment(run_every=2 if was_active else None)
    def live():
        state = benchmark.read_status(folder)
        if was_active and state['status'] not in {'queued', 'running'} and state.get('decision_status') != 'pending':
            st.rerun()
        review = json.loads((Path(folder) / 'plan.json').read_text())
        benchmark.validate_review(review)
        active = next((r for r in review['candidates'] if r['id'] == state.get('active_combination')), None)
        title = state['status'].capitalize()
        progress = benchmark.combination_progress(review, state)
        st.markdown('**' + benchmark.combination_caption(review) + '**')
        st.caption(f'{title} · {progress["completed"]} / {progress["total"]} planned evaluations finished')
        if state.get('reuse'):
            st.caption(f"Reused {state['reuse']['hits']} stage results · archived {state['reuse']['stored']}")
        if review.get('excluded'):
            with st.expander(f"Unavailable combinations ({len(review['excluded'])})"):
                st.dataframe([{'Chain': benchmark.chain_label(r), 'Reason': r['reason']} for r in review['excluded']], hide_index=True)
        st.progress(progress['completed'] / max(1, progress['total']))
        st.caption(f"Evaluation set: {len(review['records']):,} measurements + crystal validation")
        if active:
            labels = [('preparation', 'Protein'), ('pocket', active.get('pocket', benchmark.pockets.REFERENCE)['tool']), ('docking', active['docking']['tool']), ('scoring', active['scoring']['tool']), ('validation', 'Evaluate')]
            nodes = ' <span class="chain-arrow">→</span> '.join(f'<span class="chain-node {"active" if stage == state.get("stage") and state["status"] == "running" else ""}">{html.escape(label)}</span>' for stage, label in labels)
            st.markdown('<div class="evaluation-chain">' + nodes + '</div>', unsafe_allow_html=True)
            st.caption(state.get('active_compound', ''))
        st.dataframe([{'Combination': benchmark.chain_label(r), 'Status': progress['statuses'][r['id']]} for r in review['candidates']], hide_index=True, width='stretch')
        if state.get('error'):
            st.error(state['error'])
        if state['rows']:
            table = [{'Compound': r['compound_id'], 'Pocket': r.get('pocket', 'reference_box'), 'Assay group': r.get('cohort'), 'Docking': r['docking'], 'Scoring': r['scoring'], 'Score': r.get('score'), 'Units': r.get('score_unit'),
                      'Experimental endpoint': 'p' + r['endpoint'] if r.get('endpoint') else 'Crystal reference', 'Experimental pX': r.get('p_value'),
                      'Crystal-pose score': r.get('reference_score'), 'RMSD (Å)': r.get('rmsd_angstrom'), 'Physical pass': r.get('physical_pass'), 'Seconds': r.get('runtime_seconds'),
                      'Status': r['status'], 'Failure': r.get('error')} for r in state['rows']]
            st.dataframe(table, hide_index=True, width='stretch')
            from .pipeline_decisions import compare
            comparisons = {benchmark.chain_label(r['pair']): r for r in compare(review, state)}
            metrics = [{'Combination': label, 'Scored': r['completed'], 'Skipped': r['skipped'],
                        'Shared set': r['shared_measurements'],
                        'Correlation measurements': f"{r['ranked_measurements']}/{r['expected']}",
                        'Spearman ρ': r['rho'], 'Kendall τ': r['tau']}
                       for label, r in comparisons.items()]
            st.dataframe(metrics, hide_index=True, width='stretch')
            st.caption('Exploratory ranking on shared successful measurements within each assay. Failed ligands are skipped, never assigned a score.')
            available = [r for r in state['rows'] if r.get('pose_path')]
            if available:
                index = st.selectbox('Inspect docked complex', range(len(available)), key='studio_eval_inspect',
                                     format_func=lambda i: available[i]['compound_id'] + ' · ' + available[i].get('pocket', 'reference_box') + ' → ' + available[i]['docking'] + ' → ' + available[i]['scoring'])
                row = available[index]
                st.caption(f"Docked score: {row['score']:.3f} {row['score_unit']} · Experimental: " + (f"p{row['endpoint']} {row['p_value']:.3f}" if row.get('p_value') is not None else 'no matched activity measurement'))
                if row.get('reference_score') is not None:
                    st.caption(f"Crystal coordinates scored by the same tool: {row['reference_score']:.3f} {row['score_unit']} · same ligand, receptor frame")
                import streamlit.components.v1 as components
                try:
                    components.html(_scene(folder, review, row), height=525, scrolling=False)
                except (OSError, ValueError) as exc:
                    st.warning(str(exc))
            st.download_button('Download results', json.dumps(state, indent=2), 'evaluation-results.json', 'application/json', key='studio_eval_download')
        with st.expander('Live tool log'):
            for event in state.get('events', [])[-10:]:
                st.caption(f"{event.get('tool')}: {event.get('message')}")
    live()
