"""Small, source-linked evidence panels for the discovery screen."""
from __future__ import annotations

import csv
from concurrent.futures import ThreadPoolExecutor
import html
import io
import json
import os
from pathlib import Path
from urllib.parse import quote

import httpx

from . import discovery_evidence as evidence
from .retrieval import RetrievalError


def _root():
    return str(Path(os.getenv('BIOCENTRIC_HOME', '.biocentric')).resolve())


def _identity(st, draft):
    if 'studio_resolution' not in st.session_state:
        with st.spinner('Resolving target · UniProt…'):
            st.session_state['studio_resolution'] = evidence.resolve_identity(draft['target'], draft['organism'], _root())
    resolved = st.session_state['studio_resolution']
    candidates = resolved.get('candidates', [])
    if not candidates:
        st.info('No exact target match. Try a UniProt accession or refine the organism.')
        with st.expander('Source details'):
            st.write(resolved.get('limitations', []))
        if st.button('Retry target lookup'):
            st.session_state.pop('studio_resolution', None)
            st.rerun()
        return None
    if len(candidates) == 1 and resolved.get('status') == 'resolved':
        st.session_state['studio_identity'] = candidates[0]
    elif 'studio_identity' not in st.session_state:
        mapped = {r['accession']: r for r in candidates if r.get('accession')}
        if not mapped:
            st.info('Provide a UniProt accession to explore public structures.')
            return None
        accession = st.selectbox('Which target?', list(mapped), format_func=lambda k: f"{k} · {mapped[k]['name']}")
        if st.button('Use this target', type='primary'):
            st.session_state['studio_identity'] = mapped[accession]
            st.rerun()
        return None
    return st.session_state.get('studio_identity')


def _paper_url(citation):
    if citation.get('pdbx_database_id_PubMed'):
        return f"https://pubmed.ncbi.nlm.nih.gov/{quote(str(citation['pdbx_database_id_PubMed']), safe='')}/"
    if citation.get('pdbx_database_id_DOI'):
        return 'https://doi.org/' + quote(citation['pdbx_database_id_DOI'], safe='/')
    return None


def _gallery_order(rows):
    return sorted(rows, key=lambda r: (not (any('X-RAY' in m for m in r['methods']) and r['compounds']),
                                       not bool(r['compounds']), r['id']))


def _gallery_step(st, identity, direction):
    data = st.session_state['studio_evidence']
    current = data['structures']
    batch = st.session_state.get('studio_gallery_batch', 0)
    source_page = st.session_state.get('studio_structure_page', 0)
    last_batch = max(0, (len(current['rows']) - 1) // 6)
    next_batch = batch + direction
    if 0 <= next_batch <= last_batch:
        st.session_state['studio_gallery_batch'] = next_batch
        return
    new_page = source_page + direction
    if not 0 <= new_page < max(1, (current.get('total', 0) + evidence.PAGE_SIZE - 1) // evidence.PAGE_SIZE):
        return
    try:
        result = evidence.structure_page(identity, _root(), new_page)
        data['structures'] = result
        st.session_state['studio_structure_page'] = new_page
        st.session_state['studio_gallery_batch'] = 0 if direction > 0 else max(0, (len(result['rows']) - 1) // 6)
        st.session_state.pop('studio_gallery_error', None)
    except RetrievalError as exc:
        st.session_state['studio_gallery_error'] = str(exc)


def _structures(st, identity, data):
    from .structure_gallery import gallery_html
    current = data['structures']
    if current.get('error'):
        st.warning('Structure source unavailable. Other evidence may still be shown.')
        return
    rows = _gallery_order(current['rows'])
    if identity.get('requested_pdb'):
        rows.sort(key=lambda row: row['id'] != identity['requested_pdb'])
    if not rows:
        st.info('No experimental structures found for this accession.')
        return
    batch = min(st.session_state.get('studio_gallery_batch', 0), (len(rows) - 1) // 6)
    visible = rows[batch * 6:(batch + 1) * 6]
    cache_key = tuple((r['id'], r['models']) for r in visible)
    cached = st.session_state.get('studio_gallery_loaded')
    if not cached or cached['key'] != cache_key:
        items, errors = [], []
        with st.spinner('Loading structures…'):
            with ThreadPoolExecutor(max_workers=3) as pool:
                futures = [(r, pool.submit(evidence.coordinates, r['id'], _root(), 'pdb' if r['models'] > 1 else 'cif')) for r in visible]
                for record, future in futures:
                    try:
                        items.append((record, future.result()))
                    except (RetrievalError, httpx.HTTPError, OSError, ValueError) as exc:
                        errors.append(f"{record['id']}: {exc}")
        cached = {'key': cache_key, 'items': items, 'errors': errors}
        st.session_state['studio_gallery_loaded'] = cached
    with st.container(key='studio_structure_view'):
        if cached['items']:
            import streamlit.components.v1 as components
            gene = (identity.get('gene_names') or [identity['accession']])[0]
            components.html(gallery_html(cached['items'], gene), height=610, scrolling=False)
        else:
            st.warning('Experimental coordinates could not be loaded.')
        with st.container(key='studio_gallery_navigation'):
            previous, spacer, following = st.columns([1, 12, 1])
            source_page = st.session_state.get('studio_structure_page', 0)
            last_source_page = max(0, (current.get('total', 0) - 1) // evidence.PAGE_SIZE)
            previous.button('←', help='Previous structures', key='studio_gallery_previous',
                            disabled=batch == 0 and source_page == 0, on_click=_gallery_step, args=(st, identity, -1))
            following.button('→', help='Next structures', key='studio_gallery_next',
                             disabled=batch * 6 + 6 >= len(rows) and source_page >= last_source_page,
                             on_click=_gallery_step, args=(st, identity, 1))
    if cached['errors'] or st.session_state.get('studio_gallery_error'):
        with st.expander('Unavailable previews'):
            st.caption('\n'.join(cached['errors'] + [st.session_state.get('studio_gallery_error', '')]))


def _bindings(st, identity, data):
    current = data['binding']
    if current.get('error'):
        st.warning('Binding data source unavailable.')
        st.caption(current['error'])
        return
    page, total = st.session_state.get('studio_binding_page', 0), current.get('total', 0)
    st.caption(f"{len(current['rows'])} measurements on this page · {total:,} source records")
    st.caption('ChEMBL binding assays · Kd, Ki and IC50 kept separate. Measured compounds are not all confirmed binders.')
    requested = int(st.number_input('Binding data page', min_value=1, max_value=max(1, (total + 49) // 50), value=page + 1, step=1, key=f'studio_binding_page_choice_{identity["accession"]}')) - 1
    if requested != page:
        with st.spinner('Loading binding measurements…'):
            try:
                data['binding'] = evidence.binding_page(identity, _root(), requested)
                st.session_state['studio_binding_page'] = requested
                st.rerun()
            except RetrievalError as exc:
                st.warning(str(exc))
    rows = current['rows']
    if not rows:
        st.info('No matching binding measurements found.')
        return
    table = [{'Compound': r['compound'], 'Endpoint': r['endpoint'], 'Relation': r['relation'], 'Value': r['value'], 'Units': r['units'],
              'Assay': r['assay'], 'Year': r['year'], 'Quality flag': r['validity'], 'Paper': r['publication_url'], 'Source': r['source_url']} for r in rows]
    st.dataframe(table, hide_index=True, width='stretch', column_config={
        'Paper': st.column_config.LinkColumn('Paper', display_text='Paper'), 'Source': st.column_config.LinkColumn('Source', display_text='ChEMBL')})
    index = st.selectbox('Inspect measurement', range(len(rows)), format_func=lambda i: f"{rows[i]['compound']} · {rows[i]['endpoint']} {rows[i]['relation'] or ''} {rows[i]['value'] or '—'} {rows[i]['units'] or ''}")
    record = rows[index]
    st.write(record['assay_description'])
    st.markdown(f"[{record['publication'] or 'Source publication'}]({record['publication_url']})")
    if record.get('smiles'):
        with st.expander('Ligand identity & source record'):
            st.code(record['smiles'], language=None)
            st.json(record['raw'], expanded=False)
    # Preserve full source fields for later test-set curation. No endpoint pooling.
    stream = io.StringIO()
    fields = [k for k in rows[0] if k != 'raw']
    writer = csv.DictWriter(stream, fieldnames=fields)
    writer.writeheader()
    writer.writerows({k: r[k] for k in fields} for r in rows)
    st.download_button('Export this page', stream.getvalue(), f"{identity['accession']}-binding-page-{page+1}.csv", 'text/csv')


def _literature(st, data):
    # Structure-linked primary citations are stronger than generic title matches.
    seen = set()
    for row in data['structures'].get('rows', []):
        citation = row['citation']
        url = _paper_url(citation)
        if url and url not in seen:
            seen.add(url)
            st.markdown(f"**[{citation.get('title') or row['id']}]({url})**")
            st.caption(f"{citation.get('year') or ''} · Experimental structure {row['id']}")
    retrieved = data['literature']
    if retrieved.get('error'):
        st.caption('Additional literature search unavailable.')
    elif retrieved.get('evidence'):
        with st.expander('Additional literature matches'):
            st.caption('Title/abstract search matches require target and species review.')
            for row in retrieved['evidence']:
                st.markdown(f"[{row['text'].splitlines()[0]}]({row['url']})")
    if not seen and not retrieved.get('evidence'):
        st.info('No publication links retrieved yet.')


def render_experiments(st, draft):
    from .discovery_studio import combinations, selection
    from .pipeline_results_ui import render_map, render_results, render_final_visualizations, render_final_chain
    if not draft.get('request'):
        structures, mapping = st.tabs(['Structures & ligands', 'Pipeline map'])
        with structures:
            st.markdown('<div class="empty-canvas" aria-label="Empty experiments workspace"><div class="canvas-orbit"></div></div>', unsafe_allow_html=True)
        with mapping:
            render_map(st)
        return
    if draft['status'] == 'scope':
        st.markdown(f'<div class="target-heading"><span>TARGET</span><h2>{html.escape(draft.get("target") or "New target")}</h2></div>', unsafe_allow_html=True)
        st.caption('Enter a target in chat.')
        return
    identity = _identity(st, draft)
    if not identity:
        return
    gene = (identity.get('gene_names') or [identity['accession']])[0]
    if 'studio_evidence' not in st.session_state:
        with st.spinner('Finding structures, binding measurements and papers…'):
            st.session_state['studio_evidence'] = evidence.overview(identity, _root())
    data = st.session_state['studio_evidence']
    plan = combinations(st.session_state)
    run_folder = st.session_state.get('studio_eval_run')
    default_tab = 'Final visualizations' if run_folder and (Path(run_folder) / 'decision.json').exists() else 'Pipeline map' if run_folder else 'Evaluation data' if st.session_state.get('studio_eval_requested') else 'Structures & ligands'
    structures, mapping, binding, papers, pipeline, evaluation, running, results, final_chain, final_view = st.tabs(
        ['Structures & ligands', 'Pipeline map', 'Binding data', 'Literature', 'Combinations', 'Evaluation data', 'Evaluation run', 'Final results', 'Final chain', 'Final visualizations'], default=default_tab)
    with mapping:
        render_map(st)
    with results:
        render_results(st)
    with final_chain:
        render_final_chain(st)
    with final_view:
        render_final_visualizations(st)
    with structures:
        _structures(st, identity, data)
    with binding:
        _bindings(st, identity, data)
    with evaluation:
        from .discovery_evaluation_ui import render_data
        render_data(st, identity, data, _root())
    with running:
        if st.session_state.get('studio_eval_run'):
            from .discovery_evaluation_ui import render_live
            render_live(st, st.session_state['studio_eval_run'])
        else:
            st.caption('Click Start in chat to prepare the evaluation data and run the selected tools.')
    with papers:
        _literature(st, data)
    with pipeline:
        st.caption('Selected pocket × docking × scoring tools. Compatibility and literature test-set curation precede benchmarking.')
        st.dataframe(plan, hide_index=True, width='stretch')
        st.info('Predictive validation unavailable: no curated, held-out test set or benchmark run yet.')
        with st.expander('Benchmark requirements'):
            st.write('Match target species, construct, assay endpoint and units. Review censored values and quality flags; separate model training overlap from the held-out test set. Compare pose recovery, ranking accuracy and runtime on compatible protocols.')
        st.download_button('Export candidate plan', json.dumps({'target': identity['accession'], 'organism': identity['organism'], 'tools': selection(st.session_state),
                             'combinations': plan, 'validation': 'unavailable', 'test_set': 'not curated'}, indent=2), f'{gene}-candidate-plan.json', 'application/json')
    with papers:
        st.markdown(f'[UniProt target]({identity["url"]})')
        with st.expander('Sources & coverage'):
            st.caption('Browse all matched structure and binding-data pages. Literature search is a bounded sample; this is not a systematic review of every binder or conformation.')
            st.json({name: {k: value.get(k) for k in ('total', 'page', 'retrieved_at', 'error') if k in value} for name, value in data.items()}, expanded=False)
            if st.button('Retry source queries'):
                st.session_state.pop('studio_evidence', None)
                st.session_state.pop('studio_structure_page', None)
                st.session_state.pop('studio_gallery_batch', None)
                st.session_state.pop('studio_gallery_loaded', None)
                st.session_state.pop('studio_binding_page', None)
                for key in list(st.session_state):
                    if key.startswith(('studio_structure_page_choice_', 'studio_binding_page_choice_')):
                        del st.session_state[key]
                st.rerun()
