"""Landing page and independent experiment navigation."""
from pathlib import Path
import html

from . import experiment_library as library


def save_current(st, root):
    current = st.session_state.get('experiments_current')
    if not current:
        return
    if library.is_deleted(root, current):
        library.clear_workspace(st.session_state)
        st.session_state.pop('experiments_current', None)
        st.query_params.pop('experiment', None)
        st.query_params.pop('evaluation', None)
        return
    current = library.save(root, current, st.session_state)
    st.session_state['experiments_current'] = current
    st.query_params['experiment'] = current
    if st.session_state.get('studio_eval_run'):
        st.query_params['evaluation'] = Path(st.session_state['studio_eval_run']).name
    else:
        st.query_params.pop('evaluation', None)


def go_home(st, root):
    save_current(st, root)
    library.clear_workspace(st.session_state)
    st.session_state.pop('experiments_current', None)
    st.session_state['experiments_epoch'] = st.session_state.get('experiments_epoch', 0) + 1
    st.query_params.pop('experiment', None)
    st.query_params.pop('evaluation', None)


def open_experiment(st, root, experiment_id):
    if st.session_state.get('experiments_current'):
        save_current(st, root)
    saved = library.restore(root, experiment_id, st.session_state)
    st.query_params['experiment'] = experiment_id
    if saved.get('run_id'):
        st.query_params['evaluation'] = saved['run_id']
    else:
        st.query_params.pop('evaluation', None)


def new_experiment(st, root):
    if st.session_state.get('experiments_current'):
        save_current(st, root)
    open_experiment(st, root, library.create(root))


def route(st, root):
    """True renders a workspace; false means the library handled this page."""
    requested = st.query_params.get('experiment')
    legacy_run = st.query_params.get('evaluation')
    if not requested and library.valid_id(legacy_run):
        library.adopt_runs(root)
        matches = [r for r in library.summaries(root) if library.load(root, r['id']).get('run_id') == legacy_run]
        if matches:
            requested = matches[0]['id']
            st.query_params['experiment'] = requested
    if requested:
        try:
            library.load(root, requested)
            if st.session_state.get('experiments_current') != requested:
                library.restore(root, requested, st.session_state)
            return True
        except (OSError, ValueError, KeyError) as exc:
            st.error(f'Could not open experiment: {exc}')
            if st.button('All experiments'):
                st.session_state.pop('experiments_current', None)
                st.query_params.clear()
                st.rerun()
            return False
    if st.session_state.get('experiments_current'):
        # Direct browser navigation to the root also preserves the outgoing draft.
        current = st.session_state.pop('experiments_current')
        if not library.is_deleted(root, current):
            library.save(root, current, st.session_state)
        library.clear_workspace(st.session_state)
    render_home(st, root)
    return False


STATUS = {'running': '◌ Running', 'queued': '◷ Queued', 'completed': '✓ Completed', 'partial': '◐ Partial',
          'cancelled': '■ Stopped', 'interrupted': '■ Interrupted', 'failed': '⚠ Failed', 'draft': '○ Draft', 'unavailable': '⚠ Unavailable'}


def render_home(st, root):
    library.adopt_runs(root)
    with st.container(key='experiment_home'):
        st.markdown('<div class="library-brand"><span>✣</span> Bio Orchestrator</div><h1>Your experiments</h1>', unsafe_allow_html=True)
        st.button('Create new experiment', type='primary', key='experiments_new', on_click=new_experiment, args=(st, root))
        if st.session_state.get('experiments_delete_error'):
            st.error(st.session_state.pop('experiments_delete_error'))
        if st.session_state.get('experiments_deleted'):
            if st.button('Undo delete', key='experiments_undo_delete'):
                deleted = st.session_state['experiments_deleted']
                for experiment_id in deleted if isinstance(deleted, list) else [deleted]:
                    library.undo_delete(root, experiment_id)
                st.session_state.pop('experiments_deleted', None)
                st.rerun()
        @st.fragment(run_every=3)
        def table():
            if st.session_state.pop('experiments_clear_selection', False):
                st.session_state['experiments_bulk_selection'] = []
            rows = library.summaries(root)
            if not rows:
                st.caption('Your experiments will appear here.')
                return
            running = sum(r['status'] in {'queued', 'running'} for r in rows)
            if running:
                st.markdown(f'<div class="library-running"><i></i>{running} running in the background</div>', unsafe_allow_html=True)
            names = {r['id']: f"{r['target']} · {r['id'][:6]}" for r in rows}
            selected_key = 'experiments_bulk_selection'
            # Store IDs, never row positions: background refreshes/new runs must
            # not redirect a pending deletion to another experiment.
            old_selection = st.session_state.get(selected_key, [])
            if any(key not in names for key in old_selection):
                st.session_state[selected_key] = [key for key in old_selection if key in names]
            choices, all_button, delete_button = st.columns([6, 1.3, 2], vertical_alignment='bottom')
            with choices:
                selected = st.multiselect('Select experiments', list(names), format_func=names.get,
                                         key=selected_key, placeholder='Select experiments to delete', label_visibility='collapsed')
            def select_all():
                st.session_state[selected_key] = [] if len(st.session_state.get(selected_key, [])) == len(names) else list(names)
            all_button.button('Clear selection' if len(selected) == len(names) else 'Select all',
                              key='experiments_select_all', on_click=select_all, width='stretch')
            if delete_button.button(f'Delete selected ({len(selected)})', key='experiments_delete_selected',
                                    disabled=not selected, width='stretch'):
                st.session_state['experiments_delete_pending'] = list(selected)
                st.rerun()
            body = []
            esc = lambda value: html.escape(str(value), quote=True)
            for row in rows:
                name = f"{row['target']} · {row['id'][:6]}"
                cells = [f'<a href="?experiment={row["id"]}" target="_self">{esc(name)}</a>',
                         f'<span class="experiment-status {esc(row["status"])}">{esc(STATUS.get(row["status"], row["status"]))}</span>',
                         esc(row['progress']), str(row['measurements']) if row['measurements'] is not None else '—',
                         esc(row['result']), esc(row['chain']),
                         (f"{row['rho']:.3f}" + (f" · n={row['ranked_measurements']}/{row['measurements']}" if 'ranked_measurements' in row else '')
                          + (' · small sample' if row.get('small_sample') else '')) if row['rho'] is not None else '—',
                         esc(row['created'][:16].replace('T', ' ')),
                         f'<a class="experiment-delete" href="?delete_experiment={row["id"]}" target="_self" aria-label="Delete {esc(name)}" title="Delete experiment">×</a>']
                cells.insert(1, cells.pop())
                body.append('<tr>' + ''.join('<td>' + value + '</td>' for value in cells) + '</tr>')
            headings = ['Experiment', '', 'Status', 'Combinations', 'Measurements', 'Result', 'Final / best chain', 'Spearman ρ', 'Created (UTC)']
            st.markdown('<div class="experiment-table-wrap"><table class="experiment-table"><thead><tr>' +
                        ''.join('<th>' + esc(value) + '</th>' for value in headings) + '</tr></thead><tbody>' +
                        ''.join(body) + '</tbody></table></div>', unsafe_allow_html=True)
        table()
    requested_delete = st.query_params.get('delete_experiment')
    pending = [requested_delete] if requested_delete else st.session_state.get('experiments_delete_pending', [])
    if pending:
        def clear_pending():
            st.query_params.pop('delete_experiment', None)
            st.session_state.pop('experiments_delete_pending', None)
        @st.dialog('Delete experiment?' if len(pending) == 1 else f'Delete {len(pending)} experiments?')
        def confirm_delete():
            try:
                saved = [library.load(root, key) for key in pending]
            except (OSError, ValueError) as exc:
                st.error(str(exc))
                if st.button('Close'):
                    clear_pending()
                    st.rerun()
                return
            for item in saved:
                st.write((item['state'].get('studio_draft', {}).get('target') or 'New experiment') + ' · ' + item['id'][:6])
            st.caption('Removes the selected experiments and stops any that are running. Undo restores the entries; stopped runs remain stopped.')
            cancel, remove = st.columns(2)
            if cancel.button('Keep experiment' if len(pending) == 1 else 'Keep experiments', width='stretch'):
                clear_pending()
                st.rerun()
            if remove.button('Delete experiment' if len(pending) == 1 else f'Delete {len(pending)} experiments', type='primary', width='stretch'):
                try:
                    deleted, errors = library.delete_many(root, pending)
                    if deleted:
                        st.session_state['experiments_deleted'] = deleted
                        st.session_state['experiments_clear_selection'] = True
                    if errors:
                        st.session_state['experiments_delete_error'] = 'Some experiments could not be deleted: ' + '; '.join(f'{key[:6]}: {value}' for key, value in errors.items())
                    clear_pending()
                    st.rerun()
                except (OSError, ValueError) as exc:
                    st.error(str(exc))
        confirm_delete()
