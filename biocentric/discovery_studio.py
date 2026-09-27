"""Discovery candidate selection and public evidence workspace; no job execution."""
from __future__ import annotations

from copy import deepcopy
import hashlib
from itertools import product
import json
from pathlib import Path
import re

PHASES = {"structure_preparation": "3D Protein Structure Preparation", "pocket_identification": "Pocket Identification",
          "docking": "Docking", "scoring": "Scoring"}
UI_ONLY_PHASES = {"structure_preparation"}
DOCK_DEFAULTS = {"exhaustiveness": 8, "num_modes": 3, "seed": 42}
TOOLS = {
    "structure_preparation": {
        "pdbfixer": ("PDBFixer", {}),
        "pdb2pqr": ("PDB2PQR", {}),
        "reduce": ("Reduce", {}),
        "colabfold": ("ColabFold · AI folding", {}),
        "esmfold": ("ESMFold · AI folding", {}),
    },
    "pocket_identification": {
        "reference_box": ("Crystal ligand pocket", {"padding": 4, "min_size": 24}),
        "fpocket": ("fpocket", {"pocket_rank": 1, "padding": 4, "min_size": 20}),
        "p2rank": ("P2Rank · ML", {}),
        "deeppocket": ("DeepPocket · AI", {}),
    },
    "docking": {
        "vina": ("AutoDock Vina", {**DOCK_DEFAULTS, "cpu": 4, "max_evals": 0}),
        "gnina": ("GNINA", {**DOCK_DEFAULTS, "cpu": 4, "device": "gpu", "max_mc_steps": 0}),
        "unidock": ("Uni-Dock", {"search_mode": "balance", "num_modes": 3, "seed": 42}),
        "autodock4": ("AutoDock4", {"ga_runs": 10, "ga_num_evals": 2500000, "ga_pop_size": 150}),
        "smina": ("smina", DOCK_DEFAULTS),
        "qvina2": ("QuickVina 2", DOCK_DEFAULTS),
    },
    "scoring": {
        "vina_score": ("Vina score", {"cpu": 1}),
        "ad4_score": ("AutoDock4 score", {"cpu": 1}),
        "gnina_cnn": ("GNINA CNN pose", {"cnn_scoring": "rescore", "device": "gpu", "cpu": 4}),
        "gnina_affinity": ("GNINA affinity", {"device": "gpu", "cpu": 4}),
        "biocentric_affinity": ("BioCentric affinity", {"device": "cuda", "cpu": 2, "batch_size": 16, "timeout_seconds": 300}),
        "planet": ("PLANET", {"device": "gpu"}),
        "mmpbsa": ("MM/PBSA", {"start_frame": 1, "end_frame": 50, "interval": 1, "salt_molar": .15}),
        "mmgbsa": ("MM/GBSA", {"start_frame": 1, "end_frame": 50, "interval": 1, "salt_molar": .15, "igb": 5}),
    },
}
DEFAULT_SELECTION = {"structure_preparation": ["pdbfixer"], "pocket_identification": ["reference_box"], "docking": ["vina", "gnina", "unidock"], "scoring": ["vina_score", "gnina_cnn", "biocentric_affinity"]}
DOCS = {"pdbfixer": "https://github.com/openmm/pdbfixer", "pdb2pqr": "https://github.com/Electrostatics/pdb2pqr",
        "reduce": "https://github.com/rlabduke/reduce", "colabfold": "https://github.com/sokrypton/ColabFold",
        "esmfold": "https://github.com/facebookresearch/esm", "fpocket": "https://github.com/Discngine/fpocket",
        "p2rank": "https://github.com/rdk/p2rank", "deeppocket": "https://github.com/devalab/DeepPocket",
        "vina": "https://autodock-vina.readthedocs.io/", "vina_score": "https://autodock-vina.readthedocs.io/",
        "ad4_score": "https://autodock-vina.readthedocs.io/", "autodock4": "https://autodock.scripps.edu/",
        "gnina": "https://github.com/gnina/gnina", "gnina_cnn": "https://github.com/gnina/gnina", "gnina_affinity": "https://github.com/gnina/gnina",
        "unidock": "https://github.com/dptech-corp/Uni-Dock", "planet": "https://github.com/ComputArtCMCG/PLANET",
        "smina": "https://sourceforge.net/projects/smina/", "qvina2": "https://qvina.github.io/",
        "mmpbsa": "https://ambermd.org/tutorials/advanced/tutorial3/py_script/section1.htm", "mmgbsa": "https://ambermd.org/tutorials/advanced/tutorial3/py_script/section1.htm"}
LABELS = {"max_evals": "Search evaluation limit (0 = automatic)", "max_mc_steps": "Search step limit (0 = automatic)", "pocket_rank": "Pocket rank", "padding": "Box padding (Å)", "min_size": "Minimum box size (Å)", "exhaustiveness": "Exhaustiveness", "num_modes": "Poses", "seed": "Seed", "search_mode": "Search mode",
          "cpu": "CPU threads", "cnn_scoring": "CNN scoring", "ga_runs": "GA runs", "ga_num_evals": "GA evaluations",
          "ga_pop_size": "Population", "device": "Device", "start_frame": "First frame", "end_frame": "Last frame",
          "interval": "Frame interval", "salt_molar": "Salt (M)", "igb": "GB model",
          "batch_size": "Compounds per model batch", "timeout_seconds": "Batch timeout (seconds)"}
OPTIONS = {"search_mode": ["fast", "balance", "detail"], "cnn_scoring": ["rescore", "refinement"], "device": ["cpu", "gpu"], "igb": [2, 5, 7, 8]}
BOUNDS = {"max_evals": (0, 10000000), "max_mc_steps": (0, 100000), "pocket_rank": (1, 20), "padding": (1, 12), "min_size": (10, 40), "exhaustiveness": (1, 128), "num_modes": (1, 20), "seed": (1, 2147483647), "cpu": (1, 16),
          "ga_runs": (1, 100), "ga_num_evals": (1000, 25000000), "ga_pop_size": (20, 300),
          "start_frame": (1, 10000), "end_frame": (2, 10000), "interval": (1, 1000),
          "batch_size": (1, 32), "timeout_seconds": (30, 3600)}


def catalog(state, phase):
    result = dict(TOOLS[phase])
    for tool_id, plugin in state.get("studio_custom_plugins", {}).items():
        if plugin["phase"] == phase:
            result[tool_id] = (plugin["name"], plugin["parameters"])
    return result


DOCKING_PRESETS = {
    'Fast': {'vina': {'exhaustiveness': 2, 'num_modes': 1, 'cpu': 4, 'max_evals': 5000},
             'gnina': {'exhaustiveness': 2, 'num_modes': 1, 'cpu': 4, 'max_mc_steps': 20},
             'unidock': {'search_mode': 'fast', 'num_modes': 1}},
    'Balanced': {'vina': {'exhaustiveness': 8, 'num_modes': 9, 'cpu': 4, 'max_evals': 0},
                 'gnina': {'exhaustiveness': 8, 'num_modes': 3, 'cpu': 4, 'max_mc_steps': 0},
                 'unidock': {'search_mode': 'balance', 'num_modes': 3}},
}


def selected_parameters(state, phase, tool):
    values = {**catalog(state, phase)[tool][1], **state.get('studio_parameters', {}).get(f'{phase}:{tool}', {})}
    if phase == 'docking':
        values.update(DOCKING_PRESETS.get(state.get('studio_docking_preset', 'Balanced'), {}).get(tool, {}))
    return deepcopy(values)


def selection(state):
    return {phase: [{'tool': tool, 'parameters': selected_parameters(state, phase, tool)}
                   for tool in state.get(f'studio_tools_{phase}', DEFAULT_SELECTION[phase])] for phase in PHASES}


def _apply_docking_preset(state):
    for tool in catalog(state, 'docking'):
        state.setdefault('studio_parameters', {})['docking:' + tool] = selected_parameters(state, 'docking', tool)
    for key in list(state):
        if key.startswith('studio_parameter_docking_'):
            del state[key]


def selection_hash(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()


def selection_issues(state):
    issues = []
    for phase, rows in selection(state).items():
        if not rows and phase not in UI_ONLY_PHASES:
            issues.append(f"Select at least one {phase} tool.")
        for row in rows:
            if f"{phase}:{row['tool']}" in state.get('studio_parameter_errors', {}):
                issues.append('Fix the selected tool’s parameters.')
            p = row['parameters']
            if row['tool'] in {'mmpbsa', 'mmgbsa'} and not 2 <= (p['end_frame'] - p['start_frame']) // p['interval'] + 1 <= 1000:
                issues.append('MM/PBSA and MM/GBSA require 2–1000 trajectory frames.')
    return list(dict.fromkeys(issues))


def combinations(state):
    chosen = selection(state)
    result = []
    for pocket, docking, scoring in product(chosen['pocket_identification'], chosen['docking'], chosen['scoring']):
        score = scoring['tool']
        requirement = ('Reviewed MD trajectory and force fields' if score in {'mmpbsa', 'mmgbsa'} else
                       'Pocket/ligand predictor; docking-independent' if score == 'planet' else
                       'Pose recovery benchmark; not an affinity metric' if score == 'gnina_cnn' else
                       'Matched experimental endpoint and held-out test set')
        result.append({'Pocket': catalog(state, 'pocket_identification')[pocket['tool']][0], 'Docking': catalog(state, 'docking')[docking['tool']][0],
                       'Scoring': catalog(state, 'scoring')[score][0], 'Requires': requirement,
                       'Status': 'Not validated'})
    return result


def _save_parameter(state, storage, parameter, widget):
    state['studio_parameters'][storage][parameter] = state[widget]
    if storage.startswith('docking:'):
        state['studio_docking_preset'] = 'Custom'


def _save_json(state, storage, widget):
    try:
        value = json.loads(state[widget])
        if not isinstance(value, dict):
            raise ValueError('Use a JSON object.')
        state['studio_parameters'][storage] = value
        state.setdefault('studio_parameter_errors', {}).pop(storage, None)
    except (ValueError, TypeError):
        state.setdefault('studio_parameter_errors', {})[storage] = 'Enter valid JSON parameters.'


def _add_plugin(state):
    name, phase = state.get('studio_plugin_name', '').strip(), state['studio_plugin_phase']
    if not name:
        state['studio_plugin_error'] = 'Enter a plugin name.'
        return
    if name.casefold() in {label.casefold() for label, _ in catalog(state, phase).values()}:
        state['studio_plugin_error'] = 'This tool already exists in that phase.'
        return
    tool_id = 'custom_' + hashlib.sha256(f'{phase}:{name}'.encode()).hexdigest()[:12]
    state.setdefault('studio_custom_plugins', {})[tool_id] = {'name': name, 'phase': phase, 'parameters': {}}
    state[f'studio_tools_{phase}'] = [*state.get(f'studio_tools_{phase}', DEFAULT_SELECTION[phase]), tool_id]
    state.pop('studio_plugin_error', None)
    state['studio_plugin_added'] = name


def _parameters(st, phase, tool, defaults):
    storage = f'{phase}:{tool}'
    values = selected_parameters(st.session_state, phase, tool)
    st.session_state.setdefault('studio_parameters', {})[storage] = values
    for parameter, default in defaults.items():
        values.setdefault(parameter, deepcopy(default))
    if tool.startswith('custom_'):
        key = f'studio_json_{phase}_{tool}'
        st.text_area('Parameters (JSON)', value=json.dumps(values, indent=2), key=key,
                     on_change=_save_json, args=(st.session_state, storage, key))
        if st.session_state.get('studio_parameter_errors', {}).get(storage):
            st.error(st.session_state['studio_parameter_errors'][storage])
        return
    if phase in UI_ONLY_PHASES:
        st.caption(f'UI only · not connected. [Tool documentation]({DOCS[tool]})')
        return
    if phase == 'pocket_identification':
        from .pocket_identification import blocker
        reason = blocker({'tool': tool, 'parameters': values})
        if reason:
            st.caption(reason)
    for parameter, default in defaults.items():
        key = f'studio_parameter_{phase}_{tool}_{parameter}'
        kwargs = dict(key=key, on_change=_save_parameter, args=(st.session_state, storage, parameter, key))
        if parameter in OPTIONS:
            options = ['cpu', 'cuda'] if tool == 'biocentric_affinity' and parameter == 'device' else OPTIONS[parameter]
            st.selectbox(LABELS[parameter], options, index=options.index(values[parameter]), **kwargs)
        elif isinstance(default, int):
            st.number_input(LABELS[parameter], min_value=BOUNDS[parameter][0], max_value=BOUNDS[parameter][1], value=values[parameter], **kwargs)
        else:
            st.number_input(LABELS[parameter], min_value=0.0, max_value=1.0, step=.05, value=values[parameter], **kwargs)
    if tool == 'reference_box':
        st.caption('Experimental ligand coordinates define the pocket center and bounds.')
    elif tool == 'biocentric_affinity':
        st.caption('Local ensemble · predicted pK, higher is better · availability checked before execution')
    else:
        st.caption(f'[Tool documentation]({DOCS[tool]}) · Availability checked before execution')


def _header(st):
    with st.container(key='studio_header'):
        st.markdown('<div class="studio-brand"><span class="brand-symbol">✣</span> Bio <span class="brand-light">Orchestrator</span></div>', unsafe_allow_html=True)


def _tools(st):
    if 'studio_docking_preset' not in st.session_state:
        st.session_state['studio_docking_preset'] = 'Balanced'
        _apply_docking_preset(st.session_state)
    with st.container(key='studio_toolbar'):
        for index, phase in enumerate(PHASES, start=1):
            if index > 1:
                st.markdown('<div class="phase-connector" aria-hidden="true">'
                            '<svg width="16" height="22" viewBox="0 0 16 22" fill="none">'
                            '<path d="M8 0v17m-4-4 4 4 4-4" stroke="currentColor" '
                            'stroke-width="1.6" stroke-linecap="round" stroke-linejoin="round"/>'
                            '</svg></div>', unsafe_allow_html=True)
            with st.container(key='studio_phase_' + phase):
                choices = catalog(st.session_state, phase)
                widget_key = f'studio_tools_{phase}'
                chosen = st.session_state.get(widget_key, DEFAULT_SELECTION[phase])
                heading, settings = st.columns([8, 1], gap='small', vertical_alignment='center')
                with heading:
                    preset = st.session_state.get('studio_docking_preset', 'Balanced')
                    badge = (f'<span class="phase-preset">{preset}</span>'
                             if phase == 'docking' and preset in ('Fast', 'Balanced', 'Custom') else '')
                    st.markdown(f'<div class="phase-title"><span class="phase-number">{index:02}</span>'
                                f'<span>{PHASES[phase]}</span>{badge}</div>', unsafe_allow_html=True)
                with settings:
                    with st.popover(f'{PHASES[phase]} settings', icon=':material/tune:', width='content'):
                        if phase == 'docking':
                            st.selectbox('Search preset', ['Fast', 'Balanced', 'Custom'], key='studio_docking_preset',
                                         on_change=_apply_docking_preset, args=(st.session_state,))
                            st.caption('Fast: shorter search, one pose. Lower pose accuracy is possible.')
                            if st.session_state.get('studio_eval_run'):
                                st.caption('Applies to the next Start. Existing run settings stay fixed.')
                        if not chosen:
                            st.caption('Select tools first.')
                        for tool in chosen:
                            with st.expander(choices[tool][0]):
                                _parameters(st, phase, tool, choices[tool][1])
                defaults = {} if widget_key in st.session_state else {'default': DEFAULT_SELECTION[phase]}
                st.multiselect(PHASES[phase], list(choices), **defaults,
                               format_func=lambda item, c=choices: c[item][0], key=widget_key,
                               placeholder='Select tools', label_visibility='collapsed')
        with st.popover('Add plugin', icon=':material/add:', width='stretch'):
            with st.form('studio_plugin_form', border=False):
                st.text_input('Plugin name', placeholder='e.g. Schrödinger Glide', key='studio_plugin_name', max_chars=80)
                st.selectbox('Phase', list(PHASES), format_func=PHASES.get, key='studio_plugin_phase')
                st.form_submit_button('Add', on_click=_add_plugin, args=(st.session_state,))
            if st.session_state.get('studio_plugin_error'):
                st.error(st.session_state['studio_plugin_error'])
            elif st.session_state.get('studio_plugin_added'):
                st.caption(f"{st.session_state['studio_plugin_added']} added.")
            st.caption('Configuration only. Plugin execution comes later.')


def submit_message(state, message):
    message = message.strip()
    if not message:
        return
    if state.get('studio_eval_run'):
        from .discovery_benchmark import read_status
        if read_status(state['studio_eval_run'])['status'] in {'queued', 'running'}:
            state['studio_draft']['messages'].append({'role': 'assistant', 'text': 'Stop the evaluation before changing the target.'})
            return
    draft = state.setdefault('studio_draft', {'messages': [], 'status': 'empty'})
    draft['messages'].append({'role': 'user', 'text': message})
    draft['request'] = message
    if re.fullmatch(r'[A-Za-z0-9_.+\- ]{1,35}', message):
        draft['target'] = message
        state['studio_target'] = message
        _submit_scope(state)
    else:
        draft['messages'].append({'role': 'assistant', 'text': 'Enter a protein name or UniProt ID.'})
        draft['status'] = 'scope'
    draft.pop('accepted_selection', None)


def _submit_chat(state):
    previous_run = state.get('studio_eval_run')
    if previous_run:
        from .discovery_benchmark import read_status
        if read_status(previous_run)['status'] in {'queued', 'running'}:
            return
    target = (state.get('studio_chat_input') or '').strip()
    draft = state.setdefault('studio_draft', {'messages': [], 'status': 'empty'})
    if not target:
        draft['error'] = 'Enter a protein name or UniProt ID.'
        return
    draft['messages'].append({'role': 'user', 'text': target})
    draft['request'] = target
    state['studio_target'] = target
    _submit_scope(state)
    if previous_run and not state.get('studio_eval_run'):
        import streamlit as st
        st.query_params.pop('evaluation', None)


def _submit_scope(state):
    draft = state['studio_draft']
    target, organism = state['studio_target'].strip(), 'Human'
    if not target:
        draft['error'] = 'Enter a target.'
        return
    if (target, organism) != draft.get('evidence_query'):
        previous_run = state.get('studio_eval_run')
        for key in ('studio_identity', 'studio_resolution', 'studio_evidence', 'studio_structure_page', 'studio_binding_page', 'studio_pdb', 'studio_gallery_batch', 'studio_gallery_loaded'):
            state.pop(key, None)
        for key in list(state):
            if key.startswith('studio_eval_'):
                state.pop(key, None)
        if previous_run:
            state['studio_eval_detached'] = Path(previous_run).name
    draft.update(target=target, organism=organism, status='review', evidence_query=(target, organism), selection_sha256=selection_hash(selection(state)))
    draft.pop('error', None)


def _accept(state):
    if not selection_issues(state):
        state['studio_draft'].update(status='accepted', accepted_selection=selection(state), selection_sha256=selection_hash(selection(state)))


def _edit(state):
    if state.get('studio_eval_run'):
        from .discovery_benchmark import read_status
        if read_status(state['studio_eval_run'])['status'] in {'queued', 'running'}:
            return
    state['studio_draft']['status'] = 'scope'
    state['studio_draft'].pop('accepted_selection', None)


def render_studio(st):
    from .discovery_panels import render_experiments
    st.markdown('<style>' + Path(__file__).with_name('discovery_studio.css').read_text() + '</style>', unsafe_allow_html=True)
    from .discovery_panels import _root
    from .experiment_library_ui import route, go_home, save_current
    if not route(st, _root()):
        return
    if st.session_state.get('studio_version') != 3:
        # Migrate the session UI, retaining custom plugins for these two phases.
        for key in ('studio_draft', 'studio_target', 'studio_organism', 'studio_parameters', 'studio_tools_docking', 'studio_tools_scoring', 'studio_plugin_phase'):
            st.session_state.pop(key, None)
        st.session_state['studio_version'] = 3
    from .discovery_evaluation_ui import restore_run
    from .discovery_panels import _root
    restore_run(st, _root())
    draft = st.session_state.setdefault('studio_draft', {'messages': [], 'status': 'empty'})
    if draft['status'] in {'accepted', 'review'} and (selection_issues(st.session_state) or draft.get('selection_sha256') != selection_hash(selection(st.session_state))):
        draft['status'] = 'review'
        draft.pop('accepted_selection', None)
    with st.container(key='studio_split'):
        controls, canvas = st.columns([1, 2], gap='medium')
        with controls:
            with st.container(key='studio_controls', height=700, border=False):
                _header(st)
                st.button('← All experiments', key='experiments_back', on_click=go_home, args=(st, _root()))
                _tools(st)
                with st.container(key='studio_protein'):
                    if 'studio_chat_input' not in st.session_state:
                        st.session_state['studio_chat_input'] = draft.get('target', '')
                    run = st.session_state.get('studio_eval_run')
                    from .discovery_benchmark import read_status
                    active = bool(run and read_status(run)['status'] in {'queued', 'running'})
                    st.text_input('Protein target', key='studio_chat_input', placeholder='Protein name, UniProt or PDB ID',
                                  max_chars=160, disabled=active, on_change=_submit_chat, args=(st.session_state,))
                    if draft.get('error'):
                        st.error(draft['error'])
                    if draft['status'] in {'review', 'accepted'}:
                        st.caption(f'{len(combinations(st.session_state))} candidate combinations')
                        for issue in selection_issues(st.session_state):
                            st.warning(issue)
                        from .discovery_evaluation_ui import chat_review
                        chat_review(st, _root())
                    else:
                        st.button('Start', key='studio_empty_start', type='primary', disabled=True, width='stretch')
        with canvas:
            with st.container(key='studio_results', height=700, border=False):
                render_experiments(st, draft)
    save_current(st, _root())
