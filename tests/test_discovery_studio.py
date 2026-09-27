"""The discovery intake is an interactive brief, without research side effects."""
from __future__ import annotations

from copy import deepcopy
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import pytest

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def studio(tmp_path, monkeypatch):
    import streamlit as st
    from streamlit.testing.v1 import AppTest
    from biocentric.service import Service

    monkeypatch.setenv("BIOCENTRIC_HOME", str(tmp_path))
    monkeypatch.delenv("BIOCENTRIC_UI_MODE", raising=False)
    monkeypatch.delenv("OPENAI_API_KEY", raising=False)
    service = Service(tmp_path)
    saved_id = service.create_research_project("Existing ABL experiment")
    original = deepcopy(service.snapshot(saved_id))
    from biocentric import discovery_evidence as evidence
    identity = {"accession": "P00519", "name": "Tyrosine-protein kinase ABL1", "organism": "Homo sapiens", "taxon_id": 9606,
                "gene_names": ["ABL1"], "url": "https://www.uniprot.org/uniprotkb/P00519/entry"}
    monkeypatch.setattr(evidence, "resolve_identity", lambda *args: {"status": "resolved", "candidates": [identity]})
    monkeypatch.setattr(evidence, "overview", lambda *args: {"structures": {"rows": [], "total": 0}, "binding": {"rows": [], "total": 0}, "literature": {"evidence": []}})
    def workspace_app():
        from biocentric.experiment_library import create
        app = AppTest.from_file(str(ROOT / 'app.py'), default_timeout=20)
        app.query_params['experiment'] = create(tmp_path)
        return app
    st.cache_resource.clear()
    with patch("biocentric.service.Service", return_value=service) as constructor, \
            patch.object(service, "create_research_project") as create, \
            patch.object(service, "start_research") as start, \
            patch("biocentric.research_agent.research_literature") as research, \
            patch("biocentric.research_agent.recommend_tools") as recommend:
        yield SimpleNamespace(
            app=workspace_app,
            service=service, saved_id=saved_id, constructor=constructor, root=tmp_path,
        )
        # Every tested interaction, including acceptance, leaves the existing
        # persisted project untouched and never invokes planning or a worker.
        create.assert_not_called()
        start.assert_not_called()
        research.assert_not_called()
        recommend.assert_not_called()
        assert service.snapshot(saved_id) == original
        assert len(service.list_projects()) == 1
    st.cache_resource.clear()


def enter_scope(app):
    app.text_input(key="studio_chat_input").set_value("ABL1").run()
    assert not app.exception
    return app


def test_initial_screen_has_two_multiple_choices_and_empty_workspace(studio):
    app = studio.app().run()
    assert not app.exception
    assert app.multiselect(key="studio_tools_docking").value == ["vina", "gnina", "unidock"]
    assert len(app.multiselect(key="studio_tools_scoring").options) == 8
    assert len(app.multiselect(key="studio_tools_docking").options) == 6
    assert not any('preparation' in (w.key or '') for w in app.selectbox)
    assert not app.dataframe
    assert any('aria-label="Empty experiments workspace"' in item.value for item in app.markdown)
    studio.constructor.assert_not_called()


def test_exact_cartesian_candidate_pool_and_evidence_are_shown(studio):
    from biocentric.discovery_studio import combinations
    app = enter_scope(studio.app().run())
    assert {'Pipeline map', 'Final results', 'Final visualizations'} <= {tab.label for tab in app.tabs}
    app.multiselect(key="studio_tools_docking").set_value(["vina", "autodock4"]).run()
    app.multiselect(key="studio_tools_scoring").set_value(["vina_score", "planet", "mmgbsa"]).run()
    assert len(combinations(app.session_state.filtered_state)) == 6
    assert 'preparation' not in app.session_state
    assert app.session_state['studio_identity']['accession'] == 'P00519'
    assert not app.exception
    draft = app.session_state["studio_draft"]
    assert draft["status"] == "review"
    assert app.button(key="studio_eval_start").label == "Start"
    assert not any(b.label == "Accept candidate pool" for b in app.button)
    assert any('Predictive validation unavailable' in m.value for m in app.info)
    studio.constructor.assert_not_called()


def test_changing_pool_keeps_start_available(studio):
    app = enter_scope(studio.app().run())
    app.multiselect(key="studio_tools_docking").set_value(["unidock"]).run()
    assert not app.exception
    assert app.session_state["studio_draft"]["status"] == "review"
    assert "accepted_selection" not in app.session_state["studio_draft"]


def test_hyperparameters_persist_when_tool_is_removed_and_restored(studio):
    app = enter_scope(studio.app().run())
    app.number_input(key="studio_parameter_docking_vina_exhaustiveness").set_value(24).run()
    assert app.session_state["studio_draft"]["status"] == "review"
    app.multiselect(key="studio_tools_docking").set_value(["gnina"]).run()
    assert app.number_input(key="studio_parameter_docking_gnina_exhaustiveness").value == 8
    app.multiselect(key="studio_tools_docking").set_value(["vina", "gnina"]).run()
    assert app.number_input(key="studio_parameter_docking_vina_exhaustiveness").value == 24


def add_custom(app):
    app.text_input(key="studio_plugin_name").set_value("Schrödinger Glide")
    app.selectbox(key="studio_plugin_phase").set_value("docking")
    next(b for b in app.button if b.label == "Add").click().run()
    assert not app.exception
    return app.multiselect(key="studio_tools_docking").value[-1]


def test_custom_plugin_joins_pool_without_replacing_existing_tools(studio):
    app = studio.app().run()
    tool = add_custom(app)
    assert app.multiselect(key="studio_tools_docking").value == ["vina", "gnina", "unidock", tool]
    assert "Schrödinger Glide" not in app.multiselect(key="studio_tools_scoring").options
    app.text_area(key=f"studio_json_docking_{tool}").set_value('{"precision": "SP"}').run()
    enter_scope(app)
    from biocentric.discovery_studio import selection
    assert selection(app.session_state.filtered_state)["docking"][-1] == {"tool": tool, "parameters": {"precision": "SP"}}


def test_invalid_parameters_and_empty_choices_block_start(studio):
    app = studio.app().run()
    tool = add_custom(app)
    enter_scope(app)
    app.text_area(key=f"studio_json_docking_{tool}").set_value('not json').run()
    assert app.error
    assert app.button(key="studio_eval_start").disabled
    assert "accepted_selection" not in app.session_state["studio_draft"]
    app.multiselect(key="studio_tools_docking").set_value([]).run()
    assert app.button(key="studio_eval_start").disabled


def test_wrong_md_frame_range_blocks_acceptance(studio):
    app = enter_scope(studio.app().run())
    app.multiselect(key="studio_tools_scoring").set_value(["mmgbsa"]).run()
    app.number_input(key="studio_parameter_scoring_mmgbsa_start_frame").set_value(100).run()
    assert app.button(key="studio_eval_start").disabled


def test_new_target_does_not_keep_old_evidence(studio, monkeypatch):
    from biocentric import discovery_evidence as evidence
    app = enter_scope(studio.app().run())
    assert app.session_state["studio_identity"]["accession"] == "P00519"
    monkeypatch.setattr(evidence, 'resolve_identity', lambda *args: {'status': 'not_found', 'candidates': []})
    app.text_input(key="studio_chat_input").set_value("UNKNOWN").run()
    assert app.session_state['studio_draft']['target'] == "UNKNOWN"
    assert not app.exception
    assert 'studio_identity' not in app.session_state
    assert not any('EXPERIMENTAL EVIDENCE' in m.value for m in app.markdown)


def test_target_opens_directly_with_human_default(studio):
    app = studio.app().run()
    app.text_input(key="studio_chat_input").set_value("ABL1").run()
    assert not app.exception
    assert app.session_state['studio_draft']['evidence_query'] == ('ABL1', 'Human')
    assert app.session_state['studio_identity']['accession'] == 'P00519'
    assert not any(w.key == 'studio_organism' for w in app.text_input)
    assert not any(b.label == 'Find evidence' for b in app.button)


def test_explicit_project_url_preserves_saved_workbench_entry(studio):
    app = studio.app()
    app.query_params["project"] = studio.saved_id
    app.run()
    assert not app.exception
    assert app.selectbox(key="research_project_picker").value == studio.saved_id


def test_structure_pagination_replaces_structure_and_keeps_correct_coordinate_format(studio, monkeypatch):
    from biocentric import discovery_evidence as evidence
    def page(number):
        return {'total': 48, 'rows': [{'id': '1AB2' if number else '1OPL', 'title': 'Experimental structure',
                 'compounds': [], 'target_chains': ['A'], 'mutations': [], 'citation': {}, 'partners': [],
                 'resolution': None if number else 2.0, 'models': 20 if number else 1,
                 'methods': ['SOLUTION NMR' if number else 'X-RAY DIFFRACTION'], 'url': 'https://www.rcsb.org/'}]}
    monkeypatch.setattr(evidence, 'overview', lambda *a: {'structures': page(0), 'binding': {'rows': [], 'total': 0}, 'literature': {'evidence': []}})
    monkeypatch.setattr(evidence, 'structure_page', lambda identity, root, number: page(number))
    formats = []
    def coords(pdb, root, format_name):
        formats.append((pdb, format_name))
        return {}
    monkeypatch.setattr(evidence, 'coordinates', coords)
    from biocentric import structure_gallery
    monkeypatch.setattr(structure_gallery, 'gallery_html', lambda *a: '<p>Experimental geometry gallery</p>')
    app = enter_scope(studio.app().run())
    assert formats[-1] == ('1OPL', 'cif')
    app.button(key='studio_gallery_next').click().run()
    assert not app.exception
    assert app.session_state['studio_gallery_loaded']['items'][0][0]['id'] == '1AB2'
    assert formats[-1] == ('1AB2', 'pdb')
    app.multiselect(key='studio_tools_docking').set_value(['vina']).run()
    assert not app.exception
    assert app.session_state['studio_gallery_loaded']['items'][0][0]['id'] == '1AB2'


def test_target_then_one_start_click_launches_once(studio, monkeypatch, tmp_path):
    from biocentric import discovery_evaluation_ui as ui, discovery_panels, discovery_benchmark as bench
    from biocentric.discovery_studio import selection
    calls = []
    def prepare(state, root, progress):
        progress('Preparing test evaluation')
        state['studio_eval_review'] = {'identity': state['studio_identity'], 'selection': selection(state)}
        return state['studio_eval_review']
    monkeypatch.setattr(ui, 'prepare_start', prepare)
    monkeypatch.setattr(bench, 'start_run', lambda *a: calls.append(a) or str(tmp_path / ('a'*32)))
    monkeypatch.setattr(bench, 'read_status', lambda *a: {'status': 'running'})
    original = discovery_panels.render_experiments
    monkeypatch.setattr(discovery_panels, 'render_experiments', lambda st, draft: st.caption('Running') if st.session_state.get('studio_eval_run') else original(st, draft))
    app = enter_scope(studio.app().run())
    assert not calls
    app.button(key='studio_eval_start').click().run()
    assert not app.exception and len(calls) == 1
    assert app.button(key='studio_eval_stop').label == 'Stop'
    app.run()
    assert not app.exception and len(calls) == 1


def test_start_failure_is_visible_and_retryable(studio, monkeypatch):
    from biocentric import discovery_evaluation_ui as ui
    def unavailable(*args):
        raise ValueError('No compatible reference found')
    monkeypatch.setattr(ui, 'prepare_start', unavailable)
    app = enter_scope(studio.app().run())
    app.button(key='studio_eval_start').click().run()
    assert not app.exception
    assert any('No compatible reference' in error.value for error in app.error)
    assert not app.button(key='studio_eval_start').disabled
    assert not app.session_state.filtered_state.get('studio_eval_start_requested')


def test_landing_create_back_and_restore_draft(studio):
    from streamlit.testing.v1 import AppTest
    from biocentric import experiment_library as library
    app=AppTest.from_file(str(ROOT/'app.py'),default_timeout=20).run()
    assert not app.exception and app.button(key='experiments_new').label=='Create new experiment'
    app.button(key='experiments_new').click().run()
    first=app.session_state['experiments_current']
    enter_scope(app)
    app.number_input(key='studio_parameter_docking_vina_exhaustiveness').set_value(24).run()
    app.button(key='experiments_back').click().run()
    assert not app.exception and any('experiment-table' in item.value for item in app.markdown)
    assert library.load(studio.root,first)['state']['studio_draft']['target']=='ABL1'
    app.button(key='experiments_new').click().run()
    assert not app.exception
    assert app.session_state['experiments_current']!=first
    assert not app.session_state['studio_draft']['messages']
    assert app.number_input(key='studio_parameter_docking_vina_exhaustiveness').value==8
    app.query_params['experiment']=first
    app.run()
    assert not app.exception and app.session_state['studio_draft']['target']=='ABL1'
    assert app.number_input(key='studio_parameter_docking_vina_exhaustiveness').value==24


def test_restored_tool_widgets_do_not_duplicate_defaults(studio, monkeypatch):
    from streamlit.testing.v1 import AppTest
    from streamlit.elements.lib import policies
    app = enter_scope(studio.app().run())
    experiment_id = app.session_state['experiments_current']
    app.multiselect(key='studio_tools_structure_preparation').set_value(['pdb2pqr', 'reduce']).run()
    app.button(key='experiments_back').click().run()
    monkeypatch.setattr(policies, '_shown_default_value_warning', False)
    reopened = AppTest.from_file(str(ROOT/'app.py'), default_timeout=20)
    reopened.query_params['experiment'] = experiment_id
    reopened.run()
    assert not reopened.exception
    assert reopened.multiselect(key='studio_tools_structure_preparation').value == ['pdb2pqr','reduce']
    assert not any('Session State API' in item.value for item in reopened.warning)


def test_orchestrator_uses_one_line_target_and_no_chat_history(studio):
    app = enter_scope(studio.app().run())
    assert not app.chat_input and not app.chat_message
    assert app.text_input(key='studio_chat_input').value == 'ABL1'
    assert any('Orchestrator' in m.value for m in app.markdown)
    assert len(app.multiselect) == 4
    assert app.button(key='studio_eval_start').label == 'Start'


def test_library_delete_confirmation_and_undo(studio):
    app = studio.app().run()
    experiment_id = app.session_state['experiments_current']
    app.button(key='experiments_back').click().run()
    app.query_params['delete_experiment'] = experiment_id
    app.run()
    assert not app.exception
    next(b for b in app.button if b.label == 'Keep experiment').click().run()
    from biocentric import experiment_library as library
    assert library.load(studio.root, experiment_id)
    app.query_params['delete_experiment'] = experiment_id
    app.run()
    next(b for b in app.button if b.label == 'Delete experiment').click().run()
    assert not app.exception and library.is_deleted(studio.root, experiment_id)
    app.button(key='experiments_undo_delete').click().run()
    assert not app.exception and not library.is_deleted(studio.root, experiment_id)


def test_biocentric_selection_has_bounded_cuda_and_batch_settings(studio):
    app = studio.app().run()
    app.multiselect(key="studio_tools_scoring").set_value(["biocentric_affinity"]).run()
    assert not app.exception
    device = app.selectbox(key="studio_parameter_scoring_biocentric_affinity_device")
    assert device.value == "cuda" and device.options == ["cpu", "cuda"]
    app.number_input(key="studio_parameter_scoring_biocentric_affinity_batch_size").set_value(8).run()
    assert app.session_state["studio_parameters"]["scoring:biocentric_affinity"]["batch_size"] == 8
    assert not app.exception


def test_fast_balanced_and_custom_presets_match_frozen_selection(studio):
    from biocentric.discovery_studio import selection
    app = enter_scope(studio.app().run())
    assert app.selectbox(key='studio_docking_preset').value == 'Balanced'
    assert app.number_input(key='studio_parameter_docking_vina_max_evals').value == 0
    assert app.number_input(key='studio_parameter_docking_vina_num_modes').value == 9
    app.selectbox(key='studio_docking_preset').set_value('Balanced').run()
    assert not app.exception
    assert app.number_input(key='studio_parameter_docking_vina_exhaustiveness').value == 8
    assert app.number_input(key='studio_parameter_docking_vina_max_evals').value == 0
    app.number_input(key='studio_parameter_docking_vina_max_evals').set_value(10000).run()
    assert not app.exception and app.selectbox(key='studio_docking_preset').value == 'Custom'
    pair = next(x for x in selection(app.session_state.filtered_state)['docking'] if x['tool'] == 'vina')
    assert pair['parameters']['max_evals'] == 10000
    app.selectbox(key='studio_docking_preset').set_value('Fast').run()
    assert not app.exception
    assert app.number_input(key='studio_parameter_docking_vina_max_evals').value == 5000
    assert app.selectbox(key='studio_parameter_docking_unidock_search_mode').value == 'fast'
