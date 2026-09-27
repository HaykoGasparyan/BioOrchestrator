"""Conversation UI behavior is tested with persisted-looking offline records only."""
from __future__ import annotations

from copy import deepcopy
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import pytest

from biocentric.research_ui import (conversation_messages, pending_questions,
                                    research_summary, source_rows)


PID = "f" * 32
ROOT = Path(__file__).resolve().parents[1]


def snapshot():
    return {
        "project": {"id": PID, "state": "blocked", "spec": {
            "name": "Public docking research", "workflow_mode": "research", "model": "gpt-4.1-mini",
            "max_model_calls": 6, "allow_model": True, "local_only": False,
        }},
        "data": {"research": {"mode": "research", "goal": "Compare suitable docking tools on a small test set.",
            "questions": [{"id": "protocol_review", "prompt": "Use this reviewed receptor and binding site?",
                           "status": "pending", "options": [{"label": "Approve this protocol", "value": "approve"},
                                                              {"label": "Revise the inputs", "value": "revise"}]}]}},
        "steps": [{"id": "05", "label": "Choose a docking method", "status": "blocked",
                   "reason": "The literature and runtime support a small Vina pilot.",
                   "inputs": {"ligands": 3}, "preprocessing": "Verify receptor coordinates.",
                   "postprocessing": "Keep engine score names and units.",
                   "outputs": {"public_explanation": "Vina can run on the local CPU. Pose recovery is checked before screening.",
                               "members": [{"compound_id": "imatinib", "role": "reference"}],
                               "sources": [{"title": "Vina documentation", "url": "https://vina.scripps.edu/", "id": "vina-docs"}],
                               "limitations": ["A redocking pilot does not validate binding affinity."]}}],
        "events": [{"id": 1, "kind": "research_message", "stage": "05", "role": "assistant",
                    "message": "I found a CPU-compatible option. Please review its small test protocol."},
                   {"id": 2, "kind": "model_call", "stage": "05", "name": "Method selection", "status": "running", "inputs": {"model": "gpt-4.1-mini"}},
                   {"id": 3, "kind": "model_call", "stage": "05", "name": "Method selection", "status": "completed",
                    "outputs": {"response_id": "response-1", "usage": {"input_tokens": 80, "output_tokens": 15}, "text": "Use a Vina pilot."}}],
        "jobs": [{"id": "job-1", "status": "blocked"}], "artifacts": [], "evidence": [], "decisions": [], "candidates": [],
    }


class FakeService:
    def __init__(self, records=None):
        self.records = records or {}
        self.created = []
        self.started = []
        self.answers = []
        self.cancelled = []
        self.imported_measurements = []

    def list_projects(self):
        return [{**record["project"], "data": record["data"]} for record in self.records.values()]

    def snapshot(self, pid):
        return deepcopy(self.records[pid])

    def create_research_project(self, goal, **kwargs):
        self.created.append({"goal": goal, **kwargs})
        record = snapshot()
        record["project"]["state"] = "draft"
        record["data"]["research"]["goal"] = goal
        record["data"]["research"]["questions"] = []
        record["steps"] = []
        record["events"] = []
        record["jobs"] = []
        self.records[PID] = record
        return PID

    def start_research(self, pid):
        self.started.append(pid)
        self.records[pid]["jobs"] = [{"id": "job-1", "status": "running"}]
        self.records[pid]["project"]["state"] = "running"
        return "job-1"

    def answer_research(self, pid, qid, answer):
        self.answers.append((pid, qid, answer))
        question = next(item for item in self.records[pid]["data"]["research"]["questions"] if item["id"] == qid)
        question.update(status="answered", answer=answer)
        return self.start_research(pid)

    def cancel(self, pid):
        self.cancelled.append(pid)
        self.records[pid]["jobs"][0]["status"] = "cancelled"
        self.records[pid]["project"]["state"] = "cancelled"
        return True

    def import_measurements(self, pid, content):
        self.imported_measurements.append((pid, content))
        return {"count": 1}


@pytest.fixture
def ui(tmp_path, monkeypatch):
    import streamlit as st
    from streamlit.testing.v1 import AppTest

    monkeypatch.setenv("BIOCENTRIC_HOME", str(tmp_path))
    monkeypatch.delenv("BIOCENTRIC_UI_MODE", raising=False)
    monkeypatch.delenv("OPENAI_API_KEY", raising=False)
    st.cache_resource.clear()

    def build(service):
        app = AppTest.from_file(str(ROOT / "app.py"), default_timeout=20)
        # These tests preserve the executable saved-workbench contract. The
        # fresh-root discovery brief is exercised independently.
        app.query_params["view"] = "workbench"
        return app

    yield build
    st.cache_resource.clear()


def test_only_public_messages_are_rendered_and_secret_strings_are_redacted():
    record = snapshot()
    record["events"].extend([
        {"kind": "tool_call", "message": "This is a tool log, not an assistant statement."},
        {"id": 4, "kind": "research_message", "payload": {"role": "user", "text": "sk-abcdefghijklmnopqrstuvwxyz123456"}},
    ])
    actual = conversation_messages(record)
    assert len(actual) == 2
    assert actual[1]["role"] == "user"
    assert actual[1]["message"] == "[redacted]"
    assert all("tool log" not in item["message"] for item in actual)


def test_usage_counts_recorded_requests_and_unique_responses_only():
    record = snapshot()
    record["events"].append(deepcopy(record["events"][-1]))
    result = research_summary(record)
    assert result["calls"] == 1
    assert result["input_tokens"] == 80
    assert result["output_tokens"] == 15
    assert result["current"]["id"] == "05"
    assert not result["active"]
    assert research_summary({})["input_tokens"] is None


def test_answered_questions_and_unsafe_source_links_are_not_actionable():
    record = snapshot()
    record["data"]["research"]["questions"].append({"id": "old", "prompt": "Old question", "status": "answered", "answer": "yes"})
    assert [item["id"] for item in pending_questions(record)] == ["protocol_review"]
    rows = source_rows([{"title": "unsafe", "url": "javascript:alert(1)"},
                        {"title": "credentials", "url": "https://user:secret@example.com/paper"},
                        {"title": "valid", "url": "https://example.com/paper"}])
    assert rows[0]["url"] is None
    assert rows[1]["url"] is None
    assert rows[2]["url"] == "https://example.com/paper"


def test_workbench_is_conversation_and_public_example_requires_explicit_start(ui):
    service = FakeService()
    with patch("biocentric.service.Service", return_value=service):
        app = ui(service).run()
        assert not app.exception
        assert app.radio(key="workspace_mode").value == "Research conversation"
        assert any("Make the next experiment clear" in item.value for item in app.title)
        next(button for button in app.button if button.label == "Prepare public example").click().run()
        assert not app.exception
        assert service.created[0]["use_public_example"] is True
        assert service.created[0]["max_model_calls"] == 6
        assert service.started == []
        next(button for button in app.button if button.label == "Start research").click().run()
        assert not app.exception
        assert service.started == [PID]
        assert any(item.value.startswith("Working:") for item in app.info)


def test_chat_submission_creates_and_starts_a_bounded_research_project(ui):
    service = FakeService()
    with patch("biocentric.service.Service", return_value=service):
        app = ui(service).run()
        app.chat_input(key="research_goal").set_value("Which docking method fits human ABL1?").run()
        assert not app.exception
        assert service.created == [{"goal": "Which docking method fits human ABL1?", "use_public_example": False,
                                    "model": "gpt-4.1-mini", "max_model_calls": 6}]
        assert service.started == [PID]


def test_stage_cards_stay_visible_and_allow_result_and_evidence_navigation(ui):
    service = FakeService({PID: snapshot()})
    with patch("biocentric.service.Service", return_value=service):
        app = ui(service).run()
        assert not app.exception
        assert app.session_state[f"research_view_{PID}"] == "05"
        assert "### Research steps" in [item.value for item in app.markdown]
        assert any("Pose recovery is checked" in item.value for item in app.markdown)
        assert "**AI decision**" in [item.value for item in app.markdown]
        calls = next(item for item in app.expander if item.label.startswith("Exact recorded steps & calls"))
        assert not calls.proto.expanded
        assert any("Vina documentation" in item.value for item in app.markdown)
        assert len(app.dataframe) >= 1
        app.button(key=f"research_results_view_{PID}").click().run()
        assert not app.exception
        assert app.session_state[f"research_view_{PID}"] == "results"
        # The step and its reason remain visible while inspecting the ligands.
        assert any("Pose recovery is checked" in item.value for item in app.markdown)
        assert not any(item.label.startswith("Exact recorded steps & calls") for item in app.expander)
        app.button(key=f"research_toggle_{PID}_05").click().run()
        assert not app.exception
        assert app.session_state[f"research_view_{PID}"] == "05"
        app.run()
        assert app.session_state[f"research_view_{PID}"] == "05"
        assert not service.started


def test_in_chat_answer_passes_option_value_and_automatically_continues(ui):
    service = FakeService({PID: snapshot()})
    with patch("biocentric.service.Service", return_value=service):
        app = ui(service).run()
        app.radio(key=f"research_choice_{PID}_protocol_review").set_value(0)
        next(button for button in app.button if button.label == "Send & continue").click().run()
        assert not app.exception
        assert service.answers == [(PID, "protocol_review", "approve")]
        assert service.started == [PID]
        assert not any(button.label == "Send & continue" for button in app.button)
        assert any(button.label == "Cancel current work" for button in app.button)


def test_empty_question_does_not_continue_and_free_text_is_supported(ui):
    record = snapshot()
    record["data"]["research"]["questions"] = [{"id": "target_identity", "prompt": "Which accession and species?", "status": "pending"}]
    service = FakeService({PID: record})
    with patch("biocentric.service.Service", return_value=service):
        app = ui(service).run()
        next(button for button in app.button if button.label == "Send & continue").click().run()
        assert not app.exception
        assert not service.answers
        app.text_area(key=f"research_answer_text_{PID}_target_identity").set_value("P00519; Homo sapiens")
        next(button for button in app.button if button.label == "Send & continue").click().run()
        assert not app.exception
        assert service.answers == [(PID, "target_identity", "P00519; Homo sapiens")]


def test_assay_upload_uses_measurement_import_and_offers_csv_template(ui):
    service = FakeService({PID: snapshot()})
    csv = b"compound_id,canonical_smiles,endpoint,value,units,relation,assay_id,organism\nc1,c1ccccc1,Kd,10,nM,=,a1,Homo sapiens\n"
    upload = SimpleNamespace(name="measured.csv", getvalue=lambda: csv)
    with patch("biocentric.service.Service", return_value=service), patch("streamlit.file_uploader", return_value=upload):
        app = ui(service).run()
        app.selectbox(key=f"research_upload_kind_{PID}_False").set_value("measurements").run()
        assert not app.exception
        assert any(button.label == "Download assay CSV headers" for button in app.get("download_button"))
        app.button(key=f"research_save_upload_{PID}").click().run()
        assert not app.exception
        assert service.imported_measurements == [(PID, csv)]
        assert not service.started


def test_input_revision_has_explicit_rebuild_action_and_upload_instructions(ui):
    record = snapshot()
    record["data"]["research"]["questions"] = [{
        "id": "input_revision_1", "kind": "research_input_revision", "status": "pending",
        "prompt": "Upload revised inputs. Rebuilding uses 1 GPT call.",
        "options": [{"label": "Rebuild test set and proposal · 1 GPT call", "value": "rebuild"}],
    }]
    service = FakeService({PID: record})
    with patch("biocentric.service.Service", return_value=service):
        app = ui(service).run()
        assert not app.exception
        assert any("Add scientific input files below the research workspace" in item.value for item in app.caption)
        assert not any(button.label == "Send & continue" for button in app.button)
        app.radio(key=f"research_choice_{PID}_input_revision_1").set_value(0)
        next(button for button in app.button if button.label == "Rebuild test set & proposal").click().run()
        assert not app.exception
        assert service.answers == [(PID, "input_revision_1", "rebuild")]


def test_public_stream_is_live_only_and_private_fields_never_become_chat_text():
    from biocentric.research_ui import live_model_output

    record = {"events": [
        {"kind": "model_call", "status": "running", "inputs": {"purpose": "tool_selection"}},
        {"kind": "model_delta", "outputs": {"text_delta": "private text", "public_output": False}},
        {"kind": "model_delta", "outputs": {"text_delta": '{"docking_engine":"vina","reason":"Runs on the local CPU', "public_output": True}},
        {"kind": "model_delta", "outputs": {"text_delta": ' and has a documented reference protocol.","citations":[]}', "public_output": True}},
    ]}
    draft = live_model_output(record)
    assert draft["text"] == "Runs on the local CPU and has a documented reference protocol."
    assert "private text" not in draft["text"]
    record["events"].append({"kind": "model_call", "status": "completed", "outputs": {}})
    assert live_model_output(record) is None


def test_citations_join_actual_records_and_validated_spans_only():
    from biocentric.research_ui import citation_markdown, resolve_sources

    record = snapshot()
    record["data"]["research"]["sources"] = [{"id": "vina-doc", "title": "Vina docs", "url": "https://vina.scripps.edu/"}]
    assert resolve_sources(record, ["unknown", "vina-doc"])[0]["title"] == "Vina docs"
    citation = {"evidence_id": "vina-doc", "url": "https://vina.scripps.edu/", "start_index": 8, "end_index": 14}
    assert resolve_sources(record, [citation])[0]["title"] == "Vina docs"
    assert citation_markdown("A claim [cite]", [citation]) == "A claim [source](https://vina.scripps.edu/)"
    assert citation_markdown("A claim", [{**citation, "start_index": 50}]) == "A claim"
    assert citation_markdown("A claim [cite]", [{**citation, "url": "javascript:alert(1)"}]) == "A claim [cite]"


def test_real_service_public_example_stays_draft_and_preserves_typed_inputs(ui, tmp_path):
    from biocentric.service import Service

    service = Service(tmp_path)
    with patch("biocentric.service.Service", return_value=service), patch.object(service, "start_research") as start:
        app = ui(service).run()
        next(button for button in app.button if button.label == "Prepare public example").click().run()
        assert not app.exception
        assert not app.error
        project = service.list_projects()[0]
        assert project["spec"]["workflow_mode"] == "research"
        assert project["spec"]["max_model_calls"] == 6
        state = service.snapshot(project["id"])
        assert state["jobs"] == []
        assert {"reference", "receptor", "structure", "library"}.issubset({item["kind"] for item in state["artifacts"]})
        start.assert_not_called()
        # The persisted user message is rendered once, rather than duplicated
        # by the independent goal preview.
        messages = [item.value for item in app.markdown]
        assert messages.count(state["data"]["research"]["goal"]) == 1


def test_old_classic_selection_does_not_return_after_new_conversation_and_prepare(ui, tmp_path):
    from biocentric.service import Service

    service = Service(tmp_path)
    classic_id = service.create_project({"name": "Saved classic project", "target": "P00519", "organism": "Homo sapiens"})
    with patch("biocentric.service.Service", return_value=service), patch.object(service, "start_research") as start:
        app = ui(service)
        app.session_state["project_id"] = classic_id
        app.session_state["research_project_id"] = classic_id
        app.session_state["research_project_picker"] = classic_id
        app.run()
        assert not app.exception
        assert app.selectbox(key="research_project_picker").value == classic_id
        assert any("original workbench" in item.value for item in app.info)
        app.button(key="new_research_conversation").click().run()
        assert not app.exception
        assert app.selectbox(key="research_project_picker").value is None
        assert any("Make the next experiment clear" in item.value for item in app.title)
        app.button(key="research_public_example").click().run()
        assert not app.exception
        assert not app.error
        projects = service.list_projects()
        assert len(projects) == 2
        research_id = next(project["id"] for project in projects if project["id"] != classic_id)
        assert app.selectbox(key="research_project_picker").value == research_id
        assert not app.warning
        assert app.session_state["research_project_id"] == research_id
        assert any(button.label == "Start research" for button in app.button)
        assert not any("original workbench" in item.value for item in app.info)
        app.run()
        assert not app.exception
        assert app.selectbox(key="research_project_picker").value == research_id
        assert service.snapshot(research_id)["jobs"] == []
        start.assert_not_called()


def test_duplicate_project_names_have_distinct_labels_and_url_selection_survives_reload(ui, tmp_path):
    from biocentric.service import Service

    service = Service(tmp_path)
    classic_id = service.create_example(allow_model=True)
    older_id = service.create_research_project("First research question", use_public_example=True)
    newer_id = service.create_research_project("Second research question", use_public_example=True)
    assert len({project["spec"]["name"] for project in service.list_projects()}) == 1
    with patch("biocentric.service.Service", return_value=service):
        app = ui(service)
        app.query_params["project"] = older_id
        app.run()
        assert not app.exception
        assert not app.warning
        picker = app.selectbox(key="research_project_picker")
        assert picker.value == older_id
        assert len(picker.options) == len(set(picker.options))
        assert any(f"classic · {classic_id[:6]}" in option for option in picker.options)
        assert any(f"research · {older_id[:6]}" in option for option in picker.options)
        assert any(f"Project {older_id[:6]}" in caption.value for caption in app.caption)
        picker.set_value(newer_id).run()
        assert not app.exception
        assert app.selectbox(key="research_project_picker").value == newer_id
        assert app.query_params["project"] == [newer_id]
        app.run()
        assert app.selectbox(key="research_project_picker").value == newer_id
        # A fresh session models a browser reload; only the URL selection survives.
        reloaded = ui(service)
        reloaded.query_params.update(app.query_params)
        reloaded.run()
        assert not reloaded.exception
        assert not reloaded.warning
        assert reloaded.selectbox(key="research_project_picker").value == newer_id
        reloaded.button(key="new_research_conversation").click().run()
        assert "project" not in reloaded.query_params
        assert reloaded.selectbox(key="research_project_picker").value is None


def test_compact_summary_is_a_short_extract_and_testset_limits_stay_explicit():
    from biocentric.research_ui import compact_explanation, poc_testset_summary

    text = "Vina runs on the configured CPU. Published pose-recovery protocols support a small pilot. A third long discussion stays in details."
    assert compact_explanation(text) == "Vina runs on the configured CPU. Published pose-recovery protocols support a small pilot."
    assert len(compact_explanation("one " * 300, limit=100)) <= 100
    record = snapshot()
    record["data"]["research"]["testset"] = {"records": [], "pose_records": [{"id": "reference"}], "evaluation_ids": []}
    brief = poc_testset_summary(record)
    assert "no independent assay benchmark" in brief
    assert "not binding-affinity accuracy" in brief
    record["data"]["research"]["testset"]["evaluation_ids"] = ["a", "b"]
    assert "not itself a validation result" in poc_testset_summary(record)


def test_completed_workspace_keeps_results_prominent_alongside_visible_stage_reasons(ui):
    import sys
    from biocentric.research_ui import compact_explanation

    record = snapshot()
    record["project"]["state"] = "completed"
    record["jobs"][0]["status"] = "completed"
    record["steps"][0]["status"] = "completed"
    record["data"]["research"]["questions"] = []
    reason = "Vina runs on the configured CPU. The source supports testing pose recovery first. " + "Further technical discussion. " * 30
    record["data"]["research"]["tool_plan"] = {"docking_engine": "vina", "reason": reason}
    record["data"]["research"]["testset"] = {"pose_records": [{"id": "r1"}], "records": [], "evaluation_ids": []}
    record["events"].append({"id": 4, "kind": "research_message", "role": "assistant", "message": reason})
    service = FakeService({PID: record})
    received = []

    def results(st, state, root, pid):
        received.append((state, root, pid))
        st.markdown("### Ranked ligand results")
        st.download_button("Download selected ligands", b"actual test bytes", file_name="test.sdf")

    with patch("biocentric.service.Service", return_value=service), patch.dict(sys.modules, {"biocentric.ligand_results": SimpleNamespace(render_ligand_results=results)}):
        app = ui(service).run()
        assert not app.exception
        assert received[0][2] == PID and received[0][0] == record
        visible_texts = [item.value for item in app.markdown]
        assert compact_explanation(reason) in visible_texts
        assert any("no independent assay benchmark" in item.value for item in app.caption)
        history = next(item for item in app.expander if item.label == "Conversation history & full model responses")
        assert not history.proto.expanded
        assert app.session_state[f"research_view_{PID}"] == "results"
        assert "### Research steps" in visible_texts
        assert any("Pose recovery is checked" in text for text in visible_texts)
        assert app.button(key=f"research_toggle_{PID}_05")
        assert not any(item.label == "Research steps & technical details" for item in app.expander)
        assert reason.strip() in [item.value for item in history.markdown]
        assert visible_texts.index("### Ranked ligand results") < visible_texts.index(compact_explanation(reason))
        assert visible_texts.index("### Ranked ligand results") < visible_texts.index(reason.strip())
        assert len(received) == 1
        assert any(button.label == "Download selected ligands" for button in app.get("download_button"))
        assert not any(button.label in {"Start research", "Continue research"} for button in app.button)
        route = next(text for text in visible_texts if 'class="poc-route"' in text)
        assert all(label in route for label in ["Research", "Test set", "Method", "Pilot", "Screen", "Results"])


def test_live_update_remains_short_while_full_public_draft_is_optional(ui):
    record = snapshot()
    record["project"]["state"] = "running"
    record["jobs"][0]["status"] = "running"
    record["data"]["research"]["questions"] = []
    draft = "I am checking the local docking prerequisites. The pilot will use the reference pose. " + "Detailed source discussion. " * 30
    record["events"].extend([
        {"id": 4, "kind": "model_call", "status": "running", "inputs": {"purpose": "research"}},
        {"id": 5, "kind": "model_delta", "outputs": {"public_output": True, "text_delta": draft}},
    ])
    service = FakeService({PID: record})
    with patch("biocentric.service.Service", return_value=service):
        app = ui(service).run()
        assert not app.exception
        assert any("Live GPT update" in item.value for item in app.caption)
        assert "I am checking the local docking prerequisites. The pilot will use the reference pose." in [item.value for item in app.markdown]
        history = next(item for item in app.expander if item.label == "Conversation history & full model responses")
        assert draft.strip() in [item.value for item in history.markdown]
        assert not history.proto.expanded
        assert any(button.label == "Cancel current work" for button in app.button)


def test_archived_files_are_read_only_after_selection(ui):
    record = snapshot()
    record["artifacts"] = [
        {"id": "a" * 32, "name": "result.csv", "kind": "output", "sha256": "1" * 64},
        {"id": "b" * 32, "name": "result.csv", "kind": "output", "sha256": "2" * 64},
    ]
    service = FakeService({PID: record})
    with patch("biocentric.service.Service", return_value=service), patch("biocentric.research_ui.artifact_bytes", return_value=b"score\n-7\n") as read:
        app = ui(service).run()
        assert not app.exception
        read.assert_not_called()
        picker = app.selectbox(key=f"research_artifact_picker_{PID}")
        assert len(picker.options) == len(set(picker.options))
        picker.set_value("b" * 32).run()
        assert not app.exception
        assert read.call_count == 1
        assert read.call_args.args[0]["id"] == "b" * 32
        assert any(button.label == "Download result.csv" for button in app.get("download_button"))


def test_method_source_links_have_readable_labels_when_source_title_is_a_url():
    from biocentric.research_ui import method_source_label

    vina = "https://autodock-vina.readthedocs.io/en/latest/docking_basic.html"
    assert method_source_label({"title": vina, "url": vina}, 1) == "AutoDock Vina documentation"
    unknown = "https://example.org/very/long/article/path"
    assert method_source_label({"title": unknown, "url": unknown}, 2) == "Source 2"
    assert method_source_label({"title": "A documented docking comparison", "url": unknown}, 1) == "A documented docking comparison"


def test_navigation_follows_current_step_until_user_selects_then_preserves_selection(ui):
    record = snapshot()
    record["data"]["research"]["questions"] = []
    record["project"]["state"] = "running"
    record["jobs"][0]["status"] = "running"
    record["steps"][0]["status"] = "completed"
    record["steps"].append({"id": "07", "label": "Review the measured pilot", "status": "running",
                            "reason": "Checking the measured pilot before screening.", "outputs": {}})
    service = FakeService({PID: record})
    with patch("biocentric.service.Service", return_value=service):
        app = ui(service).run()
        assert not app.exception
        assert app.session_state[f"research_view_{PID}"] == "07"
        app.button(key=f"research_toggle_{PID}_05").click().run()
        assert not app.exception
        assert app.session_state[f"research_view_{PID}"] == "05"
        service.records[PID]["steps"][-1]["status"] = "completed"
        service.records[PID]["steps"].append({"id": "08", "label": "Dock prepared ligands", "status": "running",
                                             "reason": "The approved Vina command is processing prepared ligands.", "outputs": {}})
        app.run()
        assert not app.exception
        assert app.session_state[f"research_view_{PID}"] == "05"
        app.button(key=f"research_follow_current_{PID}").click().run()
        assert not app.exception
        assert app.session_state[f"research_view_{PID}"] == "08"
        assert "Recorded workflow actions" in " ".join(item.value for item in app.caption)
        assert "**AI decision**" not in [item.value for item in app.markdown]
        assert "The approved Vina command is processing prepared ligands." in [item.value for item in app.markdown]
        assert not service.started


def test_pending_question_selects_its_stage_and_ai_labels_need_actual_model_call():
    from biocentric.research_ui import _default_view, stage_has_model_call, stage_public_summary

    record = snapshot()
    record["steps"].insert(0, {"id": "04", "label": "Freeze test set", "status": "completed",
                               "reason": "Build fixed pose controls.", "outputs": {"pose_records": [{"id": "r1"}], "records": []}})
    record["data"]["research"]["questions"][0]["stage_id"] = "04"
    assert _default_view(record) == "04"
    assert not stage_has_model_call(record, "04")
    assert stage_has_model_call(record, "05")
    assert stage_public_summary(record, record["steps"][0]) == "Frozen 0 measured records and 1 pose references before choosing a tool."


def test_local_validation_uses_only_actual_steps_and_never_offers_model_start(ui):
    record = snapshot()
    record["project"]["state"] = "completed"
    record["project"]["spec"].update(allow_model=False, max_model_calls=0, compute_budget_seconds=86400)
    record["jobs"][0]["status"] = "completed"
    record["events"] = []
    record["data"]["research"].update(validation_run={"kind": "local_vina_validation"}, questions=[])
    record["steps"] = [{"id": sid, "label": label, "status": "completed", "reason": reason, "outputs": {}}
                        for sid, label, reason in [
                            ("01", "Dataset selection", "Registered the 128-compound local input library."),
                            ("04", "Preparation", "Prepared the recorded molecular structures."),
                            ("08", "Actual docking", "Ran Vina on the prepared local library."),
                            ("09", "Export checks", "Verified compound and pose downloads."),
                        ]]
    record["steps"][0]["outputs"]["validation_report"] = {
        "resources": {"wall_seconds": 543.21, "cpu_per_vina_process": 1, "llm_calls": 0},
        "protocols": {"baseline": {"duration_seconds": 450.25}},
    }
    service = FakeService({PID: record})
    with patch("biocentric.service.Service", return_value=service):
        app = ui(service).run()
        assert not app.exception
        assert any("Local validation · no planner calls" in item.value for item in app.caption)
        assert any("No literature review or independent assay benchmark" in item.value for item in app.caption)
        texts = [item.value for item in app.markdown]
        assert "### Validation steps" in texts
        route = next(text for text in texts if 'class="poc-route"' in text)
        assert all(label in route for label in ["Dataset selection", "Preparation", "Actual docking", "Export checks"])
        assert "Literature" not in route and "Pilot" not in route
        stages = [button.label for button in app.button if (button.key or "").startswith(f"research_toggle_{PID}_")]
        assert stages == ["01 · Dataset selection", "04 · Preparation", "08 · Actual docking", "09 · Export checks"]
        assert not any(button.label in {"Start research", "Continue research"} for button in app.button)
        assert not any(item.label == "Add scientific input files" for item in app.expander)
        assert not app.get("file_uploader")
        assert "Recorded validation-suite runtime: 543.2 seconds." in texts
        assert not any("86,400 second run budget" in text for text in texts)
        assert not service.started


def test_navigator_is_bounded_and_card_summaries_are_short_while_selected_detail_is_readable(ui):
    import streamlit as st
    from biocentric.research_ui import compact_explanation

    record = snapshot()
    explanation = "Vina uses the configured CPU and accepts the reviewed receptor inputs. The published protocol supports checking pose recovery before any broader compound screen."
    record["steps"][0]["outputs"]["public_explanation"] = explanation
    service = FakeService({PID: record})
    with patch("biocentric.service.Service", return_value=service), patch("streamlit.container", wraps=st.container) as container:
        app = ui(service).run()
        assert not app.exception
        assert any(call.kwargs.get("key") == f"research_navigator_{PID}" and call.kwargs.get("height") == 650
                   for call in container.call_args_list)
        texts = [item.value for item in app.markdown]
        short = compact_explanation(explanation, limit=95)
        assert short in texts and len(short) <= 95
        assert explanation in texts
        app.button(key=f"research_results_view_{PID}").click().run()
        assert not app.exception
        assert short in [item.value for item in app.markdown]
        app.button(key=f"research_toggle_{PID}_05").click().run()
        assert not app.exception
        assert explanation in [item.value for item in app.markdown]
        assert not service.started
