"""Durable research orchestration tests: no model, network or subprocess execution."""
from copy import deepcopy
from pathlib import Path

import pytest

from biocentric.research_worker import execute
from biocentric.research_workflow import AwaitingAnswer, ResearchWorkflow, pilot_checks
from biocentric.service import Service
from biocentric.testsets import build_testset, validate_testset
from biocentric.worker import Cancelled


@pytest.fixture(autouse=True)
def prohibit_network(monkeypatch):
    def fail(*args, **kwargs):
        raise AssertionError("Network calls are forbidden in workflow tests")
    monkeypatch.setattr("httpx.Client.request", fail)


def project(tmp_path, **options):
    service = Service(tmp_path)
    pid = service.create_research_project("Compare docking methods for this protein", **options)
    jid = service.start_research(pid, spawn=False)
    return service, pid, jid


def test_inline_answers_autorequeue_same_job_and_never_repeat_completed_model_step(tmp_path, monkeypatch):
    service, pid, jid = project(tmp_path)
    calls, launches = [], []
    monkeypatch.setattr("biocentric.engines.engine_catalog", lambda: [{"id": "vina", "available": True}])
    monkeypatch.setattr("biocentric.service.subprocess.Popen", lambda *args, **kwargs: launches.append(args[0]) if "biocentric.research_worker" in args[0] else None)
    def run(self):
        self.stage_run("01", "Understand the research question", self.brief)
        def model_step():
            self.emit(kind="model_call", name="Fake advisor", status="running")
            calls.append("model request")
            self.emit(kind="model_call", name="Fake advisor", status="completed", outputs={"usage": {"input_tokens": 2, "output_tokens": 3}})
            return {"summary": "Recorded sourced result"}
        self.stage_run("02", "Mock literature", model_step, model=True)
        def review():
            answer = self.ask("review", "Approve this concrete pilot?", [{"label": "Approve", "value": "approve"}])
            return {"summary": "Answer recorded", "answer": answer}
        self.stage_run("03", "Review", review)
    monkeypatch.setattr(ResearchWorkflow, "run", run)
    assert execute(tmp_path, pid, jid)
    snap = service.snapshot(pid)
    assert snap["project"]["state"] == "blocked" and snap["steps"][0]["status"] == "waiting"
    first_question = snap["data"]["research"]["questions"][0]
    with pytest.raises(ValueError, match="target accession"):
        service.answer_research(pid, first_question["id"], "missing species", spawn=False)
    assert not calls
    assert service.answer_research(pid, first_question["id"], "P00520; Mus musculus", spawn=False) == jid
    assert service.snapshot(pid)["jobs"][0]["status"] == "queued"
    assert execute(tmp_path, pid, jid)
    snap = service.snapshot(pid)
    assert len(calls) == 1
    model_artifact = next(step for step in snap["steps"] if step["id"] == "02")["outputs"]["result_artifact"]
    pending = next(q for q in snap["data"]["research"]["questions"] if q["status"] == "pending")
    with pytest.raises(ValueError, match="Answer the question"):
        service.start_research(pid, spawn=False)
    assert service.answer_research(pid, pending["id"], "approve") == jid
    assert len(launches) == 1 and "biocentric.research_worker" in launches[0]
    assert execute(tmp_path, pid, jid)
    final = service.snapshot(pid)
    assert final["project"]["state"] == "completed" and final["jobs"][0]["attempts"] == 3
    assert all(q["status"] == "answered" for q in final["data"]["research"]["questions"])
    assert len(calls) == 1
    assert next(step for step in final["steps"] if step["id"] == "02")["outputs"]["result_artifact"] == model_artifact
    assert service.start_research(pid) == jid and len(launches) == 1
    with pytest.raises(ValueError, match="no longer awaiting"):
        service.answer_research(pid, pending["id"], "approve")


def test_global_model_budget_counts_failed_attempts_across_stages_and_worker_instances(tmp_path):
    service, pid, jid = project(tmp_path, max_model_calls=2)
    first = ResearchWorkflow(service.store, pid, jid)
    first.stage = "02"
    first.emit(kind="model_call", name="one", status="running")
    first.emit(kind="model_call", name="one", status="failed")
    second = ResearchWorkflow(service.store, pid, jid)
    second.stage = "05"
    second.emit(kind="model_call", name="two", status="running")
    second.emit(kind="model_call", name="two", status="completed")
    second.stage = "07"
    with pytest.raises(ValueError, match="budget is exhausted"):
        second.emit(kind="model_call", name="never sent", status="running")
    assert sum(event.get("kind") == "model_call" and event.get("status") == "running" for event in service.store.events(pid)) == 2


def test_completed_stage_requires_intact_registered_artifact_and_does_not_rerun(tmp_path):
    service, pid, jid = project(tmp_path)
    workflow = ResearchWorkflow(service.store, pid, jid)
    result = workflow.stage_run("02", "Saved result", lambda: {"summary": "source result"})
    def must_not_run():
        raise AssertionError("Completed step was rerun")
    later = ResearchWorkflow(service.store, pid, jid)
    assert later.stage_run("02", "Saved result", must_not_run)["result_artifact"] == result["result_artifact"]
    artifact = next(row for row in service.store.artifacts(pid) if row["id"] == result["result_artifact"])
    Path(artifact["path"]).write_text("corrupt output")
    with pytest.raises(ValueError, match="integrity"):
        later.stage_run("02", "Saved result", must_not_run)


def test_failed_model_step_requires_inline_retry_answer_before_another_attempt(tmp_path):
    service, pid, jid = project(tmp_path)
    workflow = ResearchWorkflow(service.store, pid, jid)
    def failed_attempt():
        workflow.emit(kind="model_call", name="Fake failed request", status="running")
        raise ValueError("Fake incomplete response")
    with pytest.raises(ValueError, match="incomplete"):
        workflow.stage_run("02", "Literature", failed_attempt, model=True)
    calls = []
    later = ResearchWorkflow(service.store, pid, jid)
    def successful_attempt():
        calls.append(1)
        later.emit(kind="model_call", name="Fake retry", status="running")
        return {"summary": "Second attempt completed"}
    with pytest.raises(AwaitingAnswer):
        later.stage_run("02", "Literature", successful_attempt, model=True)
    assert not calls
    question = service.snapshot(pid)["data"]["research"]["questions"][0]
    assert question["id"].startswith("retry_02")
    service.store.finish(jid, "blocked", project_id=pid)
    service.answer_research(pid, question["id"], "retry", spawn=False)
    result = later.stage_run("02", "Literature", successful_attempt, model=True)
    assert result["summary"] == "Second attempt completed" and calls == [1]
    assert sum(event.get("kind") == "model_call" and event.get("status") == "running" for event in service.store.events(pid)) == 2


@pytest.mark.parametrize("rmsd,physical,failures", [(2.01, True, []), (float("nan"), True, []), (-1, True, []), (True, True, []), (1.0, False, []), (1.0, True, [{"compound_id": "missing"}])])
def test_pilot_gate_requires_finite_pose_recovery_physics_and_no_missing_results(rmsd, physical, failures):
    pilot = {"poses": [{"compound_id": "reference", "pose_number": 1, "redocking": {"rmsd_angstrom": rmsd}, "validation": {"passed": physical}}], "failures": failures}
    result = pilot_checks(pilot, {"maximum": 2.0, "require_physical_validity": True})
    assert result["status"] == "pilot_failed"
    assert result["affinity_validation"] == result["predictive_validation"] == "not_established"


def test_frozen_criteria_reject_mutation_before_protocol_computation(tmp_path):
    service, pid, jid = project(tmp_path)
    frozen = build_testset(service.snapshot(pid))
    frozen["criteria"]["pose_recovery"]["maximum"] = 100
    service.store.research_update(pid, testset=frozen, tool_plan={"docking_engine": "vina"})
    workflow = ResearchWorkflow(service.store, pid, jid)
    with pytest.raises(ValueError, match="hash mismatch"):
        workflow.protocol()


def test_failed_pilot_blocks_production_despite_advisor_recommendation(tmp_path, monkeypatch):
    service, pid, jid = project(tmp_path)
    service.store.research_update(pid, verification={"status": "pilot_failed", "model_review": {"recommendation": "proceed_exploratory"}})
    workflow = ResearchWorkflow(service.store, pid, jid)
    def no_science():
        raise AssertionError("Scientific execution must not begin after failed pilot")
    monkeypatch.setattr(workflow, "legacy", no_science)
    with pytest.raises(AwaitingAnswer):
        workflow.production()
    question = service.snapshot(pid)["data"]["research"]["questions"][0]
    assert question["options"] == [{"label": "Keep results and stop", "value": "stop"}]
    service.store.finish(jid, "blocked", project_id=pid)
    service.answer_research(pid, question["id"], "stop", spawn=False)
    with pytest.raises(Cancelled, match="Pilot checks failed"):
        workflow.production()


def test_original_chembl_retrieval_shape_builds_compatible_frozen_cohort(tmp_path, monkeypatch):
    from biocentric.retrieval import RetrievalClient
    service, pid, jid = project(tmp_path)
    service.store.research_update(pid, target={"target": "P00520", "organism": "Mus musculus"})
    assays = []
    for number, smiles in enumerate(["c1ccccc1", "CCO"], start=1):
        assays.append(RetrievalClient.parse_chembl_activity({"activity_id": number, "molecule_chembl_id": f"CHEMBL{number}",
            "canonical_smiles": smiles, "target_chembl_id": "CHEMBL_TARGET", "assay_chembl_id": "CHEMBL_ASSAY",
            "type": "Kd", "value": "10", "units": "nM", "relation": "=", "document_chembl_id": "CHEMBL_DOC"},
            {"assay_organism": "Mus musculus", "assay_type": "B", "description": "Matched assay"}, "2026-09-26"))
    monkeypatch.setattr(RetrievalClient, "resolve_target", lambda *args: {"candidates": [{"accession": "P00520", "name": "ABL1", "organism": "Mus musculus"}]})
    monkeypatch.setattr(RetrievalClient, "dossier", lambda *args: {"evidence": [], "assays": deepcopy(assays), "compounds": [], "limitations": []})
    result = ResearchWorkflow(service.store, pid, jid).evidence()
    assert result["assay_count"] == 2
    frozen = build_testset(service.snapshot(pid), organism="Mus musculus")
    assert len(frozen["records"]) == 2 and frozen["tuning_ids"] and frozen["evaluation_ids"]
    assert {row["endpoint"] for row in frozen["records"]} == {"Kd"}
    assert all(row["original"]["raw"]["value"] == "10" for row in frozen["records"])
    assert validate_testset(frozen)


MODEL_QUESTION = {"id": "scope", "question": "Continue exploratory or provide matched assay data?",
                  "options": ["Continue as exploratory", "Provide matched assay data"]}


def frozen_model_question(tmp_path, *, max_model_calls=6, legacy=False):
    service, pid, jid = project(tmp_path, max_model_calls=max_model_calls)
    service.store.finish(jid, "blocked", project_id=pid)
    service.store.research_update(pid, target={"target": "P00520", "organism": "Mus musculus"})
    workflow = ResearchWorkflow(service.store, pid, jid)
    frozen = build_testset(service.snapshot(pid), organism="Mus musculus")
    workflow.stage_run("04", "Freeze test set", lambda: frozen)
    workflow.stage_run("05", "Propose tools", lambda: {"docking_engine": "vina", "questions": [MODEL_QUESTION]})
    service.store.research_update(pid, testset=frozen, tool_plan={"docking_engine": "vina", "questions": [MODEL_QUESTION]})
    workflow.stage = "06"
    with pytest.raises(AwaitingAnswer):
        if legacy:
            workflow.ask("tool_question_0", MODEL_QUESTION["question"], MODEL_QUESTION["options"],
                         kind="model_clarification", context={"question": MODEL_QUESTION})
        else:
            workflow.model_questions([MODEL_QUESTION], "tool_question")
    return service, pid, jid, workflow, frozen


def test_legacy_model_question_snapshot_uses_coded_controls_without_mutating_database(tmp_path):
    service, pid, jid, workflow, frozen = frozen_model_question(tmp_path, legacy=True)
    visible = service.snapshot(pid)["data"]["research"]["questions"][0]
    stored = service.store.project(pid)["data"]["research"]["questions"][0]
    assert stored["kind"] == "model_clarification" and visible["kind"] == "model_clarification_control"
    assert visible["id"] == stored["id"]
    assert [choice["value"] for choice in visible["options"]] == ["continue_frozen", "revise_inputs"]
    assert visible["context"]["original_model_question"]["options"] == MODEL_QUESTION["options"]
    with pytest.raises(ValueError, match="displayed answers"):
        service.answer_research(pid, visible["id"], "Provide matched assay data", spawn=False)
    service.answer_research(pid, visible["id"], "continue_frozen", spawn=False)
    assert workflow.model_questions([MODEL_QUESTION], "tool_question")[0]["answer"] == "continue_frozen"
    assert len(service.snapshot(pid)["data"]["research"]["questions"]) == 1


def test_revision_waits_for_changed_inputs_and_explicit_paid_rebuild_then_requeues_same_job(tmp_path, monkeypatch):
    service, pid, jid, workflow, frozen = frozen_model_question(tmp_path, legacy=True)
    before = service.snapshot(pid)
    old_artifact = next(step for step in before["steps"] if step["id"] == "04")["outputs"]["result_artifact"]
    question = before["data"]["research"]["questions"][0]
    def never_spawn(*args, **kwargs):
        raise AssertionError("Revision must remain pending without a model request or subprocess")
    monkeypatch.setattr("biocentric.service.subprocess.Popen", never_spawn)
    assert service.answer_research(pid, question["id"], "revise_inputs") == jid
    pending = next(q for q in service.snapshot(pid)["data"]["research"]["questions"] if q["status"] == "pending")
    assert pending["kind"] == "research_input_revision" and pending["context"]["remaining_model_calls"] == 6
    assert "1 GPT call" in pending["options"][0]["label"]
    assert service.snapshot(pid)["jobs"][0]["status"] == "blocked"
    with pytest.raises(ValueError, match="Upload at least one changed"):
        service.answer_research(pid, pending["id"], "rebuild", spawn=False)
    service.upload(pid, "revised-library.smi", b"CCO ethanol\n", "library")
    csv = b"compound_id,endpoint,value,units,relation,assay_id,organism,smiles\na,Kd,10,nM,=,assay1,Mus musculus,CCO\n"
    service.import_measurements(pid, csv)
    waiting = service.snapshot(pid)
    assert waiting["data"]["research"]["testset"] == frozen
    assert len(waiting["data"]["research"]["input_revision"]["changed_inputs"]) == 2
    assert not any(event.get("kind") == "model_call" for event in waiting["events"])
    assert service.answer_research(pid, pending["id"], "rebuild", spawn=False) == jid
    rebuilding = service.snapshot(pid)
    assert rebuilding["jobs"][0]["status"] == "queued"
    assert rebuilding["data"]["research"]["input_revision_version"] == 1
    assert all(step["status"] == "pending" for step in rebuilding["steps"] if step["id"] >= "04")
    assert "testset" not in rebuilding["data"]["research"]
    assert service.store.read_artifact(pid, {"id": old_artifact})
    assert rebuilding["data"]["research"]["revision_history"][0]["testset"] == frozen
    calls = []
    def run(self):
        def build():
            manifest = build_testset(service.snapshot(pid), organism="Mus musculus")
            self.store.research_update(pid, testset=manifest)
            return manifest
        self.stage_run("04", "Rebuild test set", build)
        def propose():
            self.emit(kind="model_call", name="Fake revised proposal", status="running")
            calls.append(1)
            plan = {"docking_engine": "vina", "questions": [MODEL_QUESTION]}
            self.store.research_update(pid, tool_plan=plan)
            return plan
        self.stage_run("05", "Rebuild proposal", propose, model=True)
        self.stage_run("06", "Review revised inputs", lambda: {"answers": self.model_questions([MODEL_QUESTION], "tool_question")})
    monkeypatch.setattr(ResearchWorkflow, "run", run)
    assert execute(tmp_path, pid, jid)
    after = service.snapshot(pid)
    assert calls == [1] and after["jobs"][0]["status"] == "blocked"
    new_question = next(q for q in after["data"]["research"]["questions"] if q["status"] == "pending")
    assert new_question["id"] != question["id"] and new_question["context"]["input_revision_version"] == 1
    assert new_question["context"]["testset_sha256"] == after["data"]["research"]["testset"]["manifest_sha256"]
    assert after["data"]["research"]["testset"]["dataset_sha256"] != frozen["dataset_sha256"]


def test_revision_cannot_reset_budget_or_change_completed_pilot_inputs(tmp_path):
    service, pid, jid, workflow, frozen = frozen_model_question(tmp_path, max_model_calls=1)
    workflow.emit(kind="model_call", name="Previous paid attempt", status="running")
    question = service.snapshot(pid)["data"]["research"]["questions"][0]
    service.answer_research(pid, question["id"], "revise_inputs", spawn=False)
    service.upload(pid, "new.smi", b"CCC propane\n", "library")
    pending = next(q for q in service.snapshot(pid)["data"]["research"]["questions"] if q["status"] == "pending")
    assert pending["context"]["remaining_model_calls"] == 0
    with pytest.raises(ValueError, match="No model calls remain"):
        service.answer_research(pid, pending["id"], "rebuild", spawn=False)
    assert service.snapshot(pid)["data"]["research"]["testset"] == frozen
    assert next(q for q in service.snapshot(pid)["data"]["research"]["questions"] if q["id"] == pending["id"])["status"] == "pending"
    service.store.research_update(pid, pilot={"status": "completed", "poses": [{"compound_id": "reference"}]})
    with pytest.raises(ValueError, match="frozen"):
        service.upload(pid, "another.smi", b"CCCC butane\n", "library")
    with pytest.raises(ValueError, match="Completed pilot"):
        service.answer_research(pid, pending["id"], "rebuild", spawn=False)


def test_negative_legacy_answer_never_implicitly_approves_a_model_question(tmp_path):
    service, pid, jid, workflow, frozen = frozen_model_question(tmp_path, legacy=True)
    research = service.store.project(pid)["data"]["research"]
    research["questions"][0].update(status="answered", answer="Provide matched assay data")
    service.store.research_update(pid, questions=research["questions"])
    with pytest.raises(AwaitingAnswer):
        workflow.model_questions([MODEL_QUESTION], "tool_question")
    pending = next(q for q in service.snapshot(pid)["data"]["research"]["questions"] if q["status"] == "pending")
    assert "confirm" in pending["id"] and pending["context"]["legacy_answer"] == "Provide matched assay data"
    assert pending["kind"] == "model_clarification_control"


def test_historical_scientific_reupload_is_rejected_before_any_mutation(tmp_path):
    service, pid, jid = project(tmp_path)
    service.store.finish(jid, "blocked", project_id=pid)
    first = service.upload(pid, "library.smi", b"CCO first\n", "library")
    current = service.upload(pid, "library.smi", b"CCN second\n", "library")
    before = service.snapshot(pid)
    with pytest.raises(ValueError, match="new filename"):
        service.upload(pid, "library.smi", b"CCO first\n", "library")
    assert service.snapshot(pid) == before
    assert ResearchWorkflow(service.store, pid, jid).artifact("library") == Path(current["path"])
    unchanged = service.upload(pid, "same-current-renamed.smi", b"CCN second\n", "library")
    assert unchanged["unchanged"] and unchanged["id"] == current["id"]
    assert service.snapshot(pid) == before
    restored = service.upload(pid, "restored-first.smi", b"CCO first\n", "library")
    assert restored["id"] not in {first["id"], current["id"]}
    assert ResearchWorkflow(service.store, pid, jid).artifact("library") == Path(restored["path"])


def test_historical_csv_reupload_cannot_desynchronize_measurements_and_artifacts(tmp_path):
    service, pid, jid = project(tmp_path)
    service.store.finish(jid, "blocked", project_id=pid)
    header = b"compound_id,endpoint,value,units,relation,assay_id,organism,smiles\n"
    first = header + b"a,Kd,10,nM,=,assay1,Mus musculus,CCO\n"
    second = header + b"a,Kd,20,nM,=,assay1,Mus musculus,CCO\n"
    service.import_measurements(pid, first)
    current = service.import_measurements(pid, second)
    before = service.snapshot(pid)
    with pytest.raises(ValueError, match="older assay CSV"):
        service.import_measurements(pid, first)
    assert service.snapshot(pid) == before
    assert before["data"]["measurements"][0]["value"] == 20
    assert service.import_measurements(pid, second)["unchanged"]
    assert service.snapshot(pid) == before
    assert next(row for row in reversed(before["artifacts"]) if row["kind"] == "measured_results")["id"] == current["artifact"]["id"]


def test_noop_reupload_does_not_unlock_revision_or_invalidate_frozen_inputs(tmp_path):
    service, pid, jid = project(tmp_path)
    service.store.finish(jid, "blocked", project_id=pid)
    current = service.upload(pid, "library.smi", b"CCO ethanol\n", "library")
    workflow = ResearchWorkflow(service.store, pid, jid)
    frozen = build_testset(service.snapshot(pid))
    workflow.stage_run("04", "Frozen test set", lambda: frozen)
    service.store.research_update(pid, testset=frozen)
    workflow.stage = "06"
    with pytest.raises(AwaitingAnswer):
        workflow.model_questions([MODEL_QUESTION], "tool_question")
    question = service.snapshot(pid)["data"]["research"]["questions"][0]
    service.answer_research(pid, question["id"], "revise_inputs", spawn=False)
    before = service.snapshot(pid)
    assert service.upload(pid, "library.smi", b"CCO ethanol\n", "library")["id"] == current["id"]
    assert service.snapshot(pid) == before
    pending = next(q for q in before["data"]["research"]["questions"] if q["status"] == "pending")
    with pytest.raises(ValueError, match="Upload at least one changed"):
        service.answer_research(pid, pending["id"], "rebuild", spawn=False)


def test_comparative_pilots_resume_then_lock_parameters_before_review(tmp_path, monkeypatch):
    """Real durable artifacts; synthetic adapter output only, no scientific claim."""
    from biocentric.worker import Pipeline
    from biocentric.store import digest, dump
    service = Service(tmp_path)
    pid = service.create_research_project('Synthetic comparison invariant test', use_public_example=True)
    jid = service.start_research(pid, spawn=False)
    workflow = ResearchWorkflow(service.store, pid, jid)
    testset = build_testset(service.snapshot(pid), organism='Mus musculus')
    plan = {'docking_engine':'vina','affinity_engine':None,'questions':[],'candidates':[]}
    service.store.research_update(pid,testset=testset,tool_plan=plan,capabilities=[{'id':'vina','available':True}])
    calls = []
    def prepare(pipeline):
        pipeline.update(preparation={'compounds':[{'id':'fixture'}]})
        return {'compounds':[{'id':'fixture'}]}
    def pilot(pipeline):
        protocol = pipeline.data['protocol']
        calls.append((protocol['exhaustiveness'],str(pipeline.work)))
        return {'poses':[{'compound_id':'fixture_reference','pose_number':1,
                         'redocking':{'rmsd_angstrom':1.0 if protocol['exhaustiveness']==4 else .5},
                         'validation':{'passed':True}}],'failures':[]}
    monkeypatch.setattr(Pipeline,'preparation',prepare)
    monkeypatch.setattr(Pipeline,'pilot',pilot)
    with pytest.raises(AwaitingAnswer):
        workflow.stage_run('06','Comparison fixture',workflow.pilot)
    assert not calls and 'comparison' not in workflow.research
    questions=workflow.research['questions']
    questions[-1].update(status='answered',answer='approve')
    service.store.research_update(pid,questions=questions)
    first=workflow.pilot()
    assert len(calls)==2 and calls[0][1]!=calls[1][1]
    assert first['comparison']['criteria']['mode']=='pose_software'
    assert 'selection_lock' not in workflow.research
    workflow.pilot()
    assert len(calls)==2
    locked=workflow.lock_comparison()
    assert locked['selection_lock']['selected_protocol_id']=='vina-exhaustiveness-8'
    assert workflow.data['protocol']['exhaustiveness']==8
    assert locked['heldout_evaluation']['status']=='not_evaluated'
    assert workflow.lock_comparison()==locked
    service.store.research_update(pid,verification={'status':'software_pilot_passed','model_review':{},'limitations':[]})
    changed={**workflow.data['protocol'],'exhaustiveness':16}
    service.store.update(pid,data={'protocol':changed})
    with pytest.raises(ValueError,match='locked protocol'):
        workflow.production()


def test_comparative_candidates_freeze_native_parameters_and_keep_unavailable_visible(tmp_path):
    service,pid,jid=project(tmp_path)
    service.store.research_update(pid,tool_plan={'docking_engine':'vina','candidates':[
        {'engine':'gnina','role':'docking'},{'engine':'unidock','role':'docking'},{'engine':'autodock4','role':'docking'}]},
        capabilities=[{'id':engine,'available':engine!='autodock4','blockers':['Missing AutoGrid maps'] if engine=='autodock4' else []}
                      for engine in ('vina','gnina','unidock','autodock4')])
    settings={'engine':'vina','exhaustiveness':4,'seed':42,'cpu':1,'n_poses':3,'timeout_seconds':120,
              'box':{'center':[0,0,0],'size':[20,20,20]},'input_hashes':{'receptor':'a'*64}}
    candidates,excluded=ResearchWorkflow(service.store,pid,jid).comparison_candidates(settings)
    assert len(candidates)==4
    assert {row['engine_id'] for row in candidates}=={'vina','gnina','unidock'}
    gnina=next(row for row in candidates if row['engine_id']=='gnina')
    assert gnina['parameters']['native_parameters']['device']=='gpu'
    assert gnina['parameters']['native_parameters']['exhaustiveness']==4
    assert excluded==[{'engine_id':'autodock4','reason':'Missing AutoGrid maps'}]


def test_child_pipelines_share_cumulative_attempt_budget_instead_of_resetting_each_pilot(tmp_path, monkeypatch):
    service,pid,jid=project(tmp_path)
    workflow=ResearchWorkflow(service.store,pid,jid)
    workflow.started=100.0
    workflow.spent=15.0
    monkeypatch.setattr('biocentric.research_workflow.time.monotonic',lambda:100.0+workflow.spec['compute_budget_seconds']-10)
    # Time spent in previous pilots plus persisted prior attempts exhausts the
    # common budget even though this child has just been constructed.
    child=workflow.legacy()
    assert child.started==workflow.started and child.spent==workflow.spent
    assert child.cancelled()
    from biocentric.worker import Gate
    with pytest.raises(Gate,match='compute time exhausted'):
        child.check()
    with pytest.raises(ValueError,match='budget is exhausted'):
        workflow.check()


def test_frozen_candidate_adapter_code_changes_block_execution(tmp_path):
    from biocentric.store import digest
    service,pid,jid=project(tmp_path)
    code=tmp_path/'scientific-adapter.py';code.write_text('reviewed implementation')
    candidate={'provenance':{'adapter_code_sha256':{str(code):digest(code.read_bytes())}}}
    workflow=ResearchWorkflow(service.store,pid,jid)
    workflow.check_candidate_runtime(candidate)
    code.write_text('changed implementation')
    with pytest.raises(ValueError,match='adapter changed'):
        workflow.check_candidate_runtime(candidate)
