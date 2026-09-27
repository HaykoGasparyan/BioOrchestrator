import json
import hashlib
from types import SimpleNamespace

import pytest
import httpx

from biocentric.research_agent import research_literature, recommend_tools, review_pilot, ResearchError, _fetch_source_excerpt, _enrich_sources, Selection, PilotReview
from biocentric.testsets import build_testset

SPEC = {"local_only": False, "allow_model": True, "model": "gpt-4.1-mini", "max_model_calls": 3,
        "max_model_output_tokens": 1500, "organism": "Mus musculus", "target": "ABL1"}
URL = "https://autodock-vina.readthedocs.io/en/latest/"
HTTP_CLIENT = httpx.Client


@pytest.fixture(autouse=True)
def fake_source_fetching(monkeypatch):
    def factory(**kwargs):
        return HTTP_CLIENT(**kwargs, transport=httpx.MockTransport(lambda request: httpx.Response(200,
            headers={"content-type": "text/html"}, text="<article><h1>Vina documentation</h1><p>" + "Documented method and input requirements. " * 8 + "</p></article>")))
    monkeypatch.setattr("biocentric.research_agent.httpx.Client", factory)


def fake_client(output, text=None, events=None, status="completed"):
    requests = []
    response = {"id": "resp_fake", "status": status, "usage": {"input_tokens": 10, "output_tokens": 20}, "output": output}
    if text is not None:
        response["output_text"] = text
    def create(**kwargs):
        requests.append(kwargs)
        return iter([*(events or []), {"type": "response.completed", "response": response}])
    return SimpleNamespace(responses=SimpleNamespace(create=create)), requests


def search_output(url=URL):
    return [{"id": "search1", "type": "web_search_call", "status": "completed", "action": {"type": "search", "queries": ["Vina target validation"], "sources": [{"type": "url", "url": url}]}},
            {"type": "message", "content": [{"type": "output_text", "text": "Search supports a bounded pilot.", "annotations": [{"type": "url_citation", "url": url, "title": "Vina docs", "start_index": 0, "end_index": 20}]}]}]


def literature():
    client, _ = fake_client(search_output())
    return research_literature(SPEC, [], lambda **event: None, client=client)


def test_runtime_search_is_required_bounded_and_streams_public_outputs_only():
    events = []
    client, requests = fake_client(search_output(), events=[
        {"type": "response.web_search_call.searching", "item_id": "search1"},
        {"type": "response.output_text.delta", "delta": "Reading source descriptions.\n"},
        {"type": "response.reasoning_text.delta", "delta": "PRIVATE_MUST_NOT_PERSIST"}])
    result = research_literature(SPEC, [{"name": "Vina"}], lambda **event: events.append(event), client=client)
    request = requests[0]
    assert request["tool_choice"] == "required" and request["max_tool_calls"] == 4
    assert request["tools"][0]["type"] == "web_search_preview"
    assert "filters" not in request["tools"][0]
    assert request["include"] == ["web_search_call.action.sources"]
    assert request["store"] is False and request["stream"] is True
    assert result["model_calls"] == result["search_calls"] == 1
    assert result["records"][0]["url"] == URL
    assert result["records"][0]["content_access"] == "page_text_excerpt"
    assert result["records"][0]["fetch_status"] == "completed"
    assert "Documented method" in result["records"][0]["source_excerpt"]
    assert result["records"][0]["full_text_retrieved"] is False
    assert any(event["kind"] == "model_delta" for event in events)
    assert "PRIVATE_MUST_NOT_PERSIST" not in json.dumps(events)
    assert "1–2 short sentences" in json.loads(request["input"][0]["content"])["task"]


@pytest.mark.parametrize("model", ["gpt-4.1", "gpt-4.1-mini", "gpt-4.1-mini-2025-04-14"])
def test_gpt41_compatibility_route_keeps_required_search_without_unsupported_filters(model):
    client, requests = fake_client(search_output())
    result = research_literature({**SPEC, "model": model}, [], lambda **event: None, client=client)
    assert result["status"] == "completed"
    assert requests[0]["tools"] == [{"type": "web_search_preview", "search_context_size": "medium"}]
    assert requests[0]["tool_choice"] == "required" and requests[0]["model"] == model


def test_modern_search_route_includes_hosted_domain_filter():
    client, requests = fake_client(search_output())
    research_literature({**SPEC, "model": "configured-modern-model"}, [], lambda **event: None, client=client)
    tool = requests[0]["tools"][0]
    assert tool["type"] == "web_search" and "pubmed.ncbi.nlm.nih.gov" in tool["filters"]["allowed_domains"]
    assert requests[0]["tool_choice"] == "required"


@pytest.mark.parametrize("spec", [{**SPEC, "local_only": True}, {**SPEC, "allow_model": False}, {**SPEC, "max_model_calls": 0}])
def test_privacy_and_zero_budget_stop_before_any_request(spec):
    client, requests = fake_client(search_output())
    with pytest.raises(ResearchError):
        research_literature(spec, [], lambda **event: None, client=client)
    assert not requests


@pytest.mark.parametrize("output", [[], search_output("https://example.org/invented"), search_output("https://github.com/unknown/commands")])
def test_absent_search_or_disallowed_sources_cannot_be_successful_research(output):
    client, _ = fake_client(output)
    with pytest.raises(ResearchError, match="No completed search"):
        research_literature(SPEC, [], lambda **event: None, client=client)


def proposal(lit, **extra):
    cite = lit["records"][0]["id"]
    return {"docking_engine": "vina", "affinity_engine": None, "reason": "Vina is locally available for a pose-only pilot.",
            "citations": [cite], "candidates": [{"engine": "vina", "role": "docking", "reason": "CPU availability", "citations": [cite]}],
            "testset_fit": "One reference tests pose recovery only.", "uncertainties": ["No affinity benchmark"], "questions": [], **extra}


def test_selection_requires_frozen_real_testset_and_rejects_unknown_or_unavailable_tools():
    lit = literature()
    testset = build_testset({"data": {}, "artifacts": []})
    capabilities = [{"id": "vina", "available": True}, {"id": "gnina", "available": False}]
    for changed in [{"docking_engine": "gnina"}, {"docking_engine": "shell"}, {"citations": ["invented"]}]:
        client, _ = fake_client([], json.dumps(proposal(lit, **changed)))
        with pytest.raises(ResearchError):
            recommend_tools(SPEC, lit, testset, capabilities, lambda **event: None, client=client)
    client, requests = fake_client([], json.dumps(proposal(lit)))
    result = recommend_tools(SPEC, lit, testset, capabilities, lambda **event: None, client=client)
    assert result["docking_engine"] == "vina" and result["affinity_engine"] is None
    assert result["criteria"] == testset["criteria"]
    assert requests[0]["tools"] == [] and requests[0]["tool_choice"] == "none"
    testset["criteria"]["pose_recovery"]["maximum"] = 99
    with pytest.raises(ValueError, match="hash mismatch"):
        recommend_tools(SPEC, lit, testset, capabilities, lambda **event: None, client=client)
    assert len(requests) == 1


def test_global_emit_budget_rejection_precedes_send_and_cancellation_stops_work():
    client, requests = fake_client(search_output())
    def reject(**event):
        if event["kind"] == "model_call" and event["status"] == "running":
            raise RuntimeError("global project budget exhausted")
    with pytest.raises(RuntimeError, match="global project budget"):
        research_literature(SPEC, [], reject, client=client)
    with pytest.raises(ResearchError, match="cancelled"):
        research_literature(SPEC, [], lambda **event: None, cancelled=lambda: True, client=client)
    assert not requests


def test_affinity_method_cannot_be_selected_as_docking_even_if_installed():
    lit = literature()
    testset = build_testset({"data": {}, "artifacts": []})
    client, _ = fake_client([], json.dumps(proposal(lit, docking_engine="mmpbsa")))
    with pytest.raises(ResearchError, match="docking role"):
        recommend_tools(SPEC, lit, testset, {"vina": {"available": True}, "mmpbsa": {"available": True}}, lambda **event: None, client=client)


def test_provider_error_body_is_not_persisted():
    events = []
    def create(**kwargs):
        raise RuntimeError("DO_NOT_PERSIST_SECRET_BODY")
    client = SimpleNamespace(responses=SimpleNamespace(create=create))
    with pytest.raises(ResearchError, match="request failed"):
        research_literature(SPEC, [], lambda **event: events.append(event), client=client)
    assert "DO_NOT_PERSIST_SECRET_BODY" not in json.dumps(events)


@pytest.mark.parametrize("status,param", [(400, "tools"), (400, "input"), (401, "tools")])
def test_only_tool_compatibility_error_gets_bounded_redacted_diagnostics(status, param):
    events = []
    synthetic_key = "sk-proj-" + "SYNTHETIC_TEST_TOKEN_" * 8
    class FakeAPIError(Exception):
        status_code = status
        code = "unsupported_value"
    error = FakeAPIError("This raw exception must never be persisted")
    error.param = param
    error.body = {"message": "Tool unsupported. " + synthetic_key + ". " + "x" * 700}
    def create(**kwargs):
        raise error
    client = SimpleNamespace(responses=SimpleNamespace(create=create))
    with pytest.raises(ResearchError, match=f"parameter: {param}"):
        research_literature(SPEC, [], lambda **event: events.append(event), client=client)
    failed = next(event for event in events if event["status"] == "failed")
    assert failed["outputs"]["code"] == "unsupported_value" and failed["outputs"]["param"] == param
    persisted = json.dumps(events)
    assert synthetic_key not in persisted and "SYNTHETIC_TEST_TOKEN" not in persisted
    assert "raw exception" not in persisted
    if (status, param) == (400, "tools"):
        detail = failed["outputs"]["tool_compatibility_detail"]
        assert "[redacted]" in detail and len(detail) <= 500
    else:
        assert "tool_compatibility_detail" not in failed["outputs"]


def test_diagnostic_redacts_credentials_before_truncation_boundary():
    events = []
    class FakeAPIError(Exception):
        status_code = 400
        code = "unsupported_value"
        param = "tools"
        body = {"error": {"message": "x" * 490 + "sk-proj-" + "SYNTHETIC_TEST_TOKEN_" * 8}}
    def create(**kwargs):
        raise FakeAPIError("Do not persist")
    client = SimpleNamespace(responses=SimpleNamespace(create=create))
    with pytest.raises(ResearchError):
        research_literature(SPEC, [], lambda **event: events.append(event), client=client)
    failed = next(event for event in events if event["status"] == "failed")
    assert "sk-" not in failed["outputs"]["tool_compatibility_detail"]
    assert len(failed["outputs"]["tool_compatibility_detail"]) <= 500


def test_review_cannot_return_computed_metrics_or_predictive_success_status():
    lit = literature()
    testset = build_testset({"data": {}, "artifacts": []})
    plan = {"criteria_sha256": testset["criteria_sha256"], "dataset_sha256": testset["dataset_sha256"]}
    payload = {"summary": "Pose recovery passed, but predictive validation is absent.", "recommendation": "proceed_exploratory",
               "citations": [lit["records"][0]["id"]], "limitations": ["No held-out labeled data"], "questions": []}
    client, requests = fake_client([], json.dumps(payload))
    result = review_pilot(SPEC, lit, testset, plan, {"rmsd_angstrom": 1.2}, lambda **event: None, client=client)
    assert result["status"] == "advisory" and "not_established" in result["predictive_validation"]
    assert "rmsd_angstrom" not in result
    assert requests[0]["tools"] == []
    client, _ = fake_client([], json.dumps({**payload, "computed_metrics": {"rmsd": 0}}))
    with pytest.raises(ResearchError, match="structured-output"):
        review_pilot(SPEC, lit, testset, plan, {}, lambda **event: None, client=client)


def test_source_fetch_checks_redirects_before_sending_to_another_host():
    sent = []
    def redirect(request):
        sent.append(str(request.url))
        return httpx.Response(302, headers={"location": "http://127.0.0.1:8502/secrets"})
    with HTTP_CLIENT(transport=httpx.MockTransport(redirect)) as client:
        with pytest.raises(ResearchError, match="allowlist"):
            _fetch_source_excerpt(URL, client)
    assert sent == [URL]


def test_source_excerpt_has_byte_and_text_caps_and_discards_script_markup():
    html = "<html><script>UNTRUSTED_EXECUTABLE</script><nav>MENU</nav><article>" + "Published methods and caveats. " * 100 + "</article></html>"
    with HTTP_CLIENT(transport=httpx.MockTransport(lambda request: httpx.Response(200, headers={"content-type": "text/html"}, text=html))) as client:
        result = _fetch_source_excerpt(URL, client, max_bytes=1000, max_characters=150)
    assert result["fetched_bytes"] == 1000 and result["excerpt_truncated"]
    assert len(result["source_excerpt"]) == 150 and not result["full_text_retrieved"]
    assert "UNTRUSTED_EXECUTABLE" not in result["source_excerpt"] and "MENU" not in result["source_excerpt"]


def test_inaccessible_sources_remain_honestly_search_only_and_fetch_budget_is_four():
    records = [{"id": f"web_{i}", "url": f"https://pubmed.ncbi.nlm.nih.gov/{i}/", "cited": True,
                "content_access": "search_excerpt", "full_text_retrieved": False} for i in range(6)]
    sent, events = [], []
    def denied(request):
        sent.append(request.url)
        return httpx.Response(403, headers={"content-type": "text/html"}, text="Forbidden")
    with HTTP_CLIENT(transport=httpx.MockTransport(denied)) as client:
        _enrich_sources(records, lambda **event: events.append(event), lambda: False, client)
    assert len(sent) == 4
    assert all(row["content_access"] == "search_excerpt" and not row["full_text_retrieved"] for row in records)
    assert [row["fetch_status"] for row in records] == ["unavailable"] * 4 + ["not_fetched_within_budget"] * 2
    assert len([event for event in events if event["status"] == "unavailable"]) == 4


@pytest.mark.parametrize("url", ["http://localhost/private", "https://pubmed.ncbi.nlm.nih.gov:8888/data", "https://pubmed.ncbi.nlm.nih.gov@localhost/data", "https://[malformed/data"])
def test_nonallowlisted_source_urls_never_issue_network_requests(url):
    def must_not_send(request):
        raise AssertionError("Disallowed URL was fetched")
    with HTTP_CLIENT(transport=httpx.MockTransport(must_not_send)) as client:
        with pytest.raises(ResearchError, match="allowlist"):
            _fetch_source_excerpt(url, client)


def test_brief_choice_prompt_uses_actual_reference_only_scope_and_keeps_full_details(tmp_path):
    reference = tmp_path / "observed.sdf"
    reference.write_text("registered observed-pose source")
    testset = build_testset({"data": {}, "artifacts": [{"id": "ref", "kind": "reference", "path": str(reference),
                               "sha256": hashlib.sha256(reference.read_bytes()).hexdigest()}]})
    lit = literature()
    public_reason = "The ABL literature and our reference-pose control support proposing Vina for local pose recovery. Affinity accuracy remains unvalidated."
    answer = proposal(lit, reason=public_reason, uncertainties=["No measured assay benchmark", "No multi-engine comparison", "Preparation requires review"])
    client, requests = fake_client([], json.dumps(answer))
    result = recommend_tools(SPEC, lit, testset, {"vina": {"available": True}, "gnina": {"available": False}}, lambda **event: None, client=client)
    payload = json.loads(requests[0]["input"][0]["content"])
    facts = payload["testset_summary_facts"]
    assert facts["measured_record_count"] == facts["held_out_compound_count"] == 0
    assert facts["reference_pose_control_count"] == 1
    assert payload["available_docking_engines"] == ["vina"]
    assert "sole available engine" in payload["task"] and "reference-only" in payload["task"]
    assert result["reason"] == public_reason and result["uncertainties"] == answer["uncertainties"]
    assert result["citations"] == answer["citations"] and result["candidates"] == answer["candidates"]


def test_brief_summary_guidance_preserves_existing_response_schema_compatibility():
    expected = {"docking_engine", "affinity_engine", "reason", "citations", "candidates", "testset_fit", "uncertainties", "questions"}
    schema = Selection.model_json_schema()
    assert set(schema["properties"]) == set(schema["required"]) == expected
    assert "1–2 short sentences" in schema["properties"]["reason"]["description"]
    assert "maxLength" not in schema["properties"]["reason"]
    legacy = {"docking_engine": "vina", "affinity_engine": None, "reason": "Previously saved detailed explanation. " * 30,
              "citations": ["saved_source"], "candidates": [], "testset_fit": "Pose control only", "uncertainties": [], "questions": []}
    assert Selection.model_validate(legacy).reason == legacy["reason"]
    pilot_schema = PilotReview.model_json_schema()
    assert set(pilot_schema["required"]) == {"summary", "recommendation", "citations", "limitations", "questions"}
    assert "comparative superiority" in pilot_schema["properties"]["summary"]["description"]


def test_pilot_summary_prompt_receives_real_counts_without_promoting_preview_to_validation():
    rows = [{"id": f"measurement{i}", "compound_id": f"compound{i}", "canonical_smiles": smiles,
             "endpoint": "Kd", "value": 10, "units": "nM", "relation": relation, "assay_id": "a", "organism": "Mus musculus"}
            for i, (smiles, relation) in enumerate([("CCO", "="), ("c1ccccc1", "="), ("C1CCCC1", ">")])]
    testset = build_testset({"data": {"assays": rows}}, organism="Mus musculus")
    lit = literature()
    plan = {"criteria_sha256": testset["criteria_sha256"], "dataset_sha256": testset["dataset_sha256"], "docking_engine": "vina"}
    response = {"summary": "The Vina pilot recovered the reference pose within the frozen threshold. The next library run remains exploratory because affinity accuracy was not tested.",
                "recommendation": "proceed_exploratory", "citations": [lit["records"][0]["id"]],
                "limitations": ["No comparative engine results", "No completed affinity evaluation"], "questions": []}
    client, requests = fake_client([], json.dumps(response))
    result = review_pilot(SPEC, lit, testset, plan, {"engine": "vina", "measured_rmsd_angstrom": 1.2}, lambda **event: None, client=client)
    payload = json.loads(requests[0]["input"][0]["content"])
    assert payload["testset_summary_facts"]["measured_record_count"] == 3
    assert payload["testset_summary_facts"]["exact_measurement_count"] == 2
    assert payload["testset_summary_facts"]["held_out_compound_count"] == 1
    assert "Not established" in payload["testset_summary_facts"]["predictive_validation"]
    assert "without measured results from multiple engines" in payload["task"]
    assert result["summary"] == response["summary"] and result["limitations"] == response["limitations"]


def test_method_proposal_and_pilot_advisors_never_receive_final_holdout_measurements():
    rows=[{'id':'r'+str(i),'compound_id':'c'+str(i),'smiles':smiles,'assay_id':'assay',
           'endpoint':'Kd','units':'nM','value':value,'relation':'=','organism':'Mus musculus',
           'raw_provider_table':{'hidden_marker':'DO_NOT_SEND_RAW_LABEL_TABLE'}}
          for i,(smiles,value) in enumerate([('CCO',1234.56),('c1ccccc1',876.54),('C1CCC1',321.98)])]
    testset=build_testset({'data':{'measurements':rows}},organism='Mus musculus')
    lit=literature()
    client,requests=fake_client([],json.dumps(proposal(lit)))
    plan=recommend_tools(SPEC,lit,testset,{'vina':{'available':True}},lambda **event:None,client=client)
    review={'summary':'Software scope only.','recommendation':'proceed_exploratory',
            'citations':[lit['records'][0]['id']],'limitations':['Synthetic fixture'],'questions':[]}
    review_client,review_requests=fake_client([],json.dumps(review))
    review_pilot(SPEC,lit,testset,plan,{'verification':{'status':'software_pilot_passed'}},lambda **event:None,client=review_client)
    for request in [requests[0],review_requests[0]]:
        payload=json.loads(request['input'][0]['content'])
        assert payload['frozen_testset']['manifest_sha256']==testset['manifest_sha256']
        assert 'DO_NOT_SEND_RAW_LABEL_TABLE' not in request['input'][0]['content']
        for row in payload['frozen_testset']['records']:
            if row['partition']=='evaluation':
                assert not set(row)&{'value','p_measurement','value_molar','label','original','source_ids','source_url'}
        assert payload['testset_summary_facts']['held_out_compound_count']==len(testset['evaluation_ids'])
