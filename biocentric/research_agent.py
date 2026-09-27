"""Bounded, auditable literature search and advisory selection; no job execution."""
from __future__ import annotations

import hashlib
import json
import os
from html.parser import HTMLParser
from pathlib import Path
import re
import time
from typing import Literal
from urllib.parse import urljoin, urlsplit, urlunsplit

import httpx

from pydantic import BaseModel, ConfigDict, Field

from .models import utc

PROMPT = Path(__file__).parent / "prompts" / "research-v2.md"
PRIMARY_DOMAINS = ["pubmed.ncbi.nlm.nih.gov", "pmc.ncbi.nlm.nih.gov", "europepmc.org",
                   "pubs.acs.org", "jcheminf.biomedcentral.com", "nature.com", "academic.oup.com",
                   "autodock-vina.readthedocs.io", "autodock.scripps.edu", "ambermd.org", "github.com"]
OFFICIAL_REPOSITORIES = ("/gnina/gnina", "/dptech-corp/uni-dock", "/computartcmcg/planet",
                         "/ccsb-scripps/autodock-vina", "/valdes-tresanco-ms/gmx_mmpbsa", "/dkoes/docs")
ENGINE_ROLES = {"vina": "docking", "autodock4": "docking", "gnina": "docking", "unidock": "docking",
                "mmpbsa": "affinity", "mmgbsa": "affinity", "gnina_affinity": "affinity", "planet": "affinity", "biocentric_affinity": "affinity"}


class ResearchError(ValueError):
    """No valid research/selection was produced; callers retain checkpoints."""


class Candidate(BaseModel):
    model_config = ConfigDict(extra="forbid")
    engine: str
    role: Literal["docking", "affinity"]
    reason: str = Field(description="Brief candidate-specific task fit, input requirements, or tradeoff supported by the cited sources. This detail belongs here rather than in the main public reason.")
    citations: list[str]


class Question(BaseModel):
    model_config = ConfigDict(extra="forbid")
    id: str
    question: str
    options: list[str]


class Selection(BaseModel):
    model_config = ConfigDict(extra="forbid")
    docking_engine: str | None
    affinity_engine: str | None
    reason: str = Field(description="Main public explanation: 1–2 short sentences, normally at most 60 words. Link the actual protein literature and constructed test-set scope to the proposed tool and specific reason for choosing it, then give the key caveat. Do not claim an assay benchmark when only reference-pose controls exist or a comparative win when only one engine is available.")
    citations: list[str]
    candidates: list[Candidate]
    testset_fit: str = Field(description="One short sentence naming the actual test-set scope and what it can verify, based on the supplied records and counts. A preview is not a completed benchmark.")
    uncertainties: list[str] = Field(description="Additional specific evidence gaps, preparation limitations, or untested alternatives; keep them out of the short main reason unless essential.")
    questions: list[Question]


class PilotReview(BaseModel):
    model_config = ConfigDict(extra="forbid")
    summary: str = Field(description="Main public explanation: 1–2 short sentences, normally at most 60 words. State the measured pilot outcome against the frozen criterion and the recommended next action, with the key limitation. A single-engine reference redock cannot establish comparative superiority or binding-affinity accuracy.")
    recommendation: Literal["proceed_exploratory", "revise", "need_data"]
    citations: list[str]
    limitations: list[str] = Field(description="Specific limitations retained for inspection, including missing comparative or independent affinity validation; do not expand the short public summary into a long list.")
    questions: list[Question]


def _dict(value):
    if isinstance(value, dict):
        return value
    if hasattr(value, "model_dump"):
        return value.model_dump(mode="json")
    return vars(value)


def _get(value, key, default=None):
    return value.get(key, default) if isinstance(value, dict) else getattr(value, key, default)


def _client(spec, client):
    if spec.get("local_only", True) or not spec.get("allow_model", False):
        raise ResearchError("External model and web search are disabled by project policy")
    if spec.get("max_model_calls", 3) < 1:
        raise ResearchError("Model request budget is exhausted")
    if client is None:
        if not os.getenv("OPENAI_API_KEY"):
            raise ResearchError("OPENAI_API_KEY is not configured")
        from openai import OpenAI
        return OpenAI(base_url="https://api.openai.com/v1", timeout=90, max_retries=0)
    return client


def _check(cancelled):
    if cancelled():
        raise ResearchError("Research cancelled; completed checkpoints retained")


def _request(spec, request, emit, cancelled, client, purpose):
    """One request, streamed public text; the caller's emit enforces global budget."""
    _check(cancelled)
    trace = {key: value for key, value in request.items() if key != "instructions"}
    trace.update(purpose=purpose, prompt_version="research-v2",
                 prompt_sha256=hashlib.sha256(request["instructions"].encode()).hexdigest())
    emit(kind="model_call", name=f"OpenAI Responses · {spec['model']}", status="running", inputs=trace)
    response = None
    pending = ""
    stream = None
    try:
        stream = client.responses.create(**request, stream=True)
        for event in stream:
            _check(cancelled)
            event_type = _get(event, "type", "")
            if event_type == "response.output_text.delta":
                pending += _get(event, "delta", "")
                if len(pending) >= 180 or "\n" in pending:
                    emit(kind="model_delta", name=purpose, status="running",
                         message=pending, outputs={"text_delta": pending, "public_output": True})
                    pending = ""
            elif event_type.startswith("response.web_search_call."):
                emit(kind="tool_call", name="web_search", status="completed" if event_type.endswith("completed") else "running",
                     outputs={"event": event_type, "item_id": _get(event, "item_id"), "purpose": purpose})
            elif event_type in {"response.completed", "response.incomplete", "response.failed"}:
                response = _get(event, "response")
        if pending:
            emit(kind="model_delta", name=purpose, status="running",
                 message=pending, outputs={"text_delta": pending, "public_output": True})
        if response is None:
            raise ResearchError("Stream ended without a completed model response")
    except Exception as exc:
        # Error bodies may contain credentials or request contents. Do not persist them.
        safe = str(exc) if isinstance(exc, ResearchError) else f"{type(exc).__name__}; model request failed"
        diagnostic={"http_status":getattr(exc,"status_code",None)}
        # Parameter names and API error codes identify compatibility failures without
        # persisting an error body that might echo user input or credentials.
        for name in ("code","param"):
            value=getattr(exc,name,None)
            if isinstance(value,str) and re.fullmatch(r"[A-Za-z0-9_.\[\]-]{1,100}",value) and not value.startswith('sk-'):
                diagnostic[name]=value
        if diagnostic.get('param'):
            safe+=f" (parameter: {diagnostic['param']}; code: {diagnostic.get('code','unspecified')})"
        if diagnostic.get('http_status')==400 and diagnostic.get('param')=='tools':
            from .store import safe as redact
            body=getattr(exc,'body',None)
            if isinstance(body,dict):
                detail=body.get('message') or (body.get('error') or {}).get('message')
                if isinstance(detail,str): diagnostic['tool_compatibility_detail']=redact(detail)[:500]
        emit(kind="model_call", name=f"OpenAI Responses · {spec['model']}", status="failed", message=safe,
             outputs={"purpose": purpose, **diagnostic})
        raise ResearchError(safe) from None
    finally:
        if stream is not None and callable(getattr(stream, "close", None)):
            stream.close()
    public = [_dict(item) for item in _get(response, "output", [])
              if _get(item, "type") in {"message", "web_search_call"}]
    usage = _get(response, "usage")
    emit(kind="model_call", name=f"OpenAI Responses · {spec['model']}", status="completed",
         outputs={"response_id": _get(response, "id"), "status": _get(response, "status"),
                  "output": public, "usage": _dict(usage) if usage else None, "purpose": purpose})
    if _get(response, "status") != "completed":
        raise ResearchError(f"Model response incomplete: {_get(response, 'status')}")
    if any(_get(item, "type") == "function_call" for item in _get(response, "output", [])):
        raise ResearchError("Unregistered function call rejected; no jobs were executed")
    return response


def _text(response):
    if _get(response, "output_text") is not None:
        return _get(response, "output_text")
    return "".join(_get(part, "text", "") for item in _get(response, "output", [])
                   if _get(item, "type") == "message" for part in _get(item, "content", [])
                   if _get(part, "type") == "output_text")


def _url(value):
    if not isinstance(value, str):
        return None
    try:
        parsed = urlsplit(value)
    except ValueError:
        return None
    if parsed.scheme not in {"https", "http"} or not parsed.hostname or parsed.username or parsed.password:
        return None
    try:
        if parsed.port is not None and parsed.port != {"https": 443, "http": 80}[parsed.scheme]:
            return None
    except ValueError:
        return None
    host = parsed.hostname.lower()
    if not any(host == domain or host.endswith("." + domain) for domain in PRIMARY_DOMAINS):
        return None
    if host == "github.com" and not any(parsed.path.lower() == prefix or parsed.path.lower().startswith(prefix + "/") for prefix in OFFICIAL_REPOSITORIES):
        return None
    return urlunsplit((parsed.scheme, parsed.netloc, parsed.path, parsed.query, ""))


def _sources(response):
    found, annotations, actions = {}, [], []
    def retain(source, cited=False):
        source = _dict(source)
        url = _url(source.get("url"))
        if not url:
            return
        record = found.setdefault(url, {"id": "web_" + hashlib.sha256(url.encode()).hexdigest()[:16],
            "source": "OpenAI web_search", "source_id": url, "url": url,
            "title": source.get("title") or url, "retrieved_at": utc(),
            "source_kind": "official_documentation" if any(x in url for x in ["github.com/", "readthedocs.io", "ambermd.org", "autodock.scripps.edu"]) else "research_publication",
            "evidence_type": "web_search_source", "content_access": "search_excerpt",
            "full_text_retrieved": False, "instruction_trust": "untrusted_data",
            "text": "Source returned by web search; source body was not independently downloaded.",
            "location": "Hosted search context; source identity and citation linkage only", "cited": False})
        record["cited"] = record["cited"] or cited
        if cited:
            annotations.append({"evidence_id": record["id"], "url": url,
                                "start_index": source.get("start_index"), "end_index": source.get("end_index")})
    for item in _get(response, "output", []):
        if _get(item, "type") == "web_search_call":
            action = _dict(_get(item, "action", {}) or {})
            actions.append({"id": _get(item, "id"), "status": _get(item, "status"), **action})
            for source in action.get("sources") or []:
                retain(source)
        if _get(item, "type") == "message":
            for part in _get(item, "content", []):
                for annotation in _get(part, "annotations", []):
                    if _get(annotation, "type") == "url_citation":
                        retain(annotation, cited=True)
    return list(found.values()), annotations, actions


class _SourceText(HTMLParser):
    """Discard executable/navigation markup; retained text remains untrusted data."""
    skipped = {"script", "style", "noscript", "svg", "nav", "header", "footer", "form"}

    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.blocked = []
        self.parts = []

    def handle_starttag(self, tag, attrs):
        if tag in self.skipped:
            self.blocked.append(tag)

    def handle_endtag(self, tag):
        if tag in self.blocked:
            reverse_index = self.blocked[::-1].index(tag)
            del self.blocked[len(self.blocked) - reverse_index - 1:]

    def handle_data(self, data):
        if not self.blocked and data.strip():
            self.parts.append(data.strip())


def _fetch_source_excerpt(url, client, cancelled=lambda: False, *, max_bytes=500_000, max_characters=6000):
    """Read a bounded excerpt with every redirect independently allowlisted."""
    current = _url(url)
    if current is None:
        raise ResearchError("Source URL is outside the primary-source allowlist")
    deadline = time.monotonic() + 12
    for _ in range(4):
        _check(cancelled)
        if time.monotonic() > deadline:
            raise ResearchError("Source exceeded its bounded fetch time")
        with client.stream("GET", current, follow_redirects=False) as response:
            if response.status_code in {301, 302, 303, 307, 308}:
                following = _url(urljoin(current, response.headers.get("location", "")))
                if following is None:
                    raise ResearchError("Source redirect left the primary-source allowlist")
                current = following
                continue
            response.raise_for_status()
            content_type = response.headers.get("content-type", "").split(";", 1)[0].lower()
            if content_type not in {"text/html", "application/xhtml+xml", "text/plain"}:
                raise ResearchError("Source is not a supported HTML or text page; no PDF/full-text claim was made")
            body = bytearray()
            truncated = False
            for chunk in response.iter_bytes():
                _check(cancelled)
                if time.monotonic() > deadline:
                    raise ResearchError("Source exceeded its bounded fetch time")
                room = max_bytes - len(body)
                if len(chunk) > room:
                    body.extend(chunk[:room])
                    truncated = True
                    break
                body.extend(chunk)
            raw = bytes(body)
            decoded = raw.decode(response.encoding or "utf-8", errors="replace")
            if content_type in {"text/html", "application/xhtml+xml"}:
                parser = _SourceText()
                parser.feed(decoded)
                extracted = "\n".join(parser.parts)
            else:
                extracted = decoded
            extracted = re.sub(r"[ \t]+", " ", extracted).strip()
            if any(marker in extracted[:1000].lower() for marker in ("checking your browser", "access denied", "verify you are human", "enable javascript and cookies to continue")):
                raise ResearchError("Source returned an access challenge; article text was not retrieved")
            if len(extracted) < 80:
                raise ResearchError("Source returned insufficient readable text")
            excerpt = extracted[:max_characters]
            return {"source_excerpt": excerpt, "text": excerpt, "fetch_status": "completed", "fetch_url": current,
                    "content_access": "abstract_page_excerpt" if urlsplit(current).hostname == "pubmed.ncbi.nlm.nih.gov" else "page_text_excerpt",
                    "full_text_retrieved": False, "excerpt_truncated": truncated or len(extracted) > max_characters,
                    "fetched_bytes": len(raw), "fetched_bytes_sha256": hashlib.sha256(raw).hexdigest(),
                    "location": "Bounded text excerpt fetched from the cited page; not a verified complete article"}
    raise ResearchError("Source exceeded the bounded redirect limit")


def _enrich_sources(records, emit, cancelled, source_client=None):
    owned = source_client is None
    source_client = source_client or httpx.Client(timeout=httpx.Timeout(10.0, connect=5.0), follow_redirects=False,
                                                trust_env=False, headers={"User-Agent": "BioCentricGPT/0.1 research-source-reader", "Accept": "text/html,text/plain"})
    selected = sorted(records, key=lambda record: not record.get("cited", False))[:4]
    try:
        for record in selected:
            _check(cancelled)
            emit(kind="tool_call", name="read_primary_source", status="running",
                 inputs={"evidence_id": record["id"], "url": record["url"], "maximum_bytes": 500_000,
                         "maximum_excerpt_characters": 6000, "timeout_seconds": 10})
            try:
                result = _fetch_source_excerpt(record["url"], source_client, cancelled)
            except Exception as exc:
                _check(cancelled)
                reason = str(exc) if isinstance(exc, ResearchError) else f"{type(exc).__name__}; source text unavailable"
                result = {"fetch_status": "unavailable", "fetch_error": reason, "source_excerpt": None,
                          "excerpt_truncated": None, "full_text_retrieved": False}
            record.update(result)
            emit(kind="tool_call", name="read_primary_source", status="completed" if result["fetch_status"] == "completed" else "unavailable",
                 outputs={"evidence_id": record["id"], "url": record["url"], **result})
    finally:
        if owned:
            source_client.close()
    for record in records:
        record.setdefault("fetch_status", "not_fetched_within_budget")
        record.setdefault("source_excerpt", None)
        record.setdefault("excerpt_truncated", None)


def research_literature(spec, tool_catalog, emit, *, cancelled=lambda: False, client=None, source_client=None):
    """Use at most four hosted searches, then read up to four primary-page excerpts."""
    owned = client is None
    client = _client(spec, client)
    inputs = [{"role": "user", "content": json.dumps({
        "task": "Search current primary literature and official documentation before proposing any engine. Investigate target-specific validation, docking alternatives, affinity-method prerequisites and reported failure modes. Use at most four well-focused search actions. Return a public literature summary in 1–2 short sentences, normally at most 60 words: state the relevant protein/task finding and the test-set evidence that must be assembled next, with clickable citations and the key evidence gap. Detailed methods remain in the cited source records. Do not select a final engine or claim a constructed test set yet; source reading and test-set construction follow this review.",
        "project": {key: spec.get(key) for key in ("target", "organism", "context", "intended_effect", "workflow_kind", "cpu", "has_gpu", "compute_budget_seconds")},
        "untrusted_tool_catalog": tool_catalog})}]
    # The legacy GPT-4.1 deployment rejected the newer web_search tool definition
    # in the live compatibility check. The documented preview variant retains
    # mandatory search; source allowlisting is enforced on returned citations.
    legacy_search=spec['model'].startswith('gpt-4.1')
    search_tool={"type":"web_search_preview" if legacy_search else "web_search","search_context_size":"medium"}
    if not legacy_search: search_tool['filters']={"allowed_domains":PRIMARY_DOMAINS}
    request = {"model": spec["model"], "instructions": PROMPT.read_text(), "input": inputs,
        "tools": [search_tool],
        "tool_choice": "required", "max_tool_calls": 4, "include": ["web_search_call.action.sources"],
        "max_output_tokens": spec.get("max_model_output_tokens", 1500), "store": False}
    try:
        response = _request(spec, request, emit, cancelled, client, "literature_search")
    finally:
        if owned and callable(getattr(client, "close", None)):
            client.close()
    records, citations, actions = _sources(response)
    completed = [action for action in actions if action.get("status") == "completed"]
    if not completed or not records or not citations:
        emit(kind="model_validation", name="Literature source validation", status="failed",
             message="No completed search with allowed source citations; tool selection remains blocked")
        raise ResearchError("No completed search with allowed source citations; tool selection remains blocked")
    _enrich_sources(records, emit, cancelled, source_client)
    result = {"status": "completed", "records": records, "summary": _text(response), "citations": citations,
              "search_actions": actions, "model_calls": 1, "search_calls": len(actions), "model": spec["model"],
              "prompt_version": "research-v2", "response_id": _get(response, "id"),
              "limitations": ["Bounded literature review; not an exhaustive systematic review.",
                  "Up to four bounded primary-page excerpts were fetched when accessible; no complete full-text reading is claimed.",
                  "Sources and citation linkage were recorded; scientific claim entailment requires review."]}
    emit(kind="research_sources", name="Literature evidence", status="completed", outputs=result)
    return result


def _capabilities(capabilities):
    if isinstance(capabilities, dict):
        if isinstance(capabilities.get("engines"), list):
            capabilities = capabilities["engines"]
        else:
            capabilities = [{"id": key, **value} for key, value in capabilities.items() if isinstance(value, dict)]
    return {str(row.get("id") or row.get("engine") or row.get("name")): row for row in capabilities}


def _testset_summary_facts(testset):
    """Compact counts guide wording without promoting a preview to validation."""
    records=testset.get('records',[])
    return {'scope':testset.get('scope'),'status':testset.get('status'),
            'measured_record_count':len(records),
            'exact_measurement_count':sum(row.get('quantitative_eligible') is True for row in records),
            'reference_pose_control_count':len(testset.get('pose_records',[])),
            'tuning_compound_count':len(testset.get('tuning_ids',[])),
            'held_out_compound_count':len(testset.get('evaluation_ids',[])),
            'predictive_validation':'Not established by constructing this test-set preview'}


def recommend_tools(spec, literature, testset, capabilities, emit, *, cancelled=lambda: False, client=None):
    """Propose executable choices against frozen data; validate IDs and capability."""
    if literature.get("status") != "completed" or not literature.get("search_calls"):
        raise ResearchError("Completed literature search is required before tool selection")
    records = list(literature.get("records", []))
    local_sources = []
    from .biocentric_model import DOCUMENTATION
    local_doc = Path(DOCUMENTATION)
    if local_doc.is_file() and "biocentric_affinity" in _capabilities(capabilities):
        content = local_doc.read_bytes()
        local_sources.append({"id": "LOCAL-BIOCENTRIC-" + hashlib.sha256(content).hexdigest()[:16],
                              "title": "BioCentric inspected local model contract", "source_type": "local_file",
                              "path": str(local_doc), "sha256": hashlib.sha256(content).hexdigest(), "url": None,
                              "text": content.decode()[:16000]})
        records.extend(local_sources)
    evidence_ids = {record["id"] for record in records}
    if not evidence_ids or not testset.get("criteria_sha256") or not testset.get("dataset_sha256"):
        raise ResearchError("Citable research and a frozen test-set preview are required")
    from .testsets import validate_testset, model_selection_view
    validate_testset(testset)
    available = _capabilities(capabilities)
    if not available or not set(available) <= set(ENGINE_ROLES):
        raise ResearchError("Only registered docking and affinity engine capabilities may be proposed")
    schema = Selection.model_json_schema()
    for definition in [schema, schema["$defs"]["Candidate"]]:
        definition["properties"]["citations"]["items"]["enum"] = sorted(evidence_ids)
        definition["properties"]["citations"]["minItems"] = 1
    inputs = [{"role": "user", "content": json.dumps({"task": "Return the tool proposal now using only the actual source IDs and exact engine IDs supplied. Write reason as 1–2 short sentences, normally at most 60 words: connect the protein literature and actual test-set scope to the proposed tool, why it fits, and one key caveat. Use the supplied summary facts; do not describe a reference-only set as an assay benchmark or a sole available engine as a comparison winner. Retain other details in candidates, testset_fit, citations, uncertainties and questions. Choose only available engines. A null affinity engine is correct when prerequisites are missing. Do not alter frozen criteria or claim the test-set preview has been executed.",
        "project": {key: spec.get(key) for key in ("target", "organism", "context", "cpu", "has_gpu", "compute_budget_seconds")},
        "untrusted_literature": literature, "local_file_evidence": local_sources,
        "frozen_testset": model_selection_view(testset), "testset_summary_facts":_testset_summary_facts(testset),
        "available_docking_engines":[key for key,row in available.items() if ENGINE_ROLES[key]=='docking' and row.get('available') is True],
        "capabilities": available})}]
    request = {"model": spec["model"], "instructions": PROMPT.read_text(), "input": inputs, "tools": [], "tool_choice": "none",
               "max_output_tokens": spec.get("max_model_output_tokens", 1500), "store": False,
               "text": {"format": {"type": "json_schema", "name": "tool_selection", "strict": True, "schema": schema}}}
    owned = client is None
    client = _client(spec, client)
    try:
        response = _request(spec, request, emit, cancelled, client, "tool_selection")
    finally:
        if owned and callable(getattr(client, "close", None)):
            client.close()
    try:
        proposal = Selection.model_validate_json(_text(response))
        if not proposal.citations or not set(proposal.citations) <= evidence_ids:
            raise ResearchError("Tool proposal cited missing or unknown source IDs")
        for candidate in proposal.candidates:
            if candidate.engine not in available or not candidate.citations or not set(candidate.citations) <= evidence_ids:
                raise ResearchError("Candidate used an unknown engine or source ID")
            if ENGINE_ROLES[candidate.engine] != candidate.role:
                raise ResearchError("Candidate engine does not support the declared role")
        for selected, role in ((proposal.docking_engine, "docking"), (proposal.affinity_engine, "affinity")):
            if selected and not available.get(selected, {}).get("available", False):
                raise ResearchError(f"Selected engine {selected} is unavailable; no substitution or execution was performed")
            if selected and ENGINE_ROLES.get(selected) != role:
                raise ResearchError(f"Selected engine does not implement the {role} role")
        if not proposal.reason.strip() or not proposal.testset_fit.strip():
            raise ResearchError("Proposal must explain the tool choice and test-set fit")
    except ValueError as exc:
        message = str(exc) if isinstance(exc, ResearchError) else "Tool proposal failed structured-output validation"
        emit(kind="model_validation", name="Tool selection validation", status="failed", message=message)
        raise ResearchError(message) from None
    result = {**proposal.model_dump(), "status": "needs_input" if proposal.questions or not proposal.docking_engine else "proposed",
              "criteria": testset["criteria"], "criteria_sha256": testset["criteria_sha256"], "dataset_sha256": testset["dataset_sha256"],
              "model_calls": 1, "model": spec["model"], "prompt_version": "research-v2", "local_sources": local_sources,
              "validation": "Schema, source-ID linkage and local availability checked; scientific entailment requires review."}
    emit(kind="decision_summary", name="Proposed tools and reasons", status="completed", message=proposal.reason, outputs=result)
    return result


def review_pilot(spec, literature, testset, plan, pilot, emit, *, cancelled=lambda: False, client=None):
    """One public interpretation request; deterministic measured gates remain authoritative."""
    from .testsets import validate_testset, model_selection_view
    validate_testset(testset)
    if literature.get("status") != "completed" or not literature.get("search_calls"):
        raise ResearchError("Completed literature search is required before pilot interpretation")
    evidence_ids = {record["id"] for record in literature.get("records", [])}
    if not evidence_ids:
        raise ResearchError("Pilot interpretation needs recorded literature sources")
    if plan.get("criteria_sha256") != testset["criteria_sha256"] or plan.get("dataset_sha256") != testset["dataset_sha256"]:
        raise ResearchError("Pilot proposal does not match the frozen test set")
    schema = PilotReview.model_json_schema()
    schema["properties"]["citations"]["items"]["enum"] = sorted(evidence_ids)
    schema["properties"]["citations"]["minItems"] = 1
    inputs = [{"role": "user", "content": json.dumps({
        "task": "Write summary as 1–2 short sentences, normally at most 60 words: state the actual measured pilot outcome against the frozen criterion, the next action for the selected tool and protein/test-set scope, and the key caveat. Recommend proceed_exploratory, revise, or need_data. Retain detailed caveats in limitations and methodological sources in citations. Do not change metrics or criteria, claim independent affinity validation from redocking, or claim a comparative win without measured results from multiple engines. The application independently enforces all pass/fail gates. Cite supplied literature IDs for methodological interpretation; measured results are supplied data. Ask a focused follow-up question when a failed gate or missing data requires action.",
        "project": {key: spec.get(key) for key in ("target", "organism", "context")},
        "untrusted_literature": literature, "frozen_testset": model_selection_view(testset), "testset_summary_facts":_testset_summary_facts(testset), "approved_proposal": plan,
        "measured_pilot": pilot})}]
    request = {"model": spec["model"], "instructions": PROMPT.read_text(), "input": inputs,
               "tools": [], "tool_choice": "none", "store": False,
               "max_output_tokens": spec.get("max_model_output_tokens", 1500),
               "text": {"format": {"type": "json_schema", "name": "pilot_interpretation", "strict": True, "schema": schema}}}
    owned = client is None
    client = _client(spec, client)
    try:
        response = _request(spec, request, emit, cancelled, client, "pilot_interpretation")
    finally:
        if owned and callable(getattr(client, "close", None)):
            client.close()
    try:
        review = PilotReview.model_validate_json(_text(response))
        if not review.citations or not set(review.citations) <= evidence_ids:
            raise ResearchError("Pilot interpretation cited missing or unknown source IDs")
    except ValueError as exc:
        message = str(exc) if isinstance(exc, ResearchError) else "Pilot interpretation failed structured-output validation"
        emit(kind="model_validation", name="Pilot interpretation validation", status="failed", message=message)
        raise ResearchError(message) from None
    result = {**review.model_dump(), "status": "advisory", "model_calls": 1, "model": spec["model"],
              "criteria_sha256": testset["criteria_sha256"], "dataset_sha256": testset["dataset_sha256"],
              "predictive_validation": "not_established_by_this_advisory_response",
              "validation": "Source-ID linkage checked; deterministic measured pilot gates cannot be overridden by this response."}
    emit(kind="decision_summary", name="Pilot interpretation", status="completed", message=review.summary, outputs=result)
    return result
