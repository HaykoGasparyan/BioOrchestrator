"""Conversational research workspace over durable, inspectable research jobs.

The interface renders recorded messages and decisions. It never simulates model
thoughts, invokes a model, or changes scientific decisions while refreshing.
"""
from __future__ import annotations

import html
import json
import os
import platform
import re
from pathlib import Path
from typing import Any
from urllib.parse import urlsplit

from .store import safe
from .ui import artifact_bytes, table_rows

ACTIVE = {"queued", "running", "cancelling"}
FINISHED = {"completed", "complete", "succeeded"}
WAITING = {"blocked", "waiting", "awaiting_answer", "awaiting_approval", "needs_input"}
STATUS = {"pending": "Coming up", "queued": "Queued", "running": "Working", "completed": "Complete",
          "complete": "Complete", "succeeded": "Complete", "blocked": "Your input needed",
          "waiting": "Your input needed", "awaiting_answer": "Your input needed",
          "awaiting_approval": "Your input needed", "needs_input": "Your input needed",
          "failed": "Needs attention", "cancelled": "Cancelled", "cancelling": "Cancelling",
          "skipped": "Skipped", "draft": "Ready to start"}
PUBLIC_GOAL = ("Choose and verify a suitable local docking method for the public 1IEP ABL–imatinib "
               "example, then run a small exploratory screen and explain the limits of its affinity scores.")


def research_data(snapshot: dict) -> dict:
    return snapshot.get("data", {}).get("research") or {}


def is_research(snapshot: dict) -> bool:
    return bool(research_data(snapshot)) or snapshot.get("project", {}).get("spec", {}).get("workflow_mode") == "research"


def event_record(event: dict) -> dict:
    """Accept both persisted flat records and wrapped transport events."""
    payload = event.get("payload")
    return safe({**event, **payload}) if isinstance(payload, dict) else safe(event)


def pending_questions(snapshot: dict) -> list[dict]:
    questions = research_data(snapshot).get("questions") or []
    if isinstance(questions, dict):
        questions = [{"id": key, **value} for key, value in questions.items() if isinstance(value, dict)]
    return [q for q in questions if isinstance(q, dict) and q.get("status", "pending")
            in {"pending", "open", "waiting", "awaiting_answer"} and q.get("answer") is None]


def conversation_messages(snapshot: dict) -> list[dict]:
    """Only explicit public conversation records become assistant messages."""
    messages = []
    for raw in snapshot.get("events", []):
        event = event_record(raw)
        if event.get("kind") not in {"research_message", "assistant_message", "user_message"}:
            continue
        message = event.get("message") or event.get("text") or event.get("explanation")
        if not message:
            continue
        role = event.get("role", "user" if event.get("kind") == "user_message" else "assistant")
        messages.append({"id": event.get("id"), "role": "user" if role == "user" else "assistant",
                         "message": str(message), "stage": event.get("stage") or event.get("stage_id"),
                         "created_at": event.get("created_at"), "sources": event.get("sources") or []})
    return messages


def research_summary(snapshot: dict) -> dict:
    steps = snapshot.get("steps", [])
    events = [event_record(event) for event in snapshot.get("events", [])]
    calls = [event for event in events if event.get("kind") == "model_call" and event.get("status") in {"running", "started"}]
    # A response can appear in several events. Usage is counted once per response.
    responses = {}
    for event in events:
        if event.get("kind") != "model_call" or event.get("status") not in FINISHED:
            continue
        output = event.get("outputs") or {}
        if not isinstance(output, dict):
            continue
        usage = output.get("usage") or event.get("usage")
        if isinstance(usage, dict):
            responses[output.get("response_id") or event.get("response_id") or event.get("id", len(responses))] = usage
    input_tokens = sum(int(usage.get("input_tokens", 0) or 0) for usage in responses.values())
    output_tokens = sum(int(usage.get("output_tokens", 0) or 0) for usage in responses.values())
    active = any(job.get("status") in ACTIVE for job in snapshot.get("jobs", []))
    current = next((step for step in steps if step.get("status") == "running"), None)
    current = current or next((step for step in steps if step.get("status") in WAITING | {"failed"}), None)
    current = current or next((step for step in steps if step.get("status") not in FINISHED | {"skipped"}), None)
    return {"active": active, "current": current, "calls": len(calls), "input_tokens": input_tokens if responses else None,
            "output_tokens": output_tokens if responses else None,
            "recorded_stage_seconds": sum(float(event.get("duration_ms", 0) or 0) for event in events if event.get("kind") == "stage" and event.get("status") in FINISHED) / 1000,
            "complete": sum(step.get("status") in FINISHED for step in steps), "total": len(steps)}




def compact_explanation(value: Any, *, sentences: int = 2, limit: int = 330) -> str:
    """An extract of recorded public text; this does not invent a new rationale."""
    text = str(safe(value or ""))
    text = re.sub(r"\ue200[^\ue201]*\ue201", "", text)
    text = re.sub(r"\[([^\]]+)\]\([^)]+\)", r"\1", text)
    text = re.sub(r"(?:^|\n)\s*[#>*-]+\s*", "", text)
    text = re.sub(r"\s+", " ", text).strip()
    parts = re.split(r"(?<=[.!?])\s+", text)
    excerpt = " ".join(parts[:sentences])
    if len(excerpt) > limit:
        excerpt = excerpt[:limit - 1].rsplit(" ", 1)[0].rstrip(" ,;:") + "…"
    return excerpt


def poc_testset_summary(snapshot: dict) -> str:
    if research_data(snapshot).get("validation_run"):
        return ("Local software validation: these recorded calculations test preparation, docking and exports. "
                "No literature review or independent assay benchmark was run; docking scores do not establish binding-affinity accuracy.")
    testset = research_data(snapshot).get("testset") or {}
    if not testset:
        return "First research the methods and freeze a small test set, then choose a tool and verify it before screening."
    measured = len(testset.get("records") or [])
    references = len(testset.get("pose_records") or [])
    evaluation = len(testset.get("evaluation_ids") or [])
    if not evaluation:
        return (f"Test set: {references} pose reference(s), {measured} measured compound record(s); no independent assay benchmark. "
                "The pilot checks pose recovery, not binding-affinity accuracy.")
    return (f"Test set: {measured} measured compound records, including {evaluation} held-out compounds, and {references} pose reference(s). "
            "A frozen test set is not itself a validation result.")


def _route(st: Any, snapshot: dict) -> None:
    validation = research_data(snapshot).get("validation_run")
    groups = [("Research", {"01", "02", "03"}), ("Test set", {"04"}), ("Method", {"05"}),
              ("Pilot", {"06", "07"}), ("Screen", {"08"}), ("Results", {"09"})]
    if validation:
        groups = [(step.get("label") or f"Stage {step['id']}", {str(step["id"])})
                  for step in snapshot.get("steps", [])]
    steps = {str(step.get("id")): step.get("status", "pending") for step in snapshot.get("steps", [])}
    labels = []
    for name, ids in groups:
        statuses = [steps.get(sid, "pending") for sid in ids]
        status = "done" if all(value in FINISHED for value in statuses) else "now" if any(value in ACTIVE | WAITING | {"failed"} for value in statuses) else "later"
        labels.append(f'<span class="poc-route-step {status}">{html.escape(name)}</span>')
    st.markdown('<div class="poc-route">' + '<span class="poc-route-arrow">→</span>'.join(labels) + '</div>', unsafe_allow_html=True)



def method_source_label(source: dict, index: int) -> str:
    title = str(source.get("title") or "").strip()
    if not title or title.startswith(("http://", "https://", "www.")):
        url = str(source.get("url") or title).lower()
        if "vina.scripps.edu" in url or "autodock-vina.readthedocs.io" in url:
            title = "AutoDock Vina documentation"
        elif "github.com/ccsb-scripps/autodock-vina" in url:
            title = "AutoDock Vina repository"
        else:
            title = f"Source {index}"
    return compact_explanation(title, limit=55)


def _method_summary(st: Any, snapshot: dict) -> None:
    research = research_data(snapshot)
    plan = research.get("tool_plan") or next((step.get("outputs") or {} for step in snapshot.get("steps", []) if step.get("id") == "05"), {})
    reason = plan.get("reason") or plan.get("public_explanation") or plan.get("summary")
    if not reason:
        return
    names = {tool.get("id"): tool.get("name", tool.get("id")) for tool in research.get("capabilities", [])}
    engine = plan.get("docking_engine")
    title = f"Docking method · {names.get(engine, engine)}" if engine else "Method summary"
    st.markdown(f"**{title}**")
    st.write(compact_explanation(reason))
    sources = source_rows(resolve_sources(snapshot, plan.get("citations") or plan.get("sources")))
    if sources:
        links = []
        for index, source in enumerate(sources[:2], 1):
            label = method_source_label(source, index).replace("[", "\\[").replace("]", "\\]")
            if source["url"]:
                links.append(f"[{label}]({source['url'].replace(')', '%29')})")
        if links:
            suffix = f" · +{len(sources) - 2} more in details" if len(sources) > 2 else ""
            st.caption("Research basis: " + " · ".join(links) + suffix)


def _history(st: Any, snapshot: dict, messages: list[dict], live_output: dict | None) -> None:
    research = research_data(snapshot)
    with st.expander("Conversation history & full model responses", expanded=False):
        goal = research.get("goal") or snapshot.get("project", {}).get("spec", {}).get("context")
        if goal and not any(message["role"] == "user" and message["message"] == goal for message in messages):
            with st.chat_message("user"):
                st.write(safe(goal))
        for message in messages:
            with st.chat_message(message["role"], avatar=":material/science:" if message["role"] == "assistant" else None):
                literature = research.get("literature") or {}
                text = citation_markdown(message["message"], literature.get("citations")) if message["message"] == literature.get("summary") else message["message"]
                st.markdown(text)
                if message.get("sources"):
                    _sources(st, message["sources"])
        if live_output and live_output.get("text"):
            st.caption("Live public draft · not yet validated")
            st.markdown(live_output["text"])


def live_model_output(snapshot: dict) -> dict | None:
    """Collect only explicitly public output deltas from the current request."""
    current = None
    for raw in snapshot.get("events", []):
        event = event_record(raw)
        if event.get("kind") == "model_call":
            if event.get("status") in {"running", "started"}:
                current = {"text": "", "purpose": (event.get("inputs") or {}).get("purpose", "research"), "stage": event.get("stage")}
            elif event.get("status") in FINISHED | {"failed", "cancelled"}:
                current = None
        elif current is not None and event.get("kind") == "model_delta" and (event.get("outputs") or {}).get("public_output") is True:
            current["text"] += str((event.get("outputs") or {}).get("text_delta") or event.get("message") or "")
    if current and current["text"].lstrip().startswith("{"):
        # Structured output is visible in exact traces; use its public rationale
        # as the readable draft once that JSON field starts arriving.
        match = re.search(r'"(?:reason|summary)"\s*:\s*"((?:[^"\\]|\\.)*)', current["text"])
        if match:
            try:
                current["text"] = json.loads('"' + match.group(1) + '"')
            except ValueError:
                current["text"] = match.group(1).replace("\\n", "\n")
        else:
            current["text"] = ""
    return current


def citation_markdown(text: str, citations: Any) -> str:
    """Replace only validated citation spans with their recorded source links."""
    if not isinstance(citations, list):
        return text
    replacements = []
    for citation in citations:
        if not isinstance(citation, dict):
            continue
        start, end = citation.get("start_index"), citation.get("end_index")
        source = source_rows([citation])
        if not isinstance(start, int) or not isinstance(end, int) or not (0 <= start < end <= len(text)) or not source or not source[0]["url"]:
            continue
        label = str(citation.get("title") or "source").replace("[", "\\[").replace("]", "\\]")
        replacements.append((start, end, f"[{label}]({source[0]['url'].replace(')', '%29')})"))
    next_start = len(text)
    for start, end, link in sorted(replacements, reverse=True):
        if end <= next_start:
            text = text[:start] + link + text[end:]
            next_start = start
    return text


def source_rows(sources: Any) -> list[dict]:
    if isinstance(sources, dict):
        sources = list(sources.values())
    if not isinstance(sources, list):
        return []
    result = []
    for source in sources:
        if not isinstance(source, dict):
            continue
        url = source.get("url") or source.get("source_url") or source.get("link")
        if url:
            try:
                parsed = urlsplit(str(url))
                if parsed.scheme not in {"http", "https"} or not parsed.netloc or parsed.username or parsed.password:
                    url = None
            except ValueError:
                url = None
        result.append({"id": source.get("id") or source.get("evidence_id"),
                       "title": source.get("title") or source.get("name") or source.get("source") or "Recorded source",
                       "url": url, "source_type": source.get("source_type"), "path": source.get("path"), "sha256": source.get("sha256"),
                       "finding": source.get("finding") or source.get("summary") or source.get("snippet") or source.get("text")})
    return result



def resolve_sources(snapshot: dict, citations: Any) -> list[dict]:
    """Join model citation IDs to the actual retrieved source records."""
    if not isinstance(citations, list):
        return []
    available = {source.get("id"): source for source in research_data(snapshot).get("sources", []) if isinstance(source, dict)}
    available.update({source["id"]: source for source in research_data(snapshot).get("tool_plan", {}).get("local_sources", [])})
    resolved = []
    for item in citations:
        source = ({**available.get(item.get("evidence_id"), {}), **item} if isinstance(item, dict) else available.get(item))
        if source:
            resolved.append(source)
    return resolved


def _testset(st: Any, output: dict) -> None:
    st.markdown("**Frozen test-set preview**")
    records, poses = output.get("records", []), output.get("pose_records", [])
    st.write(f"{len(records)} measured compound records · {len(poses)} pose references · "
             f"{len(output.get('tuning_ids', []))} tuning / {len(output.get('evaluation_ids', []))} held-out compounds")
    if records:
        st.dataframe(table_rows(safe(records)), hide_index=True, width="stretch")
    if poses:
        st.dataframe(table_rows(safe(poses)), hide_index=True, width="stretch")
    st.caption(f"Scope: {str(output.get('scope', 'recorded in the manifest')).replace('_', ' ')}")
    if output.get("frozen_sha256") or output.get("manifest_sha256"):
        st.caption(f"Frozen manifest SHA-256: {output.get('frozen_sha256') or output.get('manifest_sha256')}")
    with st.expander("Criteria fixed before the pilot"):
        _record(st, output.get("pilot_criteria") or output.get("criteria"))
        if output.get("exclusions"):
            st.markdown("**Excluded measurements**")
            _record(st, output["exclusions"])


def _sources(st: Any, sources: Any) -> None:
    for source in source_rows(sources):
        title = str(source["title"]).replace("[", "\\[").replace("]", "\\]")
        citation = f" · {source['id']}" if source.get("id") else ""
        if source["url"]:
            st.markdown(f"[{title}]({str(source['url']).replace(')', '%29')}){citation}")
        else:
            st.write(f"{source['title']}{citation}")
        if source.get("source_type") == "local_file":
            st.caption(f"Local file: {source.get('path')} · SHA-256 {source.get('sha256')}")
        if source.get("finding"):
            finding=str(source["finding"])
            if len(finding)>500:
                with st.expander(f"Read source excerpt · {source.get('id') or source['title']}"):
                    st.write(finding)
            else:
                st.caption(finding)


def _record(st: Any, value: Any) -> None:
    if value is None or value == {} or value == [] or value == "":
        st.caption("Nothing recorded yet.")
    elif isinstance(value, (dict, list, tuple)):
        st.json(safe(value), expanded=False)
    else:
        st.write(safe(str(value)))


def _try(st: Any, action: Any, success: str) -> Any:
    try:
        result = action()
    except (ValueError, RuntimeError) as exc:
        st.error(str(safe(str(exc))))
        return None
    except Exception as exc:
        st.error(f"The action could not be completed ({type(exc).__name__}). The saved research is still available.")
        return None
    st.toast(success)
    return result if result is not None else True


def _styles(st: Any) -> None:
    st.markdown("""<style>
    .block-container {max-width:1420px;padding:4.5rem 1.5rem 2rem;}
    @media (max-width:700px) {.block-container {padding-left:1rem;padding-right:1rem;}}
    h1 {font-weight:720!important;letter-spacing:-.055em!important;font-size:2.65rem!important;}
    h2,h3 {letter-spacing:-.025em!important;}
    [data-testid="stSidebar"] {border-right:1px solid #dce6ea;}
    [data-testid="stChatMessage"] {background:transparent;border:0;padding:.7rem 0;}
    [data-testid="stVerticalBlockBorderWrapper"] {border-radius:15px;}
    .research-eyebrow {color:#087f8c;font-size:11px;letter-spacing:.14em;font-weight:750;margin-bottom:8px;}
    .research-lead {color:#54717b;font-size:16px;line-height:1.6;max-width:680px;}
    .research-brand {font-size:20px;font-weight:760;color:#143b43;letter-spacing:-.6px;margin:4px 0 8px;}
    .research-status {display:inline-block;font-size:11px;font-weight:650;padding:4px 9px;border-radius:12px;background:#e6edef;color:#48646f;}
    .research-status.running {background:#d9f1ee;color:#096d63;}
    .research-status.completed {background:#e4f1e8;color:#356c48;}
    .research-status.blocked,.research-status.failed {background:#fff0d5;color:#8b5b16;}
    .research-meta {color:#62808b;font-size:12px;padding:3px 0 12px;}
    .research-route {color:#54717b;font-size:13px;padding:15px 0 5px;line-height:1.9;}
    .poc-route {display:flex;align-items:center;flex-wrap:wrap;gap:9px;margin:6px 0 16px;font-size:13px;font-weight:650;}
    .poc-route-step {padding:5px 10px;border-radius:6px;background:#edf1f3;color:#6b8088;}
    .poc-route-step.done {background:#e1f1e7;color:#286647;}
    .poc-route-step.now {background:#d7eeed;color:#075e60;outline:1px solid #95c8c7;}
    .poc-route-arrow {color:#90a4aa;}

    [data-testid="stChatInput"] {border-color:#afcdce;}
    </style>""", unsafe_allow_html=True)


def _activate(st: Any, pid: str) -> None:
    st.session_state["research_project_id"] = pid
    # Apply the widget value before its next instantiation. Removing a widget
    # key lets the browser backfill the old selected project on a rerun.
    st.session_state["research_pending_selection"] = pid
    st.session_state["research_new"] = False



def project_label(project: dict) -> str:
    """Distinct browser-facing labels avoid ambiguous duplicate project names."""
    pid = str(project["id"])
    mode = "research" if project.get("spec", {}).get("workflow_mode") == "research" or project.get("data", {}).get("research") else "classic"
    name = str(project.get("spec", {}).get("name") or project.get("name") or "Untitled project")
    return f"{mode} · {pid[:6]} · {name}"


def _sidebar(st: Any, service: Any) -> tuple[str | None, str, int]:
    with st.sidebar:
        st.markdown('<div class="research-brand">◈ Bio Orchestrator</div>', unsafe_allow_html=True)
        st.caption("Your research, step by step")
        if st.button("＋ New conversation", width="stretch", key="new_research_conversation"):
            st.session_state["research_new"] = True
            st.session_state["research_project_id"] = None
            st.session_state["research_pending_selection"] = None
            st.rerun()
        projects = service.list_projects()
        choices = {project["id"]: project_label(project) for project in projects}
        if "research_pending_selection" in st.session_state:
            requested = st.session_state.pop("research_pending_selection")
            requested = requested if requested in choices else None
            st.session_state["research_project_picker"] = requested
            st.session_state["research_project_id"] = requested
            st.session_state["research_new"] = requested is None
        saved_url_pid = st.query_params.get("project")
        if "research_project_picker" not in st.session_state and saved_url_pid in choices:
            st.session_state["research_project_picker"] = saved_url_pid
            st.session_state["research_project_id"] = saved_url_pid
            st.session_state["research_new"] = False
        selected = (st.session_state["research_project_id"] if "research_project_id" in st.session_state
                    else st.session_state.get("project_id"))
        if selected not in choices:
            selected = next((project["id"] for project in projects if project.get("spec", {}).get("workflow_mode") == "research" or project.get("data", {}).get("research")), None)
        pid = None
        if choices:
            values = [None, *choices]
            initial = None if st.session_state.get("research_new", False) else selected
            # A persisted widget value is authoritative. Supplying a default
            # index at the same time triggers Streamlit’s dual-state warning.
            picker_index = None if "research_project_picker" in st.session_state else values.index(initial)
            pid = st.selectbox("Conversations & saved projects", values, index=picker_index,
                               format_func=lambda item: "New conversation" if item is None else choices[item],
                               key="research_project_picker")
            st.session_state["research_project_id"] = pid
            st.session_state["research_new"] = pid is None
            if pid is not None and st.query_params.get("project") != pid:
                st.query_params["project"] = pid
            elif pid is None and "project" in st.query_params:
                del st.query_params["project"]
        st.divider()
        st.markdown(f"**Local {platform.system()} workbench**")
        st.caption("Scientific calculations run locally. Web research and GPT use external services.")
        st.caption("Runtime detection, real execution and predictive validation are separate evidence. Inspect each stage for the recorded device and failures.")
        with st.expander("Model & resource budget"):
            model = st.text_input("Language model", value=os.getenv("OPENAI_MODEL", "gpt-4.1-mini"), key="research_model")
            calls = st.number_input("Maximum GPT calls per new project", min_value=1, max_value=10, value=6,
                                    help="A hard cap across planning and result reviews; web search has separate service usage.", key="research_budget")
            st.caption("Starting research authorizes this bounded API use. Each literature pass can use up to 4 hosted search actions and fetch up to 4 primary-source pages. Model tokens and scientific runtime are recorded separately.")
        st.caption("API key configured" if os.getenv("OPENAI_API_KEY") else "API key missing in this app process")
    return pid, model.strip() or "gpt-4.1-mini", int(calls)


def _create(st: Any, service: Any, goal: str, *, example: bool, model: str, calls: int, start: bool) -> None:
    pid = _try(st, lambda: service.create_research_project(goal, use_public_example=example,
                                                         model=model, max_model_calls=calls), "Research conversation created.")
    if not isinstance(pid, str):
        return
    _activate(st, pid)
    if start:
        _try(st, lambda: service.start_research(pid), "Research started. Progress appears here.")
    st.rerun()


def _welcome(st: Any, service: Any, model: str, calls: int) -> None:
    st.markdown('<div class="research-eyebrow">LITERATURE → TEST SET → VERIFICATION → EXPERIMENT</div>', unsafe_allow_html=True)
    st.title("Make the next experiment clear.")
    st.markdown('<div class="research-lead">Describe your docking or binding question. Follow the sources, the tool choices, and the small validation experiment before a larger run.</div>', unsafe_allow_html=True)
    with st.chat_message("assistant", avatar=":material/science:"):
        st.markdown("I’ll research the methods, freeze a small test set, and choose a tool to verify before screening. You’ll see short updates, the ligand list, and any questions needed to continue.")
        st.caption("Full explanations, sources, and model/tool calls stay available in optional stage details.")
    with st.container(border=True):
        st.markdown("**Try the public ABL–imatinib example**")
        st.write("Start with an existing receptor, bound reference ligand, and three small molecules. A pose-recovery pilot can check this setup; it cannot establish binding-affinity accuracy.")
        st.caption(f"Mouse ABL · P00520 · PDB 1IEP · local {platform.system()} · public inputs")
        if st.button("Prepare public example", type="primary", key="research_public_example"):
            _create(st, service, PUBLIC_GOAL, example=True, model=model, calls=calls, start=False)
    st.caption(f"Sending a question starts web research and permits up to {calls} calls to {model}. Scientific execution follows the in-chat review steps.")
    with st.container():
        prompt = st.chat_input("What docking or binding question should we investigate?", key="research_goal", max_chars=4000)
    if prompt:
        _create(st, service, prompt.strip(), example=False, model=model, calls=calls, start=True)


def _question(st: Any, service: Any, pid: str, question: dict) -> None:
    qid = str(question.get("id") or "question")
    with st.chat_message("assistant", avatar=":material/science:"):
        st.markdown("**A quick check before I continue**")
        st.write(question.get("prompt") or question.get("question") or "Please review the next action.")
        if question.get("context"):
            context = question["context"]
            if isinstance(context, dict):
                visible = {name: context[name] for name in ("engine", "cpu", "seed", "exhaustiveness", "n_poses", "timeout_seconds", "ph_assumption", "validation") if name in context}
                if visible:
                    st.dataframe([{"Setting": name.replace("_", " ").capitalize(), "Value": str(value)} for name, value in visible.items()], hide_index=True, width="stretch")
                with st.expander("Proposed inputs, checks & limits", expanded=bool(context.get("missing"))):
                    _record(st, context)
            else:
                st.caption(str(context))
        if question.get("details"):
            with st.expander("Review the supporting details"):
                _record(st, question["details"])
        options = question.get("options") or []
        labels, values = [], []
        for option in options:
            if isinstance(option, dict):
                labels.append(str(option.get("label") or option.get("value")))
                values.append(option.get("value", option.get("label")))
            else:
                labels.append(str(option))
                values.append(option)
        answer_slot = st.empty()
        with answer_slot.form(f"research_answer_{pid}_{qid}"):
            if labels:
                choice = st.radio("Your decision", list(range(len(labels))), index=None,
                                  format_func=lambda index: labels[index], key=f"research_choice_{pid}_{qid}")
                descriptions = [str(option.get("description")) for option in options if isinstance(option, dict) and option.get("description")]
                if descriptions:
                    st.caption(" · ".join(descriptions))
                answer = values[choice] if choice is not None else None
            else:
                answer = st.text_area("Your answer", placeholder=question.get("placeholder") or "Add the requested details…",
                                      key=f"research_answer_text_{pid}_{qid}", height=100)
            submit_label = "Rebuild test set & proposal" if question.get("kind") == "research_input_revision" else "Send & continue"
            submit = st.form_submit_button(submit_label, type="primary")
        if submit:
            if answer is None or isinstance(answer, str) and not answer.strip():
                st.warning("Choose an option or enter an answer to continue.")
            elif _try(st, lambda: service.answer_research(pid, qid, answer),
                      "Revision requested. Add your revised inputs below." if answer == "revise_inputs" else "Answer saved. Research is continuing."):
                answer_slot.empty()
                st.rerun()
        if question.get("kind") == "research_input_revision":
            st.caption("Open Add scientific input files below the research workspace, save revised inputs, then rebuild. The same job continues with the new test set and proposal.")
        else:
            st.caption("Your answer is saved with the research record. Continuing resumes the same job; revising waits for your new inputs.")


def _event_trace(st: Any, events: list[dict]) -> None:
    for event in events:
        status = STATUS.get(event.get("status"), str(event.get("status", "Recorded")))
        name = event.get("name") or event.get("kind", "Event").replace("_", " ").title()
        with st.expander(f"{name} · {status}"):
            if event.get("created_at"):
                st.caption(str(event["created_at"]))
            if isinstance(event.get("duration_ms"), (int, float)):
                st.caption(f"Recorded duration: {event['duration_ms'] / 1000:.2f} seconds")
            public = event.get("message") or event.get("explanation")
            if public:
                st.write(public)
            for label, key in [("Request / inputs", "inputs"), ("Response / outputs", "outputs"),
                               ("Usage", "usage"), ("Command", "command"), ("Error", "error")]:
                if event.get(key) is not None:
                    st.markdown(f"**{label}**")
                    _record(st, event[key])
            with st.expander("Complete recorded event"):
                st.json(safe(event), expanded=False)


STAGE_SHORT_LABELS = {
    "01": "Research brief", "02": "Literature search", "03": "Target evidence",
    "04": "Test set", "05": "Tool choice", "06": "Pilot", "07": "Pilot review",
    "08": "Screening", "09": "Results",
}


def stage_has_model_call(snapshot: dict, sid: str) -> bool:
    return any(event_record(event).get("kind") == "model_call" and
               str(event_record(event).get("stage") or event_record(event).get("stage_id")) == sid
               for event in snapshot.get("events", []))


def stage_public_summary(snapshot: dict, step: dict) -> str:
    """Describe the recorded decision, distinguishing model output from rules."""
    output = step.get("outputs") or {}
    if isinstance(output, dict):
        reason = ((output.get("model_review") or {}).get("summary") or output.get("public_explanation")
                  or output.get("explanation") or output.get("rationale") or output.get("summary") or output.get("reason"))
        if reason:
            return compact_explanation(reason, limit=185)
        if "pose_records" in output:
            return (f"Frozen {len(output.get('records') or [])} measured records and {len(output.get('pose_records') or [])} pose references before choosing a tool.")
        if "poses" in output:
            return f"Recorded {len(output.get('poses') or [])} pilot poses. Physical validity and pose recovery are checked before screening."
    return compact_explanation(step.get("reason") or step.get("explanation") or "No decision recorded yet.", limit=185)


def _stage_details(st: Any, snapshot: dict, step: dict) -> None:
    sid = str(step.get("id", ""))
    status = str(step.get("status", "pending"))
    label = str(step.get("label") or step.get("name") or step.get("title") or f"Stage {sid}")
    st.markdown(f"### {sid} · {label}")
    actor = "AI decision and recorded actions" if stage_has_model_call(snapshot, sid) else "Recorded workflow actions"
    st.caption(f"{STATUS.get(status, status.title())} · {actor}")
    if step.get("tool"):
        st.caption(f"Recorded tools: {step['tool']}")
    output = step.get("outputs") or {}
    if isinstance(output, dict):
        explanation = ((output.get("model_review") or {}).get("summary") or output.get("public_explanation")
                       or output.get("explanation") or output.get("rationale") or output.get("summary") or output.get("reason"))
        if explanation:
            st.markdown("**AI decision**" if stage_has_model_call(snapshot, sid) else "**Workflow update**")
            brief = compact_explanation(explanation, limit=450)
            st.write(brief)
            if str(explanation).strip() != brief:
                with st.expander("Full public explanation", expanded=False):
                    st.markdown(citation_markdown(str(safe(explanation)), output.get("citations")))
        else:
            st.write(stage_public_summary(snapshot, step))
        members = output.get("members") or output.get("test_set") or output.get("compounds")
        if isinstance(members, list) and members and all(isinstance(member, dict) for member in members):
            st.markdown("**Test set / compound records**")
            st.dataframe(table_rows(safe(members)), hide_index=True, width="stretch")
        if "docking_engine" in output:
            st.markdown("**Proposed tools**")
            st.write(f"Docking: {output.get('docking_engine') or 'No suitable tool selected'}")
            st.write(f"Affinity: {output.get('affinity_engine') or 'No affinity calculation selected'}")
            if output.get("testset_fit"):
                st.caption(str(output["testset_fit"]))
            if output.get("candidates"):
                st.dataframe(table_rows(safe(output["candidates"])), hide_index=True, width="stretch")
            for uncertainty in output.get("uncertainties") or []:
                st.caption(f"Uncertainty: {uncertainty}")
        if "pose_records" in output:
            _testset(st, output)
        if output.get("capabilities"):
            st.markdown("**Tools on this machine**")
            rows = [{"Tool": item.get("name", item.get("id")), "Role": item.get("role"),
                     "Status": "Available" if item.get("available") else "Unavailable",
                     "Why / requirements": "; ".join(item.get("blockers") or []) or item.get("validation_status", "")}
                    for item in output["capabilities"]]
            st.dataframe(rows, hide_index=True, width="stretch")
        if output.get("checks"):
            for name, passed in output["checks"].items():
                st.write(f"{'✓' if passed else '×'} {name.replace('_', ' ').capitalize()}")
        if output.get("comparison"):
            st.markdown("**Frozen comparative protocols**")
            comparison = output["comparison"]
            st.caption(f"Frozen before predictions · {comparison.get('comparison_sha256', '')}")
            candidate_rows = []
            for candidate in comparison.get("candidates", []):
                parameters = candidate.get("parameters", {})
                candidate_rows.append({"Protocol": candidate.get("protocol_id"), "Engine": candidate.get("engine_id"),
                                       "Effort": str(parameters.get("search_mode") or parameters.get("exhaustiveness", "recorded")),
                                       "Device": parameters.get("device", "cpu"), "CPU threads": parameters.get("cpu")})
            st.dataframe(candidate_rows, hide_index=True, width="stretch")
            with st.expander("Frozen parameters, criteria and input hashes"):
                st.json(comparison, expanded=False)
        if output.get("selection_lock"):
            lock = output["selection_lock"]
            st.markdown("**Locked method selection**")
            st.write(lock.get("reason") or lock.get("status"))
            objective_label = "Top-pose RMSD (Å)" if output.get("comparison", {}).get("criteria", {}).get("mode") == "pose_software" else "Tuning objective"
            summaries = [{"Protocol": row.get("protocol_id"), objective_label: row.get("tuning_objective"),
                          "Physical check": row.get("checks", {}).get("physical_validity"),
                          "Eligible": row.get("eligible"), "Training overlap": row.get("training_overlap")}
                         for row in lock.get("candidate_results", [])]
            if summaries:
                st.dataframe(summaries, hide_index=True, width="stretch")
            st.json({"selected_protocol": lock.get("selected_protocol"), "selection_sha256": lock.get("selection_sha256"),
                     "candidate_results": lock.get("candidate_results")}, expanded=False)
        if output.get("heldout_evaluation"):
            evaluation = output["heldout_evaluation"]
            st.markdown("**Independent evaluation**")
            st.write(evaluation.get("reason") or evaluation.get("status"))
            with st.expander("Held-out evidence and limitations"):
                st.json(evaluation)
        affinity = output.get("affinity")
        if isinstance(affinity, dict):
            st.markdown("**Binding-affinity calculation**")
            st.write(f"{str(affinity.get('status', 'unrecorded')).replace('_', ' ').capitalize()}: {affinity.get('reason', '')}")
        metrics = output.get("metrics")
        if isinstance(metrics, dict) and metrics:
            st.markdown("**Recorded checks**")
            for name, value in metrics.items():
                st.write(f"{name.replace('_', ' ').capitalize()}: {value}")
        for limitation in output.get("limitations") or []:
            st.caption(f"Limitation: {limitation}")
        sources = resolve_sources(snapshot, output.get("sources") or output.get("citations") or (output.get("model_review") or {}).get("citations"))
        if sources:
            st.markdown("**Sources used at this stage**")
            _sources(st, sources)
    events = [event for event in map(event_record, snapshot.get("events", []))
              if str(event.get("stage") or event.get("stage_id")) == sid
              and event.get("kind") not in {"research_message", "assistant_message", "user_message"}]
    if events:
        with st.expander(f"Exact recorded steps & calls · {len(events)}", expanded=False):
            _event_trace(st, events)
    elif status == "pending":
        st.caption("This stage has not run. Its results and calls will appear here.")
    with st.expander("Inputs, preprocessing & saved outputs"):
        for title, key in [("Inputs", "inputs"), ("Preprocessing", "preprocessing"),
                           ("Postprocessing", "postprocessing"), ("Saved outputs", "outputs")]:
            st.markdown(f"**{title}**")
            _record(st, step.get(key))



def _default_view(snapshot: dict) -> str:
    if snapshot.get("project", {}).get("state") in FINISHED:
        return "results"
    ids = {str(step["id"]) for step in snapshot.get("steps", [])}
    pending = next((str(question.get("stage_id")) for question in pending_questions(snapshot)
                    if str(question.get("stage_id")) in ids), None)
    current = research_summary(snapshot)["current"] or {}
    return pending or str(current.get("id") or (snapshot.get("steps") or [{}])[0].get("id") or "results")


def _interactive_workspace(st: Any, snapshot: dict, root: Path, pid: str) -> None:
    from .ligand_results import render_ligand_results

    steps = snapshot.get("steps", [])
    ids = {str(step["id"]) for step in steps}
    selection_key = f"research_view_{pid}"
    explicit_key = f"research_view_explicit_{pid}"
    if not st.session_state.get(explicit_key) or st.session_state.get(selection_key) not in ids | {"results"}:
        st.session_state[selection_key] = _default_view(snapshot)
    selected = st.session_state[selection_key]
    if not steps:
        render_ligand_results(st, snapshot, root, pid)
        return
    left, right = st.columns([1, 2.2], gap="medium")
    with left:
        validation = research_data(snapshot).get("validation_run")
        st.markdown("### Validation steps" if validation else "### Research steps")
        st.caption("Click a step to inspect its recorded actions." if validation else "Click a step to inspect the decision and its evidence.")
        with st.container(height=650, border=False, key=f"research_navigator_{pid}"):
            for index, step in enumerate(steps, 1):
                sid = str(step["id"])
                status = str(step.get("status", "pending"))
                label = (step.get("label") or f"Stage {index}") if validation else STAGE_SHORT_LABELS.get(sid, step.get("label") or f"Stage {index}")
                with st.container(border=True, key=f"research_card_{pid}_{sid}"):
                    if st.button(f"{sid} · {label}", key=f"research_toggle_{pid}_{sid}",
                                 type="primary" if selected == sid else "secondary", width="stretch"):
                        st.session_state[selection_key] = sid
                        st.session_state[explicit_key] = True
                        st.rerun()
                    actor = "AI" if stage_has_model_call(snapshot, sid) else "Workflow"
                    st.caption(f"{STATUS.get(status, status.title())} · {actor}")
                    st.write(compact_explanation(stage_public_summary(snapshot, step), limit=95))
    with right:
        controls = st.columns(2)
        with controls[0]:
            if st.button("Ligands & results", type="primary" if selected == "results" else "secondary",
                         width="stretch", key=f"research_results_view_{pid}"):
                st.session_state[selection_key] = "results"
                st.session_state[explicit_key] = True
                st.rerun()
        with controls[1]:
            if st.button("Follow current step", width="stretch", key=f"research_follow_current_{pid}"):
                st.session_state[explicit_key] = False
                st.rerun()
        if selected == "results":
            render_ligand_results(st, snapshot, root, pid)
            _method_summary(st, snapshot)
            verification = research_data(snapshot).get("verification") or {}
            pilot_reason = (verification.get("model_review") or {}).get("summary")
            if pilot_reason:
                st.markdown("**Pilot decision**")
                st.write(compact_explanation(pilot_reason))
            st.caption(poc_testset_summary(snapshot))
        else:
            step = next(step for step in steps if str(step["id"]) == selected)
            _stage_details(st, snapshot, step)
            with st.expander("Ligands & results alongside this step", expanded=False):
                render_ligand_results(st, snapshot, root, pid)


def _uploads(st: Any, service: Any, pid: str, *, active: bool) -> None:
    with st.expander("Add scientific input files"):
        st.caption("Files remain in this project’s local artifact store. Upload the matching coordinate set before approving preparation or execution.")
        native = st.checkbox("Advanced: native inputs for an optional engine", key=f"research_native_inputs_{pid}")
        kinds = {"structure": "Protein structure · PDB / mmCIF", "receptor": "Prepared receptor · PDBQT",
                 "reference": "Bound reference ligand · SDF", "library": "Compound library · SDF / SMILES",
                 "measurements": "Measured assay results · CSV"}
        if native:
            kinds = {"native_input": "Native scientific input file", "native_engine_pilot": "Pilot manifest · JSON",
                     "native_engine_screening": "Screening manifest · JSON", "native_engine_affinity": "Affinity manifest · JSON"}
            st.caption("Optional engines need their own reviewed files and installed runtime. Manifests reference registered artifact IDs, and their inputs are checked before execution. See docs/engines.md for the exact contract.")
        kind = st.selectbox("Input type", list(kinds), format_func=lambda item: kinds[item], key=f"research_upload_kind_{pid}_{native}")
        extensions = {"structure": ["pdb", "cif", "mmcif"], "receptor": ["pdbqt"], "reference": ["sdf"], "library": ["sdf", "smi", "smiles"], "measurements": ["csv"]}
        extensions.update({"native_input": ["prmtop", "parm7", "nc", "mdcrd", "map", "fld", "pdb", "pdbqt", "sdf", "json"],
                           "native_engine_pilot": ["json"], "native_engine_screening": ["json"], "native_engine_affinity": ["json"]})
        if kind == "measurements":
            st.caption("Use real measurements for the project’s target and species. Include molecular SMILES for scaffold separation and retain the assay, endpoint, original units and relation. Different assays are kept separate; missing or incompatible data are shown as exclusions.")
            st.download_button("Download assay CSV headers", "compound_id,canonical_smiles,endpoint,value,units,relation,assay_id,organism,construct,assay_format,conditions,publication\n",
                               file_name="assay-results-template.csv", mime="text/csv", key=f"research_assay_template_{pid}")
        upload = st.file_uploader("Choose a file", type=extensions[kind], disabled=active, key=f"research_upload_{pid}_{kind}")
        if st.button("Save input to this project", disabled=upload is None or active, key=f"research_save_upload_{pid}"):
            action = (lambda: service.import_measurements(pid, upload.getvalue())) if kind == "measurements" else (lambda: service.upload(pid, upload.name, upload.getvalue(), kind))
            if _try(st, action, "Scientific input saved."):
                st.rerun()
        if active:
            st.caption("Wait for an input question, or cancel the current work, before changing its inputs.")


def _artifacts(st: Any, snapshot: dict, root: Path, pid: str) -> None:
    artifacts = snapshot.get("artifacts", [])
    if not artifacts:
        return
    with st.expander(f"Research files & reproducibility · {len(artifacts)} artifacts", expanded=False):
        choices = {str(artifact["id"]): artifact for artifact in artifacts}
        selected = st.selectbox("Inspect an archived file", [None, *choices],
                                format_func=lambda aid: "Choose a file…" if aid is None else
                                f"{choices[aid].get('kind', 'output')} · {choices[aid].get('name', 'Artifact')} · {aid[:8]}",
                                key=f"research_artifact_picker_{pid}")
        if selected is not None:
            artifact = choices[selected]
            label = str(artifact.get("name") or selected)
            st.caption(f"{label} · {artifact.get('kind', 'output')} · Artifact ID {selected}")
            st.caption(f"SHA-256: {artifact.get('sha256', 'unrecorded')}")
            content = artifact_bytes(artifact, root / "projects" / pid)
            if content is not None:
                st.download_button(f"Download {label}", content, file_name=label,
                                   mime=artifact.get("media_type") or "application/octet-stream",
                                   key=f"research_download_{pid}_{selected}")
            else:
                st.warning("This file is missing or failed its integrity check.")
        st.download_button("Download complete event log", json.dumps(safe(snapshot.get("events", [])), indent=2),
                           file_name="research-events.json", mime="application/json", key=f"research_event_export_{pid}")


def _workspace(st: Any, service: Any, pid: str, snapshot: dict, root: Path) -> None:
    project = snapshot.get("project", {})
    spec = project.get("spec", {})
    research = research_data(snapshot)
    summary = research_summary(snapshot)
    state = project.get("state", "draft")
    validation = research.get("validation_run")
    eyebrow = "LOCAL VALIDATION · NO PLANNER CALLS" if validation else "DOCKING PROOF OF CONCEPT"
    st.markdown(f'<div class="research-eyebrow">{eyebrow}</div>', unsafe_allow_html=True)
    st.title(spec.get("name") or "Docking & affinity research")
    if validation:
        st.caption(f"Project {pid[:6]} · Local validation · no planner calls · Local {platform.system()}")
    else:
        st.caption(f"Project {pid[:6]} · {spec.get('model', 'gpt-4.1-mini')} · {summary['calls']} / {spec.get('max_model_calls', 6)} GPT calls · Local {platform.system()}")
    _route(st, snapshot)
    if summary["active"]:
        current = summary["current"] or {}
        st.info(f"Working: {current.get('label') or research.get('phase') or 'starting research'}", icon=":material/progress_activity:")
        if st.button("Cancel current work", key=f"research_cancel_{pid}"):
            if _try(st, lambda: service.cancel(pid), "Cancellation requested."):
                st.rerun()
    elif not validation and (state == "draft" or not snapshot.get("jobs")):
        st.caption(f"Start permits web research and up to {spec.get('max_model_calls', 6)} GPT calls. You'll review the pilot before it runs.")
        if st.button("Start research", type="primary", key=f"research_start_{pid}"):
            if _try(st, lambda: service.start_research(pid), "Research started."):
                st.rerun()
    elif not validation and state in {"failed", "cancelled", "blocked"} and not pending_questions(snapshot):
        st.warning("Research needs attention. Open the stage details to inspect the saved result.")
        if st.button("Continue research", type="primary", key=f"research_continue_{pid}"):
            if _try(st, lambda: service.start_research(pid), "Research is continuing from saved progress."):
                st.rerun()
    elif state in FINISHED:
        st.success("Run complete. Review the ranked ligands and download the results below.")
    if state not in FINISHED:
        st.caption(poc_testset_summary(snapshot))
    messages = conversation_messages(snapshot)
    live_output = live_model_output(snapshot) if summary["active"] else None
    if live_output:
        st.caption("Live GPT update · public draft")
        st.write(compact_explanation(live_output.get("text"), limit=250) or "Preparing the next explanation…")
    elif state not in FINISHED and not pending_questions(snapshot):
        latest = next((message for message in reversed(messages) if message["role"] == "assistant"), None)
        if latest:
            st.caption(compact_explanation(latest["message"], sentences=1, limit=220))
    for question in pending_questions(snapshot):
        _question(st, service, pid, question)

    _interactive_workspace(st, snapshot, root, pid)
    if not validation:
        _uploads(st, service, pid, active=summary["active"])
    if research.get("engine_preflight"):
        with st.expander("Selected engine: input preflight and execution requirements", expanded=False):
            preflight = research["engine_preflight"]
            for blocker in preflight.get("blockers") or []:
                st.warning(str(blocker))
            _record(st, preflight)
    sources = research.get("sources") or research.get("literature_sources")
    if sources:
        with st.expander("Literature used in this research", expanded=False):
            _sources(st, sources)
    _history(st, snapshot, messages, live_output)
    with st.expander("Run details & resource usage", expanded=False):
        if validation:
            st.write("Local validation · no planner calls. Scientific engine and affinity-model execution is recorded separately.")
        else:
            st.write(f"Language model: {spec.get('model', 'gpt-4.1-mini')}; {summary['calls']} of {spec.get('max_model_calls', 6)} permitted calls used.")
        if summary["input_tokens"] is not None:
            st.write(f"Recorded tokens: {summary['input_tokens']:,} input + {summary['output_tokens']:,} output.")
        if validation:
            report = next((step.get("outputs", {}).get("validation_report") for step in snapshot.get("steps", [])
                           if isinstance(step.get("outputs", {}).get("validation_report"), dict)), {})
            resources = report.get("resources") or validation.get("resources") or {}
            if isinstance(resources.get("wall_seconds"), (int, float)):
                st.write(f"Recorded validation-suite runtime: {resources['wall_seconds']:,.1f} seconds.")
                st.caption("Includes the original preparation, docking variants, and export/checkpoint checks; importing this saved view did not rerun them.")
            else:
                st.caption("The source report does not record total validation runtime.")
            protocol_times = [{"Protocol": name.replace("_", " "), "Recorded seconds": entry["duration_seconds"]}
                              for name, entry in (report.get("protocols") or {}).items()
                              if isinstance(entry, dict) and isinstance(entry.get("duration_seconds"), (int, float))]
            if protocol_times:
                st.dataframe(protocol_times, hide_index=True, width="stretch")
            if resources:
                _record(st, resources)
        else:
            st.write(f"{spec.get('cpu', 1)} CPU core(s); {spec.get('compute_budget_seconds', 1800):,} second run budget.")
        if summary["recorded_stage_seconds"]:
            st.write(f"{summary['recorded_stage_seconds']:.1f} seconds in completed stages.")
        if not validation:
            st.caption("Each literature pass allows up to 4 hosted search actions and 4 primary-source page reads.")
    _artifacts(st, snapshot, root, pid)


def run_app() -> None:
    import streamlit as st
    from .service import Service

    st.set_page_config(page_title="Bio Orchestrator · Research", page_icon="◈", layout="wide", initial_sidebar_state="expanded")
    # A fresh visit is the new design brief. Explicit saved-project links and the
    # existing workbench retain their execution/review behavior unchanged.
    if not st.query_params.get("project") and st.query_params.get("view") != "workbench" and os.getenv("BIOCENTRIC_UI_MODE") != "classic":
        from .discovery_studio import render_studio
        render_studio(st)
        return
    default_mode = "Classic workbench" if os.getenv("BIOCENTRIC_UI_MODE") == "classic" else "Research conversation"
    mode = st.sidebar.radio("Workspace", ["Research conversation", "Classic workbench"],
                            index=1 if default_mode == "Classic workbench" else 0, key="workspace_mode", label_visibility="collapsed")
    if mode == "Classic workbench":
        from .ui import run_app as classic_app
        classic_app(configure_page=False)
        return
    _styles(st)
    root = Path(os.getenv("BIOCENTRIC_HOME", ".biocentric")).expanduser().resolve()

    @st.cache_resource
    def get_service(home: str, implementation_version: int) -> Any:
        # Refresh the shared service when its implementation changes during a
        # live development session; existing jobs and credentials stay intact.
        return Service(home)

    service = get_service(str(root), Path(__file__).with_name("service.py").stat().st_mtime_ns)
    pid, model, calls = _sidebar(st, service)
    if pid is None:
        _welcome(st, service, model, calls)
        return
    try:
        snapshot = service.snapshot(pid)
    except Exception as exc:
        st.error(f"Cannot open the saved project ({type(exc).__name__}).")
        return
    if not is_research(snapshot):
        st.title(snapshot.get("project", {}).get("spec", {}).get("name", "Saved project"))
        st.info("This project belongs to the original workbench. Its calculations, model responses, and review decisions remain available there.")
        st.caption("Select Classic workbench in the sidebar to inspect this project, or start a new research conversation above.")
        st.session_state["project_id"] = pid
        return
    active = research_summary(snapshot)["active"]

    @st.fragment(run_every=2 if active else None)
    def live_conversation() -> None:
        fresh = service.snapshot(pid)
        if research_summary(fresh)["active"] != active:
            st.rerun()
        with st.empty().container():
            _workspace(st, service, pid, fresh, root)

    live_conversation()
