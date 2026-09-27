"""An inspectable workbench over the same persisted service used by the CLI.

Rendering helpers deliberately accept dictionaries: missing scientific evidence stays
missing, and new adapter-specific provenance remains visible without UI changes.
"""

from __future__ import annotations

import csv
import hashlib
import html
import io
import json
import os
from pathlib import Path
from typing import Any, Callable

from .explain import model_activity, stage_explanation


STAGES = [
    "Target dossier", "Known chemistry & controls", "Strategy selection",
    "Structure selection", "Protein & ligand preparation", "Protocol & pilot",
    "Benchmark gate", "Screening & inspection", "Prioritization",
    "Experimental handoff",
]
ACTIVE = {"running", "queued", "pending", "cancelling"}
DONE = {"completed", "complete", "succeeded", "success"}
STATUS_LABELS = {
    "pending": "Not started", "queued": "Queued", "running": "Running",
    "completed": "Complete", "complete": "Complete", "succeeded": "Complete",
    "blocked": "Needs a decision", "failed": "Failed", "cancelled": "Cancelled",
    "canceled": "Cancelled", "skipped": "Skipped", "draft": "Ready to plan",
    "waiting": "Waiting", "awaiting_approval": "Needs a decision",
}
SECRET_KEYS = {"api_key", "authorization", "access_token", "refresh_token", "secret", "password"}


def safe_data(value: Any) -> Any:
    """Do not render accidental credentials from a future adapter's diagnostics."""
    if isinstance(value, dict):
        return {str(k): ("[redacted]" if str(k).lower() in SECRET_KEYS else safe_data(v)) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [safe_data(v) for v in value]
    if isinstance(value, Path):
        return str(value)
    return value


def is_fixture(snapshot: dict[str, Any]) -> bool:
    project = snapshot.get("project", {})
    spec = project.get("spec", {})
    data = snapshot.get("data", {})
    markers = [snapshot.get("fixture"), project.get("fixture"), spec.get("fixture"),
               spec.get("demo"), data.get("fixture"), data.get("demo")]
    modes = {str(spec.get("mode", "")).lower(), str(data.get("mode", "")).lower(),
             str(project.get("mode", "")).lower()}
    return any(bool(v) for v in markers) or bool(modes & {"demo", "fixture", "offline", "offline_fixture"})


def pipeline_summary(snapshot: dict[str, Any]) -> dict[str, Any]:
    steps = snapshot.get("steps", [])
    complete = sum(str(s.get("status", "")).lower() in DONE for s in steps)
    running = [s for s in steps if s.get("status") == "running"]
    blocked = [s for s in steps if s.get("status") in {"blocked", "failed", "awaiting_approval"}]
    return {"complete": complete, "total": len(steps) or len(STAGES),
            "current": (running or blocked or [s for s in steps if s.get("status") not in DONE] or [{}])[0],
            "active": any(j.get("status") in ACTIVE for j in snapshot.get("jobs", []))}


def table_rows(records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Preserve nested provenance in human-inspectable tabular cells."""
    return [{str(k): json.dumps(safe_data(v), ensure_ascii=False, default=str) if isinstance(v, (dict, list, tuple))
             else safe_data(v) for k, v in row.items()} for row in records]


def candidate_csv(records: list[dict[str, Any]]) -> str:
    rows = table_rows(records)
    if not rows:
        return ""
    fields = list(dict.fromkeys(k for row in rows for k in row))
    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=fields)
    writer.writeheader()
    for row in rows:
        # CSVs are often opened in desktop spreadsheets; scientific text is not a formula.
        writer.writerow({k: "'" + v if isinstance(v, str) and v.startswith(("=", "+", "-", "@", "\t", "\r")) else v
                         for k, v in row.items()})
    return output.getvalue()


def artifact_bytes(artifact: dict[str, Any], root: Path) -> bytes | None:
    """Only serve regular files inside the local project artifact store."""
    raw = artifact.get("path")
    if not raw:
        return None
    path = Path(raw).expanduser()
    if not path.is_absolute():
        path = root / path
    try:
        path = path.resolve()
        path.relative_to(root.resolve())
        if path.is_file() and path.stat().st_size <= 50 * 1024 * 1024:
            content = path.read_bytes()
            if artifact.get("sha256") and hashlib.sha256(content).hexdigest() != artifact["sha256"]:
                return None
            return content
    except (OSError, ValueError):
        pass
    return None


def _json(st: Any, value: Any, empty: str = "No record yet.") -> None:
    if value is None or value == {} or value == [] or value == "":
        st.caption(empty)
    elif isinstance(value, (dict, list, tuple)):
        st.json(safe_data(value), expanded=False)
    else:
        st.text(str(value))


def _badge(value: Any) -> str:
    status = str(value or "pending").lower()
    color = "good" if status in DONE else "live" if status in ACTIVE - {"pending"} else "warn" if status in {"blocked", "failed", "awaiting_approval"} else "quiet"
    return f'<span class="status {color}">{html.escape(STATUS_LABELS.get(status, status.replace("_", " ").title()))}</span>'


def _attempt(st: Any, action: Callable[[], Any], success: str = "Saved.") -> Any:
    try:
        result = action()
    except Exception as exc:
        st.error(str(exc))
        return None
    st.toast(success)
    return result if result is not None else True


def _styles(st: Any) -> None:
    st.markdown("""<style>
    .block-container {padding-top:2rem; padding-bottom:3rem; max-width:1500px;}
    h1 {font-weight:740!important; letter-spacing:-.055em!important;}
    h2,h3 {letter-spacing:-.035em!important;}
    [data-testid="stSidebar"] {border-right:1px solid #dce6ea;}
    [data-testid="stMetric"] {background:#fff;border:1px solid #dce6ea;border-radius:12px;padding:15px 18px;}
    [data-testid="stMetricValue"] {font-size:1.6rem;}
    [data-testid="stVerticalBlockBorderWrapper"] > div {border-color:#dce6ea;}
    .brand {font-size:21px;font-weight:750;letter-spacing:-.7px;color:#13333d;}
    .brand-mark {display:inline-block;color:#087f8c;font-size:25px;margin-right:8px;}
    .eyebrow {color:#58707c;text-transform:uppercase;letter-spacing:.13em;font-size:11px;font-weight:750;margin:4px 0 10px;}
    .muted {color:#647b87;font-size:13px;}
    .status {display:inline-block;border-radius:5px;padding:4px 8px;font-size:11px;font-weight:650;letter-spacing:.025em;}
    .good {background:#def3e9;color:#21634b}.live {background:#ddf2f5;color:#066b77}
    .warn {background:#fff0d5;color:#89580d}.quiet {background:#eaf0f3;color:#647781}
    .stage-number {font-size:12px;font-weight:700;color:#087f8c;letter-spacing:.08em;margin-bottom:10px;}
    .stage-title {font-size:14px;font-weight:650;line-height:1.35;min-height:39px;margin-bottom:12px;}
    .fixture-banner {border:1px solid #e7bd71;background:#fff4db;color:#7a5212;padding:11px 16px;border-radius:9px;font-size:13px;margin-bottom:18px;}
    .intake-banner {background:#e6f3f3;border:1px solid #badbdc;border-radius:12px;padding:22px;margin:12px 0 24px;}
    .flow {font-size:13px;color:#496b78;background:#edf3f5;padding:12px 14px;border-radius:8px;margin:10px 0 18px;}
    button[kind="primary"] {background:#087f8c;}
    </style>""", unsafe_allow_html=True)


def _sidebar(st: Any, service: Any) -> str | None:
    with st.sidebar:
        name = os.environ.get("BIOCENTRIC_DISPLAY_NAME", "Bio Orchestrator")
        st.markdown(f'<div class="brand"><span class="brand-mark">◈</span>{html.escape(name)}</div>', unsafe_allow_html=True)
        st.caption("EVIDENCE → COMPUTATION → EXPERIMENT")
        st.divider()
        projects = service.list_projects()
        choices = {p["id"]: str(p.get("spec", {}).get("name") or p.get("name") or p["id"]) for p in projects}
        selected = st.session_state.get("project_id")
        if choices:
            if selected not in choices:
                selected = next(iter(choices))
            project_id = st.selectbox("Project", list(choices), index=list(choices).index(selected),
                                      format_func=lambda key: choices[key], key="project_picker")
            if project_id != selected:
                st.session_state["selected_stage"] = None
            st.session_state["project_id"] = project_id
        else:
            project_id = None
        with st.expander("＋ New research project", expanded=not choices):
            with st.form("new_project", clear_on_submit=False):
                name = st.text_input("Project name", placeholder="e.g. EGFR exploratory screen")
                target = st.text_area("Protein, accession, or sequence", placeholder="P00533 or a protein sequence", height=90)
                organism = st.text_input("Organism", value="Homo sapiens")
                context = st.text_area("Biological context", placeholder="What scientific question should this project address?", height=75)
                effect = st.selectbox("Intended effect", ["Unresolved", "Inhibition", "Activation", "Binding characterization"])
                isoform = st.text_input("Isoform / construct", placeholder="Unknown if unresolved")
                mutations = st.text_input("Mutations", placeholder="e.g. L858R; leave blank if unknown")
                shortlist = st.number_input("Experimental shortlist size", 1, 1000, 10)
                st.markdown("**Resources & data handling**")
                cpu = st.number_input("CPU cores", 1, 16, 2)
                has_gpu = st.checkbox("GPU available")
                commercial_tools = st.text_input("Available licensed tools", placeholder="Comma-separated names, or leave blank",
                                                help="Records availability only; a declared license does not install or validate an adapter.")
                budget = st.number_input("Compute budget (seconds)", 30, 86400, 300, 30)
                seed = st.number_input("Random seed", 1, 2147483647, 42)
                local_only = st.checkbox("Keep project inputs local", value=True,
                                         help="External target research and remote model calls are blocked in local-only mode.")
                allow_model = st.checkbox("Allow configured language model", value=False)
                model = st.text_input("Model identifier", value=os.environ.get("OPENAI_MODEL", "gpt-4.1-mini"))
                max_model_calls = st.number_input("Maximum model calls", 0, 10, 3,
                                                 help="Persisted authorization limit for this project.")
                max_model_output_tokens = st.number_input("Maximum output tokens per call", 200, 4000, 1500, 100)
                workflow = st.selectbox("Workflow", ["noncovalent", "covalent", "metalloprotein", "peptide", "other"],
                                        help="The first release executes only supported noncovalent small-molecule workflows.")
                create = st.form_submit_button("Create project", type="primary", width="stretch")
            if create:
                if not target.strip() or not organism.strip():
                    st.error("Provide a target and organism so identity can be resolved.")
                elif local_only and allow_model:
                    st.error("A remote model cannot be enabled while keeping inputs local.")
                else:
                    spec = {"name": name.strip() or f"{target.strip()[:28]} research", "target": target.strip(),
                            "organism": organism.strip(), "context": context.strip(), "intended_effect": None if effect == "Unresolved" else effect.lower(),
                            "isoform": isoform.strip() or None, "mutations": [x.strip() for x in mutations.split(",") if x.strip()],
                            "shortlist_size": int(shortlist), "cpu": int(cpu), "has_gpu": has_gpu,
                            "commercial_tools": [x.strip() for x in commercial_tools.split(",") if x.strip()],
                            "compute_budget_seconds": int(budget), "seed": int(seed), "local_only": local_only,
                            "allow_model": allow_model, "model": model.strip(), "max_model_calls": int(max_model_calls),
                            "max_model_output_tokens": int(max_model_output_tokens), "workflow_kind": workflow}
                    created = _attempt(st, lambda: service.create_project(spec), "Project created.")
                    if isinstance(created, str):
                        st.session_state["project_id"] = created
                        st.session_state.pop("project_picker", None)
                        st.rerun()
        st.divider()
        st.markdown("**Explore the workbench**")
        st.caption("An offline fixture demonstrates evidence, gates, events, and candidate inspection. It is not a calculation.")
        if st.button("Load offline demonstration", width="stretch"):
            created = _attempt(st, service.create_demo, "Offline fixture loaded.")
            if isinstance(created, str):
                st.session_state["project_id"] = created
                st.session_state.pop("project_picker", None)
                st.rerun()
        example_available = callable(getattr(service, "create_example", None))
        st.caption("The executable example uses curated inputs and actual local calculations. It verifies software operation; predictive performance remains unvalidated.")
        if st.button("Load executable smoke example", disabled=not example_available,
                     width="stretch", help=None if example_available else "The executable example is not available in this installation."):
            created = _attempt(st, service.create_example, "Executable example loaded.")
            if isinstance(created, str):
                st.session_state["project_id"] = created
                st.session_state.pop("project_picker", None)
                st.rerun()
        key_available = bool(os.environ.get("OPENAI_API_KEY"))
        st.caption("For GPT advice, load the public GPT example. Run pipeline then allows up to 3 billable model calls at Stage 03; scientific review gates remain enabled.")
        if st.button("Load GPT-enabled public example", disabled=not (example_available and key_available),
                     width="stretch", help="Uses public 1IEP inputs. Creating it does not make API calls."):
            created = _attempt(st, lambda: service.create_example(allow_model=True, model=os.environ.get("OPENAI_MODEL", "gpt-4.1-mini")), "GPT-enabled example loaded. Run pipeline to request advice.")
            if isinstance(created, str):
                st.session_state["project_id"] = created
                st.session_state.pop("project_picker", None)
                st.rerun()
        if not key_available:
            st.caption("To enable this option, start the app with the private API-key launcher.")
        st.divider()
        st.toggle("Refresh live progress", value=True, key="live_refresh",
                  help="Pause refresh while rotating a 3D structure so the camera view stays fixed.")
        st.caption("Scientific decisions and provenance persist locally. A computational shortlist is a set of hypotheses for testing.")
    return project_id


def _empty(st: Any) -> None:
    st.markdown('<div class="eyebrow">A scientific workflow you can inspect</div>', unsafe_allow_html=True)
    st.title("From target to testable hypotheses.")
    st.markdown("Build an evidence dossier, choose a supported protocol, and follow every transformation from input to result.")
    st.markdown('<div class="intake-banner"><b>Start with a scientific question.</b><br>Create a project in the sidebar, or open the clearly labeled offline demonstration to inspect the complete workbench.</div>', unsafe_allow_html=True)
    columns = st.columns(3)
    for col, title, text in zip(columns,
                              ["01 / Evidence with context", "02 / Visible execution", "03 / Reviewable handoff"],
                              ["Inspect source records, target identity, conflicting evidence, and coverage gaps.",
                               "Follow ten stages, preprocessing, model and tool calls, scientific gates, and durable jobs.",
                               "Compare candidate tradeoffs and export records with their assumptions and provenance."]):
        with col, st.container(border=True):
            st.markdown(f"**{title}**")
            st.write(text)
    st.subheader("The pipeline")
    st.markdown(" → ".join(STAGES))


def _controls(st: Any, service: Any, project_id: str, snapshot: dict[str, Any]) -> None:
    summary = pipeline_summary(snapshot)
    state = str(snapshot.get("project", {}).get("state", "draft"))
    cols = st.columns([1, 1, 1, 1, 3])
    if cols[0].button("▶ Run pipeline", type="primary", disabled=summary["active"] or is_fixture(snapshot), width="stretch"):
        if _attempt(st, lambda: service.start(project_id), "Worker started."):
            st.rerun()
    if cols[1].button("↻ Resume", disabled=summary["active"] or is_fixture(snapshot), width="stretch"):
        if _attempt(st, lambda: service.resume(project_id), "Resuming from persisted progress."):
            st.rerun()
    if cols[2].button("■ Cancel", disabled=not summary["active"], width="stretch"):
        if _attempt(st, lambda: service.cancel(project_id), "Cancellation requested."):
            st.rerun()
    if cols[3].button("Refresh", width="stretch"):
        st.rerun()
    with cols[4]:
        refresh_text = ("Updated automatically every 3 seconds" if summary["active"] else "Workflow paused or finished · view stays open") if st.session_state.get("live_refresh", True) else "Live refresh paused · use Refresh for latest state"
        st.markdown(f'{_badge(state)} <span class="muted">{refresh_text}</span>', unsafe_allow_html=True)


def _workflow(st: Any, snapshot: dict[str, Any]) -> None:
    steps = snapshot.get("steps", [])
    if not steps:
        steps = [{"id": str(i + 1), "label": label, "status": "pending"} for i, label in enumerate(STAGES)]
    summary = pipeline_summary(snapshot)
    st.markdown('<div class="eyebrow">Every stage. Every transformation.</div>', unsafe_allow_html=True)
    st.progress(summary["complete"] / max(summary["total"], 1), text=f'{summary["complete"]} of {summary["total"]} stages completed')
    current = summary["current"]
    selected = st.session_state.get("selected_stage") or current.get("id") or steps[0].get("id")
    for offset in range(0, len(steps), 5):
        columns = st.columns(5, gap="small")
        for i, (column, step) in enumerate(zip(columns, steps[offset:offset + 5]), start=offset + 1):
            with column, st.container(border=True):
                st.markdown(f'<div class="stage-number">STAGE {i:02d}</div><div class="stage-title">{html.escape(str(step.get("label", STAGES[min(i - 1, 9)])))}</div>{_badge(step.get("status"))}', unsafe_allow_html=True)
                if st.button("Inspect stage", key=f"stage_{step.get('id', i)}", width="stretch",
                             type="primary" if str(step.get("id")) == str(selected) else "secondary"):
                    st.session_state["selected_stage"] = step.get("id")
                    selected = step.get("id")
    detail = next((s for s in steps if str(s.get("id")) == str(selected)), steps[0])
    with st.container(border=True):
        st.subheader(detail.get("label", "Stage details"))
        explanation = stage_explanation(snapshot, detail)
        st.markdown("**Decision summary**")
        st.caption(explanation["actor"])
        st.write(explanation["reason"])
        st.write(f"Validation rule: {explanation['rule']}")
        for note in explanation["details"]:
            st.write(f"• {note}")
        if detail.get("id") == "03":
            activity = model_activity(snapshot, key_available=bool(os.environ.get("OPENAI_API_KEY")))
            st.info(f"GPT status: {activity['state']} · {activity['calls']} API requests started. {activity['reason']}")
            _model_advice(st, activity["advice"], snapshot)
            _model_advice(st, activity["rejected_advice"], snapshot, accepted=False)
        st.markdown('<div class="flow">Typed inputs → Preprocessing → Tool / model call → Parsed outputs → Postprocessing & validation</div>', unsafe_allow_html=True)
        for title, value in [
            ("1 · Input contract & artifacts", detail.get("inputs")),
            ("2 · Preprocessing & assumptions", detail.get("preprocessing")),
            ("3 · Tool or model selection", detail.get("tool")),
            ("4 · Outputs & evidence", detail.get("outputs")),
            ("5 · Postprocessing & validation", detail.get("postprocessing")),
        ]:
            with st.expander(title, expanded=title.startswith("2") or title.startswith("5")):
                _json(st, value)
        stage_events = [e for e in snapshot.get("events", []) if str(e.get("stage")) == str(detail.get("id"))]
        if stage_events:
            _events(st, stage_events, "stage_events", heading="Calls and events for this stage")


def _model_advice(st: Any, advice: dict[str, Any], snapshot: dict[str, Any], *, accepted: bool = True) -> None:
    if not advice.get("summary"):
        return
    st.markdown("**GPT explanation / advisory proposal**" if accepted else "**Rejected GPT answer — review only**")
    if not accepted:
        st.warning("This answer failed citation validation and is not accepted scientific advice. It is preserved here so you can inspect what the model actually returned.")
    st.write(advice["summary"])
    st.write(f"Suggested strategy: {advice.get('proposed_strategy', 'not recorded').replace('_', ' ')}")
    st.caption("This is the model’s public answer. It does not change the approved workflow route or authorize scientific jobs.")
    evidence = {record['id']: record for record in snapshot.get('evidence', [])}
    for claim in advice.get("claims", []):
        st.write(claim.get("text", ""))
        sources = [f"{identifier} ({evidence.get(identifier, {}).get('source', 'source not recorded')})" for identifier in claim.get('evidence_ids', [])]
        st.caption(f"{claim.get('category', 'claim').replace('_', ' ')} · Evidence: {', '.join(sources) or 'none cited'}")
    for title, field in [("Open questions", "unresolved_questions"), ("Limitations", "limitations")]:
        if advice.get(field):
            st.markdown(f"**{title}**")
            for item in advice[field]:
                st.write(f"• {item}")
    if advice.get("validation"):
        st.caption(advice["validation"])


def _model_panel(st: Any, service: Any, project_id: str, snapshot: dict[str, Any], activity: dict[str, Any]) -> None:
    st.subheader("Model activity & explanation")
    st.write(activity["reason"])
    st.caption("GPT is an optional advisor at Stage 03. Other stages use scientific tools and explicit workflow rules. Public answers, citations, and tool calls are shown here; hidden internal reasoning is not exposed.")
    st.write(f"Configured model: {activity['model']}")
    cols = st.columns(3)
    cols[0].metric("API requests started", f"{activity['calls']} / {activity['budget']}")
    cols[1].metric("Input tokens", activity["input_tokens"] if activity["input_tokens"] is not None else "—")
    cols[2].metric("Output tokens", activity["output_tokens"] if activity["output_tokens"] is not None else "—")
    st.caption("Counts belong to this project. The launcher’s CONNECTED test is separate and produces no research advice. A configured key does not enable GPT on an existing local-only project.")
    if activity["state"] == "disabled":
        st.info("For the public example, choose Load GPT-enabled public example in the sidebar, then Run pipeline. For your own target, create a project with Keep project inputs local off and Allow configured language model on.")
    if activity["state"] == "failed" and activity["calls"] < activity["budget"] and not snapshot.get("data", {}).get("model_repair"):
        st.caption("One correction can use evidence already opened, within the existing project budget. This makes one billable model call and preserves the rejected answer in the timeline.")
        if st.button("Correct rejected advice · 1 model call", disabled=pipeline_summary(snapshot)["active"] or not activity["key_available"]):
            with st.spinner("Requesting an evidence-grounded correction…"):
                result = _attempt(st, lambda: service.repair_model_advice(project_id), "Corrected advice saved.")
            if result:
                st.rerun()
    _model_advice(st, activity["advice"], snapshot)
    _model_advice(st, activity["rejected_advice"], snapshot, accepted=False)
    events = [event for event in snapshot.get("events", []) if event.get("kind") in {"model_call", "model_status", "model_validation"} or event.get("name") == "read_evidence"]
    _events(st, events, "model_events", heading="Model requests, responses & evidence reads")


def _events(st: Any, events: list[dict[str, Any]], key: str, heading: str = "Execution timeline") -> None:
    st.subheader(heading)
    if not events:
        st.caption("No calls have been made yet. Each call will include its recorded request, response, timing, and status.")
        return
    kinds = sorted({str(e.get("kind", "event")) for e in events})
    col1, col2 = st.columns([2, 1])
    kind = col1.selectbox("Show events", ["All event types", *kinds], key=f"{key}_kind")
    count = col2.selectbox("Recent events", [20, 50, 100, 500], key=f"{key}_count")
    filtered = [e for e in events if kind == "All event types" or str(e.get("kind", "event")) == kind]
    st.caption(f"{len(filtered)} matching events · newest first · expand any event to inspect the recorded data")
    for event in reversed(filtered[-count:]):
        title = f"{event.get('created_at', '')} · {event.get('kind', 'event')} / {event.get('name', event.get('stage', 'pipeline'))} · {event.get('status', 'recorded')}"
        with st.expander(title):
            if event.get("message"):
                st.write(event["message"])
            st.caption(f"Event {event.get('id', 'unknown')} · stage {event.get('stage', '—')} · duration {event.get('duration_ms', 'unknown')} ms")
            left, right = st.columns(2)
            with left:
                st.markdown("**Recorded request / inputs**")
                _json(st, event.get("inputs", event.get("request")))
            with right:
                st.markdown("**Recorded response / outputs**")
                _json(st, event.get("outputs", event.get("response")))
            st.markdown("**Full event record**")
            _json(st, event)


def _evidence(st: Any, snapshot: dict[str, Any]) -> None:
    data = snapshot.get("data", {})
    st.subheader("Target identity")
    identity = data.get("target_identity", data.get("identity", data.get("target")))
    _json(st, identity, "Target identity has not been resolved. Species, accession, isoform, and construct remain unverified.")
    evidence = snapshot.get("evidence", [])
    st.subheader("Evidence library")
    st.caption("Source observations, model predictions, and computed results retain their own provenance. An empty result is an evidence gap.")
    if not evidence:
        st.info("No evidence records yet. Run the target dossier stage or import supporting files.")
    else:
        left, right = st.columns([2, 1])
        query = left.text_input("Search evidence", placeholder="Source, title, accession, finding…")
        source_options = sorted({str(e.get("source", e.get("source_name", "Unknown"))) for e in evidence})
        source = right.selectbox("Source", ["All sources", *source_options])
        selected = [e for e in evidence if (not query or query.lower() in json.dumps(e, default=str).lower())
                    and (source == "All sources" or str(e.get("source", e.get("source_name", "Unknown"))) == source)]
        st.caption(f"{len(selected)} / {len(evidence)} evidence records")
        for record in selected:
            title = record.get("title") or record.get("claim") or record.get("summary") or record.get("source_id") or record.get("id", "Evidence")
            with st.expander(f"[{record.get('id', 'record')}] {str(title)[:180]}"):
                text = record.get("supporting_text") or record.get("summary") or record.get("text")
                if text:
                    st.write(text)
                url = record.get("url", record.get("source_url"))
                if isinstance(url, str) and url.startswith(("https://", "http://")):
                    st.link_button("Open source", url)
                _json(st, record)
    for title, keys in [("Assay measurements & controls", ["assays", "measurements", "known_chemistry"]),
                        ("Structure candidates", ["structures", "structure_candidates"]),
                        ("Search log & coverage limits", ["search_log", "coverage", "coverage_limits"])]:
        for field in keys:
            if data.get(field):
                with st.expander(title):
                    _json(st, data[field])
                break


def _decisions(st: Any, service: Any, project_id: str, snapshot: dict[str, Any]) -> None:
    st.subheader("Scientific decisions")
    st.caption("Decisions are persisted with your rationale. Approval cannot bypass missing inputs, unsupported chemistry, or failed validation.")
    data = snapshot.get("data", {})
    if isinstance(data.get("gate"), dict) and data["gate"].get("kind") == "identity":
        spec = snapshot.get("project", {}).get("spec", {})
        initial_identity = {"accession": None, "name": spec.get("target", "") if len(spec.get("target", "")) < 160 else "",
                            "organism": spec.get("organism", ""), "sequence": None,
                            "isoform": spec.get("isoform"), "mutations": spec.get("mutations", []),
                            "construct": None, "domain": None}
        with st.form("local_identity"):
            st.markdown("**Record a user-supplied target identity**")
            st.caption("For local or proprietary targets, supply the identity directly. This record is your reviewed input; it is not a remote database verification.")
            identity_text = st.text_area("TargetIdentity JSON", value=json.dumps(initial_identity, indent=2), height=270)
            identity_reason = st.text_input("Identity provenance and rationale", placeholder="Source of the sequence, construct, and target mapping")
            if st.form_submit_button("Record local identity", disabled=is_fixture(snapshot)):
                if not identity_reason.strip():
                    st.error("Record the provenance and rationale for this user-supplied target identity.")
                else:
                    try:
                        identity_value = json.loads(identity_text)
                        if not isinstance(identity_value, dict):
                            raise ValueError("TargetIdentity must be a JSON object.")
                    except (ValueError, json.JSONDecodeError) as exc:
                        st.error(f"Invalid target identity: {exc}")
                    else:
                        if _attempt(st, lambda: service.decide(project_id, "identity", identity_value, identity_reason)):
                            st.rerun()
    candidates = data.get("target_candidates", data.get("identity_candidates", []))
    if candidates:
        choices = {str(c.get("accession", c.get("id", ""))): c for c in candidates}
        with st.form("target_selection"):
            selection = st.selectbox("Resolve target identity", list(choices),
                                     format_func=lambda value: f"{value} · {choices[value].get('name', '')} · {choices[value].get('organism', '')}")
            _json(st, choices.get(selection))
            reason = st.text_input("Selection rationale", key="identity_reason")
            if st.form_submit_button("Record target selection", disabled=is_fixture(snapshot)):
                if not reason.strip():
                    st.error("Record why this target is the intended identity.")
                elif _attempt(st, lambda: service.decide(project_id, "target_selection", {"accession": selection}, reason)):
                    st.rerun()
    with st.expander("Review protocol / exploratory execution", expanded=any(s.get("status") == "blocked" for s in snapshot.get("steps", []))):
        _json(st, data.get("gate", data.get("benchmark", data.get("benchmark_result"))), "The worker will record the applicable gate here.")
        proposal = data.get("protocol", data.get("protocol_proposal", data.get("plan")))
        if proposal:
            st.markdown("**Recorded proposal**")
            _json(st, proposal)
        with st.form("gate_decision"):
            kind = st.selectbox("Decision", ["protocol", "exploratory"], format_func=lambda v: "Approve proposed protocol" if v == "protocol" else "Allow explicitly exploratory execution")
            st.caption("Exploratory execution means target-specific predictive performance is not established. Resulting reports must retain that designation.")
            approved = st.checkbox("I have reviewed the recorded inputs, limitations, and proposed configuration")
            receptor_reviewed = st.checkbox("I reviewed the receptor components, binding pocket, and supported ligand chemistry")
            ph_assumption = st.number_input("Preparation pH assumption", min_value=0.0, max_value=14.0, value=7.4, step=0.1,
                                            help="A recorded scientific assumption. Setting this value does not run a protonation calculation.")
            rationale = st.text_area("Scientific rationale", placeholder="Why is this protocol or exploratory scope appropriate?")
            if st.form_submit_button("Record decision", disabled=is_fixture(snapshot)):
                if not approved or not rationale.strip():
                    st.error("Review the proposal and record a rationale before approval.")
                elif kind == "protocol" and not receptor_reviewed:
                    st.error("Protocol approval requires explicit review of receptor components, the pocket, and supported chemistry.")
                elif _attempt(st, lambda: service.decide(project_id, kind,
                                                       {"approved": True, "receptor_reviewed": receptor_reviewed,
                                                        "ph_assumption": ph_assumption}, rationale)):
                    st.rerun()
    decisions = snapshot.get("decisions", [])
    with st.expander(f"Decision history ({len(decisions)})"):
        _json(st, decisions)


def _plan_runs(st: Any, service: Any, project_id: str, snapshot: dict[str, Any]) -> None:
    data = snapshot.get("data", {})
    st.subheader("Plan & assumptions")
    _json(st, data.get("plan", data.get("strategy")), "The plan will be recorded after target evidence and available inputs are assessed.")
    with st.expander("Project inputs and authorized resources"):
        _json(st, snapshot.get("project", {}).get("spec"))
    _decisions(st, service, project_id, snapshot)
    st.subheader("Durable jobs")
    jobs = snapshot.get("jobs", [])
    if jobs:
        st.dataframe(table_rows(jobs), width="stretch", hide_index=True)
    else:
        st.caption("No jobs submitted. Starting the pipeline creates a separately executing worker job.")
    _events(st, snapshot.get("events", []), "all_events")
    with st.expander("Tool availability & validation status"):
        try:
            @st.cache_data(ttl=30, show_spinner=False)
            def cached_adapters(_service: Any, project_key: str) -> Any:
                return _service.adapters()
            adapters = cached_adapters(service, project_id)
        except Exception as exc:
            st.error(f"Cannot read adapter status: {exc}")
        else:
            if isinstance(adapters, dict):
                _json(st, adapters)
            else:
                st.dataframe(table_rows(adapters), width="stretch", hide_index=True)


def _structure_viewer(st: Any, snapshot: dict[str, Any], root: Path, *, record: dict | None = None, key: str) -> None:
    from .viewer import geometry_artifacts, resolve_scene, viewer_html

    with st.expander("Interactive 3D geometry · local viewer"):
        st.caption("Rotate and zoom the actual registered geometry. Receptor cartoon and ligand sticks show coordinates; a plausible pose does not confirm binding. Pause live refresh in the sidebar to keep the camera fixed.")
        receptors, ligands = geometry_artifacts(snapshot)
        if not receptors:
            st.info("Register a PDB or mmCIF structure to inspect receptor geometry. PDBQT alone is not used as the canonical structural display.")
            return
        if not st.checkbox("Load local 3D geometry", key=f"{key}_enabled"):
            return
        protocol = snapshot.get("data", {}).get("protocol") or {}
        expected_structure = (protocol.get("input_hashes") or {}).get("structure")
        receptor_ids = [a["id"] for a in receptors]
        default_index = next((i for i, a in enumerate(receptors) if a["sha256"] == expected_structure), 0)
        lookup = {a["id"]: a for a in [*receptors, *ligands]}
        receptor_id = st.selectbox("Registered receptor", receptor_ids, index=default_index,
                                   format_func=lambda value: f"{lookup[value]['name']} · {value[:10]}", key=f"{key}_receptor")
        if record is not None:
            pose_ids = [a["id"] for a in ligands if a["id"] in record.get("artifact_ids", []) and a["sha256"] == record.get("pose_sha256")]
            if not pose_ids:
                st.error("No hash-matched registered SDF is linked to this candidate. A raw pose path is not sufficient for display.")
                return
            ligand_id = st.selectbox("Registered candidate pose", pose_ids,
                                     format_func=lambda value: f"{lookup[value]['name']} · {value[:10]}", key=f"{key}_pose")
            manual_overlay = False
        else:
            ligand_id = st.selectbox("Optional registered ligand geometry", [None, *[a["id"] for a in ligands]],
                                     format_func=lambda value: "Receptor only" if value is None else f"{lookup[value]['name']} · {value[:10]}", key=f"{key}_ligand")
            manual_overlay = st.checkbox("I verified that this structure and ligand share the same coordinate frame", key=f"{key}_frame_confirm") if ligand_id else False
        frame = (protocol.get("box") or {}).get("coordinate_frame")
        declared_frame = ""
        if frame:
            st.caption(f"Recorded coordinate frame: {frame}")
        else:
            declared_frame = st.text_input("Coordinate frame / source", placeholder="e.g. untransformed coordinates from the uploaded PDB assembly", key=f"{key}_frame")
            st.caption("This is a user-supplied display annotation. It does not validate a docking protocol or align coordinates.")
        try:
            scene = resolve_scene(snapshot, root, receptor_id, ligand_id, record, declared_frame, manual_overlay)
            document = viewer_html(scene)
        except (ValueError, OSError, UnicodeDecodeError) as exc:
            st.error(str(exc))
            return
        import streamlit.components.v1 as components
        components.html(document, height=550, scrolling=False)
        st.caption("No structures leave this application. The pinned viewer is bundled locally; its frame blocks network connections.")
        _json(st, {"coordinate_frame": scene["coordinate_frame"], "frame_source": scene["frame_source"],
                   "artifacts": [{k: v for k, v in model.items() if k != "text"} for model in scene["models"]]})


def _candidates(st: Any, snapshot: dict[str, Any], root: Path) -> None:
    st.subheader("Candidate hypotheses")
    st.caption("Docking scores are protocol-specific ranking signals. Physicochemical descriptors are not measured ADMET. Missing selectivity or safety evidence remains unknown.")
    candidates = snapshot.get("candidates", [])
    if not candidates:
        st.info("No candidate assessments yet. Preparation, protocol review, and scientific gates precede screening and prioritization.")
        return
    left, middle, right = st.columns([2, 1, 1])
    query = left.text_input("Find candidate", placeholder="Compound identity, SMILES, or assessment…")
    included_only = middle.checkbox("Shortlist only")
    flags_only = right.checkbox("Needs evidence / review")
    selected = []
    for candidate in candidates:
        text = json.dumps(candidate, default=str).lower()
        included = candidate.get("selected", candidate.get("included", candidate.get("shortlisted", False)))
        flags = candidate.get("flags", candidate.get("review_flags", candidate.get("alerts", []))) or candidate.get("missing_evidence")
        if (not query or query.lower() in text) and (not included_only or included) and (not flags_only or flags):
            selected.append(candidate)
    st.caption(f"{len(selected)} shown / {len(candidates)} assessed · inspect a compound below for full provenance")
    st.dataframe(table_rows(selected), width="stretch", hide_index=True)
    if selected:
        item = st.selectbox("Inspect compound", range(len(selected)),
                            format_func=lambda i: str(selected[i].get("name", selected[i].get("compound_id", selected[i].get("id", f"Compound {i + 1}")))))
        record = selected[item]
        with st.container(border=True):
            smiles = record.get("canonical_smiles", record.get("smiles", record.get("metrics", {}).get("smiles")))
            if smiles:
                try:
                    from rdkit import Chem
                    from rdkit.Chem import Draw
                    molecule = Chem.MolFromSmiles(str(smiles))
                    if molecule is not None:
                        st.image(Draw.MolToImage(molecule, size=(580, 240)), caption="2D chemical identity · no binding or pose claim", width=580)
                except ImportError:
                    st.caption("Install RDKit to preview this compound's 2D chemical identity.")
            _json(st, record)
        _structure_viewer(st, snapshot, root, record=record, key=f"candidate_{record.get('compound_id', item)}")
        st.download_button("Download filtered candidate CSV", candidate_csv(selected), file_name="candidate_hypotheses.csv", mime="text/csv")
    st.download_button("Download all candidate records", json.dumps(safe_data(candidates), indent=2, default=str), file_name="candidate_records.json", mime="application/json")


def _files(st: Any, service: Any, project_id: str, snapshot: dict[str, Any], root: Path) -> None:
    st.subheader("Project inputs")
    st.caption("Original files are preserved. The worker validates the declared format and records each preprocessing and postprocessing step.")
    labels = {"library": "Ligand library · SDF / SMILES", "receptor": "Curated receptor · PDBQT",
              "reference": "Reference ligand · SDF", "structure": "Structure · PDB / mmCIF",
              "identity": "User-supplied target identity · JSON", "evidence": "Evidence records with provenance · JSON"}
    with st.form("file_intake", clear_on_submit=True):
        kind = st.selectbox("Input role", list(labels), format_func=lambda value: labels[value])
        uploaded = st.file_uploader("Upload scientific input", type=["sdf", "smi", "smiles", "pdb", "cif", "mmcif", "pdbqt", "json"])
        st.caption("Use curated, prepared receptor inputs for execution. Raw structure import does not assert that repair or protonation has been validated.")
        if st.form_submit_button("Preserve and register input", disabled=is_fixture(snapshot)):
            if uploaded is None:
                st.error("Select a file first.")
            elif _attempt(st, lambda: service.upload(project_id, uploaded.name, uploaded.getvalue(), kind), "Original input registered."):
                st.rerun()
    _structure_viewer(st, snapshot, root, key="file_geometry")
    st.subheader("Artifacts & reports")
    if st.button("Generate reproducible report bundle", disabled=is_fixture(snapshot)):
        exported = _attempt(st, lambda: service.export(project_id), "Report bundle generated.")
        if isinstance(exported, dict):
            st.session_state[f"exports_{project_id}"] = exported
            st.rerun()
    artifacts = snapshot.get("artifacts", [])
    if not artifacts:
        st.caption("Artifacts appear here as inputs are preserved and jobs produce outputs.")
    for i, artifact in enumerate(artifacts):
        with st.expander(str(artifact.get("name", artifact.get("id", "Artifact")))):
            st.caption(f"Type: {artifact.get('media_type', 'unknown')} · SHA-256: {artifact.get('sha256', 'unknown')}")
            content = artifact_bytes(artifact, root)
            if content is not None:
                st.download_button("Download artifact", content,
                                   file_name=Path(str(artifact.get("name", "artifact"))).name,
                                   mime=artifact.get("media_type", "application/octet-stream"), key=f"artifact_{i}_{artifact.get('id', '')}")
                if artifact.get("media_type", "").startswith("text/") and len(content) < 150_000:
                    st.code(content.decode("utf-8", errors="replace"), language="text")
            else:
                st.caption("This artifact is unavailable for download or exceeds the 50 MB preview limit.")
            _json(st, artifact)
    exported = st.session_state.get(f"exports_{project_id}")
    if exported:
        with st.expander("Export manifest"):
            _json(st, exported)
            for key, value in exported.items():
                if isinstance(value, str):
                    content = artifact_bytes({"path": value}, root)
                    if content is not None:
                        st.download_button(f"Download {key}", content, file_name=Path(value).name, key=f"export_{key}")
    with st.expander("Import measured experimental outcomes"):
        st.caption("Imported measurements remain assay-specific observations. Include compound identity, endpoint, value, units, relation, and assay metadata.")
        st.code("compound_id,endpoint,value,units,relation,assay_id,organism,construct,assay_format,conditions,publication", language="text")
        st.caption("The first seven columns are required. Keep endpoint types such as Kd, Ki, IC50, and EC50 distinct; retain relation operators and assay context.")
        st.download_button("Download measurement CSV template",
                           "compound_id,endpoint,value,units,relation,assay_id,organism,construct,assay_format,conditions,publication\n",
                           file_name="measurement_import_template.csv", mime="text/csv")
        with st.form("measurements", clear_on_submit=True):
            measurement_file = st.file_uploader("Measured results CSV", type=["csv"], key="measurement_upload")
            if st.form_submit_button("Import measurements", disabled=is_fixture(snapshot)):
                if measurement_file is None:
                    st.error("Select a measurements CSV first.")
                elif _attempt(st, lambda: service.import_measurements(project_id, measurement_file.getvalue()), "Experimental observations imported."):
                    st.rerun()


def run_app(*, configure_page: bool = True) -> None:
    import streamlit as st
    from biocentric.service import Service

    if configure_page:
        st.set_page_config(page_title=os.environ.get("BIOCENTRIC_DISPLAY_NAME", "Bio Orchestrator"), page_icon="◈", layout="wide", initial_sidebar_state="expanded")
    _styles(st)
    root = Path(os.environ.get("BIOCENTRIC_HOME", ".biocentric")).expanduser().resolve()

    @st.cache_resource
    def get_service(home: str, implementation_version: int) -> Any:
        """Cache connection configuration; project state is always read fresh."""
        return Service(home)

    service = get_service(str(root), Path(__file__).with_name("service.py").stat().st_mtime_ns)
    project_id = _sidebar(st, service)
    if project_id is None:
        _empty(st)
        return

    refreshing = st.session_state.get("live_refresh", True) and pipeline_summary(service.snapshot(project_id))["active"]

    @st.fragment(run_every=3 if refreshing else None)
    def workspace() -> None:
        try:
            snapshot = service.snapshot(project_id)
        except Exception as exc:
            st.error(f"Cannot load project: {exc}")
            return
        if refreshing and not pipeline_summary(snapshot)["active"]:
            # Stop the timer once the worker reaches a gate or finishes. Otherwise
            # Streamlit 1.50 can reset its client-side tabs every three seconds.
            st.rerun()
        project = snapshot.get("project", {})
        spec = project.get("spec", {})
        fixture = is_fixture(snapshot)
        if fixture:
            st.markdown('<div class="fixture-banner"><b>OFFLINE DEMONSTRATION · FIXTURE DATA</b> — This named example illustrates the interface. No displayed fixture value represents a calculation or experimentally confirmed result.</div>', unsafe_allow_html=True)
        st.markdown('<div class="eyebrow">Discovery workspace / research use</div>', unsafe_allow_html=True)
        st.title(spec.get("name", project.get("name", "Research project")))
        target = spec.get("target", "Unresolved target")
        target = target[:65] + "…" if len(str(target)) > 65 else target
        privacy = "Local inputs only" if spec.get("local_only", True) else "External research enabled"
        st.caption(f"{target} · {spec.get('organism', 'Organism unresolved')} · {privacy} · Project {project_id}")
        activity = model_activity(snapshot, key_available=bool(os.environ.get("OPENAI_API_KEY")))
        model_message = f"LLM {activity['state']} · {activity['model']} · {activity['calls']} API requests started. {activity['reason']} Open Model activity for details."
        if activity["state"] in {"disabled", "failed", "unconfigured", "incomplete"}:
            st.warning(model_message)
        else:
            st.info(model_message)
        _controls(st, service, project_id, snapshot)
        summary = pipeline_summary(snapshot)
        columns = st.columns(4)
        columns[0].metric("Stages completed", f"{summary['complete']} / {summary['total']}")
        columns[1].metric("Evidence records", len(snapshot.get("evidence", [])))
        columns[2].metric("Recorded calls & events", len(snapshot.get("events", [])))
        columns[3].metric("Candidate assessments", len(snapshot.get("candidates", [])))
        current = summary["current"]
        if current.get("label"):
            current_text = f"Current stage: {current['label']}"
            if current.get("reason"):
                current_text += f" — {current['reason']}"
            if current.get("status") in {"blocked", "failed", "awaiting_approval"}:
                st.warning(current_text)
            elif current.get("status") == "running":
                st.info(current_text)
            else:
                st.caption(current_text)
        tabs = st.tabs(["Workflow", "Model activity", "Evidence", "Plan & Runs", "Candidates", "Files & handoff"])
        with tabs[0]:
            _workflow(st, snapshot)
        with tabs[1]:
            _model_panel(st, service, project_id, snapshot, activity)
        with tabs[2]:
            _evidence(st, snapshot)
        with tabs[3]:
            _plan_runs(st, service, project_id, snapshot)
        with tabs[4]:
            _candidates(st, snapshot, root)
        with tabs[5]:
            _files(st, service, project_id, snapshot, root)
        st.caption("Offline fixture — demonstration only." if fixture else "Computational results support hypotheses for experimental review; they do not establish biological efficacy.")

    workspace()
