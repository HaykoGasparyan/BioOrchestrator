"""Durable research worker; inline questions suspend and resume the same job."""
from __future__ import annotations

import argparse
import os
import threading
import time
import uuid

from .research_workflow import AwaitingAnswer, ResearchWorkflow
from .store import Store, safe
from .worker import Cancelled


def execute(root, pid, jid):
    store = Store(root)
    owner = uuid.uuid4().hex
    if not store.claim(jid, owner_token=owner, worker_pid=os.getpid(), project_id=pid):
        return False
    stop = threading.Event()

    def pulse():
        while not stop.wait(15):
            if not store.heartbeat(jid, owner):
                return

    heartbeat = threading.Thread(target=pulse, daemon=True)
    heartbeat.start()
    workflow = ResearchWorkflow(store, pid, jid, owner)
    store.update(pid, state="running")
    status, error = "completed", None
    try:
        workflow.run()
    except AwaitingAnswer as exc:
        status, error = "blocked", safe(str(exc))
        if store.owns_lease(jid, owner):
            store.research_update(pid, phase="awaiting_answer")
            store.step(pid, workflow.stage, "waiting", reason=error)
    except Cancelled as exc:
        status, error = "cancelled", safe(str(exc))
    except Exception as exc:
        status = "cancelled" if getattr(exc, "code", None) == "cancelled" else "failed"
        error = safe(f"{type(exc).__name__}: {exc}")
    finally:
        stop.set()
        heartbeat.join(timeout=1)
        if store.owns_lease(jid, owner):
            if status in {"failed", "cancelled"}:
                store.research_update(pid, phase=status)
                store.event(pid, kind="research_message", role="assistant", stage=workflow.stage,
                            status=status, message=f"{error}. Saved sources and results are retained.")
            if status == "failed":
                store.event(pid, kind="error", stage=workflow.stage, status=status, message=error)
                store.step(pid, workflow.stage, "failed", reason=error)
        finished = store.finish(jid, status, error, owner_token=owner, project_id=pid,
                                project_data={"compute_seconds": workflow.spent + time.monotonic() - workflow.started})
        if finished:
            from .reports import export_project
            try:
                export_project(store, pid)
            except Exception as exc:
                store.event(pid, kind="export", status="failed", message=safe(str(exc)))
    return True


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True)
    parser.add_argument("--project", required=True)
    parser.add_argument("--job", required=True)
    args = parser.parse_args()
    execute(args.root, args.project, args.job)


if __name__ == "__main__":
    main()
