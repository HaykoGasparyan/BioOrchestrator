"""SQLite durable jobs, immutable artifacts, append-only decisions and call traces."""
from __future__ import annotations
import hashlib
import json
import mimetypes
import re
import sqlite3
import time
import uuid
from contextlib import contextmanager
from pathlib import Path
from typing import Any
from .models import utc
from .host_identity import host_identity, process_identity
from .relocation import RelocationResolver

SECRET = re.compile(r"(?i)(api.?key|authorization|access.?token|refresh.?token|password|secret)")

def safe(value):
    if isinstance(value, dict):
        return {k: "[redacted]" if SECRET.search(str(k)) else safe(v) for k,v in value.items()}
    if isinstance(value, (list, tuple)):
        return [safe(x) for x in value]
    if isinstance(value, Path): return str(value)
    if isinstance(value, str): return re.sub(r"sk-[A-Za-z0-9_-]{16,}", "[redacted]", value)
    return value

def dump(value): return json.dumps(safe(value), ensure_ascii=False, sort_keys=True, allow_nan=False, default=str)
def digest(data: bytes) -> str: return hashlib.sha256(data).hexdigest()


class Store:
    def __init__(self, root: str | Path):
        self.root = Path(root).resolve()
        self.root.mkdir(parents=True, exist_ok=True)
        self.db = self.root / "metadata.sqlite"
        self.relocation = RelocationResolver.load(self.root)
        self.worker_host, self.worker_boot = host_identity()
        # Snapshot legacy database metadata before schema changes. Artifact bytes
        # are never modified by this migration; relocation also makes a full backup.
        if self.db.exists():
            with sqlite3.connect(f"file:{self.db}?mode=ro", uri=True) as source:
                columns = {row[1] for row in source.execute('PRAGMA table_info(jobs)')}
                if columns and not {'worker_host','worker_boot','worker_start'}.issubset(columns):
                    backup = self.root / 'metadata-backups' / uuid.uuid4().hex
                    backup.mkdir(parents=True)
                    with sqlite3.connect(backup / 'metadata.sqlite') as target:
                        source.backup(target)
                        target.row_factory = sqlite3.Row
                        (backup / 'artifacts.json').write_text(dump([dict(row) for row in target.execute('SELECT * FROM artifacts')]))
        with self.connect() as c:
            c.executescript('''
            PRAGMA journal_mode=WAL;
            CREATE TABLE IF NOT EXISTS projects(id TEXT PRIMARY KEY,spec TEXT NOT NULL,state TEXT NOT NULL,data TEXT NOT NULL,created_at TEXT NOT NULL);
            CREATE TABLE IF NOT EXISTS steps(project_id TEXT,id TEXT,status TEXT NOT NULL,data TEXT NOT NULL,PRIMARY KEY(project_id,id));
            CREATE TABLE IF NOT EXISTS jobs(id TEXT PRIMARY KEY,project_id TEXT UNIQUE,status TEXT NOT NULL,created_at TEXT NOT NULL,updated_at TEXT NOT NULL,attempts INTEGER DEFAULT 0,cancel_requested INTEGER DEFAULT 0,lease_until REAL,error TEXT,owner_token TEXT,worker_pid INTEGER,worker_host TEXT,worker_boot TEXT,worker_start TEXT);
            CREATE TABLE IF NOT EXISTS events(id INTEGER PRIMARY KEY AUTOINCREMENT,project_id TEXT,created_at TEXT,data TEXT);
            CREATE TABLE IF NOT EXISTS decisions(id INTEGER PRIMARY KEY AUTOINCREMENT,project_id TEXT,kind TEXT,value TEXT,reason TEXT,created_at TEXT);
            CREATE TABLE IF NOT EXISTS artifacts(id TEXT PRIMARY KEY,project_id TEXT,name TEXT,path TEXT,sha256 TEXT,media_type TEXT,kind TEXT,job_id TEXT,created_at TEXT);
            ''')
            c.execute('BEGIN IMMEDIATE')
            columns={row['name'] for row in c.execute('PRAGMA table_info(jobs)')}
            for name,kind in [('owner_token','TEXT'),('worker_pid','INTEGER'),('worker_host','TEXT'),('worker_boot','TEXT'),('worker_start','TEXT')]:
                if name not in columns: c.execute(f'ALTER TABLE jobs ADD COLUMN {name} {kind}')

    @contextmanager
    def connect(self):
        c = sqlite3.connect(self.db, timeout=30)
        c.row_factory = sqlite3.Row
        try:
            yield c
            c.commit()
        except Exception:
            c.rollback()
            raise
        finally: c.close()

    def directory(self, project_id: str) -> Path:
        if not re.fullmatch(r"[a-f0-9]{32}", project_id): raise ValueError("Invalid project identifier")
        p = self.root / "projects" / project_id
        p.mkdir(parents=True, exist_ok=True)
        return p

    def work_directory(self, project_id: str) -> Path:
        """Copied scientific checkpoints are evidence, never Linux cache entries."""
        path = self.directory(project_id) / 'work'
        if self.relocation.mappings:
            import platform
            scope = digest(f'{platform.system()}:{platform.machine()}:{self.worker_host}'.encode())[:16]
            path = path / f'host-{scope}'
        path.mkdir(parents=True, exist_ok=True)
        return path

    def create(self, spec: dict, steps: list[dict]) -> str:
        pid = uuid.uuid4().hex
        with self.connect() as c:
            c.execute("INSERT INTO projects VALUES (?,?,?,?,?)", (pid,dump(spec),"draft","{}",utc()))
            for step in steps:
                c.execute("INSERT INTO steps VALUES (?,?,?,?)", (pid,step['id'],"pending",dump(step)))
        self.artifact(pid, "project-spec.json", dump(spec).encode(), "immutable_input")
        self.event(pid, kind="project", name="Project created", status="completed", outputs=spec)
        return pid

    def project(self, pid):
        with self.connect() as c: row = c.execute("SELECT * FROM projects WHERE id=?", (pid,)).fetchone()
        if not row: raise ValueError("Unknown project")
        return {**dict(row), "spec":json.loads(row['spec']), "data":json.loads(row['data'])}

    def list_projects(self):
        with self.connect() as c: rows = c.execute("SELECT * FROM projects ORDER BY created_at DESC").fetchall()
        return [{**dict(r),"spec":json.loads(r['spec']),"data":json.loads(r['data'])} for r in rows]

    def update(self, pid, *, state=None, data=None):
        with self.connect() as c:
            c.execute("BEGIN IMMEDIATE")
            row = c.execute("SELECT data FROM projects WHERE id=?", (pid,)).fetchone()
            if not row: raise ValueError("Unknown project")
            if data is not None:
                current = json.loads(row[0]); current.update(data)
                c.execute("UPDATE projects SET data=? WHERE id=?", (dump(current),pid))
            if state: c.execute("UPDATE projects SET state=? WHERE id=?", (state,pid))

    def step(self, pid, sid, status, **data):
        with self.connect() as c:
            row = c.execute("SELECT data FROM steps WHERE project_id=? AND id=?", (pid,sid)).fetchone()
            if not row: raise ValueError("Unknown stage")
            current = json.loads(row[0]); current.update(data); current['status'] = status
            c.execute("UPDATE steps SET status=?,data=? WHERE project_id=? AND id=?", (status,dump(current),pid,sid))

    def ensure_step(self, pid, sid, label, **data):
        """Append a dynamic stage without rewriting its persisted progress."""
        record={'id':sid,'label':label,'status':'pending',**data}
        with self.connect() as c:
            c.execute('INSERT OR IGNORE INTO steps(project_id,id,status,data) VALUES(?,?,?,?)',
                      (pid,sid,'pending',dump(record)))

    def research_update(self, pid, **changes):
        with self.connect() as c:
            c.execute('BEGIN IMMEDIATE')
            row=c.execute('SELECT data FROM projects WHERE id=?',(pid,)).fetchone()
            if not row: raise ValueError('Unknown project')
            data=json.loads(row['data']); research=data.setdefault('research',{})
            research.update(changes)
            c.execute('UPDATE projects SET data=? WHERE id=?',(dump(data),pid))
        return research

    def steps(self,pid):
        with self.connect() as c: rows = c.execute("SELECT data FROM steps WHERE project_id=? ORDER BY id", (pid,)).fetchall()
        return [json.loads(r[0]) for r in rows]

    def event(self, pid, **data):
        with self.connect() as c:
            c.execute("INSERT INTO events(project_id,created_at,data) VALUES(?,?,?)", (pid,utc(),dump(data)))

    def events(self,pid):
        with self.connect() as c: rows = c.execute("SELECT * FROM events WHERE project_id=? ORDER BY id",(pid,)).fetchall()
        return [{"id":r['id'],"created_at":r['created_at'],**json.loads(r['data'])} for r in rows]

    def decision(self,pid,kind,value,reason):
        with self.connect() as c:
            c.execute("INSERT INTO decisions(project_id,kind,value,reason,created_at) VALUES(?,?,?,?,?)",(pid,kind,dump(value),reason,utc()))
        self.event(pid,kind="decision",name=kind,status="recorded",inputs=value,message=reason)

    def decisions(self,pid):
        with self.connect() as c: rows = c.execute("SELECT * FROM decisions WHERE project_id=? ORDER BY id",(pid,)).fetchall()
        return [{**dict(r),"value":json.loads(r['value'])} for r in rows]

    def artifact(self,pid,name,content:bytes,kind="output",job_id=None):
        sha = digest(content)
        clean = re.sub(r"[^A-Za-z0-9._-]", "_", Path(name).name)[:150] or "artifact"
        path = self.directory(pid)/"artifacts"/sha/clean
        path.parent.mkdir(parents=True,exist_ok=True)
        if not path.resolve().is_relative_to(self.directory(pid).resolve()): raise ValueError("Artifact path escapes project")
        if not path.exists():
            try:
                with path.open("xb") as f: f.write(content)
            except FileExistsError:
                pass  # Another writer may have registered identical content.
        if digest(path.read_bytes())!=sha: raise ValueError("Existing immutable artifact failed integrity check")
        aid = digest(f"{pid}:{kind}:{sha}:{clean}".encode())[:32]
        with self.connect() as c:
            c.execute("INSERT OR IGNORE INTO artifacts VALUES(?,?,?,?,?,?,?,?,?)",(aid,pid,clean,str(path),sha,mimetypes.guess_type(clean)[0] or "application/octet-stream",kind,job_id,utc()))
        with self.connect() as c:
            row=c.execute('SELECT * FROM artifacts WHERE project_id=? AND id=?',(pid,aid)).fetchone()
        return self._artifact_record(row)

    def _artifact_record(self, row):
        artifact = dict(row)
        resolved = str(self.relocation.resolve(artifact['path']))
        if resolved != artifact['path']:
            artifact['original_path'] = artifact['path']
            artifact['path'] = resolved
        return artifact

    def artifacts(self,pid):
        with self.connect() as c: rows = c.execute("SELECT * FROM artifacts WHERE project_id=? ORDER BY created_at",(pid,)).fetchall()
        return [self._artifact_record(row) for row in rows]

    def read_artifact(self,pid,artifact):
        with self.connect() as c:
            row=c.execute('SELECT * FROM artifacts WHERE project_id=? AND id=?',(pid,artifact['id'])).fetchone()
        known = self._artifact_record(row) if row else None
        if not known: raise ValueError("Artifact not registered to project")
        path=Path(known['path']).resolve()
        if not path.is_relative_to(self.directory(pid).resolve()): raise ValueError("Artifact escapes project")
        content=path.read_bytes()
        if digest(content)!=known['sha256']: raise ValueError("Artifact integrity check failed")
        return content

    def jobs(self,pid):
        with self.connect() as c: rows=c.execute("SELECT * FROM jobs WHERE project_id=?",(pid,)).fetchall()
        return [dict(r) for r in rows]

    def enqueue(self,pid):
        now=utc()
        with self.connect() as c:
            c.execute("BEGIN IMMEDIATE")
            row=c.execute("SELECT * FROM jobs WHERE project_id=?",(pid,)).fetchone()
            if row:
                if row['status']=="completed": return row['id']
                if row['status']=="running" and self._local_worker(row) and ((row['lease_until'] or 0)>time.time() or self._worker_alive(row)): return row['id']
                c.execute("UPDATE jobs SET status='queued',cancel_requested=0,updated_at=?,error=NULL,lease_until=NULL,owner_token=NULL,worker_pid=NULL,worker_host=NULL,worker_boot=NULL,worker_start=NULL WHERE id=?",(now,row['id']))
                return row['id']
            jid=uuid.uuid4().hex
            c.execute("INSERT INTO jobs(id,project_id,status,created_at,updated_at) VALUES(?,?,'queued',?,?)",(jid,pid,now,now))
            return jid

    def _local_worker(self, row):
        return row['worker_host'] == self.worker_host and row['worker_boot'] == self.worker_boot

    def _worker_alive(self, row):
        # Do not probe or signal a PID copied from another host/boot. PID reuse
        # on this boot is also rejected by the process start-time identity.
        return bool(self._local_worker(row) and row['worker_start'] and
                    process_identity(row['worker_pid']) == row['worker_start'])

    def claim(self,jid,owner_token=None,worker_pid=None,project_id=None):
        with self.connect() as c:
            c.execute("BEGIN IMMEDIATE")
            result=c.execute("UPDATE jobs SET status='running',attempts=attempts+1,lease_until=?,updated_at=?,owner_token=?,worker_pid=?,worker_host=?,worker_boot=?,worker_start=? WHERE id=? AND status='queued' AND cancel_requested=0 AND (? IS NULL OR project_id=?)",(time.time()+90,utc(),owner_token or uuid.uuid4().hex,worker_pid,self.worker_host,self.worker_boot,process_identity(worker_pid),jid,project_id,project_id))
            return result.rowcount==1

    def owns_lease(self,jid,owner_token):
        with self.connect() as c:
            row=c.execute("SELECT 1 FROM jobs WHERE id=? AND status='running' AND owner_token=? AND worker_host=? AND worker_boot=?",(jid,owner_token,self.worker_host,self.worker_boot)).fetchone()
        return bool(row)

    def heartbeat(self,jid,owner_token=None):
        with self.connect() as c:
            changed=c.execute("UPDATE jobs SET lease_until=?,updated_at=? WHERE id=? AND status='running' AND worker_host=? AND worker_boot=? AND (? IS NULL OR owner_token=?)",(time.time()+90,utc(),jid,self.worker_host,self.worker_boot,owner_token,owner_token))
            return changed.rowcount==1

    def cancelled(self,jid):
        with self.connect() as c: row=c.execute("SELECT cancel_requested FROM jobs WHERE id=?",(jid,)).fetchone()
        return bool(row and row[0])

    def cancel(self,pid):
        with self.connect() as c:
            c.execute("UPDATE jobs SET cancel_requested=1,status=CASE WHEN status='queued' THEN 'cancelled' ELSE status END WHERE project_id=? AND status IN ('running','queued')",(pid,))
        self.event(pid,kind="control",name="Cancellation requested",status="requested")

    def finish(self,jid,status,error=None,*,owner_token=None,project_id=None,project_data=None):
        if status not in {"completed","blocked","failed","cancelled"}: raise ValueError("Illegal terminal job status")
        with self.connect() as c:
            c.execute('BEGIN IMMEDIATE')
            row=c.execute('SELECT * FROM jobs WHERE id=?',(jid,)).fetchone()
            if not row or (owner_token is not None and (row['owner_token']!=owner_token or row['status']!='running' or not self._local_worker(row))): return False
            if project_id is not None and row['project_id']!=project_id: raise ValueError('Job belongs to a different project')
            c.execute("UPDATE jobs SET status=?,error=?,updated_at=?,lease_until=NULL,worker_pid=NULL,worker_start=NULL WHERE id=?",(status,safe(error),utc(),jid))
            if project_id is not None:
                project=c.execute('SELECT data FROM projects WHERE id=?',(project_id,)).fetchone()
                current=json.loads(project[0]);current.update(project_data or {})
                c.execute('UPDATE projects SET state=?,data=? WHERE id=?',(status,dump(current),project_id))
            return True
