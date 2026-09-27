"""Explicit, hash-checked relocation without rewriting historical evidence.

The mapping changes file lookup only. JSON documents, recorded commands, hashes,
and original SQLite artifact paths are retained verbatim.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import sqlite3
import tarfile
from datetime import datetime, timezone
from pathlib import Path


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


class RelocationResolver:
    def __init__(self, mappings=()):
        self.mappings = []
        for item in mappings:
            old, new = Path(item["old_root"]), Path(item["new_root"])
            if not old.is_absolute() or not new.is_absolute() or ".." in old.parts or ".." in new.parts:
                raise ValueError("Relocation roots must be absolute paths without traversal")
            if old == Path(old.anchor) or new == Path(new.anchor):
                raise ValueError("A filesystem root is too broad for relocation")
            if any(existing == old for existing, _ in self.mappings):
                raise ValueError("Duplicate relocation root")
            self.mappings.append((old, new.resolve()))
        self.mappings.sort(key=lambda pair: len(pair[0].parts), reverse=True)

    @classmethod
    def load(cls, store_root: Path):
        config = store_root / "relocation.json"
        if not config.exists():
            return cls()
        record = json.loads(config.read_text())
        if record.get("schema_version") != 1 or not isinstance(record.get("mappings"), list):
            raise ValueError("Unsupported relocation configuration")
        return cls(record["mappings"])

    def describe(self):
        return {"schema_version": 1, "mappings": [
            {"old_root": str(old), "new_root": str(new)} for old, new in self.mappings],
            "historical_contents_rewritten": False}

    def resolve(self, original: str | Path) -> Path:
        path = Path(original)
        if not path.is_absolute() or ".." in path.parts:
            raise ValueError("Recorded file path must be absolute and cannot contain traversal")
        for old, new in self.mappings:
            if path.is_relative_to(old):
                candidate = (new / path.relative_to(old)).resolve()
                if not candidate.is_relative_to(new):
                    raise ValueError("Relocated file escapes the reviewed destination root")
                return candidate
        return path.resolve()

    def verified(self, original: str | Path, sha256: str) -> Path:
        if not isinstance(sha256, str) or not re.fullmatch(r"[a-f0-9]{64}", sha256):
            raise ValueError("A valid SHA256 is required to resolve historical evidence")
        path = self.resolve(original)
        if file_sha256(path) != sha256:
            raise ValueError("Relocated file failed its recorded SHA256 integrity check")
        return path


def backup_history(store_root: Path, destination: Path, *, evidence_roots=()) -> dict:
    """Preserve a consistent SQLite snapshot plus every supplied evidence file."""
    store_root, destination = store_root.resolve(), destination.resolve()
    if destination.is_relative_to(store_root):
        raise ValueError("Backup directory must be outside the store")
    if any(destination.is_relative_to(Path(root).resolve()) for root in evidence_roots):
        raise ValueError("Backup directory must be outside supplied evidence directories")
    destination.mkdir(parents=True, exist_ok=False)
    database = store_root / "metadata.sqlite"
    with sqlite3.connect(f"file:{database}?mode=ro", uri=True) as source:
        with sqlite3.connect(destination / "metadata.sqlite") as target:
            source.backup(target)
            target.row_factory = sqlite3.Row
            artifacts = [dict(row) for row in target.execute("SELECT * FROM artifacts ORDER BY id")]
            jobs = [dict(row) for row in target.execute("SELECT * FROM jobs ORDER BY id")]
    (destination / "artifact-records.json").write_text(json.dumps(artifacts, indent=2) + "\n")
    (destination / "historical-jobs.json").write_text(json.dumps(jobs, indent=2) + "\n")
    inventory = []
    archive = destination / "evidence.tar.gz"
    with tarfile.open(archive, "w:gz") as bundle:
        roots = [(store_root, "store"), *[(Path(root).resolve(), f"evidence-{index}")
                                          for index, root in enumerate(evidence_roots)]]
        for root, prefix in roots:
            for path in sorted(root.rglob("*")):
                if path.is_symlink():
                    raise ValueError(f"Refusing an evidence symlink during backup: {path}")
                if not path.is_file() or (root == store_root and path.name.startswith("metadata.sqlite")):
                    continue
                archive_name = f"{prefix}/{path.relative_to(root)}"
                before = file_sha256(path)
                bundle.add(path, arcname=archive_name, recursive=False)
                if file_sha256(path) != before:
                    raise ValueError(f"Evidence changed during backup: {path}")
                inventory.append({"source": str(path), "archive_path": archive_name,
                                  "bytes": path.stat().st_size, "sha256": before})
    (destination / "file-inventory.json").write_text(json.dumps(inventory, indent=2) + "\n")
    record = {"schema_version": 1, "created_at": datetime.now(timezone.utc).isoformat(),
              "store_root": str(store_root), "destination": str(destination),
              "files": len(inventory), "artifact_records": len(artifacts),
              "checksums": {path.name: file_sha256(path) for path in sorted(destination.iterdir()) if path.is_file()}}
    (destination / "backup.json").write_text(json.dumps(record, indent=2) + "\n")
    return record


def configure(store_root: Path, mappings: list[dict], backup_directory: Path) -> dict:
    """Back up first; verify every registered artifact before enabling mapping."""
    store_root = store_root.resolve()
    if (store_root / "relocation.json").exists():
        raise ValueError("Relocation is already configured; review the existing mapping")
    resolver = RelocationResolver(mappings)
    backup_history(store_root, backup_directory)
    verified = []
    with sqlite3.connect(f"file:{store_root / 'metadata.sqlite'}?mode=ro", uri=True) as database:
        database.row_factory = sqlite3.Row
        for row in database.execute("SELECT * FROM artifacts ORDER BY id"):
            path = resolver.verified(row["path"], row["sha256"])
            expected = (store_root / "projects" / row["project_id"] / "artifacts" / row["sha256"] / row["name"]).resolve()
            if path != expected or not path.is_relative_to(store_root):
                raise ValueError("Relocated artifact is outside its registered project and content address")
            verified.append({"id": row["id"], "original_path": row["path"], "resolved_path": str(path), "sha256": row["sha256"]})
    record = {**resolver.describe(), "created_at": datetime.now(timezone.utc).isoformat(),
              "backup": str(backup_directory.resolve()), "verified_artifacts": len(verified),
              "backup_manifest_sha256": file_sha256(backup_directory / "backup.json")}
    with (store_root / "relocation.json").open("x") as stream:
        stream.write(json.dumps(record, indent=2) + "\n")
    (backup_directory / "relocation-verification.json").write_text(json.dumps(verified, indent=2) + "\n")
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--store", type=Path, default=Path(".biocentric"))
    parser.add_argument("--old-root", required=True)
    parser.add_argument("--new-root", required=True)
    parser.add_argument("--backup", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(configure(args.store, [{"old_root": args.old_root, "new_root": args.new_root}], args.backup), indent=2))


if __name__ == "__main__":
    main()
