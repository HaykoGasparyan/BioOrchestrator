"""Incremental disk exports with bounded SDF parts and atomic publication.

Project metadata is still a bounded in-memory snapshot. Molecular payloads and
serialized tables are written record by record instead of joining a library.
"""
from __future__ import annotations

import csv
import hashlib
import json
import os
from pathlib import Path
import tempfile

from .affinity_results import affinity_rows
from .ligand_results import build_ligand_results, iter_pose_sdf, MAX_EXPORT_BYTES


def _write_tables(directory, stem, rows):
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with (directory / f"{stem}.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            cleaned = {}
            for key, value in row.items():
                if isinstance(value, (list, dict, tuple)):
                    value = json.dumps(value, sort_keys=True, allow_nan=False)
                if isinstance(value, str) and value.startswith(("=", "+", "-", "@", "\t", "\r")):
                    value = "'" + value
                cleaned[key] = value
            writer.writerow(cleaned)
    with (directory / f"{stem}.json").open("w") as handle:
        handle.write("[\n")
        for index, row in enumerate(rows):
            if index:
                handle.write(",\n")
            json.dump(row, handle, allow_nan=False, sort_keys=True)
        handle.write("\n]\n")


def export_streamed(snapshot, root, pid, output_dir, *, poses_per_part=100, max_part_bytes=MAX_EXPORT_BYTES):
    if not 1 <= poses_per_part <= 1000 or not 1 <= max_part_bytes <= MAX_EXPORT_BYTES:
        raise ValueError("Use 1–1,000 poses per part and at most 50 MiB per part")
    output = Path(output_dir).resolve()
    if output.exists():
        raise ValueError("Choose a new output directory; exports never overwrite prior evidence")
    result = build_ligand_results(snapshot, root)
    if result["provisional"]:
        raise ValueError("Wait for registered final results before exporting")
    if snapshot.get("project", {}).get("id") != pid:
        raise ValueError("Project identity mismatch")
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".biocentric-export-", dir=output.parent) as staging:
        directory = Path(staging)
        rows = [{key: value for key, value in row.items() if key not in {"best_pose", "selected_pose", "group"}}
                for row in result["ranked_rows"]]
        _write_tables(directory, "ranked_ligands", rows)
        affinity = affinity_rows(snapshot)
        _write_tables(directory, "affinity_predictions", affinity)
        manifest = {"project_id": pid, "predictive_validation": "not_established", "ranked_rows": len(rows),
                    "affinity_rows": len(affinity), "poses_per_part": poses_per_part,
                    "max_part_bytes": max_part_bytes, "sdf_parts": [], "files": []}
        for label, poses in [("all_poses", result["poses"]),
                             ("shortlist", [row["selected_pose"] for row in result["ranked_rows"] if row["selected"] is True])]:
            handle = None
            count = size = part = 0
            try:
                for content in iter_pose_sdf(snapshot, root, pid, poses) if poses else []:
                    if len(content) > max_part_bytes:
                        raise ValueError("A single pose exceeds the selected part budget")
                    if handle is None or count >= poses_per_part or size + len(content) > max_part_bytes:
                        if handle:
                            handle.close()
                        part += 1
                        name = f"{label}-{part:04d}.sdf"
                        handle = (directory / name).open("xb")
                        manifest["sdf_parts"].append({"name": name, "scope": label, "records": 0})
                        count = size = 0
                    handle.write(content)
                    count += 1
                    size += len(content)
                    manifest["sdf_parts"][-1]["records"] = count
            finally:
                if handle:
                    handle.close()
        for path in sorted(directory.iterdir()):
            with path.open("rb") as handle:
                sha = hashlib.file_digest(handle, "sha256").hexdigest()
            manifest["files"].append({"name": path.name, "bytes": path.stat().st_size, "sha256": sha})
        (directory / "export-manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
        # A fresh directory is the completion marker: failures publish no partial
        # export, and the content-addressed source evidence is never changed.
        if output.exists():
            raise ValueError("Output appeared while exporting; choose a new directory")
        os.rename(directory, output)
    return manifest
