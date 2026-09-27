"""Readable ligand progress, rankings and verified molecular downloads.

Rendering reads immutable registered artifacts only. It never runs chemistry,
models or database writes, and never serves a worker's unregistered pose_path.
"""
from __future__ import annotations

import csv
import hashlib
import io
import math
import re
from collections import defaultdict
from functools import lru_cache
from pathlib import Path

MAX_ARTIFACT_BYTES = 10 * 1024 * 1024
MAX_EXPORT_BYTES = 50 * 1024 * 1024


def _known(snapshot, pid):
    if snapshot.get("project", {}).get("id") != pid:
        raise ValueError("Ligand results must belong to the selected project.")
    return {item["id"]: item for item in snapshot.get("artifacts", [])
            if item.get("id") and item.get("project_id") == pid and item.get("sha256")}


def _read_cached(path, sha, size, mtime_ns, ctime_ns):
    # Metadata can remain identical across immediate same-size writes on Linux.
    # Always authenticate bytes before serving them; chemistry parsing below is
    # still safely cached by the actual content rather than filesystem timestamps.
    content = Path(path).read_bytes()
    after = Path(path).stat()
    if (after.st_size, after.st_mtime_ns, after.st_ctime_ns) != (size, mtime_ns, ctime_ns):
        raise ValueError("Artifact changed while it was being read.")
    if hashlib.sha256(content).hexdigest() != sha:
        raise ValueError("Artifact hash does not match its registered content.")
    return content


def read_registered(snapshot: dict, root: str | Path, pid: str, artifact_id: str) -> bytes:
    """Verify project ownership, content-addressed path and hash before serving bytes."""
    artifact = _known(snapshot, pid).get(artifact_id)
    if artifact is None:
        raise ValueError("Molecular file is not registered to this project.")
    sha = artifact["sha256"]
    if not re.fullmatch(r"[a-f0-9]{64}", sha):
        raise ValueError("Artifact has no valid content hash.")
    directory = (Path(root).resolve() / "projects" / pid / "artifacts").resolve()
    path = Path(artifact.get("path", "")).resolve(strict=True)
    try:
        parts = path.relative_to(directory).parts
    except ValueError as error:
        raise ValueError("Molecular file is outside this project's artifact store.") from error
    if len(parts) != 2 or parts[0] != sha or not path.is_file():
        raise ValueError("Molecular file does not use its registered content-addressed path.")
    stat = path.stat()
    if stat.st_size > MAX_ARTIFACT_BYTES:
        raise ValueError("Molecular preview/download is limited to 10 MiB per file.")
    return _read_cached(str(path), sha, stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns)


def _finite(value):
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)


def _identity(mol):
    from rdkit import Chem
    return Chem.MolToSmiles(Chem.RemoveHs(mol), canonical=True, isomericSmiles=True)


@lru_cache(maxsize=8)
def _library_records(content: bytes, suffix: str):
    """Read the original list; no enumeration, optimization or invented activity labels."""
    from rdkit import Chem
    rows = []
    if suffix in {".sdf", ".sd"}:
        molecules = list(Chem.ForwardSDMolSupplier(io.BytesIO(content), removeHs=False))
        for index, mol in enumerate(molecules, 1):
            name = mol.GetProp("_Name") if mol is not None and mol.HasProp("_Name") else f"Input record {index}"
            metadata = {key: mol.GetProp(key) for key in ("CONTROL_ROLE", "ASSAY_RESULT") if mol.HasProp(key)} if mol is not None else {}
            rows.append((index, name, _identity(mol) if mol is not None else None, metadata))
    elif suffix in {".smi", ".smiles", ".txt"}:
        for line_index, line in enumerate(content.decode("utf-8-sig").splitlines(), 1):
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            parts = line.split(maxsplit=1)
            mol = Chem.MolFromSmiles(parts[0])
            rows.append((line_index, parts[1] if len(parts) > 1 else f"compound_{line_index}", _identity(mol) if mol is not None else None, {}))
    else:
        raise ValueError("Unsupported registered ligand-library format.")
    if len(rows) > 1000:
        raise ValueError("The compact ligand table supports at most 1,000 input records.")
    return tuple((index, name, identity, tuple(metadata.items())) for index, name, identity, metadata in rows)


def _library_artifact(snapshot):
    data = snapshot.get("data", {})
    known = _known(snapshot, snapshot["project"]["id"])
    prepared = data.get("preparation", {}).get("compounds", [])
    source_hash = next((row.get("source_sha256") for row in prepared if row.get("source_sha256")), None)
    libraries = [item for item in known.values() if item.get("kind") == "library"]
    if source_hash:
        match = next((item for item in libraries if item["sha256"] == source_hash), None)
        if match:
            return match
    frozen = data.get("research", {}).get("testset", {})
    library = next((item for item in frozen.get("input_artifacts", []) if item.get("kind") == "library"), None)
    if library and library.get("artifact_id") in known:
        return known[library["artifact_id"]]
    return libraries[-1] if libraries else None


def _screening_poses(snapshot):
    data = snapshot.get("data", {})
    screening = data.get("screening") or {}
    poses = screening.get("poses")
    if isinstance(poses, list) and poses:
        return poses, False
    # Root registers lightweight canonical rows only after each compound completes.
    provisional = {}
    for event in snapshot.get("events", []):
        if str(event.get("stage", "")) not in {"08", "screening"}:
            continue
        output = event.get("outputs") or event.get("data") or {}
        for pose in output.get("poses", []) if isinstance(output, dict) else []:
            if isinstance(pose, dict) and pose.get("compound_id"):
                provisional[(pose.get("engine"), pose.get("compound_id"), pose.get("pose_id") or pose.get("pose_number"))] = pose
    return list(provisional.values()), bool(provisional)


def _metric(pose):
    metrics = pose.get("metrics") or {}
    items = [{"name": key, **value} for key, value in metrics.items() if isinstance(value, dict)] if isinstance(metrics, dict) else metrics
    primary = next((item for item in items if isinstance(item, dict) and item.get("value") == pose.get("score")
                    and item.get("direction") == pose.get("score_direction")), None)
    if primary is None:
        name = "vina_score" if pose.get("engine") in {"AutoDock Vina", "vina"} else "docking_score"
        primary = {"name": name, "unit": pose.get("score_units"), "direction": pose.get("score_direction")}
    return {"name": primary.get("name", "docking_score"), "unit": primary.get("unit") or primary.get("units") or pose.get("score_units"),
            "direction": primary.get("direction") or pose.get("score_direction")}


def _extra_metric(pose, name):
    metrics = pose.get("metrics") or {}
    metric = metrics.get(name) if isinstance(metrics, dict) else next((item for item in metrics if isinstance(item, dict) and item.get("name") == name), None)
    value = metric.get("value") if isinstance(metric, dict) else None
    return value if _finite(value) else None


def _group(pose):
    metric = _metric(pose)
    return (pose.get("engine", "Unspecified engine"), pose.get("engine_version", "Unspecified version"),
            metric["name"], metric["unit"], metric["direction"], pose.get("receptor_sha256"), pose.get("_run_scope") or ("unknown", pose.get("compound_id"), pose.get("pose_id") or pose.get("pose_number")))


def _known_score_semantics(pose):
    metric = _metric(pose)
    if (metric["name"], metric["unit"], metric["direction"]) in {
        ("vina_score", "kcal/mol", "lower_is_better"), ("autodock4_score", "kcal/mol", "lower_is_better"),
        ("unidock_vina_score", "kcal/mol", "lower_is_better"), ("minimizedAffinity", "kcal/mol", "lower_is_better"),
        ("CNNaffinity", "pK", "higher_is_better"), ("CNNscore", "dimensionless", "higher_is_better"),
    }:
        return True
    return (metric["name"] == "docking_score" and metric["unit"] == "kcal/mol" and metric["direction"] == "lower_is_better"
            and pose.get("engine") in {"AutoDock4", "autodock4", "GNINA", "gnina", "Uni-Dock", "unidock"})


def _best(poses):
    finite = [pose for pose in poses if _finite(pose.get("score"))]
    valid = [pose for pose in finite if pose.get("validation", {}).get("passed") is True]
    values = valid or finite
    if not values:
        return None
    direction = _metric(values[0])["direction"]
    if direction not in {"lower_is_better", "higher_is_better"}:
        return None
    return sorted(values, key=lambda pose: pose["score"], reverse=direction == "higher_is_better")[0]


def build_ligand_results(snapshot: dict, root: str | Path | None = None) -> dict:
    """Derive input/progress/ranking rows without fabricating scores or selections."""
    pid = snapshot.get("project", {}).get("id")
    _known(snapshot, pid)
    data = snapshot.get("data", {})
    preparation = data.get("preparation") or {}
    library = _library_artifact(snapshot)
    warnings, inputs = [], []
    if root is not None and library:
        try:
            original = _library_records(read_registered(snapshot, root, pid, library["id"]), Path(library["name"]).suffix.lower())
            for record, name, identity, properties in original:
                metadata = dict(properties)
                inputs.append({"record": record, "name": name, "compound_id": "cmp_" + hashlib.sha256(identity.encode()).hexdigest()[:16] if identity else None,
                               "canonical_smiles": identity, "role": metadata.get("CONTROL_ROLE") or "Screening ligand; activity not assigned",
                               "status": "Awaiting preparation" if identity else "Needs input review", "reason": "" if identity else "Input chemical structure could not be parsed.",
                               "selected_for_calculation": None, "prepared": False})
        except (OSError, ValueError, ImportError) as error:
            warnings.append(str(error))
    by_record = {row["record"]: row for row in inputs}
    for compound in preparation.get("compounds", []):
        record = compound.get("source_record")
        row = by_record.get(record)
        if row is None:
            row = {"record": record, "name": compound.get("name") or compound["compound_id"], "compound_id": compound["compound_id"],
                   "canonical_smiles": compound.get("canonical_smiles"), "role": "Screening ligand; activity not assigned"}
            inputs.append(row)
        row.update(name=compound.get("name") or row["name"], compound_id=compound["compound_id"],
                   status="Prepared · awaiting docking", reason="Accepted by ligand preparation.", selected_for_calculation=True, prepared=True)
        if compound.get("metadata", {}).get("CONTROL_ROLE"):
            row["role"] = compound["metadata"]["CONTROL_ROLE"]
    native_preflight = data.get("research", {}).get("engine_preflight") or {}
    if native_preflight.get("phase") == "screening":
        planned = set(native_preflight.get("validation", {}).get("identities", []))
        for row in inputs:
            if row.get("canonical_smiles") in planned:
                row.update(status="Native inputs checked · awaiting docking", reason="Included in the native engine input manifest.", selected_for_calculation=True)
    exclusions = []
    native_binding = (data.get("screening") or {}).get("metadata", {}).get("frozen_input_binding", {})
    native_exclusions = native_binding.get("exclusions", []) or (native_preflight.get("frozen_input_binding", {}).get("exclusions", []) if native_preflight.get("phase") == "screening" else [])
    for exclusion in [*preparation.get("exclusions", []), *native_exclusions]:
        record = exclusion.get("record_index")
        row = by_record.get(record)
        if row is None:
            row = {"record": record, "name": exclusion.get("name") or f"Input record {record}", "compound_id": None,
                   "role": "Excluded input record", "canonical_smiles": None}
            inputs.append(row)
        row.update(status="Excluded", reason=exclusion.get("reason") or exclusion.get("code") or "Excluded by preparation.", selected_for_calculation=False, prepared=False)
        exclusions.append(dict(exclusion))
    poses, provisional = _screening_poses(snapshot)
    current_job = next((job.get("id") for job in snapshot.get("jobs", []) if job.get("status") in {"running", "queued"}), None)
    fingerprint = (data.get("screening") or {}).get("metadata", {}).get("fingerprint")
    def scoped(pose):
        record = dict(pose)
        record["_run_scope"] = (pose.get("protocol_fingerprint") or pose.get("protocol_sha256") or pose.get("job_id")
                                or (current_job if provisional else fingerprint))
        return record
    poses = [scoped(pose) for pose in poses]
    candidates = [scoped(row) for row in (snapshot.get("candidates") or data.get("candidates") or [])]
    candidate_keys = {(_group(row), row.get("compound_id")): row for row in candidates}
    grouped = defaultdict(list)
    for pose in poses:
        grouped[(_group(pose), pose.get("compound_id"))].append(pose)
    # Some historical snapshots retain candidate records without the full screening object.
    if not poses:
        for key, candidate in candidate_keys.items():
            grouped[key].append(candidate)
    ranked = []
    for (signature, compound_id), records in grouped.items():
        best = _best(records)
        chosen = best or next((row for row in records if row.get("validation", {}).get("passed") is True), records[0])
        rankable = best is not None and bool(chosen.get("_run_scope")) and _known_score_semantics(chosen)
        candidate = candidate_keys.get((signature, compound_id))
        passed = chosen.get("validation", {}).get("passed") is True
        metric = _metric(chosen)
        recorded_selection = candidate.get("selected") if candidate else None
        selection_issue = None
        if recorded_selection is True:
            matching = next((pose for pose in records if pose.get("pose_sha256") and pose.get("pose_sha256") == candidate.get("pose_sha256")
                             and pose.get("pose_id") == candidate.get("pose_id") and pose.get("score") == candidate.get("score")), None)
            if not poses or matching is None:
                selection_issue = "Recorded keep decision has no matching screening pose; review before exporting."
            elif candidate.get("validation", {}).get("passed") is not True or matching.get("validation", {}).get("passed") is not True:
                selection_issue = "Recorded keep decision references a pose that did not pass physical validation."
        reason = candidate.get("rationale", "") if candidate else "Selection has not been recorded yet."
        if not rankable:
            reason = "Score is missing/nonfinite, its metric/units/direction are unsupported, or it lacks a recorded run identity; shown without ranking. " + reason
        if selection_issue:
            reason += " " + selection_issue
        ranked.append({"name": chosen.get("name") or compound_id, "compound_id": compound_id, "engine": signature[0],
                       "engine_version": signature[1], "metric": metric["name"], "units": metric["unit"], "direction": metric["direction"],
                       "best_valid_score": chosen.get("score") if passed and rankable else None,
                       "best_available_score": chosen.get("score") if _finite(chosen.get("score")) else None,
                       "raw_score": str(chosen.get("score")), "rankable": rankable,
                       "pose_validity": "Passed" if passed else "Failed" if chosen.get("validation", {}).get("passed") is False else "Not verified",
                       "poses": len(records), "selected": recorded_selection if not selection_issue else False,
                       "recorded_selection": recorded_selection, "selection_issue": selection_issue, "reason": reason,
                       "provisional": provisional, "group": signature, "best_pose": chosen,
                       "selected_pose": candidate if recorded_selection is True and not selection_issue else None,
                       "selected_pose_score": candidate.get("score") if candidate and candidate.get("selected") is True else None,
                       "cnn_affinity_pK": _extra_metric(chosen, "CNNaffinity"), "cnn_pose_score": _extra_metric(chosen, "CNNscore")})
    signatures = list(dict.fromkeys(row["group"] for row in ranked))
    for signature in signatures:
        group = [row for row in ranked if row["group"] == signature]
        group.sort(key=lambda row: (not row["rankable"], row["best_valid_score"] is None,
                                   (-row["best_available_score"] if row["direction"] == "higher_is_better" else row["best_available_score"]) if row["best_available_score"] is not None else math.inf))
        for rank, row in enumerate(group, 1):
            row["rank"] = rank if row["rankable"] else None
    ranked.sort(key=lambda row: (signatures.index(row["group"]), row["rank"] or math.inf))
    by_compound = defaultdict(list)
    for row in ranked:
        by_compound[row["compound_id"]].append(row)
    active_events = {}
    for event in snapshot.get("events", []):
        if str(event.get("stage", "")) not in {"08", "screening"}:
            continue
        compound_id = event.get("compound_id") or (event.get("details") or {}).get("compound_id") or (event.get("inputs") or {}).get("compound_id")
        if compound_id:
            active_events[compound_id] = event
    failures = (data.get("screening") or {}).get("failures", [])
    failure_map = {row.get("compound_id"): row for row in failures}
    for row in inputs:
        compound_id = row.get("compound_id")
        if row.get("status") == "Excluded":
            continue
        if compound_id in by_compound:
            results = by_compound[compound_id]
            row["selected_for_calculation"] = True
            row["status"] = "Docked · score needs review" if any(not item["rankable"] for item in results) else "Docked · provisional" if provisional else "Docked"
            row["reason"] = "Pose checks recorded; final shortlist pending." if all(item["selected"] is None for item in results) else "Final ranking and shortlist decision recorded."
            row["pose_count"] = sum(item["poses"] for item in results)
        elif compound_id in failure_map:
            row["status"] = "Docking failed"
            row["reason"] = failure_map[compound_id].get("reason") or failure_map[compound_id].get("code") or "Tool failed."
        elif compound_id in active_events:
            event = active_events[compound_id]
            row["status"] = "Docking failed" if event.get("status") == "failed" or event.get("kind") == "compound_failed" else "Processing" if event.get("status") == "running" else "Calculation recorded · finalizing"
            row["reason"] = event.get("message") or event.get("name") or "Scientific tool activity recorded."
    testset = data.get("research", {}).get("testset") or {}
    test_rows = [{"compound_id": row.get("compound_id"), "role": row.get("partition", "Unassigned"),
                  "endpoint": row.get("endpoint"), "value": row.get("value"), "units": row.get("units"), "relation": row.get("relation"),
                  "activity_label": row.get("label") or "Unclassified measured activity", "assay": row.get("assay_id")}
                 for row in testset.get("records", [])]
    references = [dict(row) for row in testset.get("pose_records", [])]
    if root is not None:
        for reference in references:
            try:
                rows = _library_records(read_registered(snapshot, root, pid, reference["artifact_id"]), ".sdf")
                reference["name"] = rows[0][1] if len(rows) == 1 else "Observed bound-ligand reference"
            except (KeyError, OSError, ValueError, ImportError):
                reference["name"] = "Observed bound-ligand reference (file unavailable)"
    latest_library = next((item for item in reversed(snapshot.get("artifacts", [])) if item.get("kind") == "library" and item.get("project_id") == pid), None)
    draft = latest_library if library and latest_library and latest_library["id"] != library["id"] else None
    draft_rows = []
    if draft and root is not None:
        try:
            draft_rows = [{"Input record": index, "Ligand": name, "Chemical identity": identity,
                           "Status": "Uploaded draft; not included in this frozen run"}
                          for index, name, identity, _ in _library_records(read_registered(snapshot, root, pid, draft["id"]), Path(draft["name"]).suffix.lower())]
        except (OSError, ValueError, ImportError) as error:
            warnings.append("New draft library: " + str(error))
    frozen_library = any(row.get("kind") == "library" and library and row.get("artifact_id") == library["id"] for row in testset.get("input_artifacts", []))
    return {"input_rows": inputs, "frozen_library": frozen_library, "draft_library": draft, "draft_rows": draft_rows, "ranked_rows": ranked, "poses": poses, "provisional": provisional, "groups": signatures,
            "exclusions": exclusions, "failures": failures, "test_rows": test_rows, "reference_rows": references,
            "library_artifact": library, "warnings": warnings, "selected_count": sum(row.get("selected") is True for row in ranked)}


def _pose_artifact(snapshot, pid, pose):
    known = _known(snapshot, pid)
    ids = [*pose.get("artifact_ids", []), *([pose["pose_artifact_id"]] if pose.get("pose_artifact_id") else [])]
    artifact = next((known[aid] for aid in ids if aid in known and known[aid]["sha256"] == pose.get("pose_sha256")
                     and Path(known[aid]["name"]).suffix.lower() == ".sdf"), None)
    if artifact is None:
        raise ValueError("A docked pose has no matching registered SDF; download waits for artifact registration.")
    return artifact


def iter_pose_sdf(snapshot, root, pid, poses):
    """Yield one authenticated, complete SDF record at a time."""
    if not poses:
        raise ValueError("No recorded poses are available for this download.")
    for pose in poses:
        artifact = _pose_artifact(snapshot, pid, pose)
        content = read_registered(snapshot, root, pid, artifact["id"])
        if content.count(b"$$$$") != 1 or not content.rstrip().endswith(b"$$$$") or b"M  END" not in content:
            raise ValueError("Registered pose SDF is incomplete or contains more than one pose.")
        yield content.rstrip() + b"\n"


def export_pose_sdf(snapshot, root, pid, poses) -> bytes:
    """Combine only complete, hash-bound, single-pose records; never read pose_path."""
    chunks, size = [], 0
    for content in iter_pose_sdf(snapshot, root, pid, poses):
        chunks.append(content)
        size += len(content)
        if size > MAX_EXPORT_BYTES:
            raise ValueError("Combined pose download exceeds 50 MiB; download individual ligands instead.")
    return b"".join(chunks)


def ligand_csv(rows):
    """Human-readable results with typed scores and spreadsheet-safe textual cells."""
    fields = ["rank", "name", "compound_id", "engine", "engine_version", "metric", "best_valid_score", "best_available_score", "units", "direction", "pose_validity", "poses", "selected", "selected_pose_score", "cnn_affinity_pK", "cnn_pose_score", "reason", "provisional"]
    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=fields)
    writer.writeheader()
    for row in rows:
        values = {key: row.get(key, "") for key in fields}
        writer.writerow({key: "'" + value if isinstance(value, str) and value.startswith(("=", "+", "-", "@", "\t", "\r")) else value for key, value in values.items()})
    return output.getvalue().encode()


def _table_row(row):
    return {"Rank": row["rank"], "Ligand": row["name"], f"Best score ({row['units']})": row["best_valid_score"],
            "Pose check": row["pose_validity"], "Keep": "Review needed" if row["selection_issue"] else "Yes" if row["selected"] is True else "No" if row["selected"] is False else "Pending"}


def _details_row(row):
    return {"Rank": row["rank"], "Ligand": row["name"], "Compound ID": row["compound_id"], "Engine": row["engine"],
            "Version": row["engine_version"], "Metric": row["metric"], "Best valid score": row["best_valid_score"],
            "Best available score": row["best_available_score"], "Units": row["units"], "Direction": row["direction"],
            "Pose check": row["pose_validity"], "Poses": row["poses"], "Keep": row["selected"],
            "Kept pose score": row["selected_pose_score"], "CNN affinity (pK; higher better)": row["cnn_affinity_pK"],
            "CNN pose score (higher better)": row["cnn_pose_score"], "Recorded keep decision": row["recorded_selection"], "Selection issue": row["selection_issue"], "Raw score": row["raw_score"], "Reason": row["reason"]}


def render_ligand_results(st, snapshot, root, pid):
    """Prominent, read-only research results panel; usable before, during and after runs."""
    st.markdown("### Ligands and results")
    if snapshot.get("project", {}).get("id") != pid:
        st.error("Ligand results do not belong to the selected project.")
        return
    from .affinity_results import render_affinity_results
    from .result_paging import page_rows
    render_affinity_results(st, snapshot, root, pid)
    try:
        result = build_ligand_results(snapshot, root)
    except (ValueError, OSError) as error:
        st.error(str(error))
        return
    for warning in result["warnings"]:
        st.warning("Ligand input preview unavailable: " + warning)
    rows = result["ranked_rows"]
    if rows:
        st.write(f"{len(rows)} ligands · {len(result['poses']) or sum(row['poses'] for row in rows)} docked poses · {result['selected_count']} kept for review")
        if result["provisional"]:
            st.caption("Live results: completed compounds appear during the run. The final shortlist and downloads follow after screening finishes.")
        for signature in result["groups"]:
            group = [row for row in rows if row["group"] == signature]
            engine, version, metric, units, direction, _, run_scope = signature
            st.markdown(f"**{engine} · {metric}**")
            meaning = "Lower (more negative) is better" if direction == "lower_is_better" else "Higher is better" if direction == "higher_is_better" else "Score direction is not recorded; ranking is unavailable"
            st.caption(f"{meaning}. Best score uses physically valid poses; Keep means prioritized for review, not confirmed binding.")
            page = page_rows(st, group, f"docking_page_{pid}_{hashlib.sha256(repr(signature).encode()).hexdigest()[:12]}", f"{engine} results")
            st.dataframe([_table_row(row) for row in page], hide_index=True, width="stretch")
            if any(row["selection_issue"] for row in group):
                st.warning("A recorded keep decision conflicts with its validated pose. It is excluded from the kept-ligand download until reviewed.")
            if any(not row["rankable"] for row in group):
                st.caption("Some recorded results cannot be ranked; see score details for the missing information.")
        selected = [row for row in rows if row["selected"] is True]
        shortlist, all_poses, errors = None, None, []
        if not result["provisional"]:
            if selected:
                try:
                    shortlist = export_pose_sdf(snapshot, root, pid, [row["selected_pose"] for row in selected])
                except (ValueError, OSError) as error:
                    errors.append("Kept-ligand download: " + str(error))
            if result["poses"]:
                try:
                    all_poses = export_pose_sdf(snapshot, root, pid, result["poses"])
                except (ValueError, OSError) as error:
                    errors.append("All-pose download: " + str(error))
        downloads = st.columns(3)
        with downloads[0]:
            st.download_button("Ranked CSV", ligand_csv(rows), file_name="ranked_ligands.csv", mime="text/csv", key=f"ligands_ranked_{pid}", width="stretch")
        with downloads[1]:
            st.download_button("Kept ligands · SDF", shortlist or b"", file_name="shortlist.sdf", mime="chemical/x-mdl-sdfile", disabled=shortlist is None, key=f"ligands_shortlist_{pid}", width="stretch")
        with downloads[2]:
            st.download_button("All poses · SDF", all_poses or b"", file_name="all_docked_poses.sdf", mime="chemical/x-mdl-sdfile", disabled=all_poses is None, key=f"ligands_all_poses_{pid}", width="stretch")
        for error in errors:
            st.warning(error)
        with st.expander("Score details and individual ligand downloads", expanded=False):
            st.caption("Ranks are compared within recorded engine/version, score semantics, receptor hash and run identity; unknown runs are not pooled. All-pose export includes any rejected poses.")
            page = page_rows(st, rows, f"score_details_page_{pid}", "Score details")
            st.dataframe([_details_row(row) for row in page], hide_index=True, width="stretch")
            if any(row.get("cnn_affinity_pK") is not None for row in rows):
                st.caption("CNN affinity is a learned pK prediction: higher is better. CNN pose score describes pose quality, not measured binding.")
            if selected:
                st.download_button("Download kept ligand CSV", ligand_csv(selected), file_name="kept_ligands.csv", mime="text/csv", key=f"ligands_kept_csv_{pid}")
            if not result["provisional"]:
                selected_index = st.selectbox("Choose a ligand to download", list(range(len(rows))),
                                              format_func=lambda index: f"{rows[index]['name']} · {rows[index]['engine']} · rank {rows[index]['rank']}", key=f"ligand_download_choice_{pid}")
                chosen = rows[selected_index if isinstance(selected_index, int) else 0]
                try:
                    content = export_pose_sdf(snapshot, root, pid, [chosen["best_pose"]])
                    filename = re.sub(r"[^A-Za-z0-9_.-]", "_", str(chosen["name"]))[:100] or "ligand"
                    label = "Download best valid pose (SDF)" if chosen["best_valid_score"] is not None else "Download unvalidated pose for review (SDF)"
                    st.download_button(label, content, file_name=filename + ".sdf", mime="chemical/x-mdl-sdfile", key=f"ligand_individual_{pid}")
                except (ValueError, OSError) as error:
                    st.warning(str(error))
    else:
        st.info("No docking scores yet. The ligand list below tracks preparation and calculation status as the workflow runs.")
    if result["input_rows"]:
        with st.expander("Full ligand list for this frozen experiment" if result["frozen_library"] else "Full ligand list for this experiment", expanded=not bool(rows)):
            page = page_rows(st, result["input_rows"], f"input_page_{pid}", "Input library")
            st.dataframe([{"Input record": row.get("record"), "Ligand": row.get("name"), "Compound ID": row.get("compound_id"),
                           "Role": row.get("role"), "Status": row.get("status"), "Reason": row.get("reason"),
                           "Included for docking": "Yes" if row.get("selected_for_calculation") is True else "No" if row.get("selected_for_calculation") is False else "Awaiting preparation"}
                          for row in page], hide_index=True, width="stretch")
            library = result["library_artifact"]
            if library:
                try:
                    content = read_registered(snapshot, root, pid, library["id"])
                    st.download_button("Download experiment input library", content, file_name=library["name"], key=f"ligands_input_{pid}")
                except (ValueError, OSError) as error:
                    st.warning(str(error))
    elif not rows:
        st.caption("No registered ligand library or preparation records are available yet.")
    if result["draft_library"]:
        with st.expander("Newly uploaded library · not in this frozen run", expanded=True):
            st.caption("This uploaded draft is separate from the frozen experiment above. It enters the calculation only after the workflow reviews and freezes the revision.")
            if result["draft_rows"]:
                st.dataframe(result["draft_rows"], hide_index=True, width="stretch")
            try:
                draft = result["draft_library"]
                st.download_button("Download newly uploaded library", read_registered(snapshot, root, pid, draft["id"]), file_name=draft["name"], key=f"ligands_draft_{pid}")
            except (ValueError, OSError) as error:
                st.warning(str(error))
    with st.expander("Test compounds and pose-recovery controls", expanded=False):
        if result["test_rows"]:
            st.write("These measured records belong to the frozen tuning/evaluation set; they are separate from the screening library.")
            st.dataframe(result["test_rows"], hide_index=True, width="stretch")
        else:
            st.caption("No frozen measured activity test set is recorded. Untested screening compounds are not labeled inactive.")
        if result["reference_rows"]:
            st.dataframe([{"Reference": row.get("name") or row.get("id"), "Role": "Observed bound-ligand pose", "Scope": row.get("scope", "Pose recovery only"), "Artifact ID": row.get("artifact_id")} for row in result["reference_rows"]], hide_index=True, width="stretch")
            st.caption("The observed ligand tests pose recovery; this is not an independent affinity or screening benchmark.")
