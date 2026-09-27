"""Affinity evidence and downloads, independent of docking-pose ranks."""
from __future__ import annotations

import csv
import io
import json
import math


def affinity_runs(snapshot):
    data = snapshot.get("data", {})
    value = data.get("research", {}).get("production", {}).get("affinity", {})
    runs = list(data.get("affinity_results", []))
    if value.get("engine") or value.get("engine_id"):
        runs.append(value)
    return runs


def affinity_rows(snapshot):
    rows = []
    docked = {}
    for pose in snapshot.get("data", {}).get("screening", {}).get("poses", []):
        cid = pose.get("compound_id")
        if cid not in docked or pose.get("pose_number", 999) < docked[cid].get("pose_number", 999):
            docked[cid] = pose
    for run_index, run in enumerate(affinity_runs(snapshot), 1):
        native = run.get("native_result", run)
        provenance = native.get("provenance", {})
        for result in run.get("results", native.get("results", [])):
            status = result.get("status", "completed" if result.get("metrics") else "failed")
            metrics = result.get("metrics") or [{}]
            for metric in metrics:
                pose = docked.get(result.get("compound_id"), {})
                value = metric.get("value")
                finite = isinstance(value, (float, int)) and not isinstance(value, bool) and math.isfinite(value)
                valid = status == "completed" and finite
                rows.append({"run": run_index, "compound_id": result.get("compound_id") or result.get("name"),
                             "protein_id": result.get("protein_id") or native.get("input_review", {}).get("protein_id"),
                             "smiles": result.get("smiles"), "engine": run.get("engine", native.get("engine_id")),
                             "metric": metric.get("name"), "score": value if valid else None,
                             "units": metric.get("unit"), "direction": metric.get("direction"),
                             "status": status if valid or status != "completed" else "invalid_score",
                             "preprocessing": result.get("preprocessing_status", "recorded in input plan"),
                             "failure_reason": result.get("failure_reason") or (None if valid else "No complete finite prediction"),
                             "checkpoint_sha256": result.get("checkpoint_sha256") or provenance.get("weights_sha256"),
                             "code_sha256": result.get("code_sha256") or provenance.get("code_sha256"),
                             "model": result.get("models") or provenance.get("models"),
                             "device": result.get("device") or native.get("parameters", {}).get("device"),
                             "physical_check_passed": result.get("physical_check_passed"),
                             "input_structure_physical_passed": result.get("input_structure_physical_passed", result.get("physical_check_passed")),
                             "input_structure_provenance": result.get("input_structure_provenance"),
                             "shortlist_eligible": result.get("shortlist_eligible"),
                             "docking_engine": pose.get("engine", "Vina" if "vina_score" in result else None),
                             "docking_score": pose.get("score", result.get("vina_score")),
                             "docking_units": pose.get("score_units", "kcal/mol" if pose or "vina_score" in result else None),
                             "input_sha256": {key: item.get("sha256") for key, item in native.get("validated_inputs", {}).items() if key in {"receptor", "ligands", "ligand"}},
                             "result_artifact_id": run.get("result_artifact_id"), "predictive_validation": "not_established"})
        for excluded in run.get("metadata", {}).get("frozen_input_binding", {}).get("exclusions", []):
            rows.append({"run": run_index, "compound_id": f"input_record_{excluded['record_index']}", "smiles": excluded.get("canonical_smiles"),
                         "engine": run.get("engine", native.get("engine_id")), "status": "excluded", "score": None,
                         "failure_reason": excluded["reason"], "predictive_validation": "not_established"})
        if not native.get("results") and not run.get("results") and native.get("error"):
            rows.append({"run": run_index, "engine": native.get("engine_id"), "status": "failed", "score": None,
                         "failure_reason": native["error"].get("message"), "predictive_validation": "not_established"})
    return rows


def affinity_csv(rows):
    output = io.StringIO()
    fields = ["run", "compound_id", "protein_id", "smiles", "engine", "metric", "score", "units", "direction", "status", "preprocessing",
              "failure_reason", "checkpoint_sha256", "code_sha256", "model", "device", "physical_check_passed", "input_structure_physical_passed", "input_structure_provenance", "shortlist_eligible", "docking_engine", "docking_score", "docking_units", "input_sha256", "result_artifact_id", "predictive_validation"]
    writer = csv.DictWriter(output, fields)
    writer.writeheader()
    for row in rows:
        record = {}
        for key in fields:
            value = row.get(key)
            if isinstance(value, (dict, list)):
                value = json.dumps(value, sort_keys=True)
            if isinstance(value, str) and value.startswith(("=", "+", "-", "@")):
                value = "'" + value
            record[key] = value
        writer.writerow(record)
    return output.getvalue().encode()


def render_affinity_results(st, snapshot, root, pid):
    rows = affinity_rows(snapshot)
    if not rows:
        return
    st.markdown("### Affinity predictions")
    st.caption("Docking and affinity scores appear side by side with separate units; they are never combined into one rank. Missing predictions stay blank. Physically rejected structures remain excluded.")
    st.caption("Predictive accuracy has not been established for these results. These model predictions are not measured binding affinities.")
    st.write(" · ".join(sorted({f"{row.get('engine')} / {row.get('metric')} ({row.get('units')}, {str(row.get('direction')).replace('_', ' ')})" for row in rows if row.get('metric')})))
    columns = ["compound_id", "docking_score", "docking_units", "score", "units", "status", "input_structure_physical_passed", "shortlist_eligible"]
    from .result_paging import page_rows
    visible_rows = page_rows(st, rows, f"affinity_page_{pid}", "Affinity")
    st.dataframe([{key: row.get(key) for key in columns} for row in visible_rows], hide_index=True, width="stretch",
                 column_config={"compound_id": "Compound", "docking_score": st.column_config.NumberColumn("Docking score", format="%.3f"),
                                "docking_units": "Docking units", "score": st.column_config.NumberColumn("Affinity prediction", format="%.4f"),
                                "units": "Affinity units", "status": "Status", "input_structure_physical_passed": "Input physical check",
                                "shortlist_eligible": "Shortlist eligible"})
    left, right = st.columns(2)
    with left:
        st.download_button("Affinity CSV", affinity_csv(rows), file_name="affinity_predictions.csv", mime="text/csv", key=f"affinity_csv_{pid}")
    with right:
        st.download_button("Affinity JSON", json.dumps(rows, indent=2, allow_nan=False), file_name="affinity_predictions.json", mime="application/json", key=f"affinity_json_{pid}")
    with st.expander("Affinity provenance, preprocessing and actual structures"):
        st.json(visible_rows)
        from .ligand_results import read_registered
        for index, run in enumerate(affinity_runs(snapshot)):
            native = run.get("native_result", run)
            st.write({"engine": run.get("engine", native.get("engine_id")), "status": run.get("status"),
                      "elapsed_seconds": native.get("duration_seconds"), "accounting": native.get("accounting"),
                      "resource_usage": native.get("resource_usage"), "native_runtime": native.get("native_runtime"),
                      "commands": native.get("commands"), "checkpoint": native.get("provenance")})
            for artifact in run.get("artifacts", []):
                if artifact.get("name") == "scored.sdf":
                    try:
                        content = read_registered(snapshot, root, pid, artifact["id"])
                        st.download_button("Prediction structures · SDF", content, file_name=f"affinity_structures_{index+1}.sdf",
                                           mime="chemical/x-mdl-sdfile", key=f"affinity_sdf_{pid}_{index}_{artifact['id']}")
                    except (ValueError, OSError) as error:
                        st.warning(str(error))
