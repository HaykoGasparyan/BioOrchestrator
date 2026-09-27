"""Score an existing completed docking software test with BioCentric in bounded chunks.

Real model execution, no labels or predictive-validation claims. Original docking
poses and validity flags are preserved; failed physical poses never enter shortlist.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path
import time

from biocentric import engines, science
from biocentric.biocentric_model import installation, parse_output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--docking-run", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--model-dir", type=Path, required=True)
    parser.add_argument("--model-python", type=Path, required=True)
    parser.add_argument("--chunk-size", type=int, default=32)
    args = parser.parse_args()
    if not 1 <= args.chunk_size <= 32:
        raise ValueError("Comparison chunks must contain 1–32 compounds")
    source, out = args.docking_run.resolve(), args.out.resolve()
    summary = json.loads((source / "summary.json").read_text())
    baseline = json.loads((source / "baseline_result.json").read_text())
    if summary["status"] != "completed" or baseline["status"] != "completed" or baseline["failures"]:
        raise ValueError("Requires a completed baseline with one pose per selected compound")
    if out.exists() and any(out.iterdir()):
        raise ValueError("Output directory must be new; previous predictions are evidence")
    out.mkdir(parents=True, exist_ok=True)
    os.environ["BIOCENTRIC_MODEL_DIR"] = str(args.model_dir.resolve())
    os.environ["BIOCENTRIC_MODEL_PYTHON"] = str(args.model_python.resolve())
    os.environ["PYTHONDONTWRITEBYTECODE"] = "1"
    detected = installation()
    if detected["blockers"]:
        raise RuntimeError("; ".join(detected["blockers"]))
    poses = baseline["poses"]
    if len(poses) != len({p["compound_id"] for p in poses}):
        raise ValueError("Comparison must specify exactly one fixed pose per compound")
    if len(poses) != summary["selection"]["requested"]:
        raise ValueError("Incomplete source docking accounting")
    receptor = Path(baseline["metadata"]["config"]["receptor_pdb"])
    from rdkit import Chem
    all_mols = {}
    ledger = []
    for pose in poses:
        path = Path(pose["pose_path"])
        if science.sha256(path) != pose["pose_sha256"]:
            raise ValueError("Source pose hash changed")
        mols = list(Chem.SDMolSupplier(str(path), removeHs=False))
        if len(mols) != 1 or mols[0] is None or science._canonical(mols[0]) != pose["canonical_smiles"]:
            raise ValueError("Source docking pose identity is inconsistent")
        mol = mols[0]
        mol.SetProp("_Name", pose["compound_id"])
        mol.SetProp("BIOCENTRIC_PHYSICAL_CHECK_PASSED", str(pose["validation"].get("passed") is True).lower())
        mol.SetProp("BIOCENTRIC_SOURCE_POSE_SHA256", pose["pose_sha256"])
        all_mols[pose["compound_id"]] = mol
        ledger.append({"compound_id": pose["compound_id"], "name": pose["name"], "smiles": pose["canonical_smiles"],
                       "pose_path": str(path), "pose_sha256": pose["pose_sha256"], "vina_score": pose["score"],
                       "physical_check_passed": pose["validation"].get("passed") is True,
                       "physical_failed_checks": pose["validation"].get("failed_checks", [])})
    protocol = {"schema_version": 1, "purpose": "same-pose software comparison; no measured labels",
                "source_run": str(source), "source_result_sha256": science.sha256(source / "baseline_result.json"),
                "target": summary["target"], "protein_id": "mouse_ABL_1IEP_prepared_pocket",
                "receptor": str(receptor), "receptor_sha256": science.sha256(receptor),
                "models": detected["models"], "weights_sha256": detected["weights_sha256"], "code_sha256": detected["code_sha256"],
                "parameters": {"device": "cuda", "cpu": 2, "max_pairs": args.chunk_size, "timeout_seconds": 300},
                "chunk_size": args.chunk_size, "pairs": ledger,
                "shortlist_criterion": "physical_check_passed AND BioCentric preprocessing/inference completed",
                "selection_claim": "No scientific best-method claim; no held-out affinity metrics without measured records."}
    science._write_json(out / "frozen_protocol.json", protocol)
    started = time.monotonic()
    rows, batch_results = [], []
    for start in range(0, len(poses), args.chunk_size):
        batch_dir = out / ("batch_%03d" % (start // args.chunk_size + 1))
        batch_dir.mkdir()
        ligands = batch_dir / "input.sdf"
        writer = Chem.SDWriter(str(ligands))
        selected = ledger[start:start + args.chunk_size]
        for item in selected:
            writer.write(all_mols[item["compound_id"]])
        writer.close()
        result = engines.run_engine("biocentric_affinity", {"receptor": str(receptor), "ligands": str(ligands),
            "protein_id": protocol["protein_id"], "preparation_reviewed": True}, protocol["parameters"], batch_dir / "engine")
        if result["status"] not in {"completed", "partial", "failed"} or len(result["results"]) != len(selected):
            raise RuntimeError("Incomplete engine batch: " + json.dumps(result.get("error", result.get("blockers"))))
        native = json.loads((batch_dir / "engine/predictions.json").read_text())
        joined = {item["compound_id"]: item for item in selected}
        for row in result["results"]:
            item = joined[row["compound_id"]]
            rows.append({**row, "vina_score": item["vina_score"], "physical_check_passed": item["physical_check_passed"],
                         "physical_failed_checks": item["physical_failed_checks"], "pose_sha256": item["pose_sha256"],
                         "shortlist_eligible": item["physical_check_passed"] and row["status"] == "completed"})
        first_scores = [(row["status"], row["score"], row.get("failure_reason")) for row in result["results"]]
        science.run_process(result["commands"][0], batch_dir / "engine", batch_dir / "resume.log", 300)
        resumed = parse_output(batch_dir / "engine/predictions.json", result["expected_pairs"])
        if [(r["status"], r["score"], r.get("failure_reason")) for r in resumed] != first_scores or not all(r.get("reused") for r in resumed):
            raise RuntimeError("Resume did not retain all success/failure outcomes unchanged")
        batch_results.append({"batch": start // args.chunk_size + 1, "status": result["status"],
                              "accounting": result["accounting"], "resources": native["resources"], "runtime": native["runtime"],
                              "resumed_pairs": len(resumed), "new_predictions_on_resume": 0})
        science._write_json(out / "progress.json", {"completed_pairs": len(rows), "total_pairs": len(poses), "batches": batch_results})
        print(json.dumps(batch_results[-1]), flush=True)
    if {row["compound_id"] for row in rows} != set(all_mols) or len(rows) != len(all_mols):
        raise RuntimeError("Final compound accounting differs from frozen input")
    science._write_json(out / "comparison.json", {"protocol_sha256": science.sha256(out / "frozen_protocol.json"), "rows": rows})
    fields = ["protein_id", "compound_id", "smiles", "status", "vina_score", "score", "physical_check_passed",
              "shortlist_eligible", "failure_reason", "checkpoint_sha256", "code_sha256", "device", "pose_sha256"]
    with (out / "comparison.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    sdf = Chem.SDWriter(str(out / "comparison.sdf"))
    for row in rows:
        mol = Chem.Mol(all_mols[row["compound_id"]])
        for key in ("protein_id", "status", "vina_score", "shortlist_eligible", "checkpoint_sha256", "failure_reason"):
            mol.SetProp("BIOCENTRIC_" + key.upper(), str(row.get(key) if row.get(key) is not None else ""))
        mol.SetProp("BIOCENTRIC_AFFINITY_UNIT", "pK")
        mol.SetProp("BIOCENTRIC_AFFINITY_DIRECTION", "higher_is_better")
        mol.SetProp("VINA_UNIT", "kcal/mol")
        mol.SetProp("VINA_DIRECTION", "lower_is_better")
        if row["score"] is not None:
            mol.SetDoubleProp("BioCentric_affinity", row["score"])
        sdf.write(mol)
    sdf.close()
    report = {"status": "completed", "requested": len(rows), "predicted": sum(r["status"] == "completed" for r in rows),
              "failed": sum(r["status"] == "failed" for r in rows),
              "physical_passes": sum(r["physical_check_passed"] for r in rows),
              "shortlist_eligible": sum(r["shortlist_eligible"] for r in rows),
              "duration_seconds_including_resume": time.monotonic() - started,
              "prediction_worker_seconds": sum(batch["resources"]["duration_seconds"] for batch in batch_results),
              "peak_rss_kib": max(batch["resources"]["peak_rss_kib"] for batch in batch_results),
              "peak_cuda_allocated_bytes": max(batch["resources"]["peak_cuda_allocated_bytes"] for batch in batch_results),
              "checkpoint_reused_pairs": sum(batch["resumed_pairs"] for batch in batch_results),
              "new_predictions_on_resume": 0, "batches": batch_results,
              "predictive_validation": "not established: NCI software set has no target-matched affinity labels",
              "model_choice": "No predictive winner selected from unlabelled data; Vina and BioCentric scores remain distinct."}
    science._write_json(out / "summary.json", report)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
