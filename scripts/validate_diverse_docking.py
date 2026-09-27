"""Bounded, reproducible local software stress test, not an affinity benchmark.

Run: .venv/bin/python -m scripts.validate_diverse_docking
No network, LLM calls, activity labels, PLIP, or database writes are performed.
Existing evidence is reused only when pinned inputs/configuration still match.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import hashlib
import importlib.metadata
import io
import json
import math
from pathlib import Path
import platform
import resource
import sys
import time
from datetime import datetime, timezone

from biocentric import science

REPO = Path(__file__).resolve().parents[1]
DEFAULT_OUT = REPO / "runs/diverse-docking-verification"
SEED = 42
ADVERSARIAL = """C[C@H](O)C(=O)O assigned_left
C[C@@H](O)C(=O)O assigned_right
C[C@H](O)C(=O)O duplicate_left
CC(O)C(=O)O ambiguous_stereo
CCO.[Na+] disconnected_salt
[Fe] unsupported_metal
not_a_smiles invalid_smiles
"""
EXPECTED_ADVERSARIAL = {"exact_duplicate", "ambiguous_stereochemistry", "unsupported_salt", "unsupported_element", "invalid_chemistry"}


def now():
    return datetime.now(timezone.utc).isoformat()


def digest_json(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()).hexdigest()


def save(path, value):
    science._write_json(Path(path), value)


def immutable_text(path: Path, content: str):
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        if path.read_text() != content:
            raise ValueError(f"Existing evidence differs: {path}. Use a new output directory.")
    else:
        with path.open("x") as handle:
            handle.write(content)


def select_diverse(source: Path, count=128, seed=SEED):
    """One representative per Murcko scaffold, then seeded Morgan MaxMin.

    Bounds deliberately keep this software test affordable; they are not ADMET
    rules. Original row IDs, canonical identity and every omission are recorded.
    """
    from rdkit import Chem, rdBase
    from rdkit.Chem import Descriptors, Lipinski, rdFingerprintGenerator
    from rdkit.Chem.Scaffolds import MurckoScaffold
    from rdkit.SimDivFilters.rdSimDivPickers import MaxMinPicker
    if count < 1:
        raise ValueError("count must be positive")
    rows, canonical_seen, scaffolds = [], {}, {}
    with rdBase.BlockLogs():
        for line_number, line in enumerate(source.read_text().splitlines(), 1):
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            fields = line.split(maxsplit=1)
            row = {"source_line": line_number, "source_id": fields[1].strip() if len(fields) > 1 else str(line_number), "input_smiles": fields[0]}
            try:
                mol = Chem.MolFromSmiles(fields[0])
                if mol is None:
                    raise science.ScienceError("invalid_chemistry", "RDKit parsing/sanitization failed")
                science._supported(mol)
                canonical = science._canonical(mol)
                row.update(canonical_smiles=canonical, molecular_weight=Descriptors.MolWt(mol), heavy_atoms=mol.GetNumHeavyAtoms(), rotatable_bonds=Lipinski.NumRotatableBonds(mol))
                if canonical in canonical_seen:
                    row.update(status="excluded", code="exact_duplicate", duplicate_of=canonical_seen[canonical])
                elif not (100 <= row["molecular_weight"] <= 500 and 8 <= row["heavy_atoms"] <= 35 and row["rotatable_bonds"] <= 8):
                    row.update(status="excluded", code="bounded_test_domain", reason="Software-test bounds: MW 100–500, 8–35 heavy atoms, at most 8 rotatable bonds")
                    canonical_seen[canonical] = row["source_id"]
                else:
                    canonical_seen[canonical] = row["source_id"]
                    scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=False) or "ACYCLIC"
                    row.update(status="eligible_not_selected", scaffold=scaffold, name="NCI_" + row["source_id"], compound_id="cmp_" + hashlib.sha256(canonical.encode()).hexdigest()[:16])
                    scaffolds.setdefault(scaffold, []).append(row)
            except science.ScienceError as error:
                row.update(status="excluded", code=error.code, reason=str(error))
            rows.append(row)
    representatives = sorted((min(group, key=lambda row: (row["heavy_atoms"], row["rotatable_bonds"], row["canonical_smiles"])) for group in scaffolds.values()), key=lambda row: row["canonical_smiles"])
    if len(representatives) < count:
        raise ValueError(f"Requested {count} distinct scaffolds, only {len(representatives)} eligible")
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fingerprints = [generator.GetFingerprint(Chem.MolFromSmiles(row["canonical_smiles"])) for row in representatives]
    indices = list(MaxMinPicker().LazyBitVectorPick(fingerprints, len(fingerprints), count, seed=seed))
    selected = [representatives[index] for index in indices]
    for index, row in enumerate(selected, 1):
        row.update(status="selected", selection_order=index)
    return {"source_path": str(source.resolve()), "source_sha256": science.sha256(source), "rdkit_version": rdBase.rdkitVersion,
            "seed": seed, "requested": count, "input_records": len(rows), "eligible_records": sum(len(group) for group in scaffolds.values()),
            "eligible_unique_scaffolds": len(scaffolds), "method": "One smallest-heavy-atom representative per exact nonchiral Bemis–Murcko scaffold; seeded MaxMin on radius-2/2048-bit Morgan fingerprints; all acyclic molecules share one scaffold",
            "limits": {"molecular_weight": [100, 500], "heavy_atoms": [8, 35], "max_rotatable_bonds": 8},
            "selected": selected, "source_ledger": rows,
            "source_reference": "https://github.com/rdkit/rdkit/blob/master/Data/NCI/first_5K.smi",
            "limitations": "Structural diversity software test only. No activity labels, enrichment, affinity calibration, or new-binder claim."}


def verify_artifacts(artifacts):
    for artifact in artifacts:
        path = Path(artifact["path"])
        if not path.is_file() or science.sha256(path) != artifact["sha256"]:
            raise ValueError(f"Artifact integrity failed: {path}")


def prepare_cached(source, directory, emit, cancelled):
    manifest = directory / "preparation.json"
    if manifest.exists():
        result = json.loads(manifest.read_text())
        originals = [artifact for artifact in result["artifacts"] if artifact["role"] == "immutable_original"]
        if len(originals) != 1 or originals[0]["sha256"] != science.sha256(source):
            raise ValueError("Preparation input changed; refusing to replace existing evidence")
        verify_artifacts(result["artifacts"])
        return result
    return science.prepare_library(source, directory, seed=SEED, emit=emit, cancelled=cancelled)


def summarize_result(result, requested_ids):
    poses = result.get("poses", [])
    completed = {pose["compound_id"] for pose in poses}
    failures = {failure.get("compound_id"): failure for failure in result.get("failures", [])}
    requested = set(requested_ids)
    if completed - requested or set(failures) - requested:
        raise ValueError("Unexpected compound identity in docking results")
    finite = [float(pose["score"]) for pose in poses if isinstance(pose.get("score"), (int, float)) and math.isfinite(pose["score"])]
    return {"requested_ligands": len(requested), "ligands_with_poses": len(completed), "pose_count": len(poses),
            "valid_poses": sum(pose.get("validation", {}).get("passed") is True for pose in poses),
            "invalid_or_unverified_poses": sum(pose.get("validation", {}).get("passed") is not True for pose in poses),
            "failed_ligands": len(failures), "not_completed_ligands": len(requested - completed - set(failures)),
            "failure_codes": dict(Counter(failure.get("code", "unknown") for failure in failures.values())),
            "score_units": "kcal/mol", "score_direction": "lower_is_better", "score_range": [min(finite), max(finite)] if finite else None,
            "protocol_fingerprint": result.get("metadata", {}).get("fingerprint"),
            "interpretation": "Vina score is a model score, not measured affinity. Physical validity is a separate check; failures remain visible."}


def export_result(result, selected, directory):
    """Export only hash-verified canonical poses; failed physical checks stay marked."""
    from rdkit import Chem
    directory.mkdir(parents=True, exist_ok=True)
    by_id = {row["compound_id"]: row for row in selected}
    poses = result.get("poses", [])
    csv_rows, all_mols, valid_mols = [], [], []
    for pose in poses:
        if pose["compound_id"] not in by_id:
            raise ValueError("Export contains an unrequested compound")
        path = Path(pose["pose_path"])
        if science.sha256(path) != pose["pose_sha256"]:
            raise ValueError("Pose hash changed; export blocked")
        molecules = list(Chem.SDMolSupplier(str(path), removeHs=False))
        if len(molecules) != 1 or molecules[0] is None or science._canonical(molecules[0]) != by_id[pose["compound_id"]]["canonical_smiles"]:
            raise ValueError("Pose identity mismatch; export blocked")
        mol = molecules[0]
        valid = pose.get("validation", {}).get("passed") is True
        mol.SetProp("BIOCENTRIC_PHYSICAL_CHECK_PASSED", str(valid).lower())
        mol.SetProp("BIOCENTRIC_PROTOCOL_FINGERPRINT", pose.get("protocol_fingerprint", result.get("metadata", {}).get("fingerprint", "")))
        mol.SetProp("BIOCENTRIC_SOURCE_POSE_SHA256", pose["pose_sha256"])
        all_mols.append(mol)
        if valid:
            valid_mols.append(mol)
        csv_rows.append({"name": pose["name"], "compound_id": pose["compound_id"], "status": "pose_generated", "score_kcal_mol": pose["score"], "score_direction": pose["score_direction"], "physical_check_passed": valid, "reason": "" if valid else "PoseBusters failed or incomplete", "pose_sha256": pose["pose_sha256"], "canonical_smiles": by_id[pose["compound_id"]]["canonical_smiles"]})
    completed = {row["compound_id"] for row in csv_rows}
    failures = {failure.get("compound_id"): failure for failure in result.get("failures", [])}
    for row in selected:
        if row["compound_id"] not in completed:
            failure = failures.get(row["compound_id"], {})
            csv_rows.append({"name": row["name"], "compound_id": row["compound_id"], "status": "failed" if failure else "not_completed", "score_kcal_mol": "", "score_direction": "lower_is_better", "physical_check_passed": "", "reason": failure.get("reason", "No completed pose"), "pose_sha256": "", "canonical_smiles": row["canonical_smiles"]})
    fields = ["name", "compound_id", "status", "score_kcal_mol", "score_direction", "physical_check_passed", "reason", "pose_sha256", "canonical_smiles"]
    csv_path = directory / "ligand_results.csv"
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(csv_rows)
    integrity = {"csv_rows": len(csv_rows), "all_pose_records": len(all_mols), "valid_pose_records": len(valid_mols)}
    files = [csv_path]
    for name, molecules in [("all_docked_poses.sdf", all_mols), ("physically_valid_poses.sdf", valid_mols)]:
        path = directory / name
        with Chem.SDWriter(str(path)) as writer:
            for mol in molecules:
                writer.write(mol)
        recovered = list(Chem.SDMolSupplier(str(path), removeHs=False)) if molecules else []
        if len(recovered) != len(molecules) or any(mol is None for mol in recovered):
            raise ValueError("SDF export record count failed")
        files.append(path)
    integrity["artifacts"] = [science._artifact(path, "verification_export") for path in files]
    return integrity


def partial_result(directory):
    request = json.loads((directory / "request.json").read_text()) if (directory / "request.json").exists() else {}
    poses, artifacts = [], []
    for path in sorted(directory.glob("*/completed.json")):
        saved = json.loads(path.read_text())
        if saved["fingerprint"] == request.get("fingerprint"):
            verify_artifacts(saved["artifacts"])
            poses.extend(saved["poses"])
            artifacts.extend(saved["artifacts"])
    return {"status": "budget_stopped", "poses": poses, "failures": [], "artifacts": artifacts, "metadata": {"fingerprint": request.get("fingerprint"), "config": request.get("config")}}


def main(argv=None):
    from rdkit import Chem, RDConfig
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--count", type=int, default=128)
    parser.add_argument("--budget-seconds", type=int, default=1200)
    args = parser.parse_args(argv)
    if not 8 <= args.count <= 256 or not 60 <= args.budget_seconds <= 3600:
        parser.error("count must be 8–256; budget-seconds must be 60–3600")
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    source = Path(RDConfig.RDDataDir) / "NCI/first_5K.smi"
    receptor = REPO / "examples/1iep/originals/1iep_receptor.pdbqt"
    reference = REPO / "examples/1iep/originals/1iep_ligand.sdf"
    receptor_pdb = REPO / "examples/1iep/originals/1iep_receptorH.pdb"
    pinned = {"count": args.count, "selection_seed": SEED, "source_sha256": science.sha256(source), "science_adapter_sha256": science.sha256(Path(science.__file__)), "script_sha256": science.sha256(Path(__file__)),
              "versions": {name: importlib.metadata.version(name) for name in ("rdkit", "meeko", "posebusters")},
              "receptor_sha256": science.sha256(receptor), "reference_sha256": science.sha256(reference), "receptor_pdb_sha256": science.sha256(receptor_pdb)}
    immutable_text(output / "pinned_request.json", json.dumps(pinned, indent=2) + "\n")
    if (output / "summary.json").exists():
        previous = json.loads((output / "summary.json").read_text())
        if previous.get("status") == "completed":
            verify_artifacts(previous["evidence_inventory"])
            print(json.dumps({"status": "existing_completed_evidence_verified", "report": str(output / "summary.json"), "baseline": previous["protocols"]["baseline"]}, indent=2), flush=True)
            return 0
    started = time.monotonic()
    cpu_before = resource.getrusage(resource.RUSAGE_SELF)
    child_before = resource.getrusage(resource.RUSAGE_CHILDREN)
    counts, active = Counter(), {"protocol": "preparation"}
    def expired():
        return time.monotonic() - started >= args.budget_seconds
    def emit(event):
        event = {"timestamp": now(), "elapsed_seconds": round(time.monotonic() - started, 3), "protocol": active["protocol"], **event}
        with (output / "events.jsonl").open("a") as handle:
            handle.write(json.dumps(event, allow_nan=False, default=str) + "\n")
        if event.get("tool") == "AutoDock Vina" and event["kind"] == "tool_started":
            counts["vina_started_" + active["protocol"]] += 1
        if event.get("data", {}).get("checkpoint_reused"):
            counts["checkpoint_reused_" + active["protocol"]] += 1
        if event.get("tool") == "AutoDock Vina" and event["kind"] == "tool_completed":
            counts["completed_" + active["protocol"]] += 1
            if counts["completed_" + active["protocol"]] % 8 == 0:
                print(json.dumps({"protocol": active["protocol"], "completed": counts["completed_" + active["protocol"]], "elapsed_seconds": round(time.monotonic() - started, 1)}), flush=True)
        if event["kind"] == "compound_failed":
            print(json.dumps({"protocol": active["protocol"], "failed": event.get("compound_id"), "reason": event["message"]}), flush=True)
    selection_path = output / "selection.json"
    if selection_path.exists():
        selection = json.loads(selection_path.read_text())
        if selection["source_sha256"] != pinned["source_sha256"] or selection["requested"] != args.count:
            raise ValueError("Saved selection differs from pinned inputs")
    else:
        selection = select_diverse(source, args.count)
        save(selection_path, selection)
    selected = selection["selected"]
    smi_path = output / "selected_128.smi" if args.count == 128 else output / f"selected_{args.count}.smi"
    immutable_text(smi_path, "".join(f'{row["canonical_smiles"]}\t{row["name"]}\n' for row in selected))
    print(json.dumps({"status": "selected", "count": len(selected), "source_records": selection["input_records"], "eligible_scaffolds": selection["eligible_unique_scaffolds"], "source_sha256": selection["source_sha256"]}), flush=True)
    preparation = prepare_cached(smi_path, output / "preparation", emit, expired)
    prepared = Path(preparation["prepared_sdf"])
    adversarial_path = output / "adversarial.smi"
    immutable_text(adversarial_path, ADVERSARIAL)
    adversarial = prepare_cached(adversarial_path, output / "adversarial_preparation", emit, expired)
    adversarial_ok = len(adversarial["compounds"]) == 2 and {row["code"] for row in adversarial["exclusions"]} == EXPECTED_ADVERSARIAL
    save(output / "adversarial_verification.json", {"passed": adversarial_ok, "accepted": len(adversarial["compounds"]), "exclusions": adversarial["exclusions"], "expected_exclusion_codes": sorted(EXPECTED_ADVERSARIAL)})
    if not adversarial_ok:
        raise RuntimeError("Adversarial preparation did not match explicit expected behavior")
    molecules = list(Chem.SDMolSupplier(str(prepared), removeHs=False))
    subset_path = output / "comparison_8.sdf"
    subset_buffer = io.StringIO()
    with Chem.SDWriter(subset_buffer) as writer:
        for mol in molecules[:8]:
            writer.write(mol)
    immutable_text(subset_path, subset_buffer.getvalue())
    config = {"cpu": 1, "exhaustiveness": 1, "n_poses": 1, "seed": SEED, "timeout_seconds": 60, "receptor_pdb": str(receptor_pdb), "require_posebusters": True, "run_plip": False}
    requested_ids = [row["compound_id"] for row in preparation["compounds"]]
    protocols, exports, results, durations = {}, {}, {}, {}
    variants = [("baseline", prepared, config, requested_ids), ("seed_43", subset_path, {**config, "seed": 43}, requested_ids[:8]), ("exhaustiveness_4", subset_path, {**config, "exhaustiveness": 4}, requested_ids[:8])]
    for name, library, parameters, ids in variants:
        if expired():
            protocols[name] = {"status": "not_run_budget", "requested_ligands": len(ids)}
            continue
        active["protocol"] = name
        before = time.monotonic()
        try:
            result = science.dock_library(receptor, library, reference, output / name, config=parameters, emit=emit, cancelled=expired)
        except science.ScienceError as error:
            if error.code != "cancelled":
                raise
            result = partial_result(output / name)
        durations[name] = round(time.monotonic() - before, 3)
        results[name] = result
        save(output / f"{name}_result.json", result)
        protocols[name] = {"status": result["status"], **summarize_result(result, ids), "duration_seconds": durations[name], "config": parameters}
        chosen = [row for row in selected if row["compound_id"] in set(ids)]
        exports[name] = export_result(result, chosen, output / "exports" / name)
        print(json.dumps({"protocol_finished": name, **protocols[name]}), flush=True)
    reuse = {"status": "not_run_budget"}
    if "baseline" in results and results["baseline"]["status"] != "budget_stopped" and not expired():
        active["protocol"] = "baseline_resume"
        pose_signatures = {pose["pose_path"]: (pose["pose_sha256"], Path(pose["pose_path"]).stat().st_mtime_ns) for pose in results["baseline"]["poses"]}
        before = time.monotonic()
        try:
            resumed = science.dock_library(receptor, prepared, reference, output / "baseline", config=config, emit=emit, cancelled=expired)
            unchanged = all(science.sha256(Path(path)) == signature[0] and Path(path).stat().st_mtime_ns == signature[1] for path, signature in pose_signatures.items())
            reused = counts["checkpoint_reused_baseline_resume"]
            expected = len({pose["compound_id"] for pose in results["baseline"]["poses"]})
            reuse = {"status": "passed" if unchanged and reused == expected else "failed", "expected_successful_checkpoints": expected, "reused_checkpoints": reused, "new_vina_calls": counts["vina_started_baseline_resume"], "pose_hashes_and_mtimes_unchanged": unchanged, "duration_seconds": round(time.monotonic()-before, 3), "failed_ligands_are_retried": True}
        except science.ScienceError as error:
            if error.code != "cancelled":
                raise
            reuse = {"status": "budget_stopped", "reused_checkpoints": counts["checkpoint_reused_baseline_resume"]}
    save(output / "checkpoint_reuse.json", reuse)
    cpu_after = resource.getrusage(resource.RUSAGE_SELF)
    child_after = resource.getrusage(resource.RUSAGE_CHILDREN)
    comparisons = []
    baseline_by_id = {pose["compound_id"]: pose for pose in results.get("baseline", {}).get("poses", [])}
    for name in ("seed_43", "exhaustiveness_4"):
        for pose in results.get(name, {}).get("poses", []):
            original = baseline_by_id.get(pose["compound_id"])
            if original:
                comparisons.append({"compound_id": pose["compound_id"], "name": pose["name"], "variant": name, "baseline_score_kcal_mol": original["score"], "variant_score_kcal_mol": pose["score"], "difference_variant_minus_baseline": round(pose["score"] - original["score"], 6), "baseline_physical_check": original["validation"].get("passed"), "variant_physical_check": pose["validation"].get("passed"), "interpretation": "Parameter sensitivity observation, not accuracy or affinity validation"})
    save(output / "parameter_comparisons.json", comparisons)
    all_finished = len(protocols) == 3 and all(row.get("not_completed_ligands", 1) == 0 for row in protocols.values()) and reuse["status"] == "passed"
    summary = {"status": "completed" if all_finished else "partial", "created_at": now(), "pinned_inputs": pinned,
               "target": "Curated 1IEP mouse ABL pocket; bound ligand defines box only", "selection": {key: value for key, value in selection.items() if key not in {"selected", "source_ledger"}},
               "preparation": {"requested": len(selected), "accepted": len(preparation["compounds"]), "excluded": len(preparation["exclusions"]), "exclusions": preparation["exclusions"]},
               "adversarial_passed": adversarial_ok, "protocols": protocols, "exports": exports, "checkpoint_reuse": reuse,
               "resources": {"wall_seconds": round(time.monotonic()-started, 3), "budget_seconds": args.budget_seconds, "cpu_per_vina_process": 1, "vina_processes_started": sum(value for key, value in counts.items() if key.startswith("vina_started_")), "process_counts": dict(counts), "python_cpu_seconds": round(cpu_after.ru_utime+cpu_after.ru_stime-cpu_before.ru_utime-cpu_before.ru_stime, 3), "child_cpu_seconds": round(child_after.ru_utime+child_after.ru_stime-child_before.ru_utime-child_before.ru_stime, 3), "python_peak_rss_bytes": cpu_after.ru_maxrss if sys.platform == "darwin" else cpu_after.ru_maxrss*1024, "child_max_rss_bytes": child_after.ru_maxrss if sys.platform == "darwin" else child_after.ru_maxrss*1024, "hardware": platform.platform(), "gpu_used": False, "network_requests": 0, "llm_calls": 0, "plip_calls": 0},
               "limitations": ["This is a bounded software/diversity test, not a target-activity or affinity benchmark.", "Input charge and tautomer are retained; no pH/state enumeration or ADMET prediction.", "One pose and low exhaustiveness deliberately trade search depth for coverage; failed physical checks are not hits.", "Parameter comparisons use the same compounds and receptor but separate run fingerprints; no combined rank or superior-engine claim.", "No measured labels, MMPBSA/MMGBSA trajectory, GPU docking, or experimental validation."]}
    inventory = [science._artifact(path, "verification_evidence") for path in sorted(output.rglob("*")) if path.is_file() and path.name != "summary.json"]
    summary["evidence_inventory"] = inventory
    summary["resources"]["disk_bytes"] = sum(row["bytes"] for row in inventory)
    save(output / "summary.json", summary)
    print(json.dumps({"report": str(output / "summary.json"), "status": summary["status"], "resources": summary["resources"], "checkpoint_reuse": reuse}), flush=True)
    return 0 if all_finished else 2


if __name__ == "__main__":
    raise SystemExit(main())
