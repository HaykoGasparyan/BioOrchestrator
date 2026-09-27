"""Deterministic preparation of registered native inputs after protocol review.

The returned manifest still goes through native preflight and its durable review
gate. This creates files; it does not execute docking or infer scientific approval.
"""
from __future__ import annotations

import io
import json
from pathlib import Path

from . import science
from .store import digest, dump


def validate_input_structure_provenance(store, pid, record, prepared_artifact, receptor_artifact, resolve):
    """Bind each prepared affinity structure to an intact validated source pose.

    ``resolve`` registers source artifacts in the native plan so their integrity
    is checked again after execution. Physical validity describes the input
    structure; it does not describe affinity accuracy or a model-generated pose.
    """
    from rdkit import Chem
    if (record.get("prepared_ligands_artifact_id") != prepared_artifact["id"]
            or record.get("prepared_ligands_sha256") != prepared_artifact["sha256"]
            or record.get("receptor_artifact_id") != receptor_artifact["id"]
            or record.get("receptor_sha256") != receptor_artifact["sha256"]):
        raise ValueError("Input-structure provenance does not bind the prepared ligands and receptor.")
    molecules = list(Chem.ForwardSDMolSupplier(io.BytesIO(store.read_artifact(pid, prepared_artifact)), removeHs=False))
    entries = record.get("entries", [])
    if not molecules or len(molecules) != len(entries) or any(mol is None for mol in molecules):
        raise ValueError("Input-structure provenance requires one entry per prepared compound.")
    artifacts = {row["id"]: row for row in store.artifacts(pid)}
    data = store.project(pid).get("data", {})
    recorded_poses = data.get("screening", {}).get("poses", [])
    ids = []
    for index, (molecule, entry) in enumerate(zip(molecules, entries), 1):
        identity = science._canonical(molecule)
        cid = "cmp_" + digest(identity.encode())[:16]
        if (entry.get("prepared_record_index") != index or entry.get("compound_id") != cid
                or entry.get("canonical_smiles") != identity or not molecule.HasProp("_Name")
                or molecule.GetProp("_Name") != cid):
            raise ValueError("Input-structure provenance compound identity/order differs from prepared structures.")
        ids.append(cid)
        source_fields = {"source_pose_sha256": "pose_sha256", "source_pose_id": "pose_id",
                         "source_pose_number": "pose_number", "source_engine": "engine",
                         "source_receptor_sha256": "receptor_sha256",
                         "source_protocol_fingerprint": "protocol_fingerprint"}
        if not any(pose.get("compound_id") == cid and pose.get("canonical_smiles") == identity
                   and all(entry.get(key) == pose.get(field) for key, field in source_fields.items())
                   and pose.get("validation", {}).get("passed") is True
                   and (pose.get("validation", {}).get("artifact") or {}).get("sha256") == entry.get("source_validation_sha256")
                   for pose in recorded_poses):
            raise ValueError("Input-structure provenance does not match a recorded screening pose and its physical validation.")
        source_path = resolve(f"input_source_pose_{index}", entry["source_pose_artifact_id"])
        source = artifacts[entry["source_pose_artifact_id"]]
        original = list(Chem.SDMolSupplier(source_path, removeHs=False))
        if (source["sha256"] != entry.get("source_pose_sha256") or len(original) != 1 or original[0] is None
                or science._canonical(original[0]) != identity
                or Chem.MolToMolBlock(original[0]).splitlines()[3:] != Chem.MolToMolBlock(molecule).splitlines()[3:]):
            raise ValueError("Prepared affinity coordinates/chemistry differ from the mapped source pose.")
        report_path = resolve(f"input_source_validation_{index}", entry["source_validation_artifact_id"])
        report_artifact = artifacts[entry["source_validation_artifact_id"]]
        report = json.loads(Path(report_path).read_text())
        if (report_artifact["sha256"] != entry.get("source_validation_sha256") or report.get("passed") is not True
                or entry.get("source_physical_validation") != {key: report.get(key) for key in
                    ("status", "passed", "tool", "version", "mode", "receptor_checks_run", "failed_checks")}):
            raise ValueError("Mapped physical validation is absent, failed or differs from its registered report.")
        # A PDBQT source pose may have been checked against the corresponding
        # reviewed PDB. Preserve both source and model receptor hashes explicitly.
        protocol = data.get("protocol", {})
        reviewed_receptors = set(protocol.get("input_hashes", {}).get(key) for key in ("receptor", "structure"))
        if entry.get("source_receptor_sha256") not in reviewed_receptors:
            raise ValueError("Input pose provenance names a different reviewed receptor.")
    if len(ids) != len(set(ids)):
        raise ValueError("Duplicate compounds in input-structure provenance.")
    return record


def auto_prepare_native_manifest(store, pid, engine, phase, protocol, *, job_id=None, emit=None, cancelled=None):
    from rdkit import Chem
    from .engine_workflow import NativeManifest
    if engine not in {"gnina", "unidock", "gnina_affinity", "biocentric_affinity"}:
        raise ValueError("This engine requires separately reviewed native scientific inputs.")
    if (engine in {"gnina", "unidock"}) != (phase in {"pilot", "screening"}):
        raise ValueError("Native preparation engine and phase differ.")
    data = store.project(pid)["data"]
    frozen = data["research"]["testset"]
    from .testsets import validate_testset
    validate_testset(frozen)
    if frozen["manifest_sha256"] != protocol["testset_sha256"]:
        raise ValueError("Preparation must use the reviewed frozen test set.")
    artifacts = store.artifacts(pid)
    def source(kind):
        expected = protocol["input_hashes"][kind]
        artifact = next((a for a in reversed(artifacts) if a["sha256"] == expected and a["kind"] == kind), None)
        if artifact is None:
            raise ValueError("Missing reviewed input artifact: " + kind)
        store.read_artifact(pid, artifact)
        return artifact
    def register(name, content):
        return store.artifact(pid, name, content, "native_input", job_id)
    receptor, reference = source("receptor" if engine == "unidock" else "structure"), source("reference")
    protocol_sha = digest(dump(protocol).encode())
    work = store.work_directory(pid) / "native-preparation" / engine / phase / protocol_sha
    work.mkdir(parents=True, exist_ok=True)
    exclusions = []
    input_structure_provenance_artifact = None
    if phase in {"pilot", "screening"}:
        original = source("reference" if phase == "pilot" else "library")
        prep = science.prepare_library(Path(original["path"]), work / "preparation", protocol["seed"], emit, cancelled)
        for evidence in prep["artifacts"]:
            register(Path(evidence["path"]).name, Path(evidence["path"]).read_bytes())
        if not prep["compounds"]:
            raise ValueError("No compatible ligands survived preparation; exclusions are preserved.")
        content = Path(prep["prepared_sdf"]).read_bytes()
        mols = list(Chem.ForwardSDMolSupplier(io.BytesIO(content), removeHs=False))
        canonical = register("canonical-native-" + phase + ".sdf", content)
        if phase == "screening":
            path = Path(original["path"])
            originals = (list(Chem.SDMolSupplier(str(path), removeHs=False)) if path.suffix.lower() in {".sdf", ".sd"}
                         else [Chem.MolFromSmiles(line.split()[0]) for line in path.read_text().splitlines() if line.strip() and not line.lstrip().startswith("#")])
            selected = {science._canonical(mol) for mol in mols}
            reasons = {row["record_index"]: row["reason"] for row in prep["exclusions"]}
            for index, mol in enumerate(originals, 1):
                identity = science._canonical(mol) if mol else None
                if identity not in selected:
                    exclusions.append({"record_index": index, "canonical_smiles": identity, "reason": reasons.get(index, "Excluded by reviewed preparation; original retained.")})
    else:
        # Score one exact, registered, physically passing pose per compound.
        poses = data.get("screening", {}).get("poses", [])
        chosen = {}
        for pose in poses:
            if pose.get("validation", {}).get("passed") is not True:
                continue
            cid = pose["compound_id"]
            if cid not in chosen or pose.get("pose_number", 999) < chosen[cid].get("pose_number", 999):
                chosen[cid] = pose
        if not chosen:
            raise ValueError("Affinity preparation needs actual registered passing structures; no poses are invented.")
        mols, pose_sources = [], []
        for cid, pose in sorted(chosen.items()):
            expected_cid = "cmp_" + digest(pose["canonical_smiles"].encode())[:16]
            if cid != expected_cid:
                raise ValueError("Affinity compound ID differs from its canonical chemical identity.")
            if pose.get("receptor_sha256") not in {protocol["input_hashes"].get("receptor"), protocol["input_hashes"].get("structure")}:
                raise ValueError("Affinity pose was not calculated against the reviewed receptor.")
            artifact = next((a for a in artifacts if a["sha256"] == pose.get("pose_sha256") and a["path"].endswith(".sdf")), None)
            if artifact is None:
                raise ValueError("Affinity pose is not registered with its content hash.")
            records = list(Chem.ForwardSDMolSupplier(io.BytesIO(store.read_artifact(pid, artifact)), removeHs=False))
            if len(records) != 1 or records[0] is None or science._canonical(records[0]) != pose["canonical_smiles"]:
                raise ValueError("Affinity pose identity differs from its recorded compound.")
            validation_hash = (pose.get("validation", {}).get("artifact") or {}).get("sha256")
            validation_artifact = next((a for a in artifacts if a["sha256"] == validation_hash), None)
            if validation_artifact is None:
                raise ValueError("Affinity input requires the registered physical-validation report for its source pose.")
            validation_report = json.loads(store.read_artifact(pid, validation_artifact))
            if validation_report.get("passed") is not True:
                raise ValueError("The registered source-pose physical-validation report did not pass.")
            records[0].SetProp("_Name", cid)
            mols.append(records[0])
            pose_sources.append({"compound_id": cid, "canonical_smiles": pose["canonical_smiles"],
                                 "prepared_record_index": len(mols), "source_pose_artifact_id": artifact["id"],
                                 "source_pose_sha256": artifact["sha256"], "source_pose_id": pose.get("pose_id"),
                                 "source_pose_number": pose.get("pose_number"), "source_engine": pose.get("engine"),
                                 "source_receptor_sha256": pose["receptor_sha256"],
                                 "source_protocol_fingerprint": pose.get("protocol_fingerprint"),
                                 "source_validation_artifact_id": validation_artifact["id"],
                                 "source_validation_sha256": validation_artifact["sha256"],
                                 "source_physical_validation": {key: validation_report.get(key) for key in
                                     ("status", "passed", "tool", "version", "mode", "receptor_checks_run", "failed_checks")}})
        buffer = io.StringIO()
        with Chem.SDWriter(buffer) as writer:
            for mol in mols:
                writer.write(mol)
        canonical = register("affinity-input-structures.sdf", buffer.getvalue().encode())
        pose_provenance = {"schema_version": 1, "prepared_ligands_artifact_id": canonical["id"],
                           "prepared_ligands_sha256": canonical["sha256"], "receptor_artifact_id": receptor["id"],
                           "receptor_sha256": receptor["sha256"], "entries": pose_sources,
                           "selection_policy": "Lowest pose_number among registered physically passing poses per compound.",
                           "physical_validation_scope": "Prior input-structure checks only; not affinity validation or model-generated-pose validity."}
        input_structure_provenance_artifact = register("affinity-input-pose-provenance.json", json.dumps(pose_provenance, indent=2).encode())["id"]
        original = source("library")
        original_path = Path(original["path"])
        original_molecules = (list(Chem.SDMolSupplier(str(original_path), removeHs=False)) if original_path.suffix.lower() in {".sdf", ".sd"}
                              else [Chem.MolFromSmiles(line.split()[0]) for line in original_path.read_text().splitlines() if line.strip() and not line.lstrip().startswith("#")])
        selected = {science._canonical(mol) for mol in mols}
        exclusions = [{"record_index": index, "canonical_smiles": science._canonical(mol) if mol else None,
                       "reason": "No registered physically passing structure is available for this input; affinity was not imputed."}
                      for index, mol in enumerate(original_molecules, 1) if mol is None or science._canonical(mol) not in selected]
    files = {"receptor": receptor["id"]}
    review = {"preparation_reviewed": True}
    parameters = {"timeout_seconds": protocol["timeout_seconds"]}
    protein_identity_artifact = None
    if engine == "unidock":
        from meeko import MoleculePreparation, PDBQTWriterLegacy
        native = []
        for index, mol in enumerate(mols):
            science._cancel(cancelled)
            setups = MoleculePreparation().prepare(mol)
            if len(setups) != 1:
                raise ValueError("Uni-Dock input preparation produced ambiguous states.")
            text, valid, reason = PDBQTWriterLegacy.write_string(setups[0])
            if not valid:
                raise ValueError("Meeko preparation failed: " + reason)
            native.append(register(f"ligand_{index+1:06d}.pdbqt", text.encode())["id"])
        files["ligands"] = native
        parameters.update(box={key: protocol["box"][key] for key in ("center", "size")}, seed=protocol["seed"], num_modes=protocol["n_poses"], search_mode=protocol.get("search_mode", "balance"))
    elif engine in {"gnina", "gnina_affinity"}:
        files["ligand"] = canonical["id"]
        parameters.update(cpu=protocol["cpu"], device=protocol.get("device", "gpu"))
        if engine == "gnina":
            parameters.update(box={key: protocol["box"][key] for key in ("center", "size")}, seed=protocol["seed"], num_modes=protocol["n_poses"], exhaustiveness=protocol["exhaustiveness"])
    else:
        files["ligands"] = canonical["id"]
        identity = data.get("target_identity") or next((row["value"] for row in reversed(store.decisions(pid)) if row["kind"] == "identity"), {})
        if not identity.get("accession") or not identity.get("organism"):
            raise ValueError("BioCentric requires reviewed protein accession and species; resolve identity before inference.")
        chains = sorted({line[21:22].strip() for line in Path(receptor["path"]).read_text().splitlines() if line.startswith("ATOM  ")})
        identity_record = {"reviewed_identity": identity, "receptor_artifact_id": receptor["id"], "receptor_sha256": receptor["sha256"], "structure_chains": chains,
                           "construct_mapping": identity.get("construct", "Not independently established; exact structure hash is authoritative")}
        protein_identity_artifact = register("affinity-protein-identity.json", json.dumps(identity_record, indent=2).encode())["id"]
        review["protein_id"] = f"{identity['accession']}:{identity['organism']}:chains={','.join(chains)}:structure_sha256={receptor['sha256']}"
        parameters.update(device="cuda", cpu=protocol["cpu"], max_pairs=len(mols))
    manifest = NativeManifest(engine_id=engine, phase=phase, input_artifacts=files, parameters=parameters,
                              review=review, human_reviewed=True, exclusions=exclusions,
                              canonical_ligands_artifact=canonical["id"], validation_receptor_artifact=source("structure")["id"],
                              reference_artifact=reference["id"] if phase == "pilot" else None,
                              protein_identity_artifact=protein_identity_artifact,
                              input_structure_provenance_artifact=input_structure_provenance_artifact,
                              testset_sha256=frozen["manifest_sha256"], protocol_sha256=protocol_sha)
    return store.artifact(pid, f"{engine}-{phase}-manifest.json", manifest.model_dump_json(indent=2).encode(), "native_engine_" + phase, job_id)
