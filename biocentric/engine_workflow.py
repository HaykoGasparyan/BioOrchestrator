"""Bridge registered project artifacts to native scientific engine adapters.

This is the configured-engine route for the research workflow. A manifest can
reference only immutable artifacts belonging to its project, never shell text or
arbitrary filesystem paths. Native calculations and canonical pose validation
are separate recorded operations.
"""
from __future__ import annotations

import hashlib
import io
import json
import math
from copy import deepcopy
from pathlib import Path
from typing import Callable, Literal

from pydantic import BaseModel, ConfigDict, Field

from . import engines, science
from .native_preparation import auto_prepare_native_manifest


class NativeExclusion(BaseModel):
    model_config = ConfigDict(extra="forbid")
    record_index: int = Field(ge=1)
    canonical_smiles: str | None = None
    reason: str = Field(min_length=10, max_length=1000)


class NativeManifest(BaseModel):
    model_config = ConfigDict(extra="forbid")
    schema_version: Literal[1] = 1
    engine_id: Literal["vina", "autodock4", "gnina", "unidock", "mmpbsa", "mmgbsa", "gnina_affinity", "planet", "biocentric_affinity"]
    phase: Literal["pilot", "screening", "affinity"]
    input_artifacts: dict[str, str | list[str] | dict[str, str]]
    parameters: dict = Field(default_factory=dict)
    review: dict = Field(default_factory=dict)
    human_reviewed: Literal[True]
    exclusions: list[NativeExclusion] = Field(default_factory=list)
    canonical_ligands_artifact: str | None = None
    validation_receptor_artifact: str | None = None
    reference_artifact: str | None = None
    protein_identity_artifact: str | None = None
    input_structure_provenance_artifact: str | None = None
    testset_sha256: str = Field(pattern=r"^[a-f0-9]{64}$")
    protocol_sha256: str = Field(pattern=r"^[a-f0-9]{64}$")


def _artifact(store, pid, artifact_id):
    known = next((item for item in store.artifacts(pid) if item["id"] == artifact_id), None)
    if not known:
        raise ValueError("Native manifest references an artifact not registered to this project.")
    store.read_artifact(pid, known)
    return known


def _coordinate_frame(native: Path, validation: Path):
    def atoms(path):
        return [line for line in path.read_text().splitlines() if line.startswith(("ATOM  ", "HETATM"))]
    def coordinates(lines):
        return {tuple(round(float(line[start:start+8]), 2) for start in (30, 38, 46)) for line in lines}
    native_atoms = atoms(native)
    if native.suffix.lower() == ".pdbqt":
        native_atoms = [line for line in native_atoms if line.split()[-1] not in {"H", "HD", "HS"}]
    else:
        native_atoms = [line for line in native_atoms if line[76:78].strip() != "H"]
    native_coords, validation_coords = coordinates(native_atoms), coordinates(atoms(validation))
    if not native_coords or len(native_coords & validation_coords) / len(native_coords) < .95:
        raise ValueError("Native and validation receptors do not share >=95% heavy-atom coordinates at 0.01 Å.")


def plan_native_manifest(store, pid: str, manifest_artifact, output_dir: str | Path, *,
                         expected_engine: str, expected_phase: str,
                         expected_testset_sha256: str | None = None,
                         expected_protocol_sha256: str | None = None,
                         resume_plan: dict | None = None) -> dict:
    """Resolve and preflight an uploaded job manifest without creating/running anything."""
    from rdkit import Chem
    manifest_id = manifest_artifact["id"] if isinstance(manifest_artifact, dict) else str(manifest_artifact)
    manifest_source = _artifact(store, pid, manifest_id)
    expected_kind = "native_engine_" + expected_phase
    if manifest_source["kind"] != expected_kind:
        raise ValueError(f"Expected a registered {expected_kind} manifest.")
    manifest = NativeManifest.model_validate_json(store.read_artifact(pid, manifest_source))
    if manifest.engine_id != expected_engine or manifest.phase != expected_phase:
        raise ValueError("Native manifest engine/phase does not match the approved research plan.")
    if expected_testset_sha256 and manifest.testset_sha256 != expected_testset_sha256:
        raise ValueError("Native manifest does not bind the frozen test set.")
    if expected_protocol_sha256 and manifest.protocol_sha256 != expected_protocol_sha256:
        raise ValueError("Native manifest does not bind the reviewed protocol.")
    docking = manifest.engine_id in engines.ENGINE_IDS[:4]
    if docking != (manifest.phase in {"pilot", "screening"}):
        raise ValueError("Docking engines require pilot/screening; affinity engines require affinity phase.")
    resolved_artifacts = {"manifest": manifest_source}
    def resolve(role, value):
        if isinstance(value, str):
            artifact = _artifact(store, pid, value)
            resolved_artifacts[role] = artifact
            return artifact["path"]
        if isinstance(value, list):
            return [resolve(f"{role}_{i}", item) for i, item in enumerate(value)]
        return {key: resolve(f"{role}_{key}", item) for key, item in value.items()}
    inputs = {role: resolve(role, value) for role, value in manifest.input_artifacts.items()}
    if set(inputs) & set(manifest.review):
        raise ValueError("Review fields cannot override artifact-bound scientific inputs.")
    inputs.update(manifest.review)
    if manifest.protein_identity_artifact:
        resolve("protein_identity", manifest.protein_identity_artifact)
    if manifest.engine_id == "biocentric_affinity":
        if not manifest.protein_identity_artifact:
            raise ValueError("BioCentric research inference requires registered reviewed protein identity evidence.")
        identity_record = json.loads(store.read_artifact(pid, resolved_artifacts["protein_identity"]))
        identity = identity_record.get("reviewed_identity", {})
        receptor_record = resolved_artifacts["receptor"]
        chains = sorted({line[21:22].strip() for line in Path(inputs["receptor"]).read_text().splitlines() if line.startswith("ATOM  ")})
        expected_id = f"{identity.get('accession')}:{identity.get('organism')}:chains={','.join(chains)}:structure_sha256={receptor_record['sha256']}"
        approved_identity = store.project(pid).get("data", {}).get("target_identity") or next(
            (row["value"] for row in reversed(store.decisions(pid)) if row["kind"] == "identity"), {})
        if (not identity.get("accession") or not identity.get("organism")
                or identity_record.get("receptor_artifact_id") != receptor_record["id"]
                or identity_record.get("receptor_sha256") != receptor_record["sha256"]
                or identity_record.get("structure_chains") != chains
                or inputs.get("protein_id") != expected_id
                or identity != approved_identity):
            raise ValueError("BioCentric protein identity evidence does not bind the reviewed protein and receptor.")
    input_structure_provenance = None
    if manifest.input_structure_provenance_artifact:
        from .native_preparation import validate_input_structure_provenance
        resolve("input_structure_provenance", manifest.input_structure_provenance_artifact)
        provenance_record = json.loads(store.read_artifact(pid, resolved_artifacts["input_structure_provenance"]))
        if manifest.phase != "affinity":
            raise ValueError("Input-pose provenance is only applicable to affinity preparation.")
        ligand_role = "ligands" if manifest.engine_id in {"biocentric_affinity", "planet"} else "ligand"
        input_structure_provenance = validate_input_structure_provenance(
            store, pid, provenance_record, resolved_artifacts[ligand_role], resolved_artifacts["receptor"], resolve)
    out = Path(output_dir).resolve()
    if resume_plan is not None and (expected_engine != 'biocentric_affinity' or
                                   resume_plan.get('output_dir') != str(out)):
        raise ValueError('Only the original reviewed BioCentric output directory can resume.')
    if resume_plan is None and out.exists() and any(out.iterdir()):
        raise ValueError("Native workflow output directory must be new or empty.")
    validation = {}
    if docking or manifest.engine_id == "gnina_affinity":
        if not manifest.canonical_ligands_artifact or not manifest.validation_receptor_artifact:
            raise ValueError("Canonical ligand SDF and matching validation receptor PDB artifacts are required.")
        canonical = Path(resolve("canonical_ligands", manifest.canonical_ligands_artifact))
        receptor = Path(resolve("validation_receptor", manifest.validation_receptor_artifact))
        if canonical.suffix.lower() != ".sdf" or receptor.suffix.lower() != ".pdb":
            raise ValueError("Validation requires canonical SDF and receptor PDB formats.")
        molecules = list(Chem.SDMolSupplier(str(canonical), removeHs=False))
        if not molecules or any(mol is None for mol in molecules):
            raise ValueError("Canonical ligand SDF contains invalid molecules.")
        identities = [science._canonical(mol) for mol in molecules]
        if len(set(identities)) != len(identities):
            raise ValueError("Canonical ligand inputs must be unique chemical states.")
        for mol in molecules:
            science._supported(mol)
        _coordinate_frame(Path(inputs["receptor"]), receptor)
        inspection = science.inspect_structure(receptor)
        if inspection.get("routine_screening_blockers"):
            raise ValueError("Native receptor needs a special protocol: " + "; ".join(inspection["routine_screening_blockers"]))
        native_identities = []
        if manifest.engine_id in {"gnina", "gnina_affinity"}:
            native_molecules = list(Chem.SDMolSupplier(inputs["ligand"], removeHs=False))
        else:
            from meeko import PDBQTMolecule, RDKitMolCreate
            native_molecules = []
            native_files = inputs["ligands"] if manifest.engine_id == "unidock" else [inputs["ligand"]]
            for filename in native_files:
                sdf_text, conversion_failures = RDKitMolCreate.write_sd_string(PDBQTMolecule(Path(filename).read_text(), skip_typing=True))
                if conversion_failures:
                    raise ValueError("Native ligand PDBQT requires recoverable Meeko SMILES/atom mapping.")
                native_molecules.extend(Chem.ForwardSDMolSupplier(io.BytesIO(sdf_text.encode()), removeHs=False))
        if not native_molecules or any(mol is None for mol in native_molecules):
            raise ValueError("Native engine ligand chemistry could not be recovered before execution.")
        native_identities = [science._canonical(mol) for mol in native_molecules]
        if set(native_identities) != set(identities) or (docking and len(native_identities) != len(identities)):
            raise ValueError("Native engine ligand identities do not match the canonical input set.")
        validation = {"canonical_ligands": str(canonical), "receptor": str(receptor), "identities": identities}
        if manifest.phase == "pilot":
            if len(molecules) != 1 or not manifest.reference_artifact:
                raise ValueError("Reference pilot requires exactly one canonical ligand and its observed reference SDF.")
            reference = Path(resolve("reference", manifest.reference_artifact))
            refs = list(Chem.SDMolSupplier(str(reference), removeHs=False))
            if len(refs) != 1 or refs[0] is None or science._canonical(refs[0]) != identities[0]:
                raise ValueError("Pilot reference identity must match the canonical input ligand.")
            box = science.derive_box(reference)
            receptor_coords = [tuple(float(line[start:start+8]) for start in (30,38,46)) for line in receptor.read_text().splitlines() if line.startswith(("ATOM  ", "HETATM"))]
            if not any(sum((point[i]-box["center"][i])**2 for i in range(3)) < 12**2 for point in receptor_coords):
                raise ValueError("Pilot reference is outside the receptor coordinate frame.")
            validation["reference"] = str(reference)
    if manifest.engine_id in {"planet", "biocentric_affinity"}:
        query_molecules = list(Chem.SDMolSupplier(inputs["ligands"], removeHs=False))
        if not query_molecules or any(mol is None for mol in query_molecules):
            raise ValueError("Affinity inputs contain missing or invalid chemical identities.")
        query_identities = [science._canonical(mol) for mol in query_molecules]
        if len(set(query_identities)) != len(query_identities):
            raise ValueError("Affinity input chemical states must be unique.")
        validation = {"identities": query_identities}
    project_data = store.project(pid).get("data", {})
    frozen = project_data.get("research", {}).get("testset", {})
    frozen_binding = {"verified": False, "reason": "No project research test set is stored; caller must bind its reviewed context."}
    if frozen:
        from .testsets import validate_testset
        validate_testset(frozen)
        if frozen.get("manifest_sha256") != manifest.testset_sha256:
            raise ValueError("Native manifest does not match the project's stored frozen test set.")
        sources = {item["kind"]: item for item in frozen.get("input_artifacts", [])}
        if manifest.phase == "pilot":
            reference = sources.get("reference")
            registered = resolved_artifacts.get("reference")
            if not reference or not registered or reference["sha256"] != registered["sha256"]:
                raise ValueError("Native pilot must use the exact observed reference frozen before tool selection.")
        if manifest.phase == "screening" or (manifest.phase == "affinity" and manifest.engine_id in {"gnina_affinity", "planet", "biocentric_affinity"}):
            source = sources.get("library")
            if not source:
                raise ValueError("Frozen library artifact is missing from the research test set.")
            library = _artifact(store, pid, source["artifact_id"])
            if library["sha256"] != source["sha256"]:
                raise ValueError("Frozen library artifact hash has changed.")
            resolved_artifacts["frozen_library"] = library
            path = Path(library["path"])
            if path.suffix.lower() in {".sdf", ".sd"}:
                library_molecules = list(Chem.SDMolSupplier(str(path), removeHs=False))
            elif path.suffix.lower() in {".smi", ".smiles", ".txt"}:
                library_molecules = [Chem.MolFromSmiles(line.split()[0]) for line in path.read_text().splitlines() if line.strip() and not line.lstrip().startswith("#")]
            else:
                raise ValueError("Frozen screening library format is unsupported.")
            expected = [science._canonical(mol) if mol is not None else None for mol in library_molecules]
            selected = set(validation["identities"])
            if not selected <= set(expected):
                raise ValueError("Native screening introduces chemical states absent from the frozen library.")
            missing = {i: identity for i, identity in enumerate(expected, 1) if identity not in selected}
            supplied = {item.record_index: item for item in manifest.exclusions}
            if len(supplied) != len(manifest.exclusions) or set(supplied) != set(missing):
                raise ValueError("Every omitted frozen-library record needs one explicit exclusion and reason.")
            if any(supplied[index].canonical_smiles != identity for index, identity in missing.items()):
                raise ValueError("Exclusion chemical identities do not match the omitted frozen-library records.")
            frozen_binding = {"verified": True, "library_artifact_id": library["id"], "library_sha256": library["sha256"],
                              "selected_chemical_states": len(selected), "exclusions": [item.model_dump() for item in manifest.exclusions]}
        else:
            frozen_binding = {"verified": True, "testset_sha256": manifest.testset_sha256}
    resume_kwargs = {'resume_plan': resume_plan['engine_plan']} if resume_plan is not None else {}
    engine_plan = engines.plan_engine_run(manifest.engine_id, inputs, manifest.parameters, out / "engine", **resume_kwargs)
    protocol = project_data.get("protocol", {})
    if protocol:
        # Copying a protocol hash into a manifest cannot substitute its receptor.
        approved = protocol.get("input_hashes", {})
        native_receptor = resolved_artifacts.get("receptor")
        if native_receptor and native_receptor["sha256"] not in {approved.get("receptor"), approved.get("structure")}:
            raise ValueError("Native receptor differs from the reviewed protocol; freeze and review the changed preparation first.")
        validation_receptor = resolved_artifacts.get("validation_receptor")
        if validation_receptor and validation_receptor["sha256"] != approved.get("structure"):
            raise ValueError("Native validation receptor differs from the reviewed structure.")
    if manifest.phase == "screening":
        pilot = project_data.get("research", {}).get("pilot", {})
        native_pilot = pilot.get("native_result", {})
        if native_pilot:
            if native_pilot.get("engine_id") != manifest.engine_id:
                raise ValueError("Screening engine differs from the tested native pilot.")
            tested = {key: value for key, value in native_pilot["parameters"].items() if key != "timeout_seconds"}
            requested = {key: value for key, value in engine_plan["parameters"].items() if key != "timeout_seconds"}
            if tested != requested:
                raise ValueError("Screening search/scoring settings differ from the tested pilot; revise and repeat the pilot first.")
    plan = {"status": engine_plan["status"], "blockers": engine_plan["blockers"], "engine_id": manifest.engine_id,
            "phase": manifest.phase, "manifest": manifest.model_dump(), "manifest_artifact": manifest_source,
            "resolved_artifacts": resolved_artifacts, "engine_inputs": inputs, "engine_plan": engine_plan,
            "validation": validation, "frozen_input_binding": frozen_binding,
            "input_structure_provenance": input_structure_provenance, "output_dir": str(out)}
    if resume_plan is not None and plan != resume_plan:
        raise ValueError('Native resume differs from the immutable reviewed workflow plan; preserve results and review a new experiment.')
    return plan


def _native_molecules(engine_id, native_result):
    from rdkit import Chem
    if engine_id in {"gnina", "gnina_affinity"}:
        path = native_result["expected_outputs"][0]
        return list(Chem.SDMolSupplier(path, removeHs=False))
    from meeko import PDBQTMolecule, RDKitMolCreate
    paths = sorted((Path(native_result["output_dir"]) / "poses").glob("*.pdbqt")) if engine_id == "unidock" else [Path(native_result["expected_outputs"][0])]
    molecules = []
    for path in paths:
        pdbqt = PDBQTMolecule(path.read_text(), is_dlg=engine_id == "autodock4", skip_typing=True)
        sdf_text, failures = RDKitMolCreate.write_sd_string(pdbqt)
        if failures:
            raise ValueError("Native output lacks recoverable Meeko chemical identity/atom mapping; supply a mapped preparation.")
        molecules.extend(Chem.ForwardSDMolSupplier(io.BytesIO(sdf_text.encode()), removeHs=False))
    return molecules


def _canonical_poses(plan, native_result, out, emit, cancelled):
    from rdkit import Chem
    canonical = {science._canonical(mol): mol for mol in Chem.SDMolSupplier(plan["validation"]["canonical_ligands"], removeHs=False)}
    reference_path = plan["validation"].get("reference")
    reference = next(iter(Chem.SDMolSupplier(reference_path, removeHs=False))) if reference_path else None
    receptor = Path(plan["validation"]["receptor"])
    molecules = _native_molecules(plan["engine_id"], native_result)
    if len(molecules) != len(native_result["results"]):
        raise ValueError("Native pose count does not match score records.")
    poses, counts = [], {}
    metric_name = {"vina": "vina_score", "autodock4": "autodock4_score", "unidock": "unidock_vina_score", "gnina": "minimizedAffinity", "gnina_affinity": "minimizedAffinity"}[plan["engine_id"]]
    for mol, row in zip(molecules, native_result["results"]):
        science._cancel(cancelled)
        if mol is None or science._canonical(mol) not in canonical:
            raise ValueError("Native pose chemical identity does not match a registered canonical ligand.")
        identity = science._canonical(mol)
        original = canonical[identity]
        compound_id = "cmp_" + hashlib.sha256(identity.encode()).hexdigest()[:16]
        counts[identity] = counts.get(identity, 0) + 1
        number = counts[identity]
        directory = out / "validated_poses" / compound_id / f"pose_{number}"
        directory.mkdir(parents=True)
        pose_file = directory / "pose.sdf"
        with Chem.SDWriter(str(pose_file)) as writer:
            writer.write(mol)
        if emit:
            emit({"kind": "tool_started", "tool": "PoseBusters", "stage": "pose_validation", "message": f"Validate {plan['engine_id']} pose {number} chemical identity and physical geometry.", "compound_id": compound_id})
        validity = science.validate_pose(pose_file, receptor, directory,
                                        min(120, plan["engine_plan"]["parameters"]["timeout_seconds"]), cancelled,
                                        reference=reference_path and Path(reference_path))
        metric = next(item for item in row["metrics"] if item["name"] == metric_name)
        poses.append({"compound_id": compound_id, "state_id": compound_id + "_input_state", "pose_id": f"{compound_id}_pose_{number}",
                      "name": original.GetProp("_Name") if original.HasProp("_Name") else compound_id,
                      "pose_number": number, "score": metric["value"], "score_units": metric["unit"], "score_direction": metric["direction"],
                      "score_definition": metric["semantics"], "metrics": row["metrics"], "engine": plan["engine_id"],
                      "engine_version": native_result["provenance"]["engine_version"], "canonical_smiles": identity,
                      "pose_path": str(pose_file), "pose_sha256": science.sha256(pose_file), "receptor_sha256": science.sha256(receptor),
                      "validation": validity, "redocking": science._symmetry_rmsd(mol, reference) if reference is not None else {"rmsd_angstrom": None, "reason": "No cognate reference for this screening compound."},
                      "interactions": {"status": "not_requested", "tool": "PLIP"}, "descriptors": science._descriptors(original),
                      "status": "exploratory", "admet": None})
    if set(counts) != set(canonical):
        raise ValueError("One or more canonical ligands have no output poses; partial screening is not complete.")
    return poses


def inspect_native_resume(plan: dict) -> dict:
    """Authenticate every retained pair before offering an application resume."""
    if plan.get('engine_id') != 'biocentric_affinity':
        raise ValueError('Native checkpoint resume is limited to BioCentric.')
    engine_plan = plan['engine_plan']
    root = Path(engine_plan['output_dir'])
    request_sha = hashlib.sha256(engine_plan['generated_files']['request.json'].encode()).hexdigest()
    expected_paths = {root / 'pairs' / f"{row['record_index']:06d}.json" for row in engine_plan['expected_pairs']}
    if set((root / 'pairs').glob('*.json')) - expected_paths:
        raise ValueError('Unexpected pair checkpoints cannot be included in a reviewed resume.')
    checkpoints = []
    for expected in engine_plan['expected_pairs']:
        path = root / 'pairs' / f"{expected['record_index']:06d}.json"
        if not path.exists():
            continue
        raw, parsed = engines._biocentric_checkpoint(engine_plan, expected, request_sha)
        checkpoints.append({'record_index': expected['record_index'], 'compound_id': expected['compound_id'],
                            'status': raw['status'], 'sha256': parsed['source']['sha256']})
    return {'request_sha256': request_sha, 'checkpoints': checkpoints,
            'completed': sum(row['status'] == 'completed' for row in checkpoints),
            'failed': sum(row['status'] == 'failed' for row in checkpoints),
            'missing': len(engine_plan['expected_pairs']) - len(checkpoints)}


def _archive_native_resume(store, pid, plan, job_id, reviewed_state, retry_failed_pairs):
    state = inspect_native_resume(plan)
    if reviewed_state is not None and state != reviewed_state:
        raise ValueError('Pair checkpoints changed after the resume review; no recorded failure was retried.')
    out = Path(plan['output_dir'])
    archives = []
    # Save even a hard-killed worker's logs and partial exports before the runner
    # replaces its working files. Registered artifact bytes remain immutable.
    for path in sorted(out.rglob('*')):
        if not path.is_file():
            continue
        if not path.resolve().is_relative_to(out.resolve()) or path.stat().st_size > 100 * 1024 * 1024:
            raise ValueError('Prior native output is outside its run or exceeds the archive limit.')
        item = store.artifact(pid, path.name, path.read_bytes(), 'native_engine_prior_attempt', job_id)
        archives.append({'relative_path': str(path.relative_to(out)), 'artifact_id': item['id'], 'sha256': item['sha256']})
    receipt = {'reviewed_state': state, 'retry_failed_pairs': retry_failed_pairs, 'artifacts': archives}
    artifact = store.artifact(pid, 'native-resume-archive.json', json.dumps(receipt, indent=2).encode(), 'native_engine_resume_archive', job_id)
    progress = Path(plan['engine_plan']['output_dir']) / 'progress.json'
    if progress.exists():
        archived = next((item for item in archives if item['relative_path'] == str(progress.relative_to(out))), None)
        if not archived or archived['sha256'] != science.sha256(progress):
            raise ValueError('Prior progress changed while archiving; refusing resume.')
        progress.unlink()  # The new worker must not appear to report old progress.
    if retry_failed_pairs:
        # Explicit review authorizes only these already archived failures. A
        # valid completed checkpoint never gets deleted or recomputed here.
        for row in state['checkpoints']:
            if row['status'] != 'failed':
                continue
            path = Path(plan['engine_plan']['output_dir']) / 'pairs' / f"{row['record_index']:06d}.json"
            if science.sha256(path) != row['sha256']:
                raise ValueError('Failed checkpoint changed while archiving; refusing retry.')
            archived = next((item for item in archives if item['relative_path'] == str(path.relative_to(out))), None)
            if not archived or archived['sha256'] != row['sha256']:
                raise ValueError('Failed checkpoint was not archived intact; refusing retry.')
            path.unlink()
    return {'archive_artifact_id': artifact['id'], **state, 'retry_failed_pairs': retry_failed_pairs}


def run_native_manifest(store, pid: str, manifest_artifact, output_dir: str | Path, *,
                        expected_engine: str, expected_phase: str, expected_testset_sha256: str | None = None,
                        expected_protocol_sha256: str | None = None, job_id: str | None = None,
                        emit: Callable | None = None, cancelled: Callable | None = None,
                        resume_plan: dict | None = None, reviewed_resume_state: dict | None = None,
                        retry_failed_pairs: bool = False) -> dict:
    """Execute a registered native job and register its raw and canonical results."""
    plan = plan_native_manifest(store, pid, manifest_artifact, output_dir, expected_engine=expected_engine,
                                expected_phase=expected_phase, expected_testset_sha256=expected_testset_sha256,
                                expected_protocol_sha256=expected_protocol_sha256,
                                **({'resume_plan': resume_plan} if resume_plan is not None else {}))
    if retry_failed_pairs and (resume_plan is None or reviewed_resume_state is None):
        raise ValueError('Retrying failed pairs requires a reviewed resume state.')
    if plan["blockers"]:
        return {"status": "blocked", "engine": expected_engine, "blockers": plan["blockers"], "poses": [], "results": [], "artifacts": [], "executed": False}
    out = Path(plan["output_dir"])
    resume = (_archive_native_resume(store, pid, plan, job_id, reviewed_resume_state, retry_failed_pairs)
              if resume_plan is not None else None)
    out.mkdir(parents=True, exist_ok=True)
    science._write_json(out / "workflow_plan.json", plan)
    native = engines.run_engine(expected_engine, plan["engine_inputs"], plan["manifest"]["parameters"], out / "engine", emit=emit, cancelled=cancelled,
                                **({'resume_plan': resume_plan['engine_plan']} if resume_plan is not None else {}))
    result = {"status": native["status"], "engine": expected_engine, "executed": native["executed"], "native_result": native,
              "poses": [], "results": native["results"], "failures": [], "artifacts": [],
              "testset_sha256": plan["manifest"]["testset_sha256"], "protocol_sha256": plan["manifest"]["protocol_sha256"],
              "metadata": {"benchmark_status": "unvalidated_exploratory", "native_manifest_sha256": plan["manifest_artifact"]["sha256"],
                           "frozen_input_binding": plan["frozen_input_binding"]}}
    # This is prior validation of the exact model input structure, not a claim
    # that the affinity prediction itself passed a physical or predictive test.
    source_entries = {row["compound_id"]: row for row in
                      (plan.get("input_structure_provenance") or {}).get("entries", [])}
    result["results"] = [dict(row) for row in native["results"]]
    for row in result["results"]:
        source = source_entries.get(row.get("compound_id"))
        if source is not None:
            row["input_structure_provenance"] = deepcopy(source)
            row["input_structure_physical_passed"] = source["source_physical_validation"]["passed"]
    if resume is not None:
        result['metadata']['checkpoint_resume'] = resume
    if native["status"] == "completed":
        try:
            if plan["validation"].get("canonical_ligands"):
                result["poses"] = _canonical_poses(plan, native, out, emit, cancelled)
            if expected_engine == "planet":
                from rdkit import Chem
                expected = sorted(science._canonical(mol) for mol in Chem.SDMolSupplier(plan["engine_inputs"]["ligands"], removeHs=False))
                observed = []
                for row in result["results"]:
                    mol = Chem.MolFromSmiles(row["smiles"])
                    if mol is None:
                        raise ValueError("PLANET output contains invalid chemical identity.")
                    observed.append(science._canonical(mol))
                if sorted(observed) != expected:
                    raise ValueError("PLANET output does not match the registered query identities.")
            for artifact in plan["resolved_artifacts"].values():
                store.read_artifact(pid, artifact)
        except (ValueError, RuntimeError, science.ScienceError) as error:
            result["status"] = "cancelled" if getattr(error, "code", None) == "cancelled" else "failed"
            result["failures"] = [{"code": getattr(error, "code", "canonical_validation"), "reason": str(error)}]
            result["poses"] = []
    science._write_json(out / "native_workflow.json", result)
    for path in sorted(out.rglob("*")):
        if path.is_file():
            if path.stat().st_size > 100 * 1024 * 1024:
                raise ValueError("Native output exceeds the compact artifact registration limit (100 MiB).")
            artifact = store.artifact(pid, path.name, path.read_bytes(), "native_engine_output", job_id)
            result["artifacts"].append(artifact)
            # Canonical pose links should point at the registered immutable artifact.
            for pose in result["poses"]:
                if pose["pose_path"] == str(path):
                    pose["pose_path"] = artifact["path"]
                    pose["pose_artifact_id"] = artifact["id"]
    final = store.artifact(pid, f"{expected_engine}-{expected_phase}-result.json", json.dumps(result, indent=2).encode(), "native_engine_result", job_id)
    result["result_artifact_id"] = final["id"]
    return result
