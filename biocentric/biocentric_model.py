"""Reviewed contract for the separate, local BioCentric structure affinity model.

No predictor imports occur in the application process. The standalone bridge runs
in a dedicated interpreter and calls the original Prediction/run_pipeline.py API.
"""
from __future__ import annotations

import hashlib
import json
import math
import os
from functools import lru_cache
from pathlib import Path
from typing import Literal

from pydantic import BaseModel, ConfigDict, Field

from . import science

ENGINE_ID = "biocentric_affinity"
METRIC = {"name": "BioCentric_affinity", "unit": "pK", "direction": "higher_is_better",
          "semantics": "Learned PDBbind-scale affinity score, inverse-standardized ensemble mean; mixed training endpoints, not a measured Kd/Ki or binding probability.",
          "experimental_measurement": False}
MODELS = (("static_rad_15_15_15_cp_cut_false_predict_pairwise_pos_TRUE_reproduce_pair_false_big_radius", 127),
          ("static_rad_7_7_cp_cut_false_predict_pairwise_pos_TRUE_reproduce", 73),
          ("true_embeddings_ligand_PRETRAINING_LOSSES", 240))
SUPPORTED_WEIGHTS_SHA256 = "a806ae21f91be1a116054a6470c47a760a6499231488de527bfc8cbebbf42ed8"
DOCUMENTATION = str(Path(__file__).resolve().parents[1] / "docs" / "biocentric-model.md")


class BioCentricInputs(BaseModel):
    model_config = ConfigDict(extra="forbid")
    receptor: str
    ligands: str
    crop_reference: str | None = None
    protein_id: str = Field(min_length=1, max_length=200)
    preparation_reviewed: Literal[True]


class BioCentricParameters(BaseModel):
    model_config = ConfigDict(extra="forbid", allow_inf_nan=False)
    device: Literal["cpu", "cuda"] = "cuda"
    cpu: int = Field(default=2, ge=1, le=16)
    timeout_seconds: float = Field(default=300, gt=0, le=3600)
    max_pairs: int = Field(default=128, ge=1, le=1000)
    input_mode: Literal["unique_compounds", "pose_records"] = "unique_compounds"


def _digest(data):
    return hashlib.sha256(json.dumps(data, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


@lru_cache(maxsize=1024)
def _cached_hash(path, size, mtime_ns, ctime_ns):
    return science.sha256(Path(path))


def installation() -> dict:
    """File detection only; successful device execution is recorded separately."""
    root = Path(os.environ.get("BIOCENTRIC_MODEL_DIR", "/nonexistent/biocentric-model")).expanduser().resolve()
    if (root / "Prediction").is_dir():
        root = root / "Prediction"
    interpreter = os.environ.get("BIOCENTRIC_MODEL_PYTHON", "")
    blockers = []
    files = {}
    required = ["run_pipeline.py", "utils.py", "prediction.py", "bioObjects.py", "tmp_dataset.py", "Config.py", "config.py", "config.yaml", "data_config.yaml"]
    for name, epoch in MODELS:
        required += [f"models/{name}/config.yaml", f"models/{name}/dataset.py"]
        required += [f"models/{name}/Checkpoints/Epoch-{epoch}/{component}_model.pt" for component in ("ligand", "target", "aggregation")]
    for relative in required:
        if not (root / relative).is_file():
            blockers.append(f"BioCentric installation is missing {relative}.")
    if root.is_dir():
        source = sorted({*root.rglob("*.py"), *root.rglob("*.yaml")})
        weights = [root / relative for relative in required if relative.endswith(".pt")]
        for path in [*source, *weights]:
            if path.is_file():
                stat = path.stat()
                files[str(path)] = {"path": str(path), "sha256": _cached_hash(str(path), stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns), "bytes": stat.st_size}
    weights = {str(Path(path).relative_to(root)): item["sha256"] for path, item in files.items() if path.endswith(".pt")}
    weight_hash = _digest(weights)
    if weight_hash != SUPPORTED_WEIGHTS_SHA256:
        blockers.append("BioCentric checkpoint bytes do not match the inspected supported ensemble; review a new model contract before use.")
    if not interpreter or not Path(interpreter).is_file() or not os.access(interpreter, os.X_OK):
        blockers.append("Configure BIOCENTRIC_MODEL_PYTHON with the dedicated Linux model interpreter.")
    return {"root": str(root), "python": str(Path(interpreter).absolute()) if interpreter else None,
            "blockers": blockers, "files": files, "weights_sha256": weight_hash,
            "code_sha256": _digest({str(Path(path).relative_to(root)): item["sha256"] for path, item in files.items() if not path.endswith(".pt")}),
            "models": [f"{name}@{epoch}" for name, epoch in MODELS]}


def _input(path, suffix, role, hashes):
    path = Path(path).expanduser().resolve()
    if path.suffix.lower() != suffix or not path.is_file() or not 0 < path.stat().st_size <= 50 * 1024**2:
        raise ValueError(f"BioCentric {role} requires a nonempty {suffix} file of at most 50 MiB.")
    hashes[role] = {"path": str(path), "sha256": science.sha256(path), "bytes": path.stat().st_size}
    return path


def plan(inputs: dict, parameters: dict, output_dir: str | Path, *, resume_plan: dict | None = None) -> dict:
    from rdkit import Chem
    inp, par = BioCentricInputs.model_validate(inputs), BioCentricParameters.model_validate(parameters)
    out = Path(output_dir).expanduser().resolve()
    if resume_plan is not None and out.exists() and any(out.iterdir()):
        request_path = (out / 'request.json').resolve(strict=True)
        if not request_path.is_relative_to(out) or not request_path.is_file():
            raise ValueError('BioCentric saved request escaped its output directory.')
        with request_path.open('rb') as stream:
            saved_request = stream.read(16 * 1024 * 1024 + 1)
        reviewed_request = resume_plan.get('generated_files', {}).get('request.json')
        if len(saved_request) > 16 * 1024 * 1024 or not isinstance(reviewed_request, str) or saved_request != reviewed_request.encode():
            raise ValueError('BioCentric saved request differs from the reviewed resume plan.')
    elif resume_plan is None and out.exists() and any(out.iterdir()):
        raise ValueError("BioCentric output directory must be new or empty.")
    detected = installation()
    hashes = {}
    receptor = _input(inp.receptor, ".pdb", "receptor", hashes)
    ligands = _input(inp.ligands, ".sdf", "ligands", hashes)
    crop_reference = _input(inp.crop_reference, '.sdf', 'crop_reference', hashes) if inp.crop_reference else None
    if crop_reference:
        refs = list(Chem.SDMolSupplier(str(crop_reference), removeHs=False))
        if len(refs) != 1 or refs[0] is None or not refs[0].GetNumConformers() or not refs[0].GetConformer().Is3D():
            raise ValueError('Crop reference must contain exactly one prepared experimental 3D ligand')
        if not all(math.isfinite(float(v)) for xyz in refs[0].GetConformer().GetPositions() for v in xyz):
            raise ValueError('Nonfinite crop reference coordinates')
    atoms = [line for line in receptor.read_text().splitlines() if line.startswith("ATOM  ")]
    if not atoms:
        raise ValueError("BioCentric receptor needs protein ATOM records.")
    try:
        coordinates = [[float(line[start:start + 8]) for start in (30, 38, 46)] for line in atoms]
    except ValueError as error:
        raise ValueError("BioCentric receptor has malformed coordinates.") from error
    if not all(math.isfinite(value) for point in coordinates for value in point):
        raise ValueError("BioCentric receptor coordinates must be finite.")
    if any(line[16:17].strip() for line in atoms) or sum(line.startswith("MODEL ") for line in receptor.read_text().splitlines()) > 1:
        raise ValueError("Choose one reviewed receptor conformation before inference; alternate locations/models are unsupported.")
    if any(line.startswith("HETATM") for line in receptor.read_text().splitlines()):
        raise ValueError("Review and remove receptor HETATM records explicitly: upstream would silently discard them.")
    if len(atoms) > 20000:
        raise ValueError("BioCentric compact worker accepts at most 20,000 receptor atoms; prepare a reviewed pocket.")
    mols = list(Chem.SDMolSupplier(str(ligands), removeHs=False))
    if not 1 <= len(mols) <= par.max_pairs or any(mol is None for mol in mols):
        raise ValueError("BioCentric SDF must contain 1–max_pairs valid, prepared poses.")
    identities, names = [], set()
    for index, mol in enumerate(mols, 1):
        name = mol.GetProp("_Name").strip() if mol.HasProp("_Name") else ""
        if not name or name in names or len(name) > 200 or any(ord(c) < 32 for c in name):
            raise ValueError("BioCentric compound IDs must be present, unique, and free of control characters.")
        names.add(name)
        if not mol.GetNumConformers() or not mol.GetConformer().Is3D() or mol.GetNumHeavyAtoms() > 256:
            raise ValueError("BioCentric requires prepared 3D poses with at most 256 heavy atoms; no coordinates are invented.")
        coords = mol.GetConformer().GetPositions()
        if not all(math.isfinite(float(v)) for xyz in coords for v in xyz):
            raise ValueError("BioCentric ligand coordinates must be finite.")
        identities.append({"record_index": index, "compound_id": name, "protein_id": inp.protein_id,
                           "smiles": science._canonical(mol)})
    if par.input_mode == "unique_compounds" and len({item["smiles"] for item in identities}) != len(identities):
        raise ValueError("BioCentric pair inputs contain duplicate chemical states.")
    for index, item in enumerate(detected["files"].values()):
        hashes[f"model_file_{index}"] = item
    bridge = Path(__file__).with_name("biocentric_predictor.py").resolve()
    hashes["bridge"] = {"path": str(bridge), "sha256": science.sha256(bridge), "bytes": bridge.stat().st_size}
    if detected["python"]:
        executable = Path(detected["python"])
        if executable.is_file():
            hashes["executable"] = {"path": str(executable), "sha256": science.sha256(executable), "bytes": executable.stat().st_size}
    if any(Path(item["path"]).is_relative_to(out) for item in hashes.values()):
        raise ValueError("BioCentric inputs cannot live inside the output directory.")
    request = {"schema_version": 1, "model_root": detected["root"], "output_dir": str(out),
               "inputs": {**inp.model_dump(), "receptor": str(receptor), "ligands": str(ligands), "crop_reference": str(crop_reference) if crop_reference else None},
               "parameters": par.model_dump(), "pairs": identities, "validated_inputs": hashes,
               "models": detected["models"], "weights_sha256": detected["weights_sha256"], "code_sha256": detected["code_sha256"],
               "protein_structure": {"chains": sorted({line[21:22] for line in atoms}), "atom_records": len(atoms),
                                     "residues": len({(line[21:22], line[22:27]) for line in atoms})}}
    command = [detected["python"] or "<BioCentric model Python>", str(bridge), "--request", str(out / "request.json")]
    reviewed_plan = {"engine_id": ENGINE_ID, "status": "blocked" if detected["blockers"] else "ready", "blockers": detected["blockers"],
            "commands": [command], "output_dir": str(out), "expected_outputs": [str(out / "predictions.json"), str(out / "predictions.csv"), str(out / "scored.sdf")],
            "generated_files": {"request.json": json.dumps(request, indent=2)}, "validated_inputs": hashes,
            "parameters": par.model_dump(), "input_review": request["inputs"], "output_metrics": [METRIC],
            "documentation": DOCUMENTATION, "expected_pairs": identities,
            "provenance": {"adapter_schema": "biocentric-model-v1", "adapter_sha256": science.sha256(Path(__file__)),
                           "models": detected["models"], "weights_sha256": detected["weights_sha256"], "code_sha256": detected["code_sha256"],
                           "source_type": "local_file", "reference_entry_point": str(Path(detected["root"]) / "run_pipeline.py")}}
    if resume_plan is not None and reviewed_plan != resume_plan:
        raise ValueError('BioCentric resume inputs, parameters, runtime or adapter differ from the reviewed plan; use a new attempt.')
    return reviewed_plan


def parse_output(output: str | Path, expected_pairs: list[dict] | None = None) -> list[dict]:
    with Path(output).open('rb') as stream:
        content = stream.read(16 * 1024 * 1024 + 1)
    if len(content) > 16 * 1024 * 1024:
        raise ValueError('BioCentric prediction output exceeds the 16 MiB read limit.')
    return parse_payload(json.loads(content), expected_pairs)


def parse_payload(payload: dict, expected_pairs: list[dict] | None = None, *,
                  expected_provenance: dict | None = None, expected_device: str | None = None) -> list[dict]:
    """Validate native JSON or a hash-verified per-pair checkpoint uniformly."""
    if not isinstance(payload, dict) or payload.get("schema_version") != 1 or payload.get("engine_id") != ENGINE_ID:
        raise ValueError("Invalid BioCentric prediction schema.")
    rows = payload.get("predictions", [])
    expected = expected_pairs if expected_pairs is not None else payload.get("expected_pairs", [])
    if not isinstance(rows, list) or not isinstance(expected, list) or any(not isinstance(row, dict) for row in rows + expected):
        raise ValueError('BioCentric prediction and expected pair rows must be lists of objects.')
    identity_fields = ("record_index", "protein_id", "compound_id", "smiles")
    if any(type(row.get('record_index')) is not int or any(not isinstance(row.get(key), str) for key in identity_fields[1:])
           for row in rows + expected):
        raise ValueError('BioCentric pair identities require integer record indices and text identifiers.')
    observed = [tuple(row.get(key) for key in identity_fields) for row in rows]
    required = [tuple(row.get(key) for key in identity_fields) for row in expected]
    if not required or len(set(required)) != len(required) or observed != required:
        raise ValueError("BioCentric missing, duplicated, reordered or substituted protein/compound rows.")
    parsed = []
    for row in rows:
        if row.get("checkpoint_sha256") != SUPPORTED_WEIGHTS_SHA256:
            raise ValueError("BioCentric prediction checkpoint provenance is missing or incompatible.")
        if expected_provenance is not None and any(row.get(field) != expected_provenance[key] for field, key in
                (("checkpoint_sha256", "weights_sha256"), ("code_sha256", "code_sha256"), ("models", "models"))):
            raise ValueError('BioCentric prediction model/code provenance differs from the reviewed plan.')
        if expected_device is not None and row.get('device') != expected_device:
            raise ValueError('BioCentric prediction device differs from the reviewed plan.')
        if row.get("status") == "failed":
            if not row.get("failure_reason") or row.get("score") is not None:
                raise ValueError("BioCentric failures require a reason and no score.")
            parsed.append({**row, "metrics": [], "name": row["compound_id"]})
            continue
        if not isinstance(row.get('per_model'), dict):
            raise ValueError('BioCentric per-model predictions must be an object.')
        values = [row.get("score"), *row["per_model"].values()]
        if row.get("status") != "completed" or len(values) != 4 or any(isinstance(x, bool) or not isinstance(x, (int, float)) or not math.isfinite(x) for x in values):
            raise ValueError("BioCentric returned invalid/nonfinite or incomplete ensemble scores.")
        if set(row["per_model"]) != {f"{name}@{epoch}" for name, epoch in MODELS}:
            raise ValueError("BioCentric checkpoint identities do not match the supported ensemble.")
        if abs(row["score"] - sum(row["per_model"].values()) / 3) > 1e-6:
            raise ValueError("BioCentric ensemble mean differs from its native checkpoint predictions.")
        spread = row.get("ensemble_spread")
        if spread is not None and (isinstance(spread, bool) or not isinstance(spread, (int, float)) or not math.isfinite(spread) or spread < 0):
            raise ValueError("BioCentric ensemble spread must be finite and nonnegative.")
        parsed.append({**row, "metrics": [{**METRIC, "value": row["score"]}], "name": row["compound_id"]})
    return parsed
