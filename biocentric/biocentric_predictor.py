"""Fixed subprocess bridge; compatible with the upstream Python 3.9 ML runtime.

Run by absolute path, never as ``python -m biocentric``: the application package
and the independent predictor intentionally have separate interpreter imports.
"""
import argparse
import csv
import gc
import hashlib
import json
import math
import os
from pathlib import Path
import resource
import sys
import time

sys.dont_write_bytecode = True


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for part in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(part)
    return digest.hexdigest()


def write_json(path, payload):
    temporary = Path(str(path) + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, allow_nan=False))
    temporary.replace(path)


def row_sha256(row):
    return hashlib.sha256(json.dumps(row, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()).hexdigest()


def read_checkpoint(path, request_sha256, expected):
    """Reject changed checkpoints without replacing the original evidence."""
    previous = json.loads(Path(path).read_text())
    row = previous.get("row", {})
    if previous.get("request_sha256") != request_sha256 or previous.get("row_sha256") != row_sha256(row):
        raise ValueError("Existing pair checkpoint failed request/content hash verification; preserve it and use a new output directory")
    if row.get("status") not in {"completed", "failed"} or any(row.get(key) != expected[key] for key in ("record_index", "compound_id", "protein_id", "smiles")):
        raise ValueError("Existing pair checkpoint identity or status differs from its reviewed request")
    return row


def recover_cuda_after_failure(torch_module):
    # Called after the except suite, when Python has released the exception's
    # traceback and the failed model-forward frames. Collect cycles before asking
    # the CUDA allocator to release its now-unused blocks.
    gc.collect()
    torch_module.cuda.empty_cache()


def resource_usage(torch_module, device, started):
    return {"elapsed_seconds": time.monotonic() - started,
            "peak_rss_kib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            "cuda_allocated_bytes": torch_module.cuda.memory_allocated() if device == "cuda" else 0,
            "peak_cuda_allocated_bytes": torch_module.cuda.max_memory_allocated() if device == "cuda" else 0}


def event(kind, **kwargs):
    print(json.dumps({"event": kind, "engine_id": "biocentric_affinity", **kwargs}, allow_nan=False), flush=True)


def canonical(mol):
    from rdkit import Chem
    return Chem.MolToSmiles(Chem.RemoveHs(mol), canonical=True, isomericSmiles=True)


def validate_pdb_identity(mol, expected):
    """PDB lacks E/Z labels: verify them from coordinates, never from the answer.

    This validates a copy only. Native featurization and pose coordinates stay
    unchanged; altered charge, bonding or actual stereochemistry still fail.
    """
    from rdkit import Chem
    if mol is not None:
        if canonical(mol) == expected:
            return "pdb_metadata"
        if mol.GetNumConformers() and mol.GetConformer().Is3D():
            geometric = Chem.Mol(mol)
            Chem.AssignStereochemistryFrom3D(geometric, replaceExistingTags=True)
            if canonical(geometric) == expected:
                return "verified_3d_stereochemistry"
    raise ValueError("Reference SDF-to-PDB conversion changes stereochemistry, charge or bonding; pair excluded")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--request", required=True)
    args = parser.parse_args()
    request_path = Path(args.request).resolve()
    request = json.loads(request_path.read_text())
    request_digest = sha256(request_path)
    if request.get("schema_version") != 1:
        raise ValueError("Unsupported request schema")
    for item in request["validated_inputs"].values():
        if sha256(item["path"]) != item["sha256"]:
            raise ValueError("Reviewed input changed: " + item["path"])
    out = Path(request["output_dir"]).resolve()
    if Path(args.request).resolve().parent != out:
        raise ValueError("Request and output directories differ")
    out.mkdir(parents=True, exist_ok=True)
    stage = out / "prepared"
    stage.mkdir(exist_ok=True)
    checkpoints = out / "pairs"
    checkpoints.mkdir(exist_ok=True)
    started = time.monotonic()
    os.environ["OMP_NUM_THREADS"] = str(request["parameters"]["cpu"])
    os.environ["OPENBLAS_NUM_THREADS"] = str(request["parameters"]["cpu"])
    root = Path(request["model_root"]).resolve()
    os.chdir(root)
    sys.path.insert(0, str(root))
    import torch
    import torch_geometric
    from rdkit import Chem, rdBase
    import numpy as np
    torch.set_num_threads(request["parameters"]["cpu"])
    device = request["parameters"]["device"]
    if device == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("Requested CUDA is unavailable; no implicit CPU fallback")
    if device == "cuda":
        torch.cuda.reset_peak_memory_stats()
    # The upstream constructors read training/device even during inference. This
    # explicit runtime-only override preserves architecture/config hyperparameters.
    import utils
    original_config = utils.Config
    class DeviceConfig(original_config):
        def __init__(self, *a, **kw):
            super().__init__(*a, **kw)
            self._data.setdefault("experiment", {}).setdefault("training", {})["device"] = device
    utils.Config = DeviceConfig
    original_load = torch.load
    def load_on_device(*a, **kw):
        kw.setdefault("map_location", device)
        kw.setdefault("weights_only", True)
        return original_load(*a, **kw)
    torch.load = load_on_device
    import run_pipeline as reference
    from Config import Config
    cfg = Config(str(root / "config.yaml"))
    data_cfg = Config(str(root / "data_config.yaml"))
    event("model_loading", device=device, models=request["models"])
    models = reference.load_ensemble(cfg, device)
    if models["names"] != request["models"]:
        raise ValueError("Enabled native model ensemble differs from reviewed request")
    receptor = request["inputs"]["receptor"]
    protein = Chem.MolFromPDBFile(receptor, sanitize=True, removeHs=False)
    if protein is None:
        raise ValueError("Receptor fails strict sanitization; upstream atom-deletion repair is forbidden")
    if any(atom.GetPDBResidueInfo().GetIsHeteroAtom() for atom in protein.GetAtoms()):
        raise ValueError("Receptor hetero atoms require explicit preparation")
    receptor_positions = Chem.RemoveHs(protein).GetConformer().GetPositions()
    folder, internal_target, level = reference.split_pocket_path(receptor, data_cfg.get("data/protein_level"))
    target = reference.process_target(data_cfg, internal_target, folder, level, data_cfg.get("data/ply_format"), data_cfg.get("data/protein_graph_types"))
    fixed_crop = None
    if request['inputs'].get('crop_reference'):
        ref_path = request['inputs']['crop_reference']
        reference_mol = next(iter(Chem.SDMolSupplier(ref_path, removeHs=False)))
        ref_sdf = stage / 'crystal_crop.sdf'
        with Chem.SDWriter(str(ref_sdf)) as ref_writer:
            ref_writer.write(reference_mol)
        ref_pdb = reference.to_pdb(str(ref_sdf))
        roundtrip = Chem.MolFromPDBFile(ref_pdb, sanitize=True, removeHs=False)
        validate_pdb_identity(roundtrip, canonical(reference_mol))
        fixed_crop = reference.process_ligand(data_cfg, internal_target, 'crystal_crop', ref_pdb,
                                             0, True, data_cfg.get('data/ligand_graph_types'))
    molecules = list(Chem.SDMolSupplier(request["inputs"]["ligands"], removeHs=False))
    if len(molecules) != len(request["pairs"]):
        raise ValueError("Input record count changed")
    predictions = []
    writer = Chem.SDWriter(str(out / "scored.sdf"))
    try:
        for expected, mol in zip(request["pairs"], molecules):
            index = expected["record_index"]
            row = {**expected, "status": "failed", "preprocessing_status": "pending", "score": None,
                   "failure_reason": None, "per_model": {}, "models": request["models"],
                   "checkpoint_sha256": request["weights_sha256"], "code_sha256": request["code_sha256"],
                   "device": device, "crop_reference_sha256": request["validated_inputs"].get("crop_reference", {}).get("sha256")}
            event("pair_started", compound_id=expected["compound_id"], record_index=index, total=len(molecules))
            pair_started = time.monotonic()
            checkpoint = checkpoints / ("%06d.json" % index)
            previous = read_checkpoint(checkpoint, request_digest, expected) if checkpoint.exists() else None
            recover_gpu = False
            try:
                if mol is None or canonical(mol) != expected["smiles"] or mol.GetProp("_Name").strip() != expected["compound_id"]:
                    raise ValueError("Exact compound identity differs from reviewed input")
                heavy_positions = Chem.RemoveHs(mol).GetConformer().GetPositions()
                if not np.isfinite(heavy_positions).all() or not np.isfinite(receptor_positions).all():
                    raise ValueError("Nonfinite coordinates")
                if not any(np.min(np.sum((receptor_positions - point) ** 2, axis=1)) <= 7**2 for point in heavy_positions):
                    raise ValueError("Ligand does not overlap the receptor pocket within the tightest 7 angstrom crop")
                if previous is not None:
                    row = dict(previous)
                    row["reused"] = True
                else:
                    sdf = stage / ("pair_%06d.sdf" % index)
                    one = Chem.SDWriter(str(sdf))
                    one.write(mol)
                    one.close()
                    pdb = reference.to_pdb(str(sdf))
                    pdb_mol = Chem.MolFromPDBFile(pdb, sanitize=True, removeHs=False)
                    row["pdb_identity_validation"] = validate_pdb_identity(pdb_mol, expected["smiles"])
                    ligand = reference.process_ligand(data_cfg, internal_target, "pair_%06d" % index,
                                                      pdb, 0, True, data_cfg.get("data/ligand_graph_types"))
                    # Benchmark mode holds the crystal-defined protein context fixed
                    # across all docked poses; standalone calls retain native one-pose behavior.
                    merged = reference.score_pose(data_cfg, models, target, internal_target, ligand,
                                                  expected["compound_id"], fixed_crop if fixed_crop is not None else ligand, device)
                    row["per_model"] = {name: float(merged.loc[0, "model_%d" % i]) for i, name in enumerate(models["names"])}
                    row["score"] = float(merged.loc[0, "score"])
                    if not all(math.isfinite(value) for value in [row["score"], *row["per_model"].values()]):
                        raise ValueError("Nonfinite native model output")
                    row.update(status="completed", preprocessing_status="passed", failure_reason=None,
                               ensemble_spread=float(merged.loc[0, "ScoreSTD"]), spread_is_calibrated_uncertainty=False,
                               reused=False)
            except Exception as error:
                row.update(status="failed", score=None, per_model={}, preprocessing_status="failed",
                           failure_reason=type(error).__name__ + ": " + str(error))
                recover_gpu = device == "cuda"
            if recover_gpu:
                recover_cuda_after_failure(torch)
            row["duration_seconds"] = time.monotonic() - pair_started
            if previous is None:
                write_json(checkpoint, {"request_sha256": request_digest, "row_sha256": row_sha256(row), "row": row})
            predictions.append(row)
            if mol is not None:
                exported = Chem.Mol(mol)
                for key, value in {"BioCentric_protein_id": row["protein_id"], "BioCentric_compound_id": row["compound_id"],
                                   "BioCentric_status": row["status"], "BioCentric_checkpoint_sha256": row["checkpoint_sha256"],
                                   "BioCentric_score_unit": "pK", "BioCentric_score_direction": "higher_is_better",
                                   "BioCentric_failure_reason": row["failure_reason"] or ""}.items():
                    exported.SetProp(key, str(value))
                if row["score"] is not None:
                    exported.SetDoubleProp("BioCentric_affinity", row["score"])
                writer.write(exported)
            event("pair_completed", compound_id=row["compound_id"], record_index=index, total=len(molecules),
                  status=row["status"], score=row["score"], failure_reason=row["failure_reason"])
            write_json(out / "progress.json", {"completed_pairs": index, "total_pairs": len(molecules), "latest": previous if previous is not None else row,
                                               "resource_usage": resource_usage(torch, device, started)})
    finally:
        writer.close()
    for item in request["validated_inputs"].values():
        if sha256(item["path"]) != item["sha256"]:
            raise ValueError("Reviewed input changed during execution: " + item["path"])
    payload = {"schema_version": 1, "engine_id": "biocentric_affinity", "expected_pairs": request["pairs"],
               "predictions": predictions, "models": models["names"], "weights_sha256": request["weights_sha256"],
               "code_sha256": request["code_sha256"], "device": device,
               "runtime": {"python": sys.version, "torch": torch.__version__, "torch_geometric": torch_geometric.__version__,
                           "rdkit": rdBase.rdkitVersion, "cuda": torch.version.cuda,
                           "gpu_name": torch.cuda.get_device_name() if device == "cuda" else None},
               "resources": {"duration_seconds": time.monotonic() - started,
                             "peak_rss_kib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                             "peak_cuda_allocated_bytes": torch.cuda.max_memory_allocated() if device == "cuda" else 0,
                             "cpu_threads": torch.get_num_threads(), "batch_size": 1},
               "scientific_validation": "Software inference only; no measured affinity benchmark was evaluated."}
    write_json(out / "predictions.json", payload)
    fields = ["record_index", "protein_id", "compound_id", "smiles", "status", "score", "failure_reason",
              "preprocessing_status", "checkpoint_sha256", "code_sha256", "device", "duration_seconds", "per_model"]
    with (out / "predictions.csv").open("w", newline="") as stream:
        csv_writer = csv.DictWriter(stream, fields, extrasaction="ignore")
        csv_writer.writeheader()
        csv_writer.writerows({**row, "per_model": json.dumps(row["per_model"], sort_keys=True)} for row in predictions)
    event("run_completed", successful_pairs=sum(row["status"] == "completed" for row in predictions),
          failed_pairs=sum(row["status"] == "failed" for row in predictions), **payload["resources"])


if __name__ == "__main__":
    main()
