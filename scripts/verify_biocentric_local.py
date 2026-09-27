"""Real, bounded BioCentric native/adapter equivalence and checkpoint reuse.

Run with the app interpreter, e.g. ``python -m scripts.verify_biocentric_local
--model-dir /home/sargis/dti --model-python /home/sargis/anaconda3/bin/python
--out runs/linux-biocentric-reference``. No LLM or network calls; requires CUDA.
"""
import argparse
import csv
import json
import os
from pathlib import Path
import shutil

from biocentric import engines, science
from biocentric.biocentric_model import parse_output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model-dir", type=Path, required=True)
    parser.add_argument("--model-python", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    root, out = args.model_dir.resolve(), args.out.resolve()
    if out.exists() and any(out.iterdir()):
        raise ValueError("Verification output must be new or empty; preserve previous evidence")
    out.mkdir(parents=True, exist_ok=True)
    os.environ["BIOCENTRIC_MODEL_DIR"] = str(root)
    os.environ["BIOCENTRIC_MODEL_PYTHON"] = str(args.model_python.resolve())
    os.environ["PYTHONDONTWRITEBYTECODE"] = "1"
    source = root / "Inference/projects/Fxagtc000101A"
    receptor = source / "proteins/Fxagtc000101A_25A.pdb"
    original = source / "ligands/Fxagtc000101A_GTC000064A_conf0.pdb"
    native = out / "native"
    poses = native / "poses"
    poses.mkdir(parents=True)
    shutil.copyfile(original, poses / "GTC000064A.pdb")
    original_command = [str(args.model_python.resolve()), str(root / "Prediction/run_pipeline.py"),
                        "--pocket", str(receptor), "--poses", str(poses), "--out", str(native / "scores.csv"),
                        "--config", str(root / "Prediction/config.yaml"), "--data-config", str(root / "Prediction/data_config.yaml"),
                        "--device", "cuda", "--per-pose"]
    native_execution = science.run_process(original_command, root / "Prediction", native / "engine.log", 300)
    from rdkit import Chem
    mol = Chem.MolFromPDBFile(str(original), sanitize=True, removeHs=False)
    if mol is None:
        raise ValueError("Reference ligand cannot be read")
    mol.SetProp("_Name", "GTC000064A")
    ligands = out / "ligand.sdf"
    writer = Chem.SDWriter(str(ligands))
    writer.write(mol)
    writer.close()
    result = engines.run_engine("biocentric_affinity", {
        "receptor": str(receptor), "ligands": str(ligands), "protein_id": "Fxagtc000101A",
        "preparation_reviewed": True}, {"device": "cuda", "cpu": 2}, out / "adapter")
    if result["status"] != "completed" or result["results"][0]["status"] != "completed":
        raise RuntimeError("Adapter failed: " + json.dumps(result.get("error", result.get("results"))))
    with (native / "scores.csv").open() as stream:
        reference = next(csv.DictReader(stream))
    predicted = result["results"][0]
    differences = {name: abs(float(reference[name]) - score) for name, score in predicted["per_model"].items()}
    differences["ensemble"] = abs(float(reference["score"]) - predicted["score"])
    # Radius/scatter GPU kernels and float32 reductions can vary in last bits.
    tolerance = 1e-5
    if max(differences.values()) > tolerance:
        raise RuntimeError("Native/adapter numerical equivalence failed: " + json.dumps(differences))
    before = json.loads((out / "adapter/predictions.json").read_text())
    reuse_execution = science.run_process(result["commands"][0], out / "adapter", out / "resume.log", 300)
    reused = parse_output(out / "adapter/predictions.json", result["expected_pairs"])
    if not all(row.get("reused") for row in reused) or reused[0]["score"] != before["predictions"][0]["score"]:
        raise RuntimeError("Prediction checkpoint reuse failed")
    report = {"status": "passed", "native_execution": native_execution, "adapter_status": result["status"],
              "native_adapter_absolute_differences": differences, "tolerance": tolerance,
              "native_score": float(reference["score"]), "adapter_score": predicted["score"],
              "checkpoint_resume": {"reused_pairs": len(reused), "new_model_predictions": 0, "execution": reuse_execution},
              "resources": before["resources"], "runtime": before["runtime"],
              "scientific_validation": "Software/reference equivalence only. No measured affinity accuracy established."}
    science._write_json(out / "verification.json", report)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
