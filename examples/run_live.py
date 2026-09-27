"""Actually run the local Vina/PoseBusters vertical slice, with inspectable events.

    .venv/bin/python examples/run_live.py --output runs/live-smoke
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from biocentric.benchmark import smoke_report
from biocentric.science import dock_library, prepare_library, sha256


def make_library(destination: Path) -> None:
    from rdkit import Chem
    bound = next(iter(Chem.SDMolSupplier(str(ROOT / "examples/1iep/originals/1iep_ligand.sdf"), removeHs=False)))
    bound.SetProp("_Name", "imatinib_bound_control")
    bound.SetProp("SOURCE", "AutoDock Vina v1.2.7 tutorial; PDB 1IEP STI; Mus musculus ABL P00520")
    bound.SetProp("CONTROL_ROLE", "cognate bound ligand; redocking software control")
    bound.SetProp("ASSAY_RESULT", "No measured affinity imported")
    with Chem.SDWriter(str(destination)) as writer:
        writer.write(bound)
        for name, smiles in [("caffeine_untested_control", "Cn1c(=O)c2c(ncn2C)n(C)c1=O"),
                             ("aspirin_untested_control", "CC(=O)Oc1ccccc1C(=O)O")]:
            mol = Chem.MolFromSmiles(smiles)
            mol.SetProp("_Name", name)
            mol.SetProp("SOURCE", "Manually curated structure string; identity checked by RDKit; no assay assertion")
            mol.SetProp("INPUT_SMILES", smiles)
            mol.SetProp("CONTROL_ROLE", "untested chemical diversity control; NOT a tested inactive")
            mol.SetProp("ASSAY_RESULT", "unknown")
            writer.write(mol)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "runs/live-smoke")
    parser.add_argument("--exhaustiveness", type=int, default=4)
    arguments = parser.parse_args()
    output = arguments.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    originals = ROOT / "examples/1iep/originals"
    library = ROOT / "examples/1iep/library.sdf"
    if not library.exists():
        make_library(library)
    events = output / "events.jsonl"
    def emit(event):
        print(json.dumps(event), flush=True)
        with events.open("a") as handle:
            handle.write(json.dumps(event) + "\n")
    preparation = prepare_library(library, output / "preparation", emit=emit)
    result = dock_library(originals / "1iep_receptor.pdbqt", Path(preparation["prepared_sdf"]), originals / "1iep_ligand.sdf", output / "docking",
                          {"seed": 42, "cpu": 1, "exhaustiveness": arguments.exhaustiveness, "n_poses": 3,
                           "timeout_seconds": 180, "box_padding_angstrom": 4, "receptor_pdb": str(originals / "1iep_receptorH.pdb"),
                           "require_posebusters": True, "run_plip": True}, emit=emit)
    benchmark = smoke_report(result)
    (output / "benchmark.json").write_text(json.dumps(benchmark, indent=2))
    fields = ["compound_id", "name", "pose_number", "score", "score_units", "engine", "status", "pose_path"]
    with (output / "candidates.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(result["poses"])
    report = ["# BioCentricGPT — live software smoke report", "", "Status: **unvalidated / exploratory**.", "",
              "This report contains executed local calculations. The tutorial uses **mouse ABL (P00520), PDB 1IEP**, not human ABL1.", "",
              "Untested caffeine/aspirin structures are software controls, not experimentally demonstrated inactives. No enrichment estimate is possible.", "",
              "| Compound | Pose | Vina score (kcal/mol) | PoseBusters passed | Redock RMSD (Å) | PLIP |",
              "|---|---:|---:|---|---:|---|"]
    for row in result["poses"]:
        rmsd = row["redocking"].get("rmsd_angstrom")
        report.append(f"| {row['name']} | {row['pose_number']} | {row['score']:.3f} | {row['validation'].get('passed')} | {f'{rmsd:.3f}' if rmsd is not None else 'not a cognate redock'} | {row['interactions']['status']} |")
    report += ["", "## Interpretation", "", benchmark["reason"], "",
               "Scores are engine estimates, not measured binding affinity. A failed validity check remains visible and must be inspected before prioritization.", "",
               "RMSD is symmetry-aware in the receptor coordinate frame. A separately reported ligand-aligned shape RMSD is not used as evidence of correct pocket placement.", "",
               "## Reproducibility", "", f"- Vina: {result['metadata']['vina']}", f"- Request fingerprint: `{result['metadata']['fingerprint']}`",
               "- Seed 42; CPU 1; three poses; bound-ligand extents plus 4 Å padding.",
               f"- Exhaustiveness: {arguments.exhaustiveness} (software smoke budget, not a tuned production protocol).",
               "- Canonical input/prepared SDF, raw PDBQT, each pose SDF, raw logs, PoseBusters JSON and run hashes are retained beside this report.",
               "- Original tutorial receptor protonation/preparation is accepted as a curated input, not newly validated or repaired.", "",
               "## Source", "", "[PDB 1IEP](https://www.rcsb.org/structure/1IEP); [Vina v1.2.7 example](https://github.com/ccsb-scripps/AutoDock-Vina/tree/v1.2.7/example/basic_docking).", ""]
    if result["failures"]:
        report += ["## Explicit failures", "", "```json", json.dumps(result["failures"], indent=2), "```", ""]
    (output / "report.md").write_text("\n".join(report))
    manifest = {"kind": "actual_local_calculation", "target": {"accession": "P00520", "organism": "Mus musculus", "pdb_id": "1IEP"},
                "preparation": preparation, "docking": result, "benchmark": benchmark,
                "report": {"path": str(output / "report.md"), "sha256": sha256(output / "report.md")}}
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2))
    print(f"Report: {output / 'report.md'}")
    return 0 if len({row["compound_id"] for row in result["poses"]}) == len(preparation["compounds"]) and not result["failures"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
