"""Auditable local chemistry and docking; no network access and no generated commands.

The canonical chemistry record is SDF. PDBQT is an engine input/output only.
Optional tools are lazy imports and unavailable tools never become successful results.
"""
from __future__ import annotations

import hashlib
from functools import lru_cache
import importlib.metadata
import importlib.util
import io
import json
import math
import os
import platform
import re
import shutil
import signal
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable

from pydantic import BaseModel, ConfigDict, Field


class ScienceError(RuntimeError):
    def __init__(self, code: str, message: str, details: dict | None = None):
        super().__init__(message)
        self.code, self.details = code, details or {}


class DockConfig(BaseModel):
    model_config = ConfigDict(extra="forbid")
    seed: int = Field(default=42, ge=1, le=2147483647)
    cpu: int = Field(default=1, ge=1, le=16)
    exhaustiveness: int = Field(default=4, ge=1, le=128)
    n_poses: int = Field(default=3, ge=1, le=20)
    timeout_seconds: float = Field(default=120, gt=0, le=3600)
    box_padding_angstrom: float = Field(default=4.0, ge=2, le=12)
    receptor_pdb: str | None = None
    require_posebusters: bool = True
    run_plip: bool = True


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


@lru_cache(maxsize=1024)
def _runtime_hash(path, stamp):
    return sha256(Path(path))


def runtime_sha256(path: Path) -> str:
    """Hash unchanged installed runtimes once per worker; detect replacement/writes."""
    path = Path(path).resolve()
    def stamp():
        value = path.stat()
        return (value.st_dev, value.st_ino, value.st_size, value.st_mtime_ns, value.st_ctime_ns)
    before = stamp()
    # Avoid timestamp-granularity collisions during rapid same-size writes.
    # Installed binaries/weights are stable; recently touched files get a full read.
    result = (sha256(path) if time.time_ns() - max(before[-2:]) < 1_000_000_000
              else _runtime_hash(str(path), before))
    if stamp() != before:
        raise ValueError('Runtime changed during hashing: ' + str(path))
    return result


def _write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False, default=str))
    temporary.replace(path)


def _artifact(path: Path, role: str) -> dict:
    return {"path": str(path.resolve()), "role": role, "sha256": sha256(path), "bytes": path.stat().st_size}


def _version(name: str) -> str | None:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return None


def _require(name: str) -> None:
    if importlib.util.find_spec(name) is None:
        raise ScienceError("missing_dependency", f"{name} is not installed; no calculation was performed.")


def _emit(emit, kind: str, tool: str, stage: str, message: str, **data) -> None:
    if emit:
        emit({"kind": kind, "tool": tool, "stage": stage, "message": message, **data})


def _progress_poses(poses: list[dict]) -> list[dict]:
    """Publish measured per-compound results before the whole library finishes.

    These lightweight records are provisional until the stage registers its
    immutable artifacts. They never expose an unregistered file for download.
    """
    fields = ("compound_id", "name", "pose_id", "pose_number", "score", "score_units",
              "score_direction", "engine", "engine_version", "pose_sha256", "receptor_sha256", "protocol_fingerprint")
    return [{**{key: pose.get(key) for key in fields},
             "validation": {"passed": pose.get("validation", {}).get("passed")}}
            for pose in poses]


def _cancel(cancelled) -> None:
    if cancelled and cancelled():
        raise ScienceError("cancelled", "Scientific calculation cancelled by the user.")


def _input(path: Path, suffixes: set[str]) -> Path:
    path = Path(path).resolve(strict=True)
    if not path.is_file() or path.suffix.lower() not in suffixes:
        raise ScienceError("input_format", f"Expected one of {sorted(suffixes)}: {path.name}")
    if path.stat().st_size > 100 * 1024 * 1024:
        raise ScienceError("input_size", "Files larger than 100 MiB exceed this compact workflow's limits.")
    return path


def _original(path: Path, out_dir: Path) -> Path:
    """Content-addressed snapshot; an existing snapshot is never overwritten."""
    digest = sha256(path)
    destination = out_dir / "originals" / f"{digest[:16]}_{path.name}"
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists():
        if sha256(destination) != digest:
            raise ScienceError("input_integrity", "An immutable input snapshot has changed.")
    else:
        with destination.open("xb") as dst, path.open("rb") as src:
            shutil.copyfileobj(src, dst)
        destination.chmod(0o444)
    return destination


def _canonical(mol) -> str:
    from rdkit import Chem
    return Chem.MolToSmiles(Chem.RemoveHs(mol), canonical=True, isomericSmiles=True)


def _supported(mol) -> None:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Lipinski
    if len(Chem.GetMolFrags(mol)) != 1:
        raise ScienceError("unsupported_salt", "Disconnected components require an explicit parent/salt decision; nothing was stripped.")
    if any(a.GetAtomicNum() not in {1, 6, 7, 8, 9, 15, 16, 17, 35, 53} for a in mol.GetAtoms()):
        raise ScienceError("unsupported_element", "Metal-containing or unsupported-element chemistry needs a separate protocol.")
    if any(a.GetNumRadicalElectrons() for a in mol.GetAtoms()):
        raise ScienceError("unsupported_radical", "Radicals are outside the supported closed-shell workflow.")
    stereo = Chem.FindPotentialStereo(mol)
    if any(s.specified == Chem.StereoSpecified.Unspecified for s in stereo):
        raise ScienceError("ambiguous_stereochemistry", "Unspecified stereochemistry requires a measured assignment or explicit state enumeration.")
    if Descriptors.MolWt(mol) > 750 or Lipinski.NumRotatableBonds(mol) > 15:
        raise ScienceError("unsupported_flexibility", "Molecular weight >750 or >15 rotatable bonds exceeds the MVP domain.")
    if any(len(ring) >= 9 for ring in mol.GetRingInfo().AtomRings()):
        raise ScienceError("unsupported_macrocycle", "Macrocycles need an independently validated protocol.")
    amide = Chem.MolFromSmarts("C(=O)N")
    if len(mol.GetSubstructMatches(amide)) >= 4:
        raise ScienceError("unsupported_peptide", "Peptide-like chemistry with four or more amides needs a separate protocol.")


def _descriptors(mol) -> dict:
    from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors
    from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
    parameters = FilterCatalogParams()
    parameters.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    alerts = [entry.GetDescription() for entry in FilterCatalog(parameters).GetMatches(mol)]
    return {"molecular_weight": Descriptors.MolWt(mol), "cLogP": Crippen.MolLogP(mol),
            "tpsa_angstrom2": rdMolDescriptors.CalcTPSA(mol), "h_bond_donors": Lipinski.NumHDonors(mol),
            "h_bond_acceptors": Lipinski.NumHAcceptors(mol), "rotatable_bonds": Lipinski.NumRotatableBonds(mol),
            "formal_charge": sum(a.GetFormalCharge() for a in mol.GetAtoms()), "pains_review_flags": alerts,
            "admet": None, "admet_status": "not_predicted", "alerts_are_rejection_rules": False}


def prepare_library(input_path: Path, out_dir: Path, seed: int = 42,
                    emit: Callable | None = None, cancelled: Callable | None = None) -> dict:
    """Sanitize SDF/SMI, preserve explicit states/metadata, generate missing 3D geometry.

    Exact isomeric duplicates retain their source record and link to the first state.
    This is deliberately not a pH predictor or tautomer/stereoisomer enumerator.
    """
    _require("rdkit")
    from rdkit import Chem
    from rdkit.Chem import AllChem
    if not isinstance(seed, int) or not 1 <= seed <= 2147483647:
        raise ScienceError("invalid_seed", "Seed must be an integer in [1, 2147483647].")
    source = _input(input_path, {".sdf", ".sd", ".smi", ".smiles", ".txt"})
    out_dir = Path(out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    snapshot = _original(source, out_dir)
    _emit(emit, "tool_started", "RDKit", "ligand_preparation", "Sanitize, check chemical identity and create explicit 3D states.")
    if source.suffix.lower() in {".sdf", ".sd"}:
        molecules = [(i + 1, mol, {}) for i, mol in enumerate(Chem.SDMolSupplier(str(snapshot), removeHs=False))]
    else:
        molecules = []
        for line_no, line in enumerate(snapshot.read_text().splitlines(), 1):
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            fields = line.strip().split(maxsplit=1)
            mol = Chem.MolFromSmiles(fields[0])
            if mol is not None:
                mol.SetProp("_Name", fields[1] if len(fields) > 1 else f"compound_{line_no}")
                mol.SetProp("input_smiles", fields[0])
            molecules.append((line_no, mol, {"input_line": line}))
    if not molecules or len(molecules) > 1000:
        raise ScienceError("library_size", "The library must contain 1–1000 records.")
    compounds, exclusions, seen = [], [], {}
    sdf_path = out_dir / "prepared.sdf"
    with Chem.SDWriter(str(sdf_path)) as writer:
        for record_index, mol, raw in molecules:
            _cancel(cancelled)
            try:
                if mol is None:
                    raise ScienceError("invalid_chemistry", "RDKit could not parse/sanitize the record; original retained.")
                Chem.SanitizeMol(mol)
                Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
                _supported(mol)
                canonical = _canonical(mol)
                compound_id = "cmp_" + hashlib.sha256(canonical.encode()).hexdigest()[:16]
                name = mol.GetProp("_Name") if mol.HasProp("_Name") else f"compound_{record_index}"
                properties = {key: mol.GetProp(key) for key in mol.GetPropNames(includePrivate=True, includeComputed=False)}
                if canonical in seen:
                    exclusions.append({"record_index": record_index, "name": name, "code": "exact_duplicate", "reason": "Exact isomeric state already retained.", "duplicate_of": seen[canonical], "metadata": properties, **raw})
                    continue
                initial_3d = bool(mol.GetNumConformers() and mol.GetConformer().Is3D())
                hydrogenated = Chem.AddHs(mol, addCoords=initial_3d)
                minimization = "not_run_existing_3d_preserved"
                if not initial_3d:
                    hydrogenated.RemoveAllConformers()
                    params = AllChem.ETKDGv3()
                    params.randomSeed = seed
                    if AllChem.EmbedMolecule(hydrogenated, params) != 0:
                        raise ScienceError("embedding_failed", "ETKDGv3 could not generate a conformer.")
                    if AllChem.MMFFHasAllMoleculeParams(hydrogenated):
                        code = AllChem.MMFFOptimizeMolecule(hydrogenated, maxIters=500)
                        minimization = "MMFF94_converged" if code == 0 else "MMFF94_not_converged"
                    else:
                        minimization = "not_run_missing_MMFF_parameters"
                if _canonical(hydrogenated) != canonical:
                    raise ScienceError("identity_changed", "Adding hydrogen/coordinates altered the canonical isomeric identity.")
                coords = hydrogenated.GetConformer().GetPositions()
                if not all(math.isfinite(float(v)) for row in coords for v in row):
                    raise ScienceError("invalid_coordinates", "Non-finite molecular coordinates.")
                state_id = compound_id + "_input_state"
                for key, value in {"BIOCENTRIC_COMPOUND_ID": compound_id, "BIOCENTRIC_STATE_ID": state_id,
                                   "BIOCENTRIC_CANONICAL_SMILES": canonical, "BIOCENTRIC_SOURCE_RECORD": str(record_index)}.items():
                    hydrogenated.SetProp(key, value)
                writer.write(hydrogenated)
                compounds.append({"compound_id": compound_id, "state_id": state_id, "name": name,
                                  "canonical_smiles": canonical, "source_record": record_index, "source_sha256": sha256(snapshot),
                                  "metadata": properties, "descriptors": _descriptors(mol),
                                  "preparation": {"seed": seed, "input_3d_preserved": initial_3d, "geometry": "input" if initial_3d else "ETKDGv3",
                                                  "minimization": minimization, "protonation": "input formal charge retained; hydrogens added by valence",
                                                  "ph": None, "tautomer_enumeration": "not_performed", "stereoisomer_enumeration": "not_performed"}})
                seen[canonical] = compound_id
            except ScienceError as error:
                exclusions.append({"record_index": record_index, "code": error.code, "reason": str(error), **raw})
            except (ValueError, RuntimeError) as error:
                exclusions.append({"record_index": record_index, "code": "preparation_failed", "reason": str(error), **raw})
    result = {"status": "completed" if compounds else "failed", "compounds": compounds, "exclusions": exclusions,
              "prepared_sdf": str(sdf_path), "input_records": len(molecules), "tool": "RDKit", "version": _version("rdkit"),
              "artifacts": [_artifact(snapshot, "immutable_original"), _artifact(sdf_path, "canonical_prepared_library")], "created_at": _now()}
    manifest = out_dir / "preparation.json"
    _write_json(manifest, result)
    result["artifacts"].append(_artifact(manifest, "preparation_manifest"))
    _emit(emit, "tool_completed", "RDKit", "ligand_preparation", f"Prepared {len(compounds)} compounds; {len(exclusions)} explicit exclusions.", data={"accepted": len(compounds), "excluded": len(exclusions)})
    return result


def inspect_structure(path: Path) -> dict:
    """Inspect deposited records without deleting waters, cofactors, partners or loops."""
    path = _input(path, {".pdb", ".pdbqt", ".cif", ".mmcif"})
    result = {"path": str(path), "sha256": sha256(path), "format": path.suffix[1:], "changes": [],
              "missing_pocket_residues": "unknown", "biological_assembly_verified": False,
              "sequence_construct_mapping": "requires_target_mapping", "pocket_confidence": "not_assessed"}
    if path.suffix.lower() in {".cif", ".mmcif"}:
        _require("gemmi")
        import gemmi
        structure = gemmi.read_structure(str(path))
        result.update({"models": len(structure), "chains": sorted({chain.name for model in structure for chain in model}),
                       "residues": sum(len(chain) for model in structure for chain in model), "resolution_angstrom": structure.resolution or None,
                       "heterogens": sorted({res.name for model in structure for chain in model for res in chain if res.het_flag == "H"}),
                       "warning": "mmCIF metadata preserved in original; receptor preparation requires reviewed standard PDB/PDBQT export."})
        return result
    lines = path.read_text().splitlines()
    atoms = [line for line in lines if line.startswith(("ATOM  ", "HETATM"))]
    if not atoms:
        raise ScienceError("empty_structure", "No atom records in structure.")
    residues = {(line[21:22].strip(), line[22:27].strip(), line[17:20].strip()) for line in atoms}
    heterogens = sorted({line[17:20].strip() for line in atoms if line.startswith("HETATM")})
    metals = {"ZN", "FE", "MG", "MN", "CA", "CU", "CO", "NI", "CD", "HG"}
    atom_elements = {line.split()[-1].upper() if path.suffix.lower() == ".pdbqt" else line[76:78].strip().upper() for line in atoms}
    metal_residues = sorted((set(heterogens) | atom_elements) & metals)
    result.update({"atom_count": len(atoms), "chains": sorted({line[21:22].strip() for line in atoms}), "residues": len(residues),
                   "residue_map": [{"chain": c, "residue_number": n, "name": r} for c, n, r in sorted(residues)],
                   "heterogens": heterogens, "metals": metal_residues,
                   "waters": sum(1 for _, _, name in residues if name in {"HOH", "WAT"}),
                   "alternate_locations": sorted({line[16] for line in atoms if line[16].strip()}),
                   "missing_residue_remarks": [line for line in lines if line.startswith("REMARK 465")],
                   "covalent_links": [line for line in lines if line.startswith("LINK  ")],
                   "routine_screening_blockers": ["Metal/cofactor protocol requires review."] if metal_residues else [],
                   "status": "curated_input_required"})
    return result


def derive_box(reference_sdf: Path, padding: float = 4.0) -> dict:
    _require("rdkit")
    from rdkit import Chem
    source = _input(reference_sdf, {".sdf", ".sd"})
    molecules = list(Chem.SDMolSupplier(str(source), removeHs=True))
    if len(molecules) != 1 or molecules[0] is None or not molecules[0].GetNumConformers() or not molecules[0].GetConformer().Is3D():
        raise ScienceError("reference_ligand", "Box definition requires exactly one valid bound 3D ligand in the receptor coordinate frame.")
    coords = molecules[0].GetConformer().GetPositions()
    center = [(float(coords[:, i].min()) + float(coords[:, i].max())) / 2 for i in range(3)]
    size = [float(coords[:, i].max() - coords[:, i].min()) + 2 * padding for i in range(3)]
    if not all(math.isfinite(v) for v in center + size) or max(size) > 50:
        raise ScienceError("unsupported_box", "Invalid coordinates or derived box >50 Å; review the site.")
    return {"center": center, "size": size, "units": "angstrom", "method": "bound_ligand_heavy_atom_extents_plus_padding",
            "padding_angstrom": padding, "reference_sha256": sha256(source), "coordinate_frame": "input receptor, no alignment/transformation"}


def _vina_executable() -> str | None:
    # Only administrator environment/PATH or this fixed installation path; never LLM/user command text.
    candidate = os.environ.get("BIOCENTRIC_VINA") or shutil.which("vina")
    tools = Path(__file__).resolve().parents[1] / ".tools"
    local = tools / "linux" / "vina-1.2.7" if platform.system() == "Linux" else tools / "vina"
    if candidate and Path(candidate).is_file() and os.access(candidate, os.X_OK):
        return str(Path(candidate).resolve())
    return str(local) if local.is_file() and os.access(local, os.X_OK) else None


def _tool_environment() -> dict[str, str]:
    """Scientific executables and health probes never need application credentials."""
    return {key: value for key, value in os.environ.items()
            if not any(part in key.upper() for part in ("API_KEY", "TOKEN", "SECRET", "PASSWORD", "CREDENTIAL"))}


def run_process(argv: list[str], cwd: Path, log_path: Path, timeout: float,
                cancelled: Callable | None = None) -> dict:
    """Internal fixed-argv runner: no shell, bounded runtime, terminate entire job group."""
    _cancel(cancelled)
    started = time.monotonic()
    environment = _tool_environment()
    environment.update({"OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1"})
    with Path(log_path).open("w") as log:
        try:
            process = subprocess.Popen(argv, cwd=cwd, stdout=log, stderr=subprocess.STDOUT, shell=False,
                                       start_new_session=True, env=environment)
        except FileNotFoundError as error:
            raise ScienceError("missing_dependency", f"Executable unavailable: {Path(argv[0]).name}") from error
        try:
            while process.poll() is None:
                _cancel(cancelled)
                if time.monotonic() - started > timeout:
                    raise ScienceError("timeout", f"Tool exceeded {timeout:g}s timeout.")
                time.sleep(0.15)
        except BaseException:
            try:
                os.killpg(process.pid, signal.SIGTERM)
                process.wait(timeout=3)
            except (ProcessLookupError, subprocess.TimeoutExpired):
                if process.poll() is None:
                    os.killpg(process.pid, signal.SIGKILL)
                    process.wait()
            raise
    if process.returncode != 0:
        raise ScienceError("tool_failed", f"{Path(argv[0]).name} exited {process.returncode}; inspect the tool log.", {"log": str(log_path), "returncode": process.returncode})
    return {"argv": argv, "duration_seconds": time.monotonic() - started, "returncode": process.returncode,
            "log": str(log_path), "log_sha256": sha256(log_path)}


def _symmetry_rmsd(probe, reference) -> dict:
    """Same-receptor-frame RMSD is the redocking criterion, never ligand-only alignment."""
    from rdkit import Chem
    from rdkit.Chem import rdMolAlign
    probe, reference = Chem.RemoveHs(probe), Chem.RemoveHs(reference)
    if _canonical(probe) != _canonical(reference):
        return {"rmsd_angstrom": None, "reason": "Different ligand identity; not a redocking control."}
    try:
        # CalcRMS considers symmetry without superposing away pocket-placement errors.
        rmsd = rdMolAlign.CalcRMS(probe, reference, maxMatches=100000, symmetrizeConjugatedTerminalGroups=True)
        aligned = rdMolAlign.GetBestRMS(Chem.Mol(probe), reference, maxMatches=100000)
        return {"rmsd_angstrom": float(rmsd), "aligned_ligand_shape_rmsd_angstrom": float(aligned),
                "method": "symmetry-aware, receptor coordinates already aligned (same input frame)",
                "warning": "Ligand-only aligned shape RMSD is diagnostic and is not the redocking gate."}
    except RuntimeError as error:
        return {"rmsd_angstrom": None, "reason": str(error)}


def validate_pose(pose_path: Path, receptor_pdb: Path | None, out_dir: Path,
                  timeout: float = 120, cancelled=None, reference: Path | None = None) -> dict:
    """Run PoseBusters in a killable subprocess; missing receptor checks stay explicit."""
    if importlib.util.find_spec("posebusters") is None:
        return {"status": "unavailable", "tool": "PoseBusters", "reason": "Dependency not installed", "passed": None}
    payload = {"pose": str(pose_path), "receptor": str(receptor_pdb) if receptor_pdb else None,
               "reference": str(reference) if reference else None, "output": str(out_dir / "posebusters.json")}
    request_path = out_dir / "posebusters_request.json"
    _write_json(request_path, payload)
    run_process([sys.executable, str(Path(__file__).resolve()), "_posebusters", str(request_path)], out_dir,
                out_dir / "posebusters.log", timeout, cancelled)
    result = json.loads(Path(payload["output"]).read_text())
    result["artifact"] = _artifact(Path(payload["output"]), "physical_validity_report")
    return result


def _posebusters_worker(payload: dict) -> None:
    from posebusters import PoseBusters
    config = "dock" if payload["receptor"] else "mol"
    frame = PoseBusters(config=config).bust([payload["pose"]], payload.get("reference"), payload["receptor"])
    rows = json.loads(frame.reset_index().to_json(orient="records"))
    checks = {key: bool(value) for row in rows for key, value in row.items() if isinstance(value, bool)}
    if not checks:
        raise ScienceError("empty_validation", "PoseBusters returned no validation checks.")
    _write_json(Path(payload["output"]), {"status": "completed", "tool": "PoseBusters", "version": _version("posebusters"),
                "mode": config, "checks": checks, "passed": all(checks.values()), "failed_checks": [k for k, v in checks.items() if not v],
                "receptor_checks_run": payload["receptor"] is not None, "rows": rows})


def interaction_report(receptor_pdb: Path, pose_path: Path, out_dir: Path, timeout: float = 120, cancelled=None) -> dict:
    """PLIP local complex report, retained raw XML; unavailable OpenBabel is honest."""
    executable = shutil.which("plip")
    local = Path(sys.executable).parent / "plip"
    executable = executable or (str(local) if local.exists() else None)
    if not executable:
        return {"status": "unavailable", "tool": "PLIP", "reason": "PLIP and OpenBabel are not installed; no interaction report calculated."}
    from rdkit import Chem
    ligand = next(iter(Chem.SDMolSupplier(str(pose_path), removeHs=False)))
    if ligand is None:
        raise ScienceError("plip_ligand", "PLIP input SDF cannot be read.")
    expected_identity = _canonical(ligand)
    protein_lines = [line for line in receptor_pdb.read_text().splitlines() if not line.startswith(("END", "CONECT"))]
    protein_atoms = [line for line in protein_lines if line.startswith(("ATOM  ", "HETATM"))]
    used_chains = {line[21:22] for line in protein_atoms}
    ligand_chain = next((chain for chain in "ZYXWVUTSRQPONMLKJIHGFEDCBA" if chain not in used_chains), None)
    if not ligand_chain:
        raise ScienceError("plip_chain", "No free single-character PDB chain for the ligand.")
    serial_offset = max(int(line[6:11]) for line in protein_atoms)
    for atom in ligand.GetAtoms():
        info = Chem.AtomPDBResidueInfo()
        info.SetName(f"{atom.GetSymbol():>2}{atom.GetIdx()+1:<2}"[:4])
        info.SetResidueName("LIG")
        info.SetResidueNumber(1)
        info.SetChainId(ligand_chain)
        info.SetIsHeteroAtom(True)
        atom.SetMonomerInfo(info)
    ligand_lines, ligand_bonds = [], []
    for line in Chem.MolToPDBBlock(ligand).splitlines():
        if line.startswith("HETATM"):
            ligand_lines.append(line[:6] + f"{int(line[6:11]) + serial_offset:5d}" + line[11:])
        elif line.startswith("CONECT"):
            indices = [int(line[start:start+5]) + serial_offset for start in range(6, len(line), 5) if line[start:start+5].strip()]
            ligand_bonds.append("CONECT" + "".join(f"{number:5d}" for number in indices))
    complex_path = out_dir / "complex.pdb"
    complex_path.write_text("\n".join(protein_lines + ["TER"] + ligand_lines + ligand_bonds + ["END"]) + "\n")
    has_receptor_hydrogens = any(line[76:78].strip() == "H" for line in protein_atoms)
    argv = [executable, "-f", str(complex_path), "-o", str(out_dir / "plip"), "-x", "-t", "--maxthreads", "1", "--nofix"]
    if has_receptor_hydrogens:
        argv.append("--nohydro")
    run_process(argv, out_dir,
                out_dir / "plip.log", timeout, cancelled)
    import xml.etree.ElementTree as ET
    xml_path = out_dir / "plip" / "report.xml"
    if not xml_path.is_file():
        raise ScienceError("plip_missing_report", "PLIP produced no machine-readable report.")
    sites = ET.parse(xml_path).getroot().findall("bindingsite")
    ligand_site = next((site for site in sites if site.findtext("identifiers/hetid") == "LIG" and site.findtext("identifiers/chain") == ligand_chain), None)
    if ligand_site is None:
        raise ScienceError("plip_missing_ligand", "PLIP did not identify the intended docked ligand.")
    perceived_smiles = ligand_site.findtext("identifiers/smiles")
    perceived = Chem.MolFromSmiles(perceived_smiles) if perceived_smiles else None
    identity_verified = perceived is not None and _canonical(perceived) == expected_identity
    groups = ligand_site.find("interactions")
    counts = {group.tag: len(group) for group in groups} if groups is not None else {}
    artifacts = [_artifact(path, "plip_report") for path in sorted((out_dir / "plip").glob("**/*")) if path.is_file()]
    return {"status": "completed" if identity_verified else "completed_requires_review", "tool": "PLIP", "version": _version("plip"), "artifacts": artifacts,
            "interaction_counts": counts, "chemical_identity_verified": identity_verified, "perceived_smiles": perceived_smiles,
            "warnings": [] if identity_verified else ["PLIP/OpenBabel perceived ligand identity differs from canonical SDF. Inspect bond orders/protonation before interpreting contact types."],
            "hydrogen_handling": "preserved existing explicit hydrogens" if has_receptor_hydrogens else "PLIP added polar hydrogens; inspect protonated output",
            "assumptions": "PDB bond perception can differ from SDF. Contacts are geometric predictions, not affinity; canonical SDF remains authoritative."}


def dock_library(receptor_path: Path, prepared_sdf: Path, reference_sdf: Path, out_dir: Path,
                 config: dict | None = None, emit: Callable | None = None, cancelled: Callable | None = None) -> dict:
    _require("rdkit")
    _require("meeko")
    from rdkit import Chem
    from meeko import MoleculePreparation, PDBQTWriterLegacy, PDBQTMolecule, RDKitMolCreate
    cfg = DockConfig.model_validate(config or {})
    receptor = _input(receptor_path, {".pdbqt"})
    library = _input(prepared_sdf, {".sdf", ".sd"})
    reference_path = _input(reference_sdf, {".sdf", ".sd"})
    receptor_pdb = _input(Path(cfg.receptor_pdb), {".pdb"}) if cfg.receptor_pdb else None
    out_dir = Path(out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    if inspect_structure(receptor)["routine_screening_blockers"]:
        raise ScienceError("unsupported_receptor", "Metal-containing receptor requires a special validated protocol.")
    for path in [receptor, library, reference_path] + ([receptor_pdb] if receptor_pdb else []):
        _original(path, out_dir)
    box = derive_box(reference_path, cfg.box_padding_angstrom)
    receptor_atoms = [line for line in receptor.read_text().splitlines() if line.startswith(("ATOM  ", "HETATM"))]
    try:
        receptor_coords = [tuple(float(line[start:start+8]) for start in (30, 38, 46)) for line in receptor_atoms]
    except ValueError as error:
        raise ScienceError("receptor_coordinates", "Invalid PDBQT coordinate columns.") from error
    if not any(sum((point[axis]-box["center"][axis])**2 for axis in range(3)) < 12**2 for point in receptor_coords):
        raise ScienceError("coordinate_frame", "Bound-ligand box is not near receptor atoms; confirm the reference coordinate frame.")
    if receptor_pdb:
        validation_coords = {tuple(round(float(line[start:start+8]), 2) for start in (30, 38, 46))
                             for line in receptor_pdb.read_text().splitlines() if line.startswith(("ATOM  ", "HETATM"))}
        receptor_heavy = [point for line, point in zip(receptor_atoms, receptor_coords) if line.split()[-1] not in {"H", "HD", "HS"}]
        overlap = sum(tuple(round(v, 2) for v in point) in validation_coords for point in receptor_heavy) / max(1, len(receptor_heavy))
        if overlap < 0.95:
            raise ScienceError("receptor_validation_frame", "Validation PDB must match the docked receptor's heavy-atom coordinates (>=95% match at 0.01 Å).")
    reference = next(iter(Chem.SDMolSupplier(str(reference_path), removeHs=False)))
    executable = _vina_executable()
    if executable is None:
        _require("vina")
    if cfg.require_posebusters:
        _require("posebusters")
    version = subprocess.run([executable, "--version"], capture_output=True, text=True, timeout=10, env=_tool_environment()).stdout.strip() if executable else _version("vina")
    molecules = list(Chem.SDMolSupplier(str(library), removeHs=False))
    if not molecules or len(molecules) > 1000:
        raise ScienceError("library_size", "Docking library must contain 1–1000 molecules.")
    request_identity = {"adapter_sha256": sha256(Path(__file__)), "receptor": sha256(receptor), "library": sha256(library), "reference": sha256(reference_path),
                        "receptor_pdb": sha256(receptor_pdb) if receptor_pdb else None, "config": cfg.model_dump(), "vina": version,
                        "versions": {name: _version(name) for name in ["rdkit", "meeko", "posebusters", "plip"]}}
    fingerprint = hashlib.sha256(json.dumps(request_identity, sort_keys=True).encode()).hexdigest()
    request_manifest = out_dir / "request.json"
    if request_manifest.exists() and json.loads(request_manifest.read_text())["fingerprint"] != fingerprint:
        raise ScienceError("resume_mismatch", "Inputs/configuration changed; start a new run directory.")
    _write_json(request_manifest, {"fingerprint": fingerprint, **request_identity, "box": box})
    poses, failures, artifacts = [], [], []
    _emit(emit, "tool_completed", "RDKit", "pocket_definition", "Derived docking box from bound-ligand coordinates and recorded padding.", data=box)
    for index, mol in enumerate(molecules):
        _cancel(cancelled)
        if mol is None:
            failures.append({"record_index": index + 1, "code": "invalid_prepared_record", "reason": "Prepared SDF record could not be read."})
            continue
        identity = _canonical(mol)
        compound_id = mol.GetProp("BIOCENTRIC_COMPOUND_ID") if mol.HasProp("BIOCENTRIC_COMPOUND_ID") else "cmp_" + hashlib.sha256(identity.encode()).hexdigest()[:16]
        # IDs are chemical hashes only; imported metadata cannot create paths.
        safe_id = "cmp_" + hashlib.sha256(identity.encode()).hexdigest()[:16]
        if compound_id != safe_id:
            raise ScienceError("identity_metadata", "Prepared compound ID does not match chemical identity.")
        folder = out_dir / f"{index+1:04d}_{safe_id}"
        folder.mkdir(exist_ok=True)
        checkpoint = folder / "completed.json"
        if checkpoint.exists():
            cached = json.loads(checkpoint.read_text())
            if cached.get("fingerprint") == fingerprint and all(Path(a["path"]).is_file() and sha256(Path(a["path"])) == a["sha256"] for a in cached["artifacts"]):
                poses.extend(cached["poses"])
                artifacts.extend(cached["artifacts"])
                _emit(emit, "progress", "AutoDock Vina", "screening", "Reused validated compound checkpoint.", compound_id=compound_id,
                      progress=(index + 1) / len(molecules), data={"poses": _progress_poses(cached["poses"]), "checkpoint_reused": True})
                continue
        try:
            _supported(mol)
            _emit(emit, "tool_started", "Meeko", "ligand_engine_preparation", "Assign Gasteiger charges and Vina atom types; canonical SDF retained.", compound_id=compound_id)
            setups = MoleculePreparation(charge_model="gasteiger").prepare(mol)
            if len(setups) != 1:
                raise ScienceError("ambiguous_preparation", "More than one Meeko setup; explicit state selection required.")
            pdbqt_string, success, error = PDBQTWriterLegacy.write_string(setups[0])
            if not success:
                raise ScienceError("meeko_failed", error)
            ligand_path = folder / "ligand.pdbqt"
            ligand_path.write_text(pdbqt_string)
            _emit(emit, "tool_completed", "Meeko", "ligand_engine_preparation", "PDBQT generated with SMILES and atom mapping for lossless SDF recovery.", compound_id=compound_id)
            output = folder / "vina_poses.pdbqt"
            args = ["--receptor", str(receptor), "--ligand", str(ligand_path), "--out", str(output), "--cpu", str(cfg.cpu),
                    "--seed", str(cfg.seed), "--exhaustiveness", str(cfg.exhaustiveness), "--num_modes", str(cfg.n_poses)]
            for dimension, center, size in zip("xyz", box["center"], box["size"]):
                args.extend([f"--center_{dimension}", str(center), f"--size_{dimension}", str(size)])
            _emit(emit, "tool_started", "AutoDock Vina", "screening", "Docking in a cancellable local subprocess.", compound_id=compound_id,
                  data={"version": version, "seed": cfg.seed, "exhaustiveness": cfg.exhaustiveness, "cpu": cfg.cpu, "box": box})
            if executable:
                execution = run_process([executable] + args, folder, folder / "vina.log", cfg.timeout_seconds, cancelled)
            else:
                payload = {"receptor": str(receptor), "ligand": str(ligand_path), "output": str(output), "box": box, "config": cfg.model_dump()}
                request = folder / "vina_request.json"
                _write_json(request, payload)
                execution = run_process([sys.executable, str(Path(__file__).resolve()), "_vina", str(request)], folder, folder / "vina.log", cfg.timeout_seconds, cancelled)
            scores = [float(match) for match in re.findall(r"REMARK VINA RESULT:\s+(-?\d+(?:\.\d+)?)", output.read_text())]
            sdf_text, conversion_failures = RDKitMolCreate.write_sd_string(PDBQTMolecule(output.read_text(), skip_typing=True))
            if conversion_failures:
                raise ScienceError("pose_conversion", f"Meeko could not recover chemistry for pose groups {conversion_failures}.")
            converted = list(Chem.ForwardSDMolSupplier(io.BytesIO(sdf_text.encode()), removeHs=False))
            if not scores or len(converted) != len(scores):
                raise ScienceError("pose_count_mismatch", "Vina energies and reconstructed SDF poses do not match.")
            compound_poses = []
            for pose_number, (pose, score) in enumerate(zip(converted, scores), 1):
                _cancel(cancelled)
                if pose is None or _canonical(pose) != identity:
                    raise ScienceError("pose_identity_changed", "Recovered pose differs from canonical input chemistry.")
                for prop in mol.GetPropNames(includePrivate=True, includeComputed=False):
                    pose.SetProp(prop, mol.GetProp(prop))
                pose.SetDoubleProp("BIOCENTRIC_VINA_SCORE_KCAL_MOL", score)
                pose.SetIntProp("BIOCENTRIC_POSE_NUMBER", pose_number)
                pose_dir = folder / f"pose_{pose_number}"
                pose_dir.mkdir(exist_ok=True)
                pose_file = pose_dir / "pose.sdf"
                with Chem.SDWriter(str(pose_file)) as writer:
                    writer.write(pose)
                _emit(emit, "tool_started", "PoseBusters", "pose_validation", f"Checking pose {pose_number} for chemical/physical validity.", compound_id=compound_id)
                validity = validate_pose(pose_file, receptor_pdb, pose_dir, cfg.timeout_seconds, cancelled)
                _emit(emit, "tool_completed", "PoseBusters", "pose_validation", f"Pose {pose_number}: {validity['status']}; passed={validity.get('passed')}.", compound_id=compound_id, data=validity)
                rmsd = _symmetry_rmsd(pose, reference)
                interactions = {"status": "not_requested", "tool": "PLIP"}
                if cfg.run_plip and receptor_pdb:
                    _emit(emit, "tool_started", "PLIP", "interaction_inspection", f"Inspecting pose {pose_number} contacts.", compound_id=compound_id)
                    try:
                        interactions = interaction_report(receptor_pdb, pose_file, pose_dir, cfg.timeout_seconds, cancelled)
                    except ScienceError as error:
                        if error.code == "cancelled":
                            raise
                        interactions = {"status": "failed", "reason": str(error), "code": error.code, "tool": "PLIP"}
                    _emit(emit, "tool_completed", "PLIP", "interaction_inspection", f"Interaction adapter: {interactions['status']}.", compound_id=compound_id, data=interactions)
                compound_poses.append({"compound_id": compound_id, "state_id": compound_id + "_input_state", "pose_id": f"{compound_id}_pose_{pose_number}",
                                       "name": mol.GetProp("_Name") if mol.HasProp("_Name") else compound_id, "pose_number": pose_number,
                                       "score": score, "score_units": "kcal/mol", "score_direction": "lower_is_better", "engine": "AutoDock Vina",
                                       "engine_version": version, "score_definition": "Vina predicted binding energy; not measured affinity",
                                       "receptor_sha256": sha256(receptor), "protocol_fingerprint": fingerprint,
                                       "pose_path": str(pose_file), "pose_sha256": sha256(pose_file),
                                       "canonical_smiles": identity, "validation": validity, "redocking": rmsd, "interactions": interactions,
                                       "status": "exploratory", "descriptors": _descriptors(mol), "admet": None})
            compound_artifacts = [_artifact(path, "scientific_output") for path in sorted(folder.rglob("*")) if path.is_file() and path != checkpoint]
            _write_json(checkpoint, {"fingerprint": fingerprint, "poses": compound_poses, "artifacts": compound_artifacts, "execution": execution})
            poses.extend(compound_poses)
            artifacts.extend(compound_artifacts)
            _emit(emit, "tool_completed", "AutoDock Vina", "screening", f"Completed {len(compound_poses)} poses with raw scores and validity checks.", compound_id=compound_id,
                  progress=(index + 1) / len(molecules), data={"scores": scores, "duration_seconds": execution["duration_seconds"], "poses": _progress_poses(compound_poses)})
        except ScienceError as error:
            if error.code == "cancelled":
                raise
            failures.append({"compound_id": compound_id, "code": error.code, "reason": str(error), "details": error.details})
            _emit(emit, "compound_failed", "scientific_pipeline", "screening", str(error), compound_id=compound_id, data={"code": error.code})
        except (ValueError, RuntimeError) as error:
            failures.append({"compound_id": compound_id, "code": "scientific_error", "reason": str(error)})
            _emit(emit, "compound_failed", "scientific_pipeline", "screening", str(error), compound_id=compound_id)
    artifacts.extend([_artifact(request_manifest, "immutable_run_request")] + [_artifact(path, "immutable_original") for path in sorted((out_dir / "originals").iterdir())])
    # Failed/partial calculations are inspectable too; never lose their raw logs.
    by_path = {artifact["path"]: artifact for artifact in artifacts}
    for path in sorted(out_dir.rglob("*")):
        if path.is_file() and path != out_dir / "docking.json":
            by_path.setdefault(str(path.resolve()), _artifact(path, "scientific_output_or_diagnostic"))
    artifacts = list(by_path.values())
    result = {"status": "completed" if poses else "failed", "poses": poses, "failures": failures, "artifacts": artifacts,
              "metadata": {"box": box, "config": cfg.model_dump(), "fingerprint": fingerprint, "versions": request_identity["versions"],
                           "vina": version, "vina_executable_sha256": sha256(Path(executable)) if executable else None,
                           "hardware": {"platform": platform.platform(), "machine": platform.machine(), "cpu_count": os.cpu_count()},
                           "started_from": "curated receptor PDBQT; input protonation retained", "timestamp": _now(),
                           "benchmark_status": "unvalidated_exploratory", "limitation": "Software smoke calculation; no independent active/inactive evaluation."}}
    output_manifest = out_dir / "docking.json"
    _write_json(output_manifest, result)
    result["artifacts"].append(_artifact(output_manifest, "docking_manifest"))
    return result


def _vina_worker(payload: dict) -> None:
    from vina import Vina
    cfg, box = payload["config"], payload["box"]
    engine = Vina(sf_name="vina", cpu=cfg["cpu"], seed=cfg["seed"])
    engine.set_receptor(payload["receptor"])
    engine.set_ligand_from_file(payload["ligand"])
    engine.compute_vina_maps(center=box["center"], box_size=box["size"])
    engine.dock(exhaustiveness=cfg["exhaustiveness"], n_poses=cfg["n_poses"])
    engine.write_poses(payload["output"], n_poses=cfg["n_poses"], overwrite=True)


def prepare_receptor(input_path: Path, out_dir: Path, method: str, *, ph: float = 7.4,
                     add_missing_heavy_atoms: bool = False, reviewed: bool = False,
                     timeout: float = 120, cancelled=None) -> dict:
    """Explicit inspectable PDBFixer/PDB2PQR/Meeko paths. No automatic loop filling.

    The service must record scientific review before passing reviewed=True.
    Raw input is snapshotted; success means a tool ran, not a validated receptor.
    """
    if method not in {"pdbfixer", "pdb2pqr", "meeko"}:
        raise ScienceError("unsupported_adapter", "Allowed preparation methods: pdbfixer, pdb2pqr, meeko.")
    if not reviewed:
        raise ScienceError("review_required", "Inspect retained chains/waters/cofactors, missing atoms and protonation assumptions before preparation.")
    if not math.isfinite(ph) or not 0 <= ph <= 14:
        raise ScienceError("invalid_ph", "pH must be finite and between 0 and 14.")
    source = _input(input_path, {".pdb"})
    out_dir = Path(out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    original = _original(source, out_dir)
    report = inspect_structure(source)
    if report["metals"] or report["alternate_locations"]:
        raise ScienceError("receptor_review", "Resolve metals/alternate conformations with a dedicated protocol first.")
    # Both PDBFixer and Meeko may otherwise fetch CCD templates for unknown
    # residues. Restrict this local-only preparation route before tool dispatch.
    standard_residues = set("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split())
    unsupported_residues = sorted({row["name"] for row in report["residue_map"]} - standard_residues)
    if unsupported_residues:
        raise ScienceError("unsupported_receptor_residue", "Automatic preparation supports standard amino acids only; retained waters/cofactors/nonstandard residues require a reviewed curated-input path. No template lookup was sent.", {"residues": unsupported_residues})
    payload = {"input": str(original), "output": str(out_dir / "prepared.pdb"), "ph": ph,
               "add_missing_heavy_atoms": add_missing_heavy_atoms, "report": str(out_dir / "repair_report.json")}
    if method == "pdbfixer":
        _require("pdbfixer")
        request = out_dir / "pdbfixer_request.json"
        _write_json(request, payload)
        argv = [sys.executable, str(Path(__file__).resolve()), "_pdbfixer", str(request)]
    elif method == "pdb2pqr":
        _require("pdb2pqr")
        payload["output"] = str(out_dir / "prepared.pqr")
        executable = shutil.which("pdb2pqr30") or str(Path(sys.executable).parent / "pdb2pqr30")
        argv = [executable, "--ff=AMBER", "--keep-chain", "--titration-state-method=propka", f"--with-ph={ph}",
                f"--pdb-output={out_dir / 'protonated.pdb'}", str(original), payload["output"]]
    else:
        _require("meeko")
        executable = shutil.which("mk_prepare_receptor.py") or str(Path(sys.executable).parent / "mk_prepare_receptor.py")
        payload["output"] = str(out_dir / "receptor.pdbqt")
        argv = [executable, "--read_pdb", str(original), "-o", str(out_dir / "receptor"), "-p", "-j"]
    execution = run_process(argv, out_dir, out_dir / f"{method}.log", timeout, cancelled)
    if not Path(payload["output"]).is_file():
        raise ScienceError("missing_output", "Preparation tool exited without its required output.")
    result = {"status": "completed_requires_inspection", "method": method, "ph": ph, "input_inspection": report,
            "execution": execution, "output": payload["output"], "artifacts": [_artifact(p, "receptor_preparation") for p in out_dir.rglob("*") if p.is_file()],
            "versions": {name: _version(name) for name in ["pdbfixer", "openmm", "pdb2pqr", "propka", "meeko"]},
            "adapter_sha256": sha256(Path(__file__)), "created_at": _now(),
            "limitations": "No receptor local validation performed. pH is an assumption; PQR is not docking PDBQT. Inspect repair/protonation logs before using output."}
    manifest = out_dir / "preparation_result.json"
    _write_json(manifest, result)
    result["artifacts"].append(_artifact(manifest, "receptor_preparation_manifest"))
    return result


def _pdbfixer_worker(payload: dict) -> None:
    from pdbfixer import PDBFixer
    from openmm.app import PDBFile
    fixer = PDBFixer(filename=payload["input"])
    fixer.findMissingResidues()
    missing_residues = {str(key): values for key, values in fixer.missingResidues.items()}
    # Missing loops remain unresolved, never silently modeled.
    fixer.missingResidues = {}
    fixer.findNonstandardResidues()
    fixer.findMissingAtoms()
    report = {"missing_residues_not_filled": missing_residues,
              "nonstandard_residues_not_replaced": [str(res) for res, _ in fixer.nonstandardResidues],
              "missing_atoms": {str(res): [atom.name for atom in atoms] for res, atoms in fixer.missingAtoms.items()},
              "waters_chains_cofactors": "not removed", "ph": payload["ph"],
              "heavy_atom_repair_requested": payload["add_missing_heavy_atoms"]}
    _write_json(Path(payload["report"]), report)
    if fixer.nonstandardResidues:
        raise ScienceError("nonstandard_residues", "Nonstandard residues require explicit templates; automatic replacement disabled.")
    if payload["add_missing_heavy_atoms"]:
        fixer.addMissingAtoms(seed=42)
    elif fixer.missingAtoms:
        raise ScienceError("missing_atoms", "Missing heavy atoms detected; inspect repair_report.json and explicitly choose repair.")
    fixer.addMissingHydrogens(payload["ph"])
    with Path(payload["output"]).open("w") as handle:
        PDBFile.writeFile(fixer.topology, fixer.positions, handle, keepIds=True)


def health() -> list[dict]:
    checks = []
    for label, module, package in [("RDKit", "rdkit", "rdkit"), ("Meeko", "meeko", "meeko"), ("PoseBusters", "posebusters", "posebusters"),
                                   ("PLIP", "plip", "plip"), ("PDBFixer", "pdbfixer", "pdbfixer"), ("PDB2PQR", "pdb2pqr", "pdb2pqr"),
                                   ("PROPKA", "propka", "propka")]:
        installed = importlib.util.find_spec(module) is not None
        checks.append({"name": label, "available": installed, "status": "installed_not_probe_tested" if installed else "unavailable",
                       "version": _version(package), "implementation": "implemented", "local_only": True})
    executable = _vina_executable()
    version = None
    if executable:
        try:
            version = subprocess.run([executable, "--version"], capture_output=True, text=True, timeout=5, check=True, env=_tool_environment()).stdout.strip()
        except (OSError, subprocess.SubprocessError):
            executable = None
    python_vina = importlib.util.find_spec("vina") is not None
    checks.append({"name": "AutoDock Vina", "available": bool(executable or python_vina), "status": "executable_verified" if executable else "installed_not_probe_tested" if python_vina else "unavailable",
                   "version": version or _version("vina"), "executable": executable, "implementation": "implemented", "local_only": True})
    return checks


if __name__ == "__main__":
    # This private entry point only dispatches fixed scientific operations.
    if len(sys.argv) != 3 or sys.argv[1] not in {"_vina", "_posebusters", "_pdbfixer"}:
        raise SystemExit("Internal scientific worker: use the application CLI.")
    data = json.loads(Path(sys.argv[2]).read_text())
    {"_vina": _vina_worker, "_posebusters": _posebusters_worker, "_pdbfixer": _pdbfixer_worker}[sys.argv[1]](data)
