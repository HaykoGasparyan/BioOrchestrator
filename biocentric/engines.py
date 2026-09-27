"""Typed, local docking and affinity execution contracts.

Plans do not execute tools. Executions use fixed argv, bounded resources, immutable
input hashes and the scientific runner's cancellation/credential scrubbing.
Availability is installation readiness, never a claim of scientific validation.
"""
from __future__ import annotations

import csv
import hashlib
import io
import json
import math
import os
import platform
import re
import shutil
import sys
import time
from pathlib import Path
from typing import Callable, Literal

from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

from . import science, biocentric_model

ENGINE_IDS = ("vina", "autodock4", "gnina", "unidock", "mmpbsa", "mmgbsa", "gnina_affinity", "planet", "biocentric_affinity")
_DOCS = {
    "vina": "https://autodock-vina.readthedocs.io/en/stable/docking_basic.html",
    "autodock4": "https://autodock.scripps.edu/wp-content/uploads/sites/56/2022/04/AutoDock4.2.6_UserGuide.pdf",
    "gnina": "https://github.com/gnina/gnina",
    "unidock": "https://github.com/dptech-corp/Uni-Dock/blob/main/unidock/README.md",
    "mmpbsa": "https://ambermd.org/tutorials/advanced/tutorial3/py_script/section1.htm",
    "mmgbsa": "https://ambermd.org/tutorials/advanced/tutorial3/py_script/section1.htm",
    "gnina_affinity": "https://github.com/dkoes/docs/wiki/Gnina",
    "planet": "https://github.com/ComputArtCMCG/PLANET",
    "biocentric_affinity": biocentric_model.DOCUMENTATION,
}


class Contract(BaseModel):
    model_config = ConfigDict(extra="forbid", allow_inf_nan=False)


class Box(Contract):
    center: tuple[float, float, float]
    size: tuple[float, float, float]

    @field_validator("size")
    @classmethod
    def bounded_size(cls, values):
        if not all(0 < x <= 50 for x in values):
            raise ValueError("Each box side must be in (0, 50] angstrom for this compact workflow.")
        return values


class DockInputs(Contract):
    receptor: str
    ligand: str
    preparation_reviewed: Literal[True]


class BatchInputs(Contract):
    receptor: str
    ligands: list[str] = Field(min_length=1, max_length=1000)
    preparation_reviewed: Literal[True]


class AutoDockInputs(Contract):
    receptor: str
    ligand: str
    grid_field: str
    atom_maps: dict[str, str] = Field(min_length=1, max_length=32)
    electrostatic_map: str
    desolvation_map: str
    preparation_reviewed: Literal[True]
    maps_match_receptor_reviewed: Literal[True]


class MDReview(Contract):
    """Human-reviewed MD assumptions; these are attestations, not measured facts."""
    equilibrated_trajectory: Literal[True]
    matching_topologies_and_atom_order: Literal[True]
    force_field: str = Field(min_length=2, max_length=200)
    ligand_parameterization: str = Field(min_length=2, max_length=300)
    sampling_notes: str = Field(min_length=10, max_length=2000)
    reviewed_input_sha256: dict[str, str]


class MMInputs(Contract):
    solvated_topology: str
    complex_topology: str
    receptor_topology: str
    ligand_topology: str
    trajectory: str
    review: MDReview


class PlanetInputs(Contract):
    receptor: str
    pocket_reference: str
    ligands: str
    preparation_reviewed: Literal[True]


class RunParameters(Contract):
    timeout_seconds: float = Field(default=120, gt=0, le=3600)


class DockParameters(RunParameters):
    box: Box
    cpu: int = Field(default=1, ge=1, le=16)
    seed: int = Field(default=42, ge=1, le=2147483647)
    exhaustiveness: int = Field(default=8, ge=1, le=128)
    num_modes: int = Field(default=3, ge=1, le=20)


class VinaParameters(DockParameters):
    max_evals: int = Field(default=0, ge=0, le=10000000)


class GninaParameters(DockParameters):
    max_mc_steps: int = Field(default=0, ge=0, le=100000)
    device: Literal["cpu", "gpu"] = "cpu"


class ScoreParameters(RunParameters):
    cpu: int = Field(default=1, ge=1, le=16)
    device: Literal["cpu", "gpu"] = "cpu"


class UniDockParameters(RunParameters):
    box: Box
    seed: int = Field(default=42, ge=1, le=2147483647)
    num_modes: int = Field(default=3, ge=1, le=20)
    search_mode: Literal["fast", "balance", "detail"] = "balance"


class AutoDockParameters(RunParameters):
    seed: int = Field(default=42, ge=1, le=2147483646)
    ga_runs: int = Field(default=10, ge=1, le=100)
    ga_num_evals: int = Field(default=2500000, ge=1000, le=25000000)
    ga_pop_size: int = Field(default=150, ge=20, le=300)


class MMParameters(RunParameters):
    start_frame: int = Field(default=1, ge=1)
    end_frame: int = Field(default=50, ge=2, le=10000)
    interval: int = Field(default=1, ge=1, le=1000)
    salt_molar: float = Field(default=0.15, ge=0, le=1)
    igb: Literal[2, 5, 7, 8] = 5

    @model_validator(mode="after")
    def multiple_frames(self):
        frames = (self.end_frame - self.start_frame) // self.interval + 1
        if not 2 <= frames <= 1000:
            raise ValueError("Require 2–1000 reviewed trajectory snapshots; a docking pose is not an MD ensemble.")
        return self


def _metric(name, unit, direction, semantics):
    return {"name": name, "unit": unit, "direction": direction, "semantics": semantics,
            "experimental_measurement": False}


_METRICS = {
    "vina": [_metric("vina_score", "kcal/mol", "lower_is_better", "Empirical Vina docking score; not experimental affinity or rigorous free energy.")],
    "autodock4": [_metric("autodock4_score", "kcal/mol", "lower_is_better", "AutoDock4 estimated binding energy; not interchangeable with Vina scores.")],
    "gnina": [_metric("minimizedAffinity", "kcal/mol", "lower_is_better", "GNINA Vina empirical score."),
              _metric("CNNscore", "dimensionless", "higher_is_better", "Learned pose score, 0–1; not binding affinity."),
              _metric("CNNaffinity", "pK", "higher_is_better", "Learned affinity prediction; not an experimental Kd/Ki.")],
    "unidock": [_metric("unidock_vina_score", "kcal/mol", "lower_is_better", "Uni-Dock Vina empirical score; protocol-specific ranking.")],
    "mmpbsa": [_metric("delta_total", "kcal/mol", "lower_is_better", "MM/PBSA endpoint estimate without entropy; MD and dielectric assumptions apply.")],
    "mmgbsa": [_metric("delta_total", "kcal/mol", "lower_is_better", "MM/GBSA endpoint estimate without entropy; MD and dielectric assumptions apply.")],
    "planet": [_metric("PLANET_affinity", "pK", "higher_is_better", "Learned pocket/2D-ligand prediction; neither a docking pose nor experimental affinity.")],
}
_METRICS["gnina_affinity"] = _METRICS["gnina"]
_METRICS["biocentric_affinity"] = [biocentric_model.METRIC]


def _executable(name: str) -> str | None:
    if name == "vina":
        return science._vina_executable()
    variable = {"MMPBSA.py": "BIOCENTRIC_MMPBSA", "autodock4": "BIOCENTRIC_AUTODOCK4",
                "gnina": "BIOCENTRIC_GNINA", "unidock": "BIOCENTRIC_UNIDOCK"}[name]
    candidate = os.environ.get(variable) or shutil.which(name)
    if not candidate:
        local = Path(sys.executable).parent / name
        candidate = str(local) if local.is_file() else None
    if not candidate and platform.system() == "Linux":
        local = Path(__file__).resolve().parents[1] / '.tools' / 'linux' / name
        candidate = str(local) if local.is_file() else None
    return str(Path(candidate).resolve()) if candidate and Path(candidate).is_file() and os.access(candidate, os.X_OK) else None


def _planet_installation() -> tuple[list[str], list[str]]:
    location = os.environ.get("BIOCENTRIC_PLANET_DIR")
    interpreter = os.environ.get("BIOCENTRIC_PLANET_PYTHON")
    blockers = []
    if not location or not interpreter:
        return [], ["PLANET requires a reviewed upstream checkout, weights and dedicated Python environment; configure BIOCENTRIC_PLANET_DIR and BIOCENTRIC_PLANET_PYTHON."]
    root = Path(location).resolve()
    for file in ("PLANET_run.py", "PLANET_model.py", "PLANET.param", "layers.py", "chemutils.py", "nnutils.py"):
        if not (root / file).is_file():
            blockers.append(f"PLANET installation is missing {file}.")
    weights = root / "PLANET.param"
    expected = os.environ.get("BIOCENTRIC_PLANET_WEIGHTS_SHA256", "")
    if not re.fullmatch(r"[a-f0-9]{64}", expected) or (weights.is_file() and science.sha256(weights) != expected):
        blockers.append("PLANET weights need an explicitly reviewed matching BIOCENTRIC_PLANET_WEIGHTS_SHA256.")
    if not Path(interpreter).is_file() or not os.access(interpreter, os.X_OK):
        blockers.append("PLANET dedicated Python interpreter is unavailable.")
    return [str(Path(interpreter).resolve()), str(root / "PLANET_run.py")], blockers


def engine_catalog() -> list[dict]:
    """Read-only local discovery; no network, imports of ML runtimes, or provisioning."""
    names = {"vina": "AutoDock Vina", "autodock4": "AutoDock4", "gnina": "GNINA docking", "unidock": "Uni-Dock",
             "mmpbsa": "MM/PBSA (Amber MMPBSA.py)", "mmgbsa": "MM/GBSA (Amber MMPBSA.py)",
             "gnina_affinity": "GNINA CNN affinity", "planet": "PLANET", "biocentric_affinity": "BioCentric affinity"}
    schemas = {"vina": DockInputs, "gnina": DockInputs, "gnina_affinity": DockInputs, "autodock4": AutoDockInputs,
               "unidock": BatchInputs, "mmpbsa": MMInputs, "mmgbsa": MMInputs, "planet": PlanetInputs,
               "biocentric_affinity": biocentric_model.BioCentricInputs}
    rows = []
    for engine_id in ENGINE_IDS:
        binary = {"gnina_affinity": "gnina", "mmpbsa": "MMPBSA.py", "mmgbsa": "MMPBSA.py"}.get(engine_id, engine_id)
        blockers = []
        executable = None
        if engine_id == "biocentric_affinity":
            detected = biocentric_model.installation()
            executable, blockers = detected["python"], list(detected["blockers"])
        elif engine_id == "planet":
            command, blockers = _planet_installation()
            executable = command[0] if command else None
        else:
            executable = _executable(binary)
            if not executable:
                blockers.append(f"{binary} executable is not installed/configured.")
        if engine_id == "unidock":
            if platform.system() != "Linux":
                blockers.append("Uni-Dock requires a configured Linux NVIDIA GPU worker; Apple GPU/macOS is unsupported.")
            elif not Path("/dev/nvidiactl").exists():
                blockers.append("No local NVIDIA driver device detected; configure and verify a CUDA worker.")
        hardware = {"platform": platform.system(), "architecture": platform.machine(), "cpu_supported": engine_id != "unidock",
                    "gpu_required": engine_id == "unidock", "gpu_supported": engine_id in {"unidock", "gnina", "gnina_affinity", "planet", "biocentric_affinity"},
                    "gpu_ready": platform.system() == "Linux" and Path("/dev/nvidiactl").exists(),
                    "remote_dispatch_implemented": False}
        rows.append({"id": engine_id, "name": names[engine_id], "role": "docking" if engine_id in ENGINE_IDS[:4] else "affinity",
                     "input_contract": schemas[engine_id].model_json_schema(), "output_metrics": _METRICS[engine_id],
                     "available": not blockers, "blockers": blockers, "docs": [_DOCS[engine_id]], "hardware": hardware,
                     "executable": executable, "implementation_status": "implemented typed local adapter",
                     "documentation_provenance": "local_file" if engine_id == "biocentric_affinity" else "official_web",
                     "runtime_scope": ("compatibility_only: bundled Uni-Dock 1.1.0; boxes below 30 A, common atom types, non-macrocycles with identity metadata"
                                       if engine_id == 'unidock' and executable and Path(executable).name == 'unidock-1.1.0' else 'requires_task_specific_validation'),
                     "validation_status": "Runtime detection only; actual execution and predictive validation are separate recorded evidence."})
    return rows


def _file(value: str, suffixes: set[str], role: str, hashes: dict) -> Path:
    path = Path(value).expanduser().resolve(strict=True)
    if not path.is_file() or path.suffix.lower() not in suffixes or path.stat().st_size == 0:
        raise ValueError(f"{role} requires a nonempty {sorted(suffixes)} file.")
    if path.stat().st_size > 1024 * 1024 * 1024:
        raise ValueError(f"{role} exceeds the 1 GiB compact-run limit.")
    hashes[role] = {"path": str(path), "sha256": science.sha256(path), "bytes": path.stat().st_size}
    return path


def _molecule_file(path: Path, *, pdbqt=False, receptor=False, require_3d=True):
    if pdbqt or receptor:
        lines = path.read_text().splitlines()
        atoms = [line for line in lines if line.startswith(("ATOM  ", "HETATM"))]
        if not atoms:
            raise ValueError(f"{path.name}: no atomic coordinates.")
        for line in atoms:
            try:
                xyz = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
            except ValueError as error:
                raise ValueError(f"{path.name}: invalid coordinates.") from error
            if not all(math.isfinite(x) for x in xyz):
                raise ValueError(f"{path.name}: nonfinite coordinates.")
        if pdbqt and not receptor and not any(line.startswith("TORSDOF ") for line in lines):
            raise ValueError("Prepared ligand PDBQT must contain its torsion tree and TORSDOF.")
    else:
        from rdkit import Chem
        molecules = list(Chem.SDMolSupplier(str(path), removeHs=False))
        if not molecules or len(molecules) > 1000 or any(m is None for m in molecules):
            raise ValueError("SDF must contain 1–1000 valid molecules; invalid records are not silently skipped.")
        if require_3d and any(not mol.GetNumConformers() or not mol.GetConformer().Is3D() for mol in molecules):
            raise ValueError("Reviewed three-dimensional SDF input is required.")


def _box_arguments(box: Box) -> list[str]:
    return [item for stem, values in (("center", box.center), ("size", box.size))
            for axis, value in zip("xyz", values) for item in (f"--{stem}_{axis}", str(value))]


def _review_unidock_110(receptor: Path, ligands: list[Path], box: Box):
    """Bound the bundled legacy binary to its locally exercised chemistry scope."""
    if any(side >= 30 for side in box.size):
        raise ValueError('Uni-Dock 1.1.0 compatibility runtime requires every box side below 30 angstrom')
    allowed = {'C','A','N','NA','O','OA','S','SA','P','H','HD','F','Cl','Br','I'}
    for path in [receptor, *ligands]:
        atom_types = {line.split()[-1] for line in path.read_text().splitlines() if line.startswith(('ATOM  ','HETATM'))}
        if atom_types - allowed:
            raise ValueError('Uni-Dock 1.1.0 compatibility runtime rejects unsupported atom types or macrocycle glue atoms: ' + ', '.join(sorted(atom_types - allowed)))
    from rdkit import Chem
    for path in ligands:
        smiles = next((line[len('REMARK SMILES '):] for line in path.read_text().splitlines()
                       if line.startswith('REMARK SMILES ') and not line.startswith('REMARK SMILES IDX ')), None)
        mol = Chem.MolFromSmiles(smiles) if smiles else None
        if mol is None:
            raise ValueError('Uni-Dock 1.1.0 compatibility runtime requires Meeko SMILES identity metadata to exclude unsupported macrocycles')
        if any(len(ring) >= 8 for ring in mol.GetRingInfo().AtomRings()):
            raise ValueError('Uni-Dock 1.1.0 compatibility runtime does not accept macrocycles')


def _prmtop_field(content: str, field: str) -> list[str]:
    match = re.search(r"%FLAG " + re.escape(field) + r"\s*\n%FORMAT\((\d+)([aAiIeEfF])(\d+)(?:\.\d+)?\)\s*\n(.*?)(?=\n%FLAG |\Z)", content, re.S)
    if not match:
        raise ValueError(f"Amber topology is missing a supported {field} array.")
    width = int(match.group(3))
    if not 1 <= width <= 32:
        raise ValueError("Unsupported Amber topology field width.")
    return [line[i:i+width].strip() for line in match.group(4).splitlines()
            for i in range(0, len(line), width) if line[i:i+width].strip()]


def _validate_prmtop(path: Path) -> dict:
    content = path.read_text()
    if not content.startswith("%VERSION"):
        raise ValueError(f"{path.name}: expected an Amber topology.")
    pointers = _prmtop_field(content, "POINTERS")
    count = int(pointers[0])
    names = _prmtop_field(content, "ATOM_NAME")
    charges = [float(x.replace("D", "E")) for x in _prmtop_field(content, "CHARGE")]
    if count <= 0 or len(names) != count or len(charges) != count or not all(math.isfinite(x) for x in charges):
        raise ValueError(f"{path.name}: topology atom count, names or charge array is invalid.")
    return {"count": count, "names": names, "charges": charges}


def plan_engine_run(engine_id: str, inputs: dict, parameters: dict, output_dir: str | Path, *, resume_plan: dict | None = None) -> dict:
    """Validate a request and create a reviewable plan without creating files/running tools."""
    if engine_id == "biocentric_affinity":
        return biocentric_model.plan(inputs, parameters, output_dir, resume_plan=resume_plan)
    if resume_plan is not None:
        raise ValueError('Reviewed native checkpoint resume is only supported for BioCentric.')
    if engine_id not in ENGINE_IDS:
        raise ValueError(f"Unknown engine {engine_id!r}.")
    entry = next(row for row in engine_catalog() if row["id"] == engine_id)
    blockers = list(entry["blockers"])
    out = Path(output_dir).expanduser().resolve()
    if out.exists() and any(out.iterdir()):
        raise ValueError("Engine output directory must be new or empty; previous scientific results are immutable.")
    hashes, generated, outputs = {}, {}, []
    executable = entry["executable"]
    if engine_id in {"vina", "gnina", "gnina_affinity"}:
        inp = DockInputs.model_validate(inputs)
        par = ({"vina": VinaParameters, "gnina": GninaParameters, "gnina_affinity": ScoreParameters}[engine_id]).model_validate(parameters)
        receptor = _file(inp.receptor, {".pdbqt"} if engine_id == "vina" else {".pdb", ".pdbqt"}, "receptor", hashes)
        ligand = _file(inp.ligand, {".pdbqt"} if engine_id == "vina" else {".sdf"}, "ligand", hashes)
        _molecule_file(receptor, receptor=True)
        _molecule_file(ligand, pdbqt=engine_id == "vina")
        output = out / ("poses.pdbqt" if engine_id == "vina" else "poses.sdf")
        argv = [executable or ("vina" if engine_id == "vina" else "gnina"), "--receptor", str(receptor), "--ligand", str(ligand), "--out", str(output), "--cpu", str(par.cpu)]
        if engine_id != "gnina_affinity":
            argv += _box_arguments(par.box) + ["--seed", str(par.seed), "--exhaustiveness", str(par.exhaustiveness), "--num_modes", str(par.num_modes)]
        if engine_id == "vina":
            argv += ["--scoring", "vina"]
            if par.max_evals:
                argv += ["--max_evals", str(par.max_evals)]
        else:
            argv += ["--scoring", "vina", "--cnn_scoring", "rescore"]
            if engine_id == "gnina_affinity":
                argv += ["--score_only"]
            elif par.max_mc_steps:
                argv += ["--max_mc_steps", str(par.max_mc_steps)]
            if par.device == "cpu":
                argv += ["--no_gpu"]
            elif not entry["hardware"]["gpu_ready"]:
                blockers.append("GPU execution needs a verified local NVIDIA CUDA device; choose CPU or configure a worker.")
        outputs = [str(output)]
    elif engine_id == "unidock":
        inp, par = BatchInputs.model_validate(inputs), UniDockParameters.model_validate(parameters)
        receptor = _file(inp.receptor, {".pdbqt"}, "receptor", hashes)
        _molecule_file(receptor, receptor=True)
        ligands = [_file(value, {".pdbqt"}, f"ligand_{i}", hashes) for i, value in enumerate(inp.ligands)]
        if len({path.stem for path in ligands}) != len(ligands):
            raise ValueError("Uni-Dock batch ligand basenames must be unique to prevent output overwrites.")
        for ligand in ligands:
            _molecule_file(ligand, pdbqt=True)
        if executable and Path(executable).name == 'unidock-1.1.0':
            _review_unidock_110(receptor, ligands, par.box)
        argv = [executable or "unidock", "--receptor", str(receptor), "--gpu_batch", *map(str, ligands), "--scoring", "vina", "--search_mode", par.search_mode, "--seed", str(par.seed), "--num_modes", str(par.num_modes), "--dir", str(out / "poses"), *_box_arguments(par.box)]
        outputs = [str(out / "poses")]
    elif engine_id == "autodock4":
        inp, par = AutoDockInputs.model_validate(inputs), AutoDockParameters.model_validate(parameters)
        receptor = _file(inp.receptor, {".pdbqt"}, "receptor", hashes)
        ligand = _file(inp.ligand, {".pdbqt"}, "ligand", hashes)
        _molecule_file(receptor, receptor=True)
        _molecule_file(ligand, pdbqt=True)
        field = _file(inp.grid_field, {".fld"}, "grid_field", hashes)
        atom_types = sorted({line.split()[-1] for line in ligand.read_text().splitlines() if line.startswith(("ATOM  ", "HETATM"))})
        if not all(re.fullmatch(r"[A-Za-z][A-Za-z0-9]{0,3}", atom) for atom in atom_types) or set(atom_types) != set(inp.atom_maps):
            raise ValueError("AutoDock4 requires exactly one reviewed atom affinity map for each ligand atom type.")
        maps = [_file(inp.atom_maps[atom], {".map"}, f"atom_map_{atom}", hashes) for atom in atom_types]
        electro = _file(inp.electrostatic_map, {".map"}, "electrostatic_map", hashes)
        desolv = _file(inp.desolvation_map, {".map"}, "desolvation_map", hashes)
        if any(any(character.isspace() for character in str(path)) for path in [ligand, field, electro, desolv, *maps]):
            raise ValueError("AutoDock4 DPF paths cannot contain whitespace; stage them into a reviewed path first.")
        grid_headers = []
        for path in [*maps, electro, desolv]:
            header = {parts[0]: parts[1:] for line in path.read_text().splitlines()[:6] if (parts := line.split())}
            if not all(key in header for key in ("SPACING", "NELEMENTS", "CENTER")):
                raise ValueError("AutoGrid map is missing SPACING, NELEMENTS or CENTER metadata.")
            grid_headers.append(tuple(tuple(header[key]) for key in ("SPACING", "NELEMENTS", "CENTER")))
        if len(set(grid_headers)) != 1:
            raise ValueError("AutoDock4 maps do not share the same grid geometry.")
        atom_lines = [line for line in ligand.read_text().splitlines() if line.startswith(("ATOM  ", "HETATM"))]
        center = [sum(float(line[a:b]) for line in atom_lines) / len(atom_lines) for a, b in ((30,38),(38,46),(46,54))]
        torsion = next(line.split()[1] for line in ligand.read_text().splitlines() if line.startswith("TORSDOF "))
        if not torsion.isdigit() or int(torsion) > 32:
            raise ValueError("Unsupported AutoDock4 torsion count.")
        dpf = ["autodock_parameter_version 4.2", "outlev adt", f"seed {par.seed} {par.seed + 1}", "unbound_model bound", "ligand_types " + " ".join(atom_types), f"fld {field}", *[f"map {path}" for path in maps], f"elecmap {electro}", f"desolvmap {desolv}", f"move {ligand}", "about " + " ".join(map(str, center)), "tran0 random", "quat0 random", "dihe0 random", f"torsdof {torsion}", f"ga_pop_size {par.ga_pop_size}", f"ga_num_evals {par.ga_num_evals}", "ga_num_generations 27000", "ga_elitism 1", "ga_mutation_rate 0.02", "ga_crossover_rate 0.8", "set_ga", "sw_max_its 300", "sw_max_succ 4", "sw_max_fail 4", "sw_rho 1.0", "sw_lb_rho 0.01", "ls_search_freq 0.06", "set_psw1", f"ga_run {par.ga_runs}", "rmstol 2.0", "analysis"]
        generated["protocol.dpf"] = "\n".join(dpf) + "\n"
        argv = [executable or "autodock4", "-p", str(out / "protocol.dpf"), "-l", str(out / "docking.dlg")]
        outputs = [str(out / "docking.dlg")]
    elif engine_id in {"mmpbsa", "mmgbsa"}:
        inp, par = MMInputs.model_validate(inputs), MMParameters.model_validate(parameters)
        tops = {role: _file(getattr(inp, role), {".prmtop", ".parm7", ".top"}, role, hashes)
                for role in ("solvated_topology", "complex_topology", "receptor_topology", "ligand_topology")}
        topology_data = {role: _validate_prmtop(path) for role, path in tops.items()}
        counts = {role: data["count"] for role, data in topology_data.items()}
        if counts["complex_topology"] != counts["receptor_topology"] + counts["ligand_topology"] or counts["solvated_topology"] < counts["complex_topology"]:
            raise ValueError("Amber topology atom counts do not match complex = receptor + ligand.")
        # This compact adapter supports receptor then ligand, followed by solvent/ions.
        # Other orderings require a reviewed reordering protocol, not silent remapping.
        for attribute in ("names", "charges"):
            complex_values = topology_data["complex_topology"][attribute]
            components = topology_data["receptor_topology"][attribute] + topology_data["ligand_topology"][attribute]
            solvent_prefix = topology_data["solvated_topology"][attribute][:counts["complex_topology"]]
            matching = (complex_values == components == solvent_prefix if attribute == "names" else
                        all(abs(a-b) < 1e-5 and abs(a-c) < 1e-5 for a, b, c in zip(complex_values, components, solvent_prefix)))
            if not matching:
                raise ValueError("Amber topology atom names/charges do not match receptor+ligand and solvated-prefix ordering.")
        trajectory = _file(inp.trajectory, {".nc", ".netcdf", ".mdcrd", ".crd"}, "trajectory", hashes)
        expected = {role: data["sha256"] for role, data in hashes.items()}
        if inp.review.reviewed_input_sha256 != expected:
            raise ValueError("MD review must bind the exact SHA-256 hashes of all four topologies and the trajectory.")
        general = f"&general\n startframe={par.start_frame}, endframe={par.end_frame}, interval={par.interval}, verbose=1,\n/\n"
        method = f"&pb\n istrng={par.salt_molar},\n/\n" if engine_id == "mmpbsa" else f"&gb\n igb={par.igb}, saltcon={par.salt_molar},\n/\n"
        generated["mmpbsa.in"] = "Reviewed endpoint estimate; entropy not included\n" + general + method
        argv = [executable or "MMPBSA.py", "-i", str(out / "mmpbsa.in"), "-o", str(out / "FINAL_RESULTS.dat"), "-sp", str(tops["solvated_topology"]), "-cp", str(tops["complex_topology"]), "-rp", str(tops["receptor_topology"]), "-lp", str(tops["ligand_topology"]), "-y", str(trajectory)]
        outputs = [str(out / "FINAL_RESULTS.dat")]
    else:
        inp, par = PlanetInputs.model_validate(inputs), RunParameters.model_validate(parameters)
        receptor = _file(inp.receptor, {".pdb"}, "receptor", hashes)
        pocket = _file(inp.pocket_reference, {".sdf"}, "pocket_reference", hashes)
        ligands = _file(inp.ligands, {".sdf"}, "ligands", hashes)
        _molecule_file(receptor, receptor=True)
        _molecule_file(pocket)
        _molecule_file(ligands, require_3d=False)
        if not any(line.startswith("ATOM  ") and line[12:16].strip() == "CA" for line in receptor.read_text().splitlines()):
            raise ValueError("PLANET receptor must contain correctly prepared residue alpha-carbon atoms.")
        installation, _ = _planet_installation()
        if installation:
            for filename in ("PLANET_run.py", "PLANET_model.py", "PLANET.param", "layers.py", "chemutils.py", "nnutils.py"):
                path = Path(installation[1]).parent / filename
                if path.is_file():
                    _file(str(path), {path.suffix}, f"installation_{filename}", hashes)
        argv = (installation or ["<configured PLANET Python>", "<configured PLANET_run.py>"]) + ["-p", str(receptor), "-l", str(pocket), "-m", str(ligands), "--prefix", str(out / "result")]
        outputs = [str(out / "result.csv"), str(out / "result.sdf")]
    if any(Path(item["path"]).is_relative_to(out) for item in hashes.values()):
        raise ValueError("Inputs cannot live inside the output directory.")
    if executable and Path(executable).is_file():
        hashes["executable"] = {"path": executable, "sha256": science.runtime_sha256(Path(executable)), "bytes": Path(executable).stat().st_size}
        # The reviewed GNINA launcher sets isolated CUDA paths. Bind the actual
        # native binary (including embedded CNN weights) and libraries as well.
        local_tools = Path(__file__).resolve().parents[1] / '.tools' / 'linux'
        runtime_manifest = local_tools / 'gpu-runtime-manifest.json'
        if Path(executable) == local_tools / 'gnina' and runtime_manifest.is_file():
            runtime = json.loads(runtime_manifest.read_text())
            _file(str(runtime_manifest), {'.json'}, 'runtime_manifest', hashes)
            for index, item in enumerate(runtime['files']):
                runtime_path = Path(item['path']).resolve()
                if not runtime_path.is_relative_to(local_tools):
                    raise ValueError('GPU runtime manifest escapes the configured local installation')
                if not runtime_path.is_file() or runtime_path.stat().st_size != item['bytes']:
                    raise ValueError('GPU runtime file size differs from the installation record')
                hashes[f'runtime_{index}'] = {'path': str(runtime_path), 'bytes': runtime_path.stat().st_size,
                                             'sha256': science.runtime_sha256(runtime_path)}
                if hashes[f'runtime_{index}']['sha256'] != item['sha256']:
                    raise ValueError('GPU runtime file changed since installation verification')
    return {"engine_id": engine_id, "status": "blocked" if blockers else "ready", "blockers": blockers,
            "commands": [argv], "output_dir": str(out), "expected_outputs": outputs, "generated_files": generated,
            "validated_inputs": hashes, "parameters": par.model_dump(mode="json"), "input_review": inp.model_dump(mode="json"),
            "output_metrics": entry["output_metrics"], "documentation": entry["docs"],
            "provenance": {"adapter_schema": "engines-v1", "adapter_sha256": science.sha256(Path(__file__)), "platform": platform.platform(), "engine_version": "not probed; executable hash recorded when installed"}}


_NUMBER = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?"


def parse_engine_output(engine_id: str, output: str | Path) -> list[dict]:
    """Parse engine-native fields only; never substitute one score for another."""
    if engine_id == "biocentric_affinity":
        return biocentric_model.parse_output(output)
    if engine_id not in _METRICS:
        raise ValueError(f"Unknown engine {engine_id!r}.")
    path = Path(output)
    metrics = {item["name"]: item for item in _METRICS[engine_id]}
    rows = []
    if engine_id in {"vina", "unidock"}:
        metric = "vina_score" if engine_id == "vina" else "unidock_vina_score"
        for line in path.read_text().splitlines():
            if not line.strip().startswith("REMARK VINA RESULT:"):
                continue
            values = line.split(":", 1)[1].split()
            if len(values) != 3 or any(not re.fullmatch(_NUMBER, value) for value in values):
                raise ValueError("Vina output contains a malformed score/RMSD record; partial pose scores are not success.")
            score, lower_rmsd, upper_rmsd = map(float, values)
            if not all(math.isfinite(value) for value in (score, lower_rmsd, upper_rmsd)) or lower_rmsd < 0 or upper_rmsd < 0:
                raise ValueError("Vina output contains nonfinite scores or invalid RMSD bounds.")
            rows.append({"pose_index": len(rows) + 1, "metrics": [{**metrics[metric], "value": score}]})
    elif engine_id in {"gnina", "gnina_affinity"}:
        from rdkit import Chem
        molecules = list(Chem.SDMolSupplier(str(path), removeHs=False))
        for index, mol in enumerate(molecules, 1):
            if mol is None:
                raise ValueError("GNINA returned an invalid SDF record.")
            parsed = []
            for name in ("minimizedAffinity", "CNNscore", "CNNaffinity"):
                if not mol.HasProp(name):
                    raise ValueError(f"GNINA output is missing {name}; no substituted/fabricated score.")
                value = float(mol.GetProp(name))
                if not math.isfinite(value) or name == "CNNscore" and not 0 <= value <= 1:
                    raise ValueError(f"GNINA returned invalid {name}.")
                parsed.append({**metrics[name], "value": value})
            rows.append({"pose_index": index, "name": mol.GetProp("_Name") if mol.HasProp("_Name") else None, "metrics": parsed})
    elif engine_id == "autodock4":
        marker = r"DOCKED:\s+USER\s+Estimated Free Energy of Binding\s*="
        for line in path.read_text().splitlines():
            if not re.match(marker, line.strip()):
                continue
            match = re.fullmatch(marker + r"\s*(" + _NUMBER + r")\s+kcal/mol(?:\s+.*)?", line.strip())
            if match is None:
                raise ValueError("AutoDock4 returned an invalid binding-energy record or unsupported units; no run is silently skipped.")
            rows.append({"pose_index": len(rows) + 1, "metrics": [{**metrics["autodock4_score"], "value": float(match.group(1))}]})
    elif engine_id in {"mmpbsa", "mmgbsa"}:
        content = path.read_text()
        heading = "POISSON BOLTZMANN" if engine_id == "mmpbsa" else "GENERALIZED BORN"
        methods = list(re.finditer(heading + r"\s*:(.*?)(?=POISSON BOLTZMANN\s*:|GENERALIZED BORN\s*:|\Z)", content, re.S))
        if not methods:
            raise ValueError(f"Amber output does not contain the requested {heading} calculation.")
        if len(methods) != 1:
            raise ValueError("Amber output has multiple requested method sections; the estimate is ambiguous.")
        difference = re.search(r"Differences\s*\(Complex\s*-\s*Receptor\s*-\s*Ligand\)\s*:(.*)", methods[0].group(1), re.S)
        if not difference:
            raise ValueError("Amber output lacks the complex-minus-components summary.")
        units = re.findall(r"All units are(?: reported)? in\s+([^\r\n.]+)", content, re.I)
        if not units or any(unit.strip().lower() not in {"kcal/mol", "kcal/mole"} for unit in units):
            raise ValueError("Amber energy units must explicitly be kcal/mol or kcal/mole; no unit conversion is inferred.")
        totals = [line.strip() for line in difference.group(1).splitlines()
                  if re.match(r"(?:DELTA TOTAL|DELTA G binding\s*=)", line.strip())]
        if len(totals) > 1:
            raise ValueError("Amber output has multiple difference totals; the estimate is ambiguous.")
        if totals:
            match = re.fullmatch(r"(?:DELTA TOTAL|DELTA G binding\s*=)[ \t]+(" + _NUMBER + r")[ \t]+(?:\+/-[ \t]*)?(" + _NUMBER + r")[ \t]+(" + _NUMBER + r")", totals[0])
            if match is None:
                raise ValueError("Amber difference summary must contain complete mean, standard deviation and standard error values.")
            mean, std, sem = map(float, match.groups())
            if not all(math.isfinite(value) for value in (mean, std, sem)) or std < 0 or sem < 0:
                raise ValueError("Amber mean and uncertainty must be finite; uncertainty cannot be negative.")
            rows = [{"metrics": [{**metrics["delta_total"], "value": mean, "standard_deviation": std, "standard_error": sem}], "entropy_included": False}]
    elif engine_id == "planet":
        from rdkit import Chem
        reader = csv.DictReader(io.StringIO(path.read_text()))
        if not {"mol_name", "SMILES", "PLANET_affinity"} <= set(reader.fieldnames or []):
            raise ValueError("PLANET native CSV requires mol_name, SMILES and PLANET_affinity columns.")
        if len(reader.fieldnames) != len(set(reader.fieldnames)):
            raise ValueError("PLANET native CSV contains duplicate columns.")
        for index, row in enumerate(reader, 1):
            if None in row or any(value is None for value in row.values()):
                raise ValueError("PLANET returned an incomplete or malformed CSV record.")
            molecule = Chem.MolFromSmiles(row["SMILES"])
            if molecule is None or molecule.GetNumAtoms() == 0:
                raise ValueError("PLANET returned an invalid or missing molecular identity.")
            rows.append({"record_index": index, "name": row["mol_name"], "smiles": row["SMILES"],
                         "metrics": [{**metrics["PLANET_affinity"], "value": float(row["PLANET_affinity"])}], "generated_docking_pose": False})
    if not rows:
        raise ValueError("No complete native engine scores found; the run cannot be reported as successful.")
    if any(not math.isfinite(metric["value"]) for row in rows for metric in row["metrics"]):
        raise ValueError("Nonfinite engine scores are invalid.")
    for row in rows:
        row["source"] = {"path": str(path.resolve()), "sha256": science.sha256(path)}
    return rows


def _bounded_native_json(path: Path, root: Path, limit: int = 512 * 1024):
    """Worker progress/checkpoint data never authorizes reads outside its run."""
    resolved = path.resolve(strict=True)
    if not resolved.is_relative_to(root.resolve()) or not resolved.is_file():
        raise ValueError('Native progress/checkpoint escaped the output directory')
    with resolved.open('rb') as stream:
        content = stream.read(limit + 1)
    if len(content) > limit:
        raise ValueError('Native progress/checkpoint exceeds the bounded read limit')
    def nonfinite(value):
        raise ValueError(f'Native JSON contains nonfinite value {value}')
    payload = json.loads(content, parse_constant=nonfinite)
    if not isinstance(payload, dict):
        raise ValueError('Native progress/checkpoint must be a JSON object')
    return content, payload


def _biocentric_checkpoint(plan: dict, expected: dict, request_sha256: str):
    out = Path(plan['output_dir'])
    path = out / 'pairs' / f"{expected['record_index']:06d}.json"
    content, checkpoint = _bounded_native_json(path, out)
    row = checkpoint.get('row')
    if not isinstance(row, dict):
        raise ValueError('Pair checkpoint has no prediction row')
    row_sha = hashlib.sha256(json.dumps(row, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()
    if checkpoint.get('request_sha256') != request_sha256 or checkpoint.get('row_sha256') != row_sha:
        raise ValueError('Pair checkpoint request or row hash differs from the reviewed request')
    parsed = biocentric_model.parse_payload({'schema_version': 1, 'engine_id': 'biocentric_affinity',
                                            'predictions': [row]}, [expected], expected_provenance=plan['provenance'],
                                            expected_device=plan['parameters']['device'])[0]
    return row, {**parsed, 'source': {'path': str(path), 'sha256': hashlib.sha256(content).hexdigest()},
                 'checkpoint_verified': True}


def _biocentric_failure_row(plan: dict, expected: dict, reason: str, code: str):
    return {**expected, 'status': 'failed', 'score': None, 'per_model': {}, 'metrics': [],
            'name': expected['compound_id'], 'failure_reason': reason, 'failure_code': code,
            'preprocessing_status': 'not_completed', 'checkpoint_sha256': plan['provenance']['weights_sha256'],
            'code_sha256': plan['provenance']['code_sha256'], 'models': plan['provenance']['models'],
            'device': plan['parameters']['device'], 'checkpoint_recovered': False}


def _recover_biocentric(plan: dict, reason: str, code: str):
    """Keep verified finished pairs and explicitly account for every other pair."""
    out = Path(plan['output_dir'])
    integrity_error = None
    try:
        request_bytes, _ = _bounded_native_json(out / 'request.json', out, 16 * 1024 * 1024)
        if request_bytes != plan['generated_files']['request.json'].encode():
            raise ValueError('Native request changed after planning')
        request_sha = science.sha256(out / 'request.json')
        for item in plan['validated_inputs'].values():
            if science.sha256(Path(item['path'])) != item['sha256']:
                raise ValueError('A reviewed input or runtime file changed during execution')
    except (OSError, ValueError) as error:
        integrity_error = str(error)
    rows = []
    for expected in plan['expected_pairs']:
        checkpoint_issue = integrity_error
        if integrity_error is None:
            try:
                _, recovered = _biocentric_checkpoint(plan, expected, request_sha)
                rows.append({**recovered, 'checkpoint_recovered': True})
                continue
            except FileNotFoundError:
                pass
            except (OSError, ValueError, KeyError, TypeError) as error:
                checkpoint_issue = str(error)
        explanation = reason + (f'; checkpoint not accepted: {checkpoint_issue}' if checkpoint_issue else '')
        rows.append(_biocentric_failure_row(plan, expected, explanation, code))
    science._write_json(out / 'partial-predictions.json', {'schema_version': 1, 'engine_id': 'biocentric_affinity',
                         'expected_pairs': plan['expected_pairs'], 'predictions': rows,
                         'systemic_failure': {'code': code, 'message': reason},
                         'checkpoint_files_rewritten': False})
    return rows


def run_engine(engine_id: str, inputs: dict, parameters: dict, output_dir: str | Path,
               emit: Callable | None = None, cancelled: Callable | None = None, *, resume_plan: dict | None = None) -> dict:
    """Run a fresh validated plan. Missing tools return blockers without execution."""
    plan = plan_engine_run(engine_id, inputs, parameters, output_dir, **({'resume_plan': resume_plan} if resume_plan is not None else {}))
    if plan["blockers"]:
        return {**plan, "status": "blocked", "results": [], "executed": False}
    out = Path(plan["output_dir"])
    out.mkdir(parents=True, exist_ok=True)
    science._write_json(out / "plan.json", plan)
    for name, content in plan["generated_files"].items():
        (out / name).write_text(content)
    if engine_id == "unidock":
        (out / "poses").mkdir()
    started = time.monotonic()
    execution_attempted = False
    progress_state = {'content': None, 'completed_pairs': 0, 'resource_usage': {}, 'warnings': []}

    def poll_progress():
        if engine_id != 'biocentric_affinity':
            return
        progress_event = None
        try:
            content, progress = _bounded_native_json(out / 'progress.json', out)
            if content == progress_state['content']:
                return
            progress_state['content'] = content
            completed = progress.get('completed_pairs')
            total = len(plan['expected_pairs'])
            if (type(completed) is not int or not 1 <= completed <= total or
                    progress.get('total_pairs') != total or completed < progress_state['completed_pairs']):
                raise ValueError('Native progress count does not match the reviewed pair list')
            expected = plan['expected_pairs'][completed - 1]
            reviewed_request_sha = hashlib.sha256(plan['generated_files']['request.json'].encode()).hexdigest()
            raw, latest = _biocentric_checkpoint(plan, expected, reviewed_request_sha)
            if progress.get('latest') != raw:
                raise ValueError('Native progress latest row differs from its hash-verified checkpoint')
            progress_state['completed_pairs'] = completed
            resources = progress.get('resource_usage', progress.get('resources', {}))
            if isinstance(resources, dict):
                progress_state['resource_usage'] = resources
            progress_event = {'kind': 'progress', 'tool': engine_id, 'stage': 'engine_execution',
                      'message': f"BioCentric processed {completed}/{total} protein–compound pairs.",
                      'compound_id': latest['compound_id'], 'progress': completed / total,
                      'data': {'completed_pairs': completed, 'total_pairs': total, 'latest': latest,
                               'device': plan['parameters']['device'], 'resource_usage': progress_state['resource_usage']}}
        except FileNotFoundError:
            return
        except (OSError, ValueError, KeyError, TypeError) as error:
            warning = str(error)
            if warning not in progress_state['warnings']:
                progress_state['warnings'].append(warning)
        if emit and progress_event is not None:
            emit(progress_event)

    def check_cancelled():
        if cancelled and cancelled():
            return True
        poll_progress()
        return False

    try:
        for role, item in plan["validated_inputs"].items():
            hash_file = science.runtime_sha256 if role == "executable" or role.startswith(("runtime_", "model_file_")) else science.sha256
            if hash_file(Path(item["path"])) != item["sha256"]:
                raise science.ScienceError('input_changed', "An input changed after planning; create a new reviewed plan.")
        if emit:
            emit({"kind": "tool_started", "tool": engine_id, "stage": "engine_execution", "message": f"Run {engine_id} with the reviewed typed protocol.", "argv": plan["commands"][0]})
        execution_attempted = True
        execution = science.run_process(plan["commands"][0], out, out / "engine.log", plan["parameters"]["timeout_seconds"],
                                      check_cancelled if engine_id == 'biocentric_affinity' else cancelled)
        poll_progress()
        if engine_id == "unidock":
            outputs = sorted((out / "poses").glob("*.pdbqt"))
            expected_names = {Path(ligand).stem + "_out.pdbqt" for ligand in plan["input_review"]["ligands"]}
            if {path.name for path in outputs} != expected_names:
                raise ValueError("Uni-Dock did not produce one output file per ligand with the expected input-bound filename; partial/substituted results are not success.")
        else:
            outputs = [Path(plan["expected_outputs"][0])]
        if engine_id == 'biocentric_affinity':
            _, native_payload = _bounded_native_json(outputs[0], out, 16 * 1024 * 1024)
            rows = biocentric_model.parse_payload(native_payload, plan['expected_pairs'],
                       expected_provenance=plan['provenance'], expected_device=plan['parameters']['device'])
        else:
            rows = [row for path in outputs for row in parse_engine_output(engine_id, path)]
        for role, item in plan["validated_inputs"].items():
            hash_file = science.runtime_sha256 if role == "executable" or role.startswith(("runtime_", "model_file_")) else science.sha256
            if hash_file(Path(item["path"])) != item["sha256"]:
                raise ValueError("An input changed during execution; results do not have stable provenance.")
        if engine_id == "planet":
            from rdkit import Chem
            expected = sorted(science._canonical(mol) for mol in Chem.SDMolSupplier(plan["input_review"]["ligands"], removeHs=False))
            observed = sorted(science._canonical(Chem.MolFromSmiles(row["smiles"])) for row in rows)
            if observed != expected:
                raise ValueError("PLANET skipped, duplicated or substituted input molecular identities; incomplete output requires review.")
            native_sdf = Path(plan["expected_outputs"][1])
            _file(str(native_sdf), {".sdf"}, "PLANET output SDF", {})
            _molecule_file(native_sdf, require_3d=False)
            sdf_identities = sorted(science._canonical(mol) for mol in Chem.SDMolSupplier(str(native_sdf), removeHs=False))
            if sdf_identities != expected:
                raise ValueError("PLANET output SDF does not match the complete input molecular identities.")
        result = {**plan, "status": "completed", "executed": True, "execution": execution,
                  "results": rows, "duration_seconds": time.monotonic() - started,
                  "scientific_validation": "Not established by tool completion; inspect pose checks and the independent test set."}
        if engine_id == "biocentric_affinity":
            result['resource_usage'] = native_payload.get('resource_usage', native_payload.get('resources', {}))
            result['native_runtime'] = native_payload.get('runtime', {})
            failed = [row for row in rows if row.get("status") == "failed"]
            result["accounting"] = {"requested": len(plan["expected_pairs"]), "completed": len(rows)-len(failed), "failed": len(failed)}
            if failed:
                result["status"] = "partial" if len(failed) < len(rows) else "failed"
    except (science.ScienceError, ValueError, OSError) as error:
        result = {**plan, "status": "cancelled" if getattr(error, "code", None) == "cancelled" else "failed",
                  "executed": execution_attempted, "results": [], "duration_seconds": time.monotonic() - started,
                  "error": {"code": getattr(error, "code", "invalid_engine_output"), "message": str(error)}}
        if engine_id == 'biocentric_affinity':
            result['results'] = _recover_biocentric(plan, str(error), result['error']['code'])
            failed = sum(row['status'] == 'failed' for row in result['results'])
            result['accounting'] = {'requested': len(plan['expected_pairs']),
                                    'completed': len(result['results']) - failed, 'failed': failed}
            result['resource_usage'] = progress_state['resource_usage']
            if result['status'] != 'cancelled' and result['accounting']['completed']:
                result['status'] = 'partial'
    if engine_id == 'biocentric_affinity':
        result['progress_warnings'] = progress_state['warnings']
    science._write_json(out / "result.json", result)
    if emit and not (cancelled and cancelled()):
        emit({"kind": "tool_completed" if result["status"] == "completed" else "tool_failed", "tool": engine_id,
              "stage": "engine_execution", "message": f"{engine_id}: {result['status']}", "result_path": str(out / "result.json")})
    return result
