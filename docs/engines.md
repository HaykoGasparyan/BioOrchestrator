# Docking and affinity engine adapters

The local adapter layer is `biocentric/engines.py`. It distinguishes **implemented adapter**, **installed executable**, **completed calculation**, and **validated scientific method**. These are different states. Vina is the installed, previously live-tested local baseline. The other engines are implemented optional adapters and stay blocked until their dependencies and scientific inputs are supplied. No GPU server is configured and no software/weights are downloaded automatically.

## Contracts and availability

| ID | Native inputs and preprocessing | Native outputs and postprocessing | Resources / local readiness |
|---|---|---|---|
| `vina` | Reviewed prepared receptor and ligand PDBQT; an explicit box in Å. Use the established RDKit/Meeko pipeline to preserve ligand identity in SDF. | PDBQT poses; `REMARK VINA RESULT` parsed as empirical Vina score, kcal/mol, lower better. Existing workflow performs SDF export, identity and PoseBusters checks. | CPU, installed locally. |
| `autodock4` | Reviewed receptor and ligand PDBQT, AutoGrid field, one atom map per ligand atom type, electrostatic and desolvation maps. Adapter checks common grid geometry and generates a bounded deterministic-seed LGA DPF. Map/receptor correspondence remains an explicit scientific review. | DLG, per-pose estimated binding energies in kcal/mol, lower better. Retains raw DLG; downstream SDF conversion and pose validation must occur before use in screening. | CPU; AutoDock4 installation and prepared maps required. This is AutoDock4, not Vina with `--scoring ad4`. |
| `gnina` | Prepared receptor PDB/PDBQT and ligand SDF, explicit box; CNN rescore with Vina empirical scoring. CPU mode is explicit (`--no_gpu`); GPU requires a local Linux NVIDIA device. | SDF with separate `minimizedAffinity` (kcal/mol), `CNNscore` (0–1 pose score), and `CNNaffinity` (predicted pK). Never merges the metrics. | GNINA executable and bundled CNNs required; CPU can be slow. Current Mac has no configured GNINA. |
| `unidock` | Reviewed receptor PDBQT and unique-basename ligand PDBQT batch, box, Vina scoring and one declared search mode. | One output PDBQT per ligand; Vina-score records with provenance. Missing per-ligand files fail the run. | Officially Linux/NVIDIA GPU; blocked on this Mac. Remote dispatch is not implemented. |
| `mmpbsa` | Four Amber topologies: solvated, complex, receptor, ligand; an equilibrated MD trajectory. A human MD review binds the exact hashes and records force field, ligand parameters, atom ordering and sampling review. Atom-count/name/charge checks enforce receptor+ligand order and its solvated prefix; equilibration, full force-field correspondence and trajectory atom order also require review. Other topological orderings are rejected rather than silently remapped. | Amber PB complex-minus-receptor-minus-ligand mean, standard deviation and standard error, kcal/mol. Entropy is explicitly absent. | AmberTools `MMPBSA.py`, its solvers and trajectory-processing dependencies; bounded 2–1000-frame serial calculation. Not derivable from a docked SDF alone. |
| `mmgbsa` | Same reviewed MD contract, explicit GB model and salt concentration. | Amber GB difference estimate with uncertainties, kcal/mol; entropy absent. | Same AmberTools requirement. No MD or force-field parameterization is silently generated. |
| `gnina_affinity` | Reviewed receptor and already placed 3D ligand SDF, preserving receptor coordinates; `--score_only`, without docking. | GNINA CNN affinity pK plus its separate empirical and pose scores. Higher pK is better; this is not an experimental Kd/Ki or a rigorous free energy. | Uses the same GNINA installation; defaults to CPU. |
| `planet` | Prepared protein PDB with intact alpha carbons, a 3D pocket-defining reference ligand SDF, and query ligand SDF representing reviewed chemical states. PLANET uses the ligand's 2D graph; a docked pose is not required. | Native CSV columns `mol_name`, `SMILES`, `PLANET_affinity` and output SDF. Predicted pK, higher better. The output SDF is not a generated 3D docking pose. Silent upstream skipped molecules become a failed/incomplete result. | Reviewed upstream PLANET checkout, dedicated Python environment and matching weight SHA-256. CPU supported by upstream. Unavailable until configured. This adapter targets original PLANET, not PLANET v2.0. |

Installation readiness is discovered without executing software, importing ML frameworks, connecting to a server or provisioning hardware. A detected executable can still fail due to missing runtime libraries; that failure is recorded. GPU device detection is a prerequisite, not a successful CUDA benchmark. Optional tools have contract/format tests but have not been live-run here.

## API

```python
from biocentric.engines import engine_catalog, plan_engine_run, run_engine

catalog = engine_catalog()
inputs = {
    "receptor": "/absolute/path/receptor.pdbqt",
    "ligand": "/absolute/path/ligand.pdbqt",
    "preparation_reviewed": True,
}
parameters = {
    "box": {"center": [15.190, 53.903, 16.917], "size": [20, 20, 20]},
    "cpu": 1, "seed": 42, "exhaustiveness": 8, "num_modes": 3,
    "timeout_seconds": 120,
}
plan = plan_engine_run("vina", inputs, parameters, "/new/output/directory")
# Display plan, input hashes, resource limits and blockers before execution.
result = run_engine("vina", inputs, parameters, "/new/output/directory",
                    emit=on_event, cancelled=is_cancelled)
```

`plan_engine_run` makes no directories and runs no programs. Missing software becomes explicit blockers. Malformed requests, unreviewed chemical preparation, unsupported formats, invalid numbers, unsupported options, changed MD inputs and nonempty output directories fail closed. It returns fixed argv and generated configuration contents so the UI can display exactly what will run. No LLM-generated shell text or unbounded extra arguments are accepted.

`run_engine` creates `plan.json`, generated configuration files, `engine.log` and `result.json`; the latter includes tool argv, timing, input and executable hashes, raw output hashes, named scores, score directions and scientific limitations. Inputs are hashed again before and after execution. Previously populated output directories cannot be reused. Cancellation and timeouts terminate the scientific process group through the established scientific runner. The scientific subprocess receives an environment without application API credentials. Successful exit without complete native scores is a failed result, not success. A low score does not itself establish binding or validate a method.

## Setup guidance

Configuration is administrative environment configuration, never agent-supplied command text:

- Vina: `BIOCENTRIC_VINA`, `PATH`, or the project's existing `.tools/vina`.
- AutoDock4: `BIOCENTRIC_AUTODOCK4` or `PATH`. Prepare AutoGrid4 maps separately after receptor/protonation/box review. The adapter generates a DPF from these files, not arbitrary external DPF instructions.
- GNINA: `BIOCENTRIC_GNINA` or `PATH`; follow the official release/build requirements. GPU use is opt-in through typed `device="gpu"`; CPU uses `--no_gpu`.
- Uni-Dock: `BIOCENTRIC_UNIDOCK` or `PATH` on a verified Linux NVIDIA installation. An environment path cannot turn macOS into a CUDA worker. Do not claim this is available through Apple Metal.
- Amber methods: `BIOCENTRIC_MMPBSA` or `PATH` pointing to `MMPBSA.py` in a functional AmberTools environment. The input review's `reviewed_input_sha256` must contain exactly the five named scientific input hashes. Require an independently reviewed MD setup and sampling evidence.
- PLANET: `BIOCENTRIC_PLANET_DIR` for the upstream repository containing `PLANET_run.py`, `PLANET_model.py`, `PLANET.param`, `layers.py`, `chemutils.py`, `nnutils.py`; `BIOCENTRIC_PLANET_PYTHON` for its dedicated environment; and `BIOCENTRIC_PLANET_WEIGHTS_SHA256` for the reviewed weight digest. Configure CPU dependencies per upstream instructions; do not load unreviewed model weights. The adapter invokes the documented `-p`, `-l`, `-m`, `--prefix` interface. Code and weight files are hashed into each plan. Dependencies are validated by actual execution rather than claimed from file presence.

Each optional installation must pass a native reference example plus the application's held-out scientific test set before being selected for production screening. Do not compare raw Vina and AutoDock4 energies as if calibrated on the same scale; compare pose recovery/enrichment/error against the same held-out references. Keep training-set overlap, scaffold/target leakage, pose validity and assay endpoint compatibility explicit for neural predictors. MM endpoint estimates require convergence and uncertainty review; docking scores and ML pK are different types of evidence.

## Primary references consulted

- [AutoDock Vina official basic docking tutorial](https://autodock-vina.readthedocs.io/en/stable/docking_basic.html): native preparation, box configuration, affinity maps, CLI and score limitations.
- [AutoDock 4.2.6 official user guide](https://autodock.scripps.edu/wp-content/uploads/sites/56/2022/04/AutoDock4.2.6_UserGuide.pdf): atom maps, DPF and LGA controls.
- [GNINA official repository](https://github.com/gnina/gnina): docking/rescoring CLI, CPU flag, native output and runtime prerequisites.
- [GNINA developers' usage guide](https://github.com/dkoes/docs/wiki/Gnina) and [maintainer clarification of CNNaffinity units](https://github.com/gnina/gnina/issues/253): distinct pose, energy and affinity metrics.
- [Uni-Dock official native-engine instructions](https://github.com/dptech-corp/Uni-Dock/blob/main/unidock/README.md): Linux/NVIDIA support, batch PDBQT CLI and search modes. [Original paper](https://doi.org/10.1021/acs.jctc.2c01145).
- [Amber official MM/PB(GB)SA tutorial](https://ambermd.org/tutorials/advanced/tutorial3/index.php) and [MMPBSA.py invocation/output](https://ambermd.org/tutorials/advanced/tutorial3/py_script/section1.htm): equilibrated snapshots, matching topologies, native solver selection and output summaries.
- [PLANET official repository](https://github.com/ComputArtCMCG/PLANET), [native CLI](https://raw.githubusercontent.com/ComputArtCMCG/PLANET/main/PLANET_run.py), and [model definition](https://raw.githubusercontent.com/ComputArtCMCG/PLANET/main/PLANET_model.py): protein pocket/2D ligand inputs, actual output columns and model weights.

References inspected 2026-09-26. Executable/weight hashes in each run are the reproducibility identifiers; an unpinned upstream web page is not a runtime version guarantee.

## Configured engines in the research workflow

The interactive research workflow uses `biocentric/engine_workflow.py` to preflight and run a **registered native manifest** for configured optional tools. This route executes the real adapters; it is not a catalog-only placeholder. The default Vina path continues to use the established end-to-end preparation/docking/pose-check pipeline.

Upload native scientific files as `native_input` artifacts, then upload one JSON manifest as `native_engine_pilot`, `native_engine_screening` or `native_engine_affinity`. Every input references its registered project artifact ID, never an arbitrary path. The manifest engine and phase must match the LLM's reviewed selection. Hashes must match the frozen test set and reviewed protocol displayed in the research card.

Example shape for a GNINA pilot (replace placeholders with the current project's artifact IDs and exact frozen hashes):

```json
{
  "schema_version": 1,
  "engine_id": "gnina",
  "phase": "pilot",
  "input_artifacts": {
    "receptor": "REGISTERED_RECEPTOR_PDB_ARTIFACT_ID",
    "ligand": "REGISTERED_PREPARED_REFERENCE_SDF_ARTIFACT_ID"
  },
  "parameters": {
    "box": {"center": [15.190, 53.903, 16.917], "size": [20, 20, 20]},
    "cpu": 1,
    "seed": 42,
    "exhaustiveness": 8,
    "num_modes": 3,
    "device": "cpu",
    "timeout_seconds": 120
  },
  "review": {"preparation_reviewed": true},
  "human_reviewed": true,
  "canonical_ligands_artifact": "REGISTERED_PREPARED_REFERENCE_SDF_ARTIFACT_ID",
  "validation_receptor_artifact": "REGISTERED_MATCHING_RECEPTOR_PDB_ARTIFACT_ID",
  "reference_artifact": "REGISTERED_OBSERVED_REFERENCE_SDF_ARTIFACT_ID",
  "testset_sha256": "EXACT_64_HEX_FROZEN_TESTSET_HASH",
  "protocol_sha256": "EXACT_64_HEX_REVIEWED_PROTOCOL_HASH",
  "exclusions": []
}
```

`input_artifacts` has the same file-role keys as the chosen engine's typed input contract. For Uni-Dock `ligands` is a list of artifact IDs; for AutoDock4 `atom_maps` maps atom types to artifact IDs. Non-file values belong in `review` and cannot override artifact bindings. MM methods put their `review` object under the manifest's `review.review`, while the five file fields are artifact IDs. The inner MD review retains the original input hashes and force-field/sampling assumptions. PLANET supplies `receptor`, `pocket_reference` and `ligands` artifact IDs; canonical docking-pose validation fields are not required because PLANET does not produce docking poses.

```python
from biocentric.engine_workflow import plan_native_manifest, run_native_manifest

preflight = plan_native_manifest(
    store, project_id, manifest_artifact, new_output_directory,
    expected_engine="gnina", expected_phase="pilot",
    expected_testset_sha256=frozen_testset_hash,
    expected_protocol_sha256=reviewed_protocol_hash,
)
# Display exact argv, input bindings, resource limits and blockers; ask for review.
result = run_native_manifest(
    store, project_id, manifest_artifact, new_output_directory,
    expected_engine="gnina", expected_phase="pilot",
    expected_testset_sha256=frozen_testset_hash,
    expected_protocol_sha256=reviewed_protocol_hash,
    job_id=job_id, emit=record_event, cancelled=is_cancelled,
)
```

Before execution the helper checks native ligand chemical identities against canonical SDF and verifies the receptor coordinate frame. For a pilot it requires one cognate reference with the exact frozen reference artifact hash. For screening, selected chemical states must be a subset of the frozen library; each omitted original record needs an `exclusions` entry with `record_index` (one-based), `canonical_smiles` (or null for an unparseable record), and a concrete `reason`. Arbitrary new compounds cannot be introduced by copying a test-set hash into a manifest. Protocol-specific protonation changes require a separately reviewed frozen input set rather than a silent identity change. Native and validation receptor hashes must match the reviewed protocol; after a native pilot, screening must preserve its box, search, seed, device and scoring settings (a longer authorized timeout is permitted). A copied protocol hash cannot substitute different inputs or an untested search configuration.

After docking, GNINA native SDF or Meeko-mapped Vina/Uni-Dock PDBQT/AutoDock4 DLG is converted to canonical SDF. Missing atom mapping, changed chemical identity, mismatched score/pose counts or missing compounds fail validation. Each pose receives the existing physical-validity checks; pilot RMSD is symmetry-aware in the receptor frame. Standard pose dictionaries carry native named metrics, validity checks, RMSD, descriptors, receptor/pose hashes and immutable artifact links. The registry includes raw logs and failed outputs as well as successful canonical poses. Optional-native runs do not add PLIP automatically; its status remains `not_requested`.

Local verification: the new Vina native runner completed a one-pose CPU calculation on the bundled prepared caffeine/1IEP inputs in 0.79 seconds, producing an empirical score of −5.399 kcal/mol. Recorded under `runs/engine-adapter-vina-verification/`. This checks execution and parsing, not caffeine binding or predictive accuracy. Optional engines remain contract-tested only on this host.
