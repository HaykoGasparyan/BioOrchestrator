# Local BioCentric model contract

Inspected on Ubuntu on 2026-09-26. This document describes the model in `/home/sargis/dti`, not PLANET and not the application's Python package named `biocentric`. Source provenance is **local file**, not a web citation. The upstream checkout identifies commit `70af11ccaffdcab41f76a893a5fec9e57713d6a8` but has substantial modified/untracked code; execution therefore binds individual source/config/checkpoint SHA-256 hashes, not just that commit.

## Entry point, supported weights, and actual status

The functioning reference entry point is [Prediction/run_pipeline.py](/home/sargis/dti/Prediction/run_pipeline.py). Its `--poses` path accepts a prepared protein-pocket PDB and already positioned ligand PDB/SDF files; it runs without BioNeMo. Its `--smiles` route invokes DiffDock and requires additional BioNeMo assets; the application bridge does not invoke that route.

All nine checkpoint files were present and successfully loaded for actual inference:

| Directory under `Prediction/models` | Epoch | Target crop |
| --- | --- | --- |
| `static_rad_15_15_15_cp_cut_false_predict_pairwise_pos_TRUE_reproduce_pair_false_big_radius` | 127 | 15 Å reference-ligand crop |
| `static_rad_7_7_cp_cut_false_predict_pairwise_pos_TRUE_reproduce` | 73 | 7 Å reference-ligand crop |
| `true_embeddings_ligand_PRETRAINING_LOSSES` | 240 | 10 Å crop inside complex aggregation |

Each directory has `Checkpoints/Epoch-N/{ligand,target,aggregation}_model.pt`, a frozen `config.yaml`, `dataset.py` and `Networks/*.py`. Despite its name, the third checkpoint's config sets both `ligand_pretraining` and `protein_pretraining` to false. No ESM or molecular-language-model embedding is required by this ensemble. The supported weights aggregate hash is `a806ae21f91be1a116054a6470c47a760a6499231488de527bfc8cbebbf42ed8`: SHA-256 of compact sorted JSON mapping relative checkpoint filenames to SHA-256 hashes. Every individual hash is retained in execution `request.json` and `plan.json`.

The older [Inference/run.py](/home/sargis/dti/Inference/run.py) is not used: its `Config('config.yaml')` call does not match that directory's two-argument Config implementation, and its output assembly differs. [Prediction/inference.py](/home/sargis/dti/Prediction/inference.py) is a dataset-oriented script with historical absolute paths. [Prediction/evaluate.py](/home/sargis/dti/Prediction/evaluate.py) evaluates saved datasets from the repository root; it is not the interactive pair adapter.

## Scientific input contract

This is a **structure/graph affinity regressor**, not a sequence-only DTI probability classifier. Both protein and ligand coordinates in one receptor frame affect predictions. A real prepared pose is required for this particular model. The bridge does not invent a docking result, align a ligand, generate a conformer or assign a target sequence.

The protein input is a PDB with reviewed construct and chains. The native code does not check species, UniProt identity, mutations, completeness or experimental resolution. The application records the supplied `protein_id`, receptor hash, chain IDs, residue count and atom records. Users must retain construct/species provenance in the registered project artifacts. There is no sequence tokenizer, maximum amino-acid token length or sequence padding. The compact adapter bounds protein ATOM records to 20,000 and ligand heavy atoms to 256 for resources; larger systems require an explicitly reviewed pocket/protocol.

The default `unique_compounds` input mode accepts one prepared SDF record per compound, with unique `_Name` IDs and unique canonical isomeric SMILES. The explicit `pose_records` mode also accepts repeated chemical states with distinct pose IDs; the discovery integration uses a manifest to bind those IDs to measurements and pose indices. All records need finite three-dimensional coordinates. Charge, stereochemistry, tautomer and protonation must already be selected. The bridge does not neutralize, enumerate tautomers, strip salts or silently assign stereochemistry. A receptor PDB with HETATM, alternate locations or multiple models is rejected for explicit preparation because upstream silently drops heterogens and chooses/repairs some structures. Strict RDKit receptor sanitization occurs before upstream construction; invalid structures are rejected rather than allowing upstream `BioObject.check` to remove a problematic atom.

Native [Prediction/bioObjects.py](/home/sargis/dti/Prediction/bioObjects.py), `Protein`, removes hydrogens, removes all HETATM records, adds coordinate-bearing hydrogens and computes Gasteiger charges. It constructs atom/bond and residue graphs in RDKit/PDB order. The model's selected graph is `graph_mol`; residue graphs are generated but not selected by these configurations. The chain encoder has a fixed alphabet defined in this file, so unusual chain IDs can fail explicitly.

Native ligand `read_ref` parses a PDB, removes hydrogens, derives canonical **non-isomeric** SMILES for its atom-order permutation, and renumbers the original molecule using substructure matches; ligand construction then adds hydrogens and Gasteiger charges. This ordering is reproduced exactly, not replaced with a new featurizer. For input SDFs, native `to_pdb` serializes through RDKit. PDB cannot reliably preserve every chemical state: the bridge compares canonical **isomeric**, charge-sensitive identities before and after this conversion. A changed state becomes a visible failed pair with no score. Successful outputs retain the original SDF coordinates and chemistry, not a fabricated pose from the predictor.

Standalone input records retain the single-pose native crop behavior unless a `crop_reference` SDF is supplied. Discovery evaluation now always supplies the frozen experimental crystal ligand as that shared crop reference, matching the benchmark's conf0-defined pocket for all docked poses. Native protein/ligand preprocessing, model-specific 7/15 Å crops, weights, and output scaling are unchanged. The crop reference is hash-checked and its chemistry must survive native PDB conversion. All generated poses are scored separately, and the highest complete ensemble prediction selects the top pose. If any pose in a compound group fails, that group has no selected score, preventing selection bias from silently dropping failed poses.

## Featurization, tensors, model, and postprocessing

The definitive files are [Prediction/Datasets/MakeData.py](/home/sargis/dti/Prediction/Datasets/MakeData.py), [Prediction/Networks/Featurizer.py](/home/sargis/dti/Prediction/Networks/Featurizer.py), [Prediction/tmp_dataset.py](/home/sargis/dti/Prediction/tmp_dataset.py), each checkpoint's config/Networks, and [Prediction/prediction.py](/home/sargis/dti/Prediction/prediction.py).

Native construction produces NetworkX graphs and converts them to PyTorch Geometric Data objects. Graph nodes carry atom-number category, chirality, formal charge, valence, hybridization, Gasteiger charge, position and hydrogen flag; protein nodes also carry residue/chain metadata. The selected model configurations embed atom-number categories (64 features) and bond-type categories (32 features), use hidden node dimension 128 and hidden edge dimension 32, and remove hydrogens for encoding. Thus retaining a state in the input does not imply the learned predictor explicitly distinguishes every state feature.

`edge_index` is an integer `2 × E` tensor, categorical node/edge arrays are integer tensors, coordinates are `N × 3` floating tensors, and graph batching adds an integer `N` batch-membership tensor. Native encoders convert their working embeddings/coordinates to required floating dtypes. Dense auxiliary ligand–protein tensors use native `to_dense_batch` masks; the bridge does not insert custom padding or truncate sequence. Sparse radius operations cap each query at 500 neighbors in upstream code. That cap is an existing model behavior and can affect dense pockets. Runtime batches are deliberately one pair, evaluated sequentially while keeping all three models loaded.

The models use ligand and target graph encoders and `ComplexGraphAggregate`; configs include auxiliary pairwise-distance training. Inference calls original `StandartEvaluation` via `run_pipeline.score_pose`; auxiliary training tasks are not presented as affinity uncertainty. Native `.eval()` and `torch.no_grad()` disable training/dropout. No retraining or architecture modification is performed.

The raw regression head `y` is inverse-transformed by the actual native constants:

```
checkpoint_score = 6.3569711464280205 + 1.8564121996627636 * y
ensemble_score = arithmetic_mean(checkpoint_scores)
```

The application metric is **`BioCentric_affinity`, pK, higher is better**. Reference `run_pipeline.py` describes a pKd scale, but the training provenance is PDBbind and the local configuration/labels do not establish a homogeneous Kd-only endpoint. Therefore this is labeled a learned PDBbind-scale affinity estimate; it is not measured Kd, Ki, kcal/mol or a binding probability. No additional inverse transformation to molar Kd is made. `ScoreSTD` is the sample spread across the three different checkpoints. It is preserved as `ensemble_spread` with `spread_is_calibrated_uncertainty=false`; it is not a calibrated confidence interval or repeated-training uncertainty.

## Isolation, invocation, accounting, and exports

The engine ID is `biocentric_affinity`. [biocentric_model.py](/home/sargis/BioGPT/BioGPT/biocentric/biocentric_model.py) supplies the typed app plan and parser; [biocentric_predictor.py](/home/sargis/BioGPT/BioGPT/biocentric/biocentric_predictor.py) is a standalone, Python-3.9-compatible worker that imports the original predictor only inside its dedicated interpreter.

Configure:

```bash
export BIOCENTRIC_MODEL_DIR=/home/sargis/dti
export BIOCENTRIC_MODEL_PYTHON=/home/sargis/anaconda3/bin/python
```

This is an existing **Linux** scientific environment, separate from the copied macOS `.venv` and the new `.venv-linux` app environment. Observed runtime: Python 3.9.7, torch 2.1.0+cu118, torch-geometric 2.3.0, torch-scatter 2.1.2+pt21cu118, torch-cluster 1.6.3+pt21cu118, torch-sparse 0.6.18+pt21cu118, RDKit 2022.09.5, NumPy 1.26.4. This environment was verified in place, not rebuilt or silently changed. The original `dti/requirements.txt` is a historical dependency list with old Torch/CUDA pins and is not a tested installation lock for this observed runtime. See `runs/linux-biocentric-verification/runtime-versions.json` for additional observed dependencies.

[requirements-biocentric-observed-linux.txt](/home/sargis/BioGPT/BioGPT/requirements-biocentric-observed-linux.txt) records the relevant exact scientific/runtime pins, including CUDA-specific PyG extensions and their transitive utilities. It excludes unrelated conda applications and training-only dependencies. It is an **observed version record**, not a hash-locked clean installation that has been tested. No fresh model-environment rebuild from this file has been performed.

For a separate candidate rebuild, use Python 3.9, first install torch 2.1.0 from the official CUDA 11.8 index, then the recorded requirements. PyTorch documents the historical CUDA 11.8 binary in its [previous-version installation instructions](https://pytorch.org/get-started/previous-versions/#v210). The [PyG installation guide](https://pytorch-geometric.readthedocs.io/en/latest/install/installation.html) requires extension wheels matched to the PyTorch/CUDA combination; its [torch 2.1.0/CUDA 11.8 wheel listing](https://data.pyg.org/whl/torch-2.1.0%2Bcu118.html) includes the exact CPython 3.9 Linux extension builds recorded here. [RDKit's installation documentation](https://www.rdkit.org/docs/Install.html#cross-platform-using-pip) explains that older wheels use the historical `rdkit-pypi` distribution name; that is the actual name observed for version 2022.9.5.

```bash
python3.9 -m venv .venv-biocentric-rebuild
.venv-biocentric-rebuild/bin/python -m pip install torch==2.1.0 \
  --index-url https://download.pytorch.org/whl/cu118
.venv-biocentric-rebuild/bin/python -m pip install \
  -r requirements-biocentric-observed-linux.txt
.venv-biocentric-rebuild/bin/python -m pip check
.venv-linux/bin/python -m scripts.verify_biocentric_local \
  --model-dir /home/sargis/dti \
  --model-python .venv-biocentric-rebuild/bin/python \
  --out runs/biocentric-candidate-environment-verification
```

The final command is required to establish actual preprocessing, checkpoint and device compatibility for that candidate environment. Successful installation or import alone does not establish equivalence with the recorded working runtime; preserve the original Linux environment until the fresh candidate passes.

The parent application's fixed-argv runner strips credentials, supports cancellation by terminating the process group, and enforces a 1–3600-second timeout. Parameters bound CPU threads to 1–16 and pair count to 1–1000 (default 128); batch size stays one. The worker records peak RSS, peak allocated CUDA memory, elapsed time, framework/device versions, and per-pair time. No CUDA request silently falls back to CPU. CPU execution overrides only the historical training-device selector and checkpoint map location; it does not edit the source, weights or architecture.

Plans bind receptor, SDF, bridge code, all upstream Python/YAML source files, nine weights and interpreter hashes. Detection caches hashes by size/mtime/ctime for UI speed; execution rechecks bytes before and after calculation. A hash change invalidates execution. Inputs are immutable and output directories must be new/empty. The request is generated from typed fields and inspected source; no arbitrary executable or arguments come from the planner.

The worker writes real JSON progress events, `progress.json` with observed resource usage, immutable prepared pair inputs and per-pair JSON checkpoints. Each checkpoint binds the exact request hash and a SHA-256 of compact sorted JSON for its result row. Resume verifies both hashes and the protein/compound identity; a changed checkpoint fails without replacing that evidence. A rerun of the **identical generated request** may reuse completed matching checkpoints; it records `reused=true` and preserves scores. Matching failed checkpoints retain their original failure reason and no score; they are not silently retried with altered chemistry. At batch size one, an individual preprocessing/model/OOM failure is retained and later pairs can continue; GPU cleanup runs after the exception scope has released failed-forward traceback frames, collects unreachable objects, then calls `torch.cuda.empty_cache()`. There is no smaller batch than one, no silent CPU fallback, and no alteration of chemistry to recover an OOM. Systemic model-load failure fails the run. Cancellation leaves checkpoints for an explicit same-request resume. OOM cleanup has a synthetic object-lifetime test; a physical GPU OOM was not induced during verification.

`predictions.json` retains one status-bearing row for every exact protein/compound pair, original canonical isomeric SMILES, per-checkpoint scores, ensemble score, device and hashes. `predictions.csv` exports this accounting. `scored.sdf` contains the actual supplied structures, IDs, prediction status, pK/direction and checkpoint metadata; failed pairs have a reason and no score. The parser rejects missing/extra/duplicated/reordered/substituted rows, wrong model/checkpoint identities, NaN/Inf, booleans masquerading as scores, and inconsistent ensemble means. App execution is `partial` if some pairs fail and `failed` if all fail; failures are never replaced by zero.

## Training domain and validation limits

Application interruption recovery uses the registered, immutable approved workflow plan and the original output directory. Fresh preflight must reproduce the complete plan, request bytes and input/code/weight/runtime hashes before the nonempty directory is accepted. Prior logs, partial exports and checkpoints are archived as immutable project artifacts before retry. Completed checkpoints retain their exact bytes and modification times; only the final result view marks them reused. Missing pairs run normally. Recorded failures stay retained by default; a separate durable review can explicitly authorize retrying failed pairs after their original records are archived. Changed plans, checkpoints or reviewed artifacts block recovery.

The real application-path verification in `runs/linux-biocentric-research-resume-v1/verification.json` interrupted an RTX 4090 run after 2 of 16 actual previously validated Linux Vina structures. It resumed the same request and output directory, retained both checkpoint hashes and modification times, computed the remaining 14 pairs, and finished 16/16 with zero failures or progress warnings. Seventeen prior output artifacts were archived and verified. This is execution/recovery evidence, not affinity accuracy. `scripts/verify_biocentric_research_resume.py` reproduces this bounded test in a fresh directory; synthetic orchestration tests separately cover opt-in failed-pair retries and provenance tampering.

The local [OpenBind results README](/home/sargis/dti/results/OpenBind_EVA71/README.md) reports training on PDBbind general v2020 excluding CASF-2013/2016 core sets (18,934 complexes). It describes substantial protein overlap even after those core exclusions and reports different checkpoint behavior on affinity ranking, pose selection and binder/nonbinder screening, including negative screening results. Those are **historical local claims and tables**, not new calculations by this integration. They were not independently reproduced here. An exact checkpoint-specific training-membership manifest is not established by the integration. ABL/1IEP pretrained overlap remains unknown and cannot be assumed absent. Training configuration includes MSE affinity regression and docking/auxiliary losses; it does not establish DTI-probability calibration.

Real new checks:

| Check | Evidence / outcome |
| --- | --- |
| Device computation | RTX 4090 torch CUDA arithmetic and torch-cluster radius kernel both executed successfully. |
| Native reference | Original `Prediction/run_pipeline.py`, original weights, packaged Fxagtc000101A / GTC000064A input. |
| App integration + numerical equivalence | `runs/linux-biocentric-reference-v1/verification.json`: native 8.1943181356, app 8.1943179766; largest checkpoint difference 9.54e-7, tolerance 1e-5. |
| Durable resume | Same verifier: one reused pair, zero new model predictions, identical stored score. |
| Resource measurement | App reference inference 2.44 s worker time, 897,460 KiB peak RSS, 79,432,192 bytes peak allocated CUDA memory. |
| ABL prepared-pose smoke | `runs/linux-biocentric-verification/abl-gpu`: actual protonated imatinib/1IEP prediction 6.7495233218; all three checkpoint scores preserved; 2.68 s and 100,292,096 bytes peak CUDA allocation. |
| CPU computation | `runs/linux-biocentric-verification/bridge-cpu`: 8.1977421443; differs from CUDA by about 0.00342. CPU and CUDA must remain separately fingerprinted protocols; no cross-device numerical equivalence claim. |
| 128 new Linux Vina poses | `runs/ubuntu-biocentric-diverse/summary.json`: 128/128 GPU predictions, zero inference/preprocessing failures, four chunks of 32. All 128 checkpoints reused unchanged with zero new predictions. 121 poses remain shortlist eligible; seven physical failures remain visible, scored but excluded. |
| Bounded workload resources | 18.545 s summed prediction-worker time; 32.807 s including restart/resume checks; peak RSS 925,556 KiB and peak allocated GPU memory 87,532,032 bytes. |
| Updated integrity/recovery bridge | `runs/linux-biocentric-reference-v2/verification.json`: native/app largest checkpoint difference 9.54e-7; one row with a verified content checksum reused unchanged, zero new predictions. Successful inference/preprocessing is unchanged from the 128-pair run. |

These are software/inference checks. No matched measured ABL labels, held-out accuracy estimate, calibrated uncertainty or evidence of a new binder was produced by them. Copied poses are software inputs, not new docking calculations. Docking scores and this model's learned pK estimates must remain separate columns/rankings.

Reproduce the real native/adapter/resume check into a new output directory:

```bash
.venv-linux/bin/python -m scripts.verify_biocentric_local \
  --model-dir /home/sargis/dti \
  --model-python /home/sargis/anaconda3/bin/python \
  --out runs/linux-biocentric-reference-new
```

This command performs genuine local GPU inference and spends no API credits. In a sandbox that hides NVIDIA devices, GPU commands need the host's execution approval; `nvidia-smi` failure inside that sandbox alone is not evidence of a broken driver.

The bounded same-pose software comparison has its own reproducible runner. It freezes pose hashes, source Vina results, all physical-validity statuses, model/code identity and shortlist criteria before predictions, then verifies all checkpoint reuse:

```bash
.venv-linux/bin/python -m scripts.score_biocentric_screening \
  --docking-run runs/ubuntu-diverse-docking-verification \
  --model-dir /home/sargis/dti \
  --model-python /home/sargis/anaconda3/bin/python \
  --chunk-size 32 --out runs/ubuntu-biocentric-diverse-new
```

Its `comparison.json`, `comparison.csv` and `comparison.sdf` preserve separate Vina and BioCentric scores. This unlabelled NCI set cannot select a scientifically superior affinity predictor or estimate target-specific predictive accuracy.


## Bio Orchestrator discovery integration

`BioCentric affinity` is available in the Scoring multiselect and participates in pocket × docking × scoring combinations. Existing default selections are unchanged. The settings expose CUDA/CPU, CPU threads, compounds per model batch (default 16, maximum 32), and batch timeout (default 300 seconds). The CUDA option requires a detected local GPU; it never falls back silently. The selected device and settings are frozen into each reviewed chain.

`discovery_affinity.py` flattens each chunk of docked poses into uniquely named records and loads the three-model ensemble once per chunk. Actual GPU forward batch size remains one. Docking poses are reused across scoring tools for the same immutable pocket, docking parameters, and chemistry. The crystal pose is included in the first scoring chunk for a same-model reference comparison. Original docked coordinates are preserved; the best scored pose then passes through the existing PoseBusters and symmetry-aware receptor-frame RMSD checks. Only chemically matching crystal ligands have reference comparisons.

Native predictions, pose manifests, per-model values, uncalibrated ensemble spread, device, model/source/checkpoint hashes, failures, and input/output artifacts are retained under each evaluation run. The approved runtime fingerprint includes the local model source, weights, interpreter, and both adapters. Batch inference time is allocated equally among the chunk's measurement/control rows for chain runtime accounting; original native timing and resource usage remain in the batch result. Progress still counts completed tool combinations over the whole evaluation set, never poses or compounds. Cancellation remains available during docking and inference; missing records stay explicitly unrun.

The output is learned pK, higher is better, not an experimental measurement. Experimental endpoints remain separate, and ranking uses the existing within-assay metrics. Model training overlap is unknown; selecting a chain on these data does not establish held-out predictive accuracy.

Local integration artifacts are under `runs/biocentric-discovery-integration`: `verification.json` records real CUDA scoring of two supplied ABL docked pose records and the crystal reference, while `chain-verification.json` records a bounded Vina → BioCentric → pose-validation chain on three existing ABL measurements and a crystal control. These are software checks only, not predictive validation.


### Crystal-pocket correction (2026-09-27)

New experiments default to the crystal-ligand box (4 Å padding, minimum 24 Å per side). The minimum avoids using a tiny fragment-sized grid for larger ligands. This docking grid is distinct from BioCentric's learned 7/15 Å graph crops. Saved runs retain their original inputs and parameters.

The investigated run `8d9390d9374a4ef6adb8ddfe11fc4461` used 1O41/300 and fpocket rank 1, not 1A1B. Its 35 measurements mixed 14 SH2 binding records and 21 kinase inhibition records. SRC SH2 structures now require explicitly SH2-described assays; explicit PDB input is resolved via RCSB and cannot silently select another PDB. Peptide/polymer binders are not replaced by isolated modified residues or another ligand.

`runs/crystal-pocket-correction/verification.json` compares original native GPU inference with the adapter on identical crystal/docked inputs and one fixed crystal crop. Ensemble differences are below 1e-6. This is numerical parity, not evidence that docking reproduces affinity measurements. The initial bounded Fast Vina control had RMSD 3.45 Å: better site placement than the archived 24.10 Å, but it still fails the 2 Å pose-recovery gate.

The HDD search did not locate archived 1A1B model predictions. The local `General_Labels.csv` contains its experimental pKd 6.4 (0.4 µM); this label is not a reproduced model prediction. A matched historical-score comparison remains unavailable until those exact inputs/outputs are located.


New experiments include BioCentric by default and use Balanced docking. Vina uses exhaustiveness 8, up to 9 poses, four CPU threads and its automatic search budget. Fast remains selectable. On the same crystal control, the 24 Å box with Fast yielded RMSD 12.20 Å; Balanced yielded 3.81 Å and passed physical checks, but still failed the 2 Å recovery criterion. Neither result is represented as successful crystal recovery. Detailed artifacts are in `default-24A/control.json` and `balanced-24A/control.json` beneath the correction run folder.


### TYK2 ligand failure and score parity (2026-09-27)

The latest TYK2/5F20 evaluation rejected CHEMBL514409 because PDB conversion omitted its E/Z stereo metadata. The adapter now accepts missing metadata only if stereochemistry independently assigned from the existing 3D coordinates exactly matches the input isomeric identity. Actual stereo, charge or connectivity mismatches still fail. The validation copy does not alter native model inputs.

The copied experiment `3bfd3f1dba9a4d44891636efe03b75aa` (run `7940c17f2ceb436ba940696631f3fac6`) completed all eight measured ligands and the crystal control for all three chains. Direct original-script CUDA inference on all nine selected BioCentric poses matched adapter ensemble scores within 9.54e-7 and individual model outputs within 2.39e-6. Artifacts: `runs/tyk2-biocentric-repair/native-parity/verification.json`.

Within the one rankable assay (4 of 8 measurements; other assay sizes 1, 2 and 1), GNINA CNN rho=1.0, BioCentric rho=0.4, and Vina rho=0.0. The GNINA correlation is a valid four-point rank result (exact two-sided permutation p=0.0833), not a full-eight-measurement result or predictive-validation claim. Native parity confirms no scoring discrepancy for these poses; it does not imply BioCentric must outrank other models.

Failed ligands remain explicit and no longer prevent shared-success-set comparison or final selection. Failed chemical identity checks are not bypassed to obtain more scores. See `bayesian-pipeline-search.md` for the shared-subset selection policy.
