# Scientific contracts, scope and verification

`biocentric.science` implements local scientific tools with explicit preprocessing and postprocessing. The service validates project-scoped paths and persists scientific decisions; these adapters never invoke a shell, accept arbitrary command text, or contact external molecular services. Dependency availability is separate from successful execution and from scientific validity.

## Input and output contracts

| Operation | Accepted input | Preprocessing | Returned output / postprocessing |
|---|---|---|---|
| `prepare_library` | SDF/SD; whitespace-separated SMILES and optional name (`.smi`, `.smiles`, `.txt`) | Content-addressed untouched snapshot; sanitize; explicit stereochemistry/domain checks; exact isomeric identity deduplication; explicit H; deterministic ETKDGv3 and MMFF94 only when 3D is absent | Canonical `prepared.sdf` with original properties and compound/state IDs; descriptors; every excluded record with reason; preparation JSON and hashes |
| `inspect_structure` | PDB, mmCIF/CIF, PDBQT | Read records only; no repair, chain, water or cofactor deletion | Chains, heterogens, alternate locations, missing-residue remarks where available; unresolved mapping/quality questions |
| `dock_library` | Curated receptor PDBQT; canonical explicit-H 3D SDF; exactly one bound 3D reference SDF; optional matching receptor PDB | Box from reference heavy-atom extents plus explicit padding; frame checks; Meeko Gasteiger charges and Vina types | Vina raw PDBQT/log; Meeko SMILES/index-map reconstruction to per-pose SDF; canonical identity verification; raw scores with units/version; PoseBusters reports; PLIP reports if installed; resumable compound checkpoints |
| `prepare_receptor` | Reviewed standard PDB; explicit method, pH and heavy-atom repair decision | PDBFixer, PDB2PQR/PROPKA or Meeko through fixed arguments | Auditable prepared copy and raw log; a completed preparation still requires scientific inspection |

`prepare_library(input_path, out_dir, seed=42, emit=None, cancelled=None)` returns `compounds`, `exclusions`, `prepared_sdf`, `artifacts` and metadata. `dock_library(receptor_path, prepared_sdf, reference_sdf, out_dir, config, emit=None, cancelled=None)` returns `poses`, `failures`, `artifacts`, `metadata`. Each event is a dictionary containing `kind`, `tool`, `stage`, `message`, and optionally `compound_id`, `progress` and `data`. Events expose actual adapter calls and intermediate values, not hidden chain-of-thought.

`DockConfig` accepts only seed, CPU count, exhaustiveness, number of poses, timeout, box padding, optional receptor PDB, require-PoseBusters and run-PLIP flags. Bounds cap CPU at 16, exhaustiveness at 128, poses at 20, file size at 100 MiB and library size at 1,000. Scientific child processes have wall-time limits and cancellation terminates their process group. Generated structures retain their canonical parent-state IDs. A configuration, input hash or relevant tool version change cannot reuse an existing checkpoint directory.

## Scientific boundaries

- No unrequested salt stripping, neutralization, tautomer/stereoisomer expansion, water deletion, loop reconstruction, automatic nonstandard-residue conversion or mutation. Unassigned stereocenters/double-bond stereo, unsupported elements, radicals, large/flexible compounds, macrocycles and peptide-like compounds are excluded with reasons. Existing 3D coordinates are preserved.
- Added ligand hydrogens follow the input formal charge and valence. This is **not** a physiological pH assignment. pH is unknown unless a reviewed protein preparation method explicitly receives one. Protonation/tautomer ensembles require a separate planned validation step.
- PDBQT is not canonical chemical identity. Meeko writes SMILES and atom maps and recovers SDF bond orders; the recovered canonical isomeric SMILES must exactly match the input. SDF metadata survives the round trip.
- Descriptor values (molecular weight, cLogP, TPSA, donors, acceptors, rotatable bonds, formal charge) and PAINS review flags are not ADMET predictions. Missing ADMET values are null. Alerts do not automatically reject compounds.
- Receptor metal atoms block routine docking. Structure/mutation/construct/pocket confidence, cofactor decisions and biological assembly mapping need review. No automatic AlphaFold prediction or model repair is claimed.
- PoseBusters runs `dock` when a matching receptor PDB is provided and `mol` otherwise. Its report explicitly states whether receptor checks ran. A missing tool is unavailable; a failed check is preserved. Ligand identity is independently checked against input SDF. PLIP geometry describes predicted contacts, not affinity.
- Cognate redocking RMSD uses RDKit symmetry-aware `CalcRMS` in the shared receptor frame. Ligand-only aligned `GetBestRMS` is separately labeled shape similarity and never used to hide translation into the wrong pocket. Cross-receptor alignment is not implemented and must be supplied as reviewed aligned inputs.

## Optional receptor preparation paths

Use `prepare_receptor(..., method="pdbfixer", reviewed=True, ph=7.4, add_missing_heavy_atoms=False)` only after an explicit project decision. PDBFixer records missing residues and atoms; missing loops are disabled, nonstandard residue replacement is disabled, and water/chain deletion is disabled. Missing heavy atoms block unless explicitly enabled. Its atom-building seed is 42; pH/protonation remains a modeling assumption. Input PDB files with alternate locations or metals block the ordinary path.

`method="pdb2pqr"` runs AMBER parameters, preserves chain IDs, invokes PROPKA titration-state assignment at the requested pH and writes PQR plus a PDB copy and logs. **PQR must not be sent directly to Vina.** Inspect protonation/repair decisions, then use `method="meeko"` on an approved PDB with explicit hydrogens to obtain receptor PDBQT and the Meeko JSON record. Unknown residues/charges fail rather than being silently skipped.

The default application vertical slice uses curated receptor PDBQT. These preparation functions are library APIs for explicit review; they are not automatically executed during screening. Automatic preparation accepts only standard amino-acid residues. Waters, cofactors and unknown residues remain in the original and block this path rather than being deleted or triggering an implicit external CCD template lookup; use separately reviewed curated inputs for those cases. The tutorial already supplies a prepared receptor; the application does not claim those preparation decisions were experimentally validated.

## Installation

Install `requirements-science.txt` after the main application dependencies. Vina's Python wheels are platform-dependent. A verified official **AutoDock Vina 1.2.7** executable may be placed at `.tools/vina`, found on PATH, or configured through administrator environment variable `BIOCENTRIC_VINA`. A Python Vina installation is the fallback and runs inside a subprocess. A missing executable/library is an explicit error.

For this Apple Silicon environment, the official [mac_aarch64 v1.2.7 release executable](https://github.com/ccsb-scripps/AutoDock-Vina/releases/tag/v1.2.7) was downloaded into `.tools/vina`; its SHA-256 and source URL are in `examples/1iep/sources.json`. Nothing was globally installed. On another supported platform, `python examples/install_vina.py` downloads the matching pinned official release (Mac ARM asset additionally checked against its previously verified hash) and records the source/version/hash. It will not overwrite an existing executable.

Optional adapters have extra dependencies. PLIP requires OpenBabel; on supported machines `openbabel-wheel` supplies its Python import. PLIP's source build may require `--no-build-isolation` so its build hook can see that installed module. Use the versions recorded in `requirements-science-optional.txt` and review the upstream license terms (PLIP/OpenBabel GPL; other tools have separate licenses). See the run records below for what was actually exercised.

## Live calculation and tests

Run `.venv/bin/python examples/run_live.py --output runs/live-smoke-new` to execute the complete scientific slice. `examples/1iep/README.md` explains source/licensing and controls. It is a public **mouse** ABL/P00520 example (PDB 1IEP); substituting human ABL1 without mapping is incorrect. The verified local runs produced three compounds/eight poses. Top imatinib Vina score was -13.211 kcal/mol; receptor-frame symmetry-aware RMSD was 0.380 Å and PoseBusters passed. These numbers verify an actual calculation on this machine, not biological screening predictivity. Initial report (PLIP unavailable at that time): `runs/live-smoke/report.md`. Full rerun with PLIP executed and chemistry-review flags: `runs/live-verified/report.md`; machine evidence: `runs/live-verified/manifest.json`.

Scientific tests cover stereochemical identity, enantiomer separation, ambiguous/invalid/unsupported compound exclusion, original SDF property/coordinate retention, derived boxes, timeout/cancellation, correct-frame RMSD, held-out/scaffold leakage and an honestly failed enrichment result. Run `.venv/bin/python -m pytest tests/test_science.py -q`.

## Independent evaluation harness

`biocentric.benchmark` and `examples/evaluate_benchmark.py` are separate from the docking workflow. Dataset records require a compound ID, tuning/evaluation partition, scaffold/analog group, measured active or tested-inactive label, endpoint definition and source identifier. Artificial decoys and untested compounds cannot be passed as tested inactive. Partition IDs and scaffold groups must be disjoint. Dataset/version SHA-256 and pass criteria are frozen with exclusive file creation **before** score evaluation; scores must carry their metric, job ID and artifact SHA-256. Every evaluation sample requires a score, preventing selective failure exclusion.

The report includes class balance, top-k precision, enrichment, Wilson 95% intervals, sample-size checks, failed criteria, tie handling and unknown training overlap. Intervals state their independent-binomial assumptions and condition enrichment uncertainty on observed prevalence; they do not imply a large validation dataset. No affinity correlation is calculated from mixed or censored endpoints. The example has no independent active/inactive evaluation and therefore fails the screening-predictivity gate by missing controls. Protocol hashes can reveal file changes, but chronological pre-registration and source-label truth remain the scientist's responsibility.

## Optional adapter execution record

PDBFixer 1.12.0/OpenMM 8.3.1, PDB2PQR 3.7.1/PROPKA 3.5.1 and Meeko receptor preparation were successfully executed on PDB 1CRN under an explicit pH 7.4 test assumption. See `examples/receptor_prep/README.md` and raw `runs/adapter-checks/1crn-*` outputs. The interrupted SER438 residue in the curated tutorial 1IEP PDB triggered honest PDB2PQR/Meeko failures; no automatic reorder/repair was attempted.

PLIP 2.4.0/OpenBabel 3.1.1 was also executed on actual docked poses. Its reported ligand SMILES can disagree with the canonical SDF after PDB bond/hydrogen perception. Every report includes a chemical-identity check and such disagreement returns `completed_requires_review`, with raw XML/text retained and contact types not treated as validated chemistry. Explicit hydrogens are retained when supplied; unique atom serials and ligand CONECT records are written to the temporary complex. This remains an inspectable limitation of this adapter path.

## Official implementation references

- [RDKit getting started and SDF suppliers](https://www.rdkit.org/docs/GettingStartedInPython.html), [RDKit stereochemistry](https://www.rdkit.org/docs/RDKit_Book.html)
- [Meeko preparation and mapped export](https://meeko.readthedocs.io/en/develop/), [Meeko export implementation](https://github.com/forlilab/Meeko/blob/develop/docs/source/export_usage.rst)
- [Vina Python API](https://autodock-vina.readthedocs.io/en/latest/docking_python.html), [Vina official tutorial data](https://github.com/ccsb-scripps/AutoDock-Vina/tree/v1.2.7/example/basic_docking)
- [PoseBusters API and validation modes](https://posebusters.readthedocs.io/en/latest/api_notebook.html)
- [PLIP local structure and report CLI](https://github.com/pharmai/plip/blob/master/DOCUMENTATION.md)
- [PDBFixer manual](https://github.com/openmm/pdbfixer/blob/master/Manual.html), [PDB2PQR/PROPKA usage](https://pdb2pqr.readthedocs.io/en/latest/using/index.html)
