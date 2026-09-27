# Diverse local validation — 26 September 2026

The larger test completed **144 actual Vina calculations**: 128 different compounds on the baseline protocol, then eight compounds with each of two parameter changes. All calculations returned a pose. Physical checks passed for **121 of 128 baseline poses** and all 16 comparison poses. The seven baseline failures remain recorded and are excluded from the shortlist.

This tests software operation, chemical identity, physical checks, progress, checkpoint reuse and exports. It is not an independent target-activity benchmark, a validation of binding-affinity accuracy, or a comparison of different docking engines.

## Inputs and diversity

The input source is RDKit's bundled [NCI first-5K SMILES sample](https://github.com/rdkit/rdkit/blob/master/Data/NCI/first_5K.smi). The installed file contains 4,999 records; its hash and RDKit version are pinned in the [complete report](../runs/diverse-docking-verification/summary.json). After supported-chemistry and workload bounds, 2,536 records across 492 distinct scaffold groups were eligible. The selection takes one representative per nonchiral Bemis–Murcko scaffold, then uses seeded Morgan-fingerprint MaxMin selection, following RDKit's [diversity-picking interface](https://rdkit.org/docs/GettingStartedInPython.html#picking-diverse-molecules-using-fingerprints).

The selected 128 compounds represent 128 scaffold groups. Their molecular weights range from 112.1 to 461.5, with 0–8 rotatable bonds. The set contains 96 compounds with aromatic atoms, 62 with nonaromatic ring atoms, and one acyclic compound; these categories overlap. It includes nitrogen-, oxygen-, sulfur-, fluorine-, chlorine- and bromine-containing compounds. All acyclic molecules share one selection group, so this is deliberately limited chemical coverage rather than a universal chemistry test.

The [selection ledger](../runs/diverse-docking-verification/selection.json) preserves every original row, inclusion/exclusion reason, molecular identity and selection order. All 128 selected compounds prepared successfully. A separate seven-record adversarial library accepted two explicitly assigned stereoisomers and rejected a duplicate, ambiguous stereochemistry, a disconnected salt, unsupported metal and invalid SMILES with the expected individual codes.

## Actual calculations

The receptor is the reviewed 1IEP mouse ABL example; the bound imatinib coordinates define the search box. RDKit preparation, Meeko conversion, AutoDock Vina 1.2.7, canonical SDF recovery and PoseBusters ran locally. Input charge and tautomer are retained; no pH-state enumeration, PLIP or ADMET calculations were performed.

| Protocol | Compounds / poses | Physical checks passed | Calculation failures | Time |
|---|---:|---:|---:|---:|
| Baseline: seed 42, exhaustiveness 1 | 128 / 128 | 121 | 0 | 440.6 s |
| Seed 43, exhaustiveness 1 | 8 / 8 | 8 | 0 | 23.6 s |
| Seed 42, exhaustiveness 4 | 8 / 8 | 8 | 0 | 39.8 s |

All protocols use one Vina CPU, one requested pose and a 60-second per-compound timeout. The two variants reuse the same eight prepared coordinate sets. Their fingerprints and results remain separate; scores are not pooled across protocols. [Parameter comparisons](../runs/diverse-docking-verification/parameter_comparisons.json) record observed changes without declaring an accuracy winner. Low exhaustiveness and one pose prioritize broad software coverage over deep pose search.

The baseline scores range from −10.811 to +19.789 kcal/mol. A score does not override failed physical checks. There are no matched measured activity labels in this run, so neither a good score nor shortlist membership establishes binding.

## Reuse, files and resources

All 128 completed baseline checkpoints were reused in 0.43 seconds with **zero new Vina calls** and unchanged pose hashes and modification times. CSV/SDF exports were parsed and checked against molecular identities and source hashes. The baseline exports contain 128 total poses and 121 physically passing poses; invalid poses remain in the all-pose export with explicit labels.

- Total wall time: **513.8 seconds (8 min 34 sec)**, within a 20-minute cap.
- Parent Python CPU: 17.3 seconds; child-process CPU: 584.2 seconds. These are cumulative process times, not wall time.
- Observed parent peak resident memory: 172 MB; maximum child resident memory: 701 MB. These are separate measurements, not a measured combined peak.
- Source evidence on disk: 10.9 MB, before the separate app artifact copy.
- **No GPU, network requests or language-model calls.**

The app import verifies the complete report inventory, source-library/preparation binding, baseline protocol, per-compound accounting and pose hashes before creating a separate saved run. It preserves source logs and manifests, ranks only the baseline poses and creates a ten-pose shortlist from physically valid candidates using the existing diversity rule. The saved run is labeled **Local validation · no model calls**. The original GPT research conversation is unchanged.

Open the [128-ligand saved run](http://127.0.0.1:8502/?project=726fbd0fea4e4b7496088bc3e458fd5b) to inspect the table and download the ten shortlisted poses, all 128 baseline poses or ranked scores. Reimporting the same verified source bundle returns this existing project. The [app verification record](../runs/interactive-workspace-verification/verification.json) confirms counts, hash-verified exports and the unchanged completed original research project.

## Other engines and reproduction

The separate [engine contract matrix](engine-validation-matrix.md) adds 138 synthetic cases covering all eight requested methods. Together with 19 existing adapter tests, 157 focused tests passed. These tests found and corrected malformed-score handling, Amber units/uncertainty checks, PLANET molecular-identity matching and Uni-Dock filename binding. Optional-engine processes were mocked: AutoDock4, GNINA, Uni-Dock, MM/PBSA, MM/GBSA, GNINA CNN affinity and PLANET still require their actual runtimes and scientific inputs before live validation.

```sh
.venv/bin/python -m scripts.validate_diverse_docking
.venv/bin/python -m pytest -q tests/test_diverse_docking.py tests/test_engine_matrix.py tests/test_engines.py
```

The runner verifies and reuses existing completed evidence. To perform a new calculation, supply a new `--output` directory; it will not silently replace changed pinned inputs. The [CSV](../runs/diverse-docking-verification/exports/baseline/ligand_results.csv), [all-pose SDF](../runs/diverse-docking-verification/exports/baseline/all_docked_poses.sdf), [physically passing SDF](../runs/diverse-docking-verification/exports/baseline/physically_valid_poses.sdf) and [complete evidence manifest](../runs/diverse-docking-verification/summary.json) are available for inspection.
