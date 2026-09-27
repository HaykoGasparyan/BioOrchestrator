# Engine contract validation matrix

This matrix describes deterministic software checks in `tests/test_engine_matrix.py`. These checks use synthetic molecular/native-format fixtures and mocked process boundaries. They do **not** run AutoDock4, GNINA, Uni-Dock, AmberTools or PLANET, measure binding, compare engine accuracy, or establish a scientifically validated screening method. Vina is the installed, previously live-tested local engine; optional adapters remain unavailable until configured and independently exercised.

The matrix complements `tests/test_engines.py`, which covers basic plans, installation blockers, input hashes, the registered native-manifest workflow, frozen input/protocol binding, chemical recovery and explicit exclusions. It adds varied valid and malformed formats so an apparently successful process cannot silently turn incomplete or mislabeled outputs into a successful result.

| Adapter | Additional valid fixtures | Rejected or guarded cases | Meaning preserved |
|---|---|---|---|
| Vina | Two poses; signed decimal and scientific-notation values; per-file SHA provenance | Numeric prefixes such as `-7bad`; malformed/nonfinite second pose; missing or invalid RMSD fields | Empirical `vina_score`, kcal/mol, lower better; no measured affinity claim |
| AutoDock4 | Repeated DLG docking runs; scientific notation; positive energy; bracketed native energy explanation; unrelated summary ignored | Malformed or wrong-unit docking run; missing/extra atom map; unsupported torsion count; missing grid metadata | `autodock4_score`, kcal/mol; not pooled with Vina scores |
| GNINA docking | Multiple SDF pose records with conflicting energy and CNN orderings; CNN pose-score endpoints 0 and 1 | Any missing score field; nonfinite/non-numeric values; CNN pose score outside [0,1]; invalid later SDF record; unsupported local GPU request | Empirical energy, CNN pose quality and CNN predicted affinity remain three separate metrics |
| Uni-Dock | Mocked two-ligand batch, two poses each; separate source files and hashes | Incomplete batch; correct file count but substituted filename; malformed pose score; unsupported Mac GPU execution | `unidock_vina_score`, kcal/mol; expected `<input-stem>_out.pdbqt` binding; chemical mapping still belongs to native workflow validation |
| MM/PBSA | PB difference mean/SD/SEM; scientific notation; both supported native total-label variants | Wrong/missing units; GB-only file; duplicate requested method or total; truncated summary; nonfinite or negative uncertainty; invalid frame count | PB complex-minus-receptor-minus-ligand endpoint estimate without entropy; uncertainties remain explicit |
| MM/GBSA | GB difference mean/SD/SEM; scientific notation; both supported native total-label variants | Wrong/missing units; PB-only file; duplicate requested method or total; truncated summary; nonfinite or negative uncertainty; invalid frame count | GB endpoint estimate without entropy, distinct from PB and measured affinity |
| GNINA CNN affinity | Same strict multi-record score parser; explicit `--score_only` and CPU flag | Docking-only parameters are excluded by its typed request; missing/invalid metrics; unsupported GPU request | Predicted pK, higher better, kept separate from empirical energy and pose score; no docking performed by this mode |
| PLANET | Native CSV with optional index and quoted name; multiple molecules; equivalent SMILES and reordered output; matching 2D output SDF | Invalid/blank SMILES; malformed CSV row; nonfinite score; skipped/duplicated/substituted molecules; changed stereochemistry or charge; missing/mismatched native SDF; invalid query record; 2D pocket reference; protein lacking alpha carbons | Learned pK from protein pocket and ligand graph; returned SDF is not claimed to be a generated docking pose |

Every adapter additionally receives a complete typed input fixture and is checked for a reviewable bounded plan, immutable input hashes, unavailable-tool blockers, missing preparation/MD review, empty/wrong-format scientific files, arbitrary extra arguments and invalid timeout limits. Plans create no output directory and execute no program. The MD topology and trajectory files are deliberately minimal contract fixtures; their review flags are synthetic test inputs, not evidence of equilibrated simulations.

## Regression corrections

The new fixtures exposed and corrected these failure modes in `biocentric/engines.py`:

- Vina-family and AutoDock4 parsing previously skipped malformed score records if another valid score remained. A malformed marked record now fails the complete output rather than silently dropping a pose/run.
- Amber parsing previously accepted a numeric prefix or split a two-number line into three fields, ignored declared units and failed to check SD/SEM finiteness. It now requires explicit kcal/mol or kcal/mole, an unambiguous requested-method difference summary, and complete finite mean/SD/SEM with nonnegative uncertainties.
- PLANET parsing previously accepted missing or invalid molecular identities. The direct runner checked counts only, so same-count substitutions could appear complete before downstream workflow checks. Native CSV and output SDF now both match the complete input identity multiset, preserving charge and stereochemistry.
- Uni-Dock previously checked output-file count only. It now requires exactly the expected input-derived PDBQT filenames; downstream molecular-identity checks remain necessary and unchanged.

These changes reject incomplete evidence; they do not invent missing scores, convert energy units, impute labels or select a different chemistry automatically. Unsupported native format variants remain explicit errors until an adapter for that variant is reviewed.

## Running and interpreting the checks

```sh
.venv/bin/python -m pytest -q tests/test_engine_matrix.py tests/test_engines.py
```

The matrix has 138 parametrized cases; together with the 19 existing engine tests, the focused run passes 157 tests. Execution is mocked throughout these files, including tests whose fake result has `executed: true`: that field verifies the process-boundary contract, not a real optional-engine calculation. An autouse matrix fixture rejects any subprocess call unless that specific test replaces it with a fake runner.

Actual runtime verification is recorded separately under the Vina run artifacts described in [engine documentation](engines.md). Each optional engine still needs installation/runtime validation, a real native reference example, chemical and physical checks, and an appropriate held-out measured benchmark before claims about predictive performance. A one-pose redocking control checks software and pose recovery only; this matrix establishes neither affinity accuracy nor enrichment.

## Native-format references

The existing [engine documentation and primary references](engines.md#primary-references-consulted) describe the supported versions and file contracts. Uni-Dock filename matching follows upstream [`default_output`](https://github.com/dptech-corp/Uni-Dock/blob/main/unidock/src/lib/utils.cpp) and the GPU-batch use of that function in [`main.cpp`](https://github.com/dptech-corp/Uni-Dock/blob/main/unidock/src/main/main.cpp). The fixtures are locally generated examples of those contracts, not copied benchmark datasets or evidence that optional runtimes are installed.
