# Frozen comparative protocols and independent evaluation

The research conversation now freezes a bounded set of **operational docking protocols** before its first pilot. It includes the literature proposal's available Vina, GNINA and Uni-Dock candidates; the Vina pilot compares exhaustiveness 4 and 8. Other unsupported or unselected methods remain explicit exclusions. Native GNINA/Uni-Dock inputs are prepared from registered files and still pass native preflight and a concrete review. Candidate records include exact input hashes, parameters, executable hash, dependency versions, source IDs and unresolved training overlap. Separate protocol directories prevent checkpoint reuse between candidates.

Stage 06 records every calculation or failure in an immutable artifact. Stage 07 checks each protocol's top-ranked symmetry-aware RMSD **in the receptor frame**, physical validity, and complete reference accounting. Candidates that fail any criterion are ineligible. The predeclared selection rule chooses the passing candidate with lowest mean tuning RMSD, breaking ties by protocol ID. This is selection for the supplied reference pose; it does not establish that one engine is predictively best. A single remaining candidate is identified explicitly. Execution failures retain reasons and never receive zero scores.

The selected protocol is locked before the workflow records the held-out evaluation status or requests the advisor's interpretation. Production verifies that its engine and parameters exactly match that lock. The advisor cannot rewrite failed checks or select a different protocol. Per-protocol artifacts and same-job question continuation avoid repeating completed pilots. A changed executable, parameter set, test set or final prediction table requires a new comparison.

**When no suitable measured protocol exists**, the conversation records `heldout_evaluation.status = not_evaluated`, `predictive_validation = not_established`, and `selection_scope = software_pose_only`. Even a project with a measured-record preview needs a separately reviewed affinity protocol; the application does not invent predictive cutoffs. The 1IEP reference and unlabeled NCI sample support software and pose tests only. No binding, affinity-accuracy or screening-enrichment claim follows from these results.

## Endpoint-matched affinity evaluation API

`biocentric.comparative` also provides a real, deterministic evaluation harness for an explicitly supplied measured cohort:

- `freeze_comparison(testset, candidates, criteria, destination)` exclusively creates the candidate/protocol record before predictions.
- `select_and_lock(testset, tuning_results, frozen, destination)` accepts every tuning pair from every candidate, never held-out rows.
- `evaluate_locked(testset, heldout_results, frozen, lock, destination)` evaluates only the selected protocol and exclusively writes one final result.
- `run_comparison(testset, candidates, criteria, output_dir, runner, emit=..., cancelled=...)` orchestrates these calls. Its real adapter callback is invoked on tuning data for each candidate, then on held-out data **only after** the lock exists. Saved prediction tables are reused on resume, with hash verification and no new runner calls.

This API is available for integration with a scientifically compatible predictor; it does not imply that arbitrary installed models satisfy its contract. The ordinary interactive workflow currently executes pose comparisons and reports predictive evaluation as unavailable. A measured-data affinity run must supply a real adapter callback and a reviewed endpoint-specific protocol through the API; there is no automatic conversion of BioCentric/GNINA/PLANET learned scores to experimental pKd or pKi.

A candidate declares an engine ID, parameters, input/code/weight provenance, output metric, units, direction and known training overlap. The regression path requires identical target, species, construct, assay, endpoint and recorded conditions, exact uncensored Kd or Ki measurements, separate tuning/evaluation scaffolds and exact canonical ligand identities. Missing target or construct information blocks regression. IC50/EC50, binding probabilities, unspecified pK scores and energies cannot be relabeled as pKd or pKi. Known training compound/scaffold overlap blocks the final holdout; unknown overlap stays visible.

Every requested protein/compound pair must have exactly one finite result or an explicit failure/exclusion, with job/artifact provenance and the comparison hash. Held-out rows must additionally carry the prior selection-lock hash. Failures disqualify a tuning candidate; a failed final pair fails evaluation without computing a selectively filtered error. No uncertainty or score is invented for missing results.

Affinity selection uses tuning MAE. Final evaluation reports MAE, RMSE, signed error, Pearson correlation, tie-aware Spearman correlation and deterministic scaffold-cluster bootstrap percentile intervals. Minimum sample/scaffold counts, MAE and correlation cutoffs must be explicitly approved before freezing. Few scaffolds and remaining analog dependence limit interval coverage. These intervals describe evaluation uncertainty; they are not model epistemic uncertainty. The selection lock does not prove that an external caller has never viewed labels; it enforces sequencing and provenance inside this harness.

## Local verification

Run a fresh Linux integration with the installed Linux Vina path:

```bash
BIOCENTRIC_VINA="$PWD/.tools/linux/vina" .venv-linux/bin/python -m scripts.verify_research_local \
  --output runs/linux-research-comparison-new
```

`--engines vina gnina unidock` includes those proposed available engines, with real native preparation and durable review continuation. `--affinity biocentric_affinity` or `--affinity gnina_affinity` requests an actual optional affinity integration after screening. Each requested runtime still needs its configured executable/weights and required inputs. Choose a new output directory to preserve older evidence. The runner explicitly mocks the literature/advisor/retrieval, spends no API credits, and executes real registered scientific calculations. Its output reports actual selected protocols, raw metrics, evaluation limitations and checkpoint behavior.

`tests/test_comparative.py` uses clearly labeled synthetic fixtures solely to test sequencing, identity, endpoint, leakage, failure and immutability invariants. Those tests are not biological or predictive validation.

For a real browser continuation, add `--pause-before-screening`. The runner stops with stages 01–07 complete and an actual unanswered exploratory-execution question; `browser-review-ready.json` records the same project/job. Launch the app against that store and answer the question to continue ordinary scientific workers. No completed model stage is rerun.

## Advisor and runtime integrity

Method-proposal and pilot-review model requests receive a redacted view of the frozen test set: final-heldout measurement values, activity labels and source-record payloads are withheld. Tuning values remain available; the complete immutable dataset remains in the project for deterministic evaluation and scientist review. This prevents direct table leakage through these advisor requests, but does not prove that public literature or pretrained models never contained related information.

New candidate freezes include scientific adapter source hashes as well as runtime fingerprints. A changed adapter, executable, GPU runtime manifest or preprocessing dependency blocks continuation of that newly frozen candidate. Child pipelines share the parent worker clock and persisted usage, so multiple pilots cannot each reset the project compute-time budget. Historical artifacts retain their original provenance fields.

New automatically prepared affinity requests retain an immutable mapping from each canonical compound to the exact registered source pose, receptor and passing physical-validation report. Native preflight checks molecule order, chemistry, coordinates and artifact hashes; results carry this mapping separately as `input_structure_provenance` and `input_structure_physical_passed`. These fields describe prior checks on the model input, including on partial model results. They do not establish affinity accuracy or shortlist eligibility, and the raw model output remains unchanged. BioCentric identity evidence must exactly match all supplied reviewed identity fields, including sequence, isoform, mutations and domain when present; absent optional identity fields are not invented.
