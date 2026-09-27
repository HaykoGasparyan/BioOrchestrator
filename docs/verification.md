# Ubuntu verification — 2026-09-26

The current Linux suite passed **572 tests**. Actual Vina, GNINA, Uni-Dock and BioCentric runs are recorded in [Ubuntu verification](ubuntu-verification.md), with [setup and launch instructions](ubuntu-setup.md). The records below are preserved Mac history, including their test counts and runtime limitations; they are not new Linux results.

# Verification — 2026-09-26

## Interactive stages and diverse ligand validation

The research workspace now keeps clickable stage cards visible beside the selected decision or ligand results. Cards show short recorded explanations and distinguish model activity from workflow rules. The selected pane exposes the public decision, test set, sources and expandable exact calls. Selection survives periodic refresh; **Follow current step** resumes automatic following. Inline review questions continue the same job. Browser checks opened the saved tool-choice stage and its completed GPT response, then inspected the results view. [UI screenshot](../runs/interactive-workspace-verification/ai-decision.png).

The [larger validation](diverse-validation.md) completed **144 actual Vina calculations**: 128 chemically diverse baseline compounds and two eight-compound parameter checks. All returned poses; physical checks passed for 121 baseline poses and all 16 comparison poses. Seven baseline poses failed and remain visible, excluded from the shortlist. All 128 checkpoints reused unchanged files with zero new docking calls; all 1,197 source evidence hashes and exported molecular identities verified. Runtime was **513.8 seconds**, with no GPU, network or model calls.

The [eight-engine contract matrix](engine-validation-matrix.md) adds 138 synthetic cases and fixes score parsing, units, uncertainty and output-identity bugs. Optional engines were not executed; their runtimes still need configuration and independent live verification. The larger Vina run tests software behavior and parameter sensitivity on one receptor, not binding accuracy or screening enrichment. Historical verification checkpoints below are retained as records of earlier application states.

The final full suite passed **404 tests** in 8.08 seconds, with the three existing nonfatal Pydantic/RDKit warnings. The separate saved [128-ligand app run](http://127.0.0.1:8502/?project=726fbd0fea4e4b7496088bc3e458fd5b) displays all baseline results and a ten-pose physically valid shortlist. Its import is idempotent, preserves original source provenance and requires complete hash-verified evidence before creating app state. [App/export verification](../runs/interactive-workspace-verification/verification.json) records the counts and export-payload hashes.

## Ligand results and simplified demo update

The user's research project `6a2b1de5ee184b1b8c03b35249ed2057` has now completed all nine stages. The main screen shows its ligand ranking immediately below completion, followed by direct CSV/SDF downloads. Long explanations, source records, exact calls and resource usage are optional details. The [demo walkthrough](demo.md) explains the visible research → test set → method → pilot → screen → results sequence.

The final automated suite passed **234 tests** in 6.86 seconds with the three existing nonfatal warnings. Tests cover ranked/live ligand views, score direction and units, unknown/unranked scores, valid shortlist membership, native inputs, frozen-versus-draft libraries, same-project artifact ownership, content hashes, path restrictions, CSV exports and completed-run display order.

[Verified exports](../runs/ligand-ui-verification/results-verification.json) contain three shortlisted ligand poses and eight total docked poses; every exported SDF record parsed successfully with coordinates. The browser's **Kept ligands · SDF** download matched the verified shortlist hash. A [separate small real Vina run](../runs/ligand-ui-verification/progress-verification.json) confirmed per-compound progress and checkpoint reuse without any model call. The completed user's calculation was not rerun or changed.

The current example ranks imatinib, aspirin and caffeine at −13.211, −7.143 and −5.888 kcal/mol respectively, all with recorded passing pose checks. These are demonstration controls and exploratory docking scores. No independent affinity or screening benchmark is established, and no new-binder claim is made.

## Research conversation update

Final automated suite: **213 passed, 0 failed** in 6.40 seconds. The same three nonfatal Pydantic/RDKit warnings from the baseline remain. Coverage includes model/web compatibility, strict engine contracts, assay uploads, input revision, unchanged/historical reuploads, frozen hashes, budgets, same-job continuation, and Streamlit interaction.

The GPT-enabled conversation is running at [localhost:8502](http://127.0.0.1:8502/?project=6a2b1de5ee184b1b8c03b35249ed2057). Browser verification confirmed expandable stage cards, the cited Vina proposal, model usage, and the pending inline review. Public model explanations and actual request/response events are visible; private internal reasoning is not exposed.

The [live GPT record](../runs/research-conversation-verification/live-gpt-verification.json) contains two successful GPT-4.1-mini responses: hosted web research and tool selection. Two earlier rejected API requests remain in its history. Recorded usage is 34,489 input and 1,427 output tokens, with 17 literature records. Of four direct page reads, one returned an accessible excerpt and three were unavailable. Target/assay retrieval timed out, so no independent matched assay benchmark is claimed. The live project is waiting for the user's review before its scientific pilot.

The [separate local integration run](../runs/research-conversation-verification/local-science-verification.json) used mocked advisors and actual RDKit, Meeko, Vina, PoseBusters and PLIP. All nine stages completed, the same job resumed through two review gates, and 132 immutable artifacts were hash-checked. The top reference pose recovered the crystal pose at 0.3804 Å RMSD; screening produced eight poses. This validates software operation and reference pose recovery, not binding-affinity accuracy or screening enrichment.

Vina is installed and live-tested on this Mac. AutoDock4, GNINA, Uni-Dock, MM/PBSA, MM/GBSA, GNINA CNN affinity and PLANET have typed adapters and native input contracts; their execution is unavailable until the required runtime, hardware and scientific inputs are configured. They have not been live-validated here. Remote GPU dispatch is not implemented.

## Initial classic-workbench verification (earlier checkpoint)

The remaining sections preserve the earlier baseline results; the research-conversation update above supersedes their model and interface availability statements.

The implemented application is running locally at http://127.0.0.1:8501. Its health endpoint returned `ok` after the final changes. The UI, CLI and background worker share the same SQLite-backed service.

## Automated checks actually run

`.venv/bin/python -m pytest -q`: **78 passed, 0 failed** in 3.45 seconds.

Coverage includes species/target ambiguity, requested mutation and sequence identity, endpoint/unit/relation fidelity, stereochemical identity and SDF metadata, unsupported compounds, offline privacy, template-fetch prevention, strict model tools, unread/unknown citations, untrusted retrieval content, timeouts/cancellation, atomic job claiming, stale-worker ownership, idempotent resume, input/checkpoint integrity, explicit failed/exploratory gates, source/result provenance, actual Streamlit control interactions, viewer artifact/frame checks and injected-script escaping.

Three nonfatal warnings were observed: two Pydantic warnings about the scientific field name `construct` shadowing its legacy method, and an RDKit duplicate converter registration warning. `.venv/bin/python -m pip check` reported no broken requirements. Python compile checks passed. The generated viewer JavaScript passed Node syntax checking.

Streamlit AppTest rendered and exercised intake, executable example loading, fixture labeling, all ten stage inspectors, local identity decisions and 3D component generation. **Actual WebGL/browser visual rendering remains unverified**: no browser connector was available and the native browser was locked. No authentication boundary was bypassed.

## Full application integration run

[Report](../runs/integration-verified/report.md), [verification JSON](../runs/integration-verified/verification.json), [manifest](../runs/integration-verified/manifest.json), [shortlist SDF](../runs/integration-verified/shortlist.sdf).

- Project: `86e3f91a7c0f4c9dbefee105695a3055`.
- Durable job: `ffaf912a1df5478c81f6a7e685076d2e`.
- All ten stages completed using a **separate background worker**, with a recorded protocol gate and explicit exploratory test-run decision.
- Three prepared compounds; eight actual Vina poses; three distinct shortlisted hypotheses; 98 persisted events; 128 registered artifacts at verification.
- Every pose passed the configured PoseBusters checks. The reference is a curated mouse ABL/P00520 example, PDB 1IEP. The verified full-length UniProt sequence is shown separately from the experimental kinase-domain construct.
- PLIP generated raw XML/contact reports, but its perceived ligand chemistry differed from canonical SDF. These reports are marked **completed_requires_review** and do not establish correct interaction chemistry.
- Completed-stage checksums, canonical pose SDF, raw PDBQT, logs, configurations, CPU, seeds, executable hashes, tool versions, original inputs and decisions are retained.
- Actual compute time recorded: approximately 95.4 seconds. This is a software smoke workload, not a general performance estimate.

The independent scientific-adapter runner also succeeded: [live-verified report](../runs/live-verified/report.md). Its cognate imatinib top pose has Vina score −13.211 kcal/mol and symmetry-aware receptor-frame RMSD 0.380 Å. These are computed outputs, **not measured affinity or evidence of screening enrichment**. Caffeine and aspirin are untested chemical controls, not tested inactives. No independent active/inactive dataset was supplied. The screening-predictivity gate therefore remains unvalidated/exploratory.

## Live integrations

All six public-data clients were actually exercised: UniProt, Europe PMC, Open Targets, RCSB PDB, ChEMBL and PubChem. The public SRC dossier returned 40 original assay records and four structure records. Ambiguous Open Targets gene mapping was blocked; an explicitly verified Ensembl query was separately exercised. See [retrieval verification](../runs/retrieval-verification.json) and [coverage details](retrieval.md). These searches are bounded and are not exhaustive evidence reviews.

RDKit, Meeko ligand preparation, Vina, PoseBusters and PLIP ran on the real 1IEP workflow. PDBFixer, PDB2PQR/PROPKA and Meeko receptor preparation were additionally exercised on standard protein PDB 1CRN. The tutorial 1IEP PDB failed automatic preparation due to its interrupted residue records; those failures are retained. Curated 1IEP receptor PDBQT succeeded. Protein preparation outputs still require scientific review; no universal repair capability is claimed.

## Not live-verified / deferred

- OpenAI Responses: implemented strict schema, allowlisted read-only tools, bounded calls, request/result tracing and provenance checks. Offline adapter tests pass. **No live call was made: `OPENAI_API_KEY` is absent.** [Model evaluation harness](../evals/README.md) is available for explicit credentialed testing; no live model scientific-accuracy score is claimed.
- No independent target-specific screening benchmark; no automatic promotion to a validated result.
- Automatic target-to-construct/isoform alignment, predicted-pocket qualification and general receptor repair remain scientist-reviewed or deferred.
- Ligand-based strategy is advisory; no validated QSAR/pharmacophore execution is claimed.
- Optional ML rescoring, ADMET, structure prediction, dynamics, free energies and generation are not implemented. Registry entries show them as deferred.

The current application is a single-user local research tool, not a clinically validated system or a remotely secured multiuser deployment.
