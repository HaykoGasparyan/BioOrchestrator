# Bio Orchestrator

A local scientific workspace for comparing protein pocket, docking and scoring pipelines against experimental literature data.

## Quick start on Ubuntu

```sh
git clone git@github.com:HaykoGasparyan/BioOrchestrator.git
cd BioOrchestrator
bash scripts/setup_ubuntu.sh
bash scripts/launch_ubuntu.sh
```

Open **http://127.0.0.1:8502/**. Create an experiment, select tools, enter a protein target, and start the evaluation. Saved experiments, live pipeline chains, final metrics and complex visualizations appear in the workspace. See [discovery workflow](docs/discovery-workspace.md), [evaluation](docs/discovery-evaluation.md), and [pipeline search and selection](docs/bayesian-pipeline-search.md).

The application setup installs the Python environment and Vina. GNINA, Uni-Dock and the BioCentric predictor require separate runtimes; follow [Ubuntu setup](docs/ubuntu-setup.md) and [GPU setup](docs/ubuntu-gpu-engines.md). BioCentric model source and trained weights are external dependencies; see its [model contract](docs/biocentric-model.md).

This repository contains source, tests, documentation and small public example fixtures. Machine-local experiments (`.biocentric/`), generated verification artifacts (`runs/`), downloaded executables, model weights, environments and credentials are not committed. Historical reports may reference those local artifacts; cloning does not restore them or reproduce their results automatically. The initial commits organize the existing implementation by subsystem; they are not backdated development snapshots.

Run offline tests with `.venv-linux/bin/python -m pytest -m 'not live'`. GPU inference and live docking verification are separate, explicitly launched checks.

A local, interactive research workbench for early, reversible noncovalent protein–small-molecule discovery. Python + Streamlit, a shared CLI/service, SQLite, immutable local artifacts, and a separate worker process.

**Every pipeline stage is visible.** Select a stage to see its input contract, preprocessing, selected tool, parsed output, and validation. Inspect actual model/tool requests and results in the event timeline. Review recorded scientific gates, cancel running work, resume from checkpoints, and download results with source citations and hashes.

The bundled example performs actual RDKit → Meeko → AutoDock Vina → PoseBusters calculations and creates reports. Offline UI fixtures are separately labeled and cannot start real jobs. The example is a software smoke test, **not a validated virtual-screening protocol or evidence of a new binder**.

## Research conversation — legacy workbench

The app is served at **http://127.0.0.1:8502/** by the Ubuntu launcher. New GPT research requires a privately configured API key; this continuation made no paid requests. Open the legacy workbench with `?view=workbench`, then choose **New conversation → Prepare public example → Start research**, or type a question. Preparing the example creates a draft. Starting research authorizes the displayed model budget. The default is **GPT-4.1 mini**; existing projects remain in **Classic workbench**.

The demo follows **Research → Frozen test set → Candidate methods → Comparative pilots → Locked protocol → Independent evaluation → Review → Experiment → Results**. The model uses literature and the frozen test set to propose a method; a measured pilot checks the proposal before screening. A visible column of stage cards shows progress and short recorded explanations. Click any card to inspect its decision, evidence, test set and exact calls alongside the journey. **Follow current step** tracks live progress; **Ligands & results** returns to rankings and downloads. Explicit stage selections stay selected during automatic refresh. Full public responses and raw inputs/outputs remain expandable.

**Ligands & results** shows the submitted screening library and updates its status as preparation and docking progress. The completed view shows names, docking scores, pose validity and the recorded shortlist. Download the shortlist SDF, ranked CSV, all registered docked poses, or an individual ligand's poses. Scores retain their method and units: Vina's more negative score ranks better within the same protocol, while learned CNN affinity uses different semantics. These scores are predictions, not measured binding. The bundled imatinib/aspirin/caffeine library is a software demonstration, not a discovery of new binders. See the [POC walkthrough](docs/demo.md).

**Send & continue** saves an inline answer and automatically requeues the same job. Completed stages are reused. Failed requests consume their reserved call; a retry requires an explicit answer within the saved budget. A successful path uses three model requests (literature, selection, pilot interpretation), within the default six-request project cap. Each literature pass allows at most four hosted search actions and four bounded primary-page reads. Actual model tokens, call counts and stage times are visible. The Vina example uses one local CPU. GNINA, Uni-Dock and BioCentric have separately verified local GPU runtimes; no remote worker is provisioned.

Before a pilot has run, **Revise inputs before continuing** pauses the same job for new files. Open **Add scientific input files** to upload revised structures, ligands, or a measured-assay CSV with molecular SMILES. **Rebuild test set & proposal** then freezes the revised inputs and requests one new tool proposal within the remaining model budget. Earlier manifests and decisions remain in the record. Once a pilot has completed, changed scientific inputs require a new conversation.

Search prefers original papers and official documentation and records accessible source excerpts. It is bounded research, not an exhaustive review. The deterministic test-set builder preserves assay, endpoint, species, construct, censoring, real labels and scaffold-separated partitions. The model explains tool suitability against that frozen manifest. A reference redock checks pose recovery and physical validity only; it does not validate affinity or screening enrichment.

| Method | Ubuntu execution status | Required inputs / result meaning |
|---|---|---|
| AutoDock Vina | Actual Linux baseline and application runs passed | Reviewed receptor/ligand PDBQT; canonical poses and empirical kcal/mol scores |
| AutoDock4 | Adapter implemented; blocked by runtime/maps | AD4 PDBQT and matching AutoGrid maps/field |
| GNINA docking | Actual CUDA docking and application checks passed | Reviewed receptor/3D SDF; empirical energy, CNN pose score and CNN affinity stay separate |
| Uni-Dock | Actual CUDA native/application checks passed with bounded 1.1.0 compatibility guards | Reviewed PDBQT; see supported chemistry/grid limits in GPU documentation |
| MM/PBSA / MM/GBSA | Blocked by Amber and simulation inputs | Matching topologies and sampled trajectory; docked SDF is insufficient |
| GNINA CNN affinity | Actual CUDA native/application checks passed | Actual poses; learned pK prediction |
| PLANET | Adapter implemented; no configured environment/weights | Protein pocket/reference and ligand graph; distinct from BioCentric |
| BioCentric affinity | Actual CUDA inference, original-reference equivalence and 128 predictions passed | Protein structure and ligand graph; learned PDBbind-scale pK, not measured affinity |

Vina preparation runs automatically after review. GNINA/Uni-Dock pilots and screening, plus BioCentric/GNINA affinity preparation from compatible registered inputs, now produce typed native manifests automatically. Their exact commands, hashes and resources still pass durable review. Advanced native inputs remain available for other contracts. See [engine contracts](docs/engines.md), [GPU verification](docs/ubuntu-gpu-engines.md) and [research planning](docs/research-planning.md).


The [local integration script](scripts/verify_research_local.py) uses mocked advisors and actual chemistry/docking tools. Its [recorded result](runs/research-conversation-verification/local-science-verification.json) completed all nine stages with automatic continuation, top-reference RMSD 0.3804 Å, and eight screening poses. This is software validation, not live-model or predictive validation.

The [diverse-library runner](scripts/validate_diverse_docking.py) exercises 128 distinct scaffold representatives from RDKit's bundled NCI sample with actual local Vina, plus two eight-compound parameter checks, adversarial preparation and checkpoint/export verification. [Engine contract tests](docs/engine-validation-matrix.md) cover synthetic contracts separately from the actual Ubuntu runtime evidence. Predictive performance remains unavailable without independent measured data.

## Historical Mac launch (preserved context)

The copied snapshot contained an Apple Silicon environment and binary. These historical launch commands apply to that original Mac; on Ubuntu use the launcher above:

```sh
cd /Users/hayk/Documents/BioGPT
.venv/bin/python -m streamlit run app.py --server.address 127.0.0.1
```

Open [the local workbench](http://127.0.0.1:8501). The default is the research conversation above. For local-only calculations, switch to **Classic workbench** and use **Load executable smoke example**, then **Run pipeline**. Under **Plan & Runs**, review the protocol, resume, inspect the pilot and record the exploratory choice when no held-out controls exist. The scientist controls these decisions.

Use **Load offline demonstration** to inspect the interface without scientific computation. Every fixture view/report is visibly marked. A fixture has no real docking or model results.

## Enable GPT with a private local key prompt

The launcher below creates a classic example for compatibility. After launch, select **Research conversation → New conversation** for the new workflow. The three-call Stage 03 descriptions below refer to **Classic workbench** only.

Create an OpenAI project API key with API billing/credits and access to the selected model. Follow the [official API setup guide](https://developers.openai.com/api/docs/quickstart). Do not paste the key into chat, a project file or a command-line argument.

Run this in your own interactive Terminal:

```sh
cd /Users/hayk/Documents/BioGPT
.venv/bin/python -m biocentric.launch --with-openai --example
```

The launcher asks for the key with terminal echo disabled (or uses an existing `OPENAI_API_KEY` environment variable). It does not save the key to disk or print it. It makes one small paid Responses connection check, then opens a separate loopback instance at **http://127.0.0.1:8502** and starts a new **1IEP · GPT-enabled public research** project. The model defaults to `gpt-4.1-mini`; `--model` changes it. The existing port-8501 local-only server is left running.

The public example enables external evidence retrieval and at most **three advisory model calls**, each capped at 1,500 output tokens. Input tokens are billed separately. The initial connection check adds one request capped at 16 output tokens. Scientific protocol/experimental-status review gates remain in place. Inspect **Model activity** for readable advice, evidence citations, request counts, recorded token usage, and request/response events.

An API key does not change an existing project's saved model permission. **Load executable smoke example** creates a local-only project with GPT disabled. To use GPT from the interface, choose **Load GPT-enabled public example**, then **Run pipeline**. Loading creates a separate project without starting a job or making a paid call. GPT is an optional advisor at **Stage 03: Strategy selection**; all other stages run scientific tools and workflow rules. Every stage inspector shows a **Decision summary**, explicitly distinguished from the model's public explanation.

Keep that terminal running. The key is held in the app and opted-in worker process environments; it is never stored in project metadata. Use **Cancel** in the app to stop a running worker; stopping the web server alone does not cancel an already detached scientific job. Starting again after these processes exit requires entering the key again.

## Installation elsewhere

Use Python 3.11–3.13 (3.12 tested). No system-wide install, Redis or GPU is required.

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -e '.[science,test]'
.venv/bin/python -m streamlit run app.py
```

Install the matching official [AutoDock Vina 1.2.7 release](https://github.com/ccsb-scripps/AutoDock-Vina/releases/tag/v1.2.7) at `.tools/vina`, put `vina` on PATH, or set `BIOCENTRIC_VINA` to an administrator-controlled executable. Vina's Python distribution is a supported fallback, but source compilation may require a compiler and Boost. The Mac binary shipped into this workspace is **not portable to Linux/Windows**; its source and SHA-256 are recorded in [example provenance](examples/1iep/sources.json).

Pinned direct scientific dependencies: [requirements-science.txt](requirements-science.txt). Optional PLIP/OpenBabel and receptor-preparation dependencies: [requirements-science-optional.txt](requirements-science-optional.txt). The complete tested environment is recorded in `requirements-lock-macos-arm64.txt`; it includes platform-specific packages. Review tool/data licenses separately before redistribution or commercial use.

## Classic workflow and file contracts

| Stage | Inputs | Processing and outputs |
|---|---|---|
| Target dossier | Protein name/accession or sequence + organism | Species-safe mapping; ambiguity blocks; source-linked records, retrieval dates and coverage gaps |
| Known chemistry | Resolved identity | Original assays retain Kd/Ki/IC50/EC50, units, relation, construct and source lineage |
| Strategy | Evidence + uploaded inputs | Deterministic strategy proposal; optional bounded language-model advice |
| Structure review | Original PDB/mmCIF, curated receptor PDBQT, observed site SDF | Inspect chains/components/uncertainty; derive site box from bound-ligand geometry; explicit protocol review |
| Preparation | `.sdf`, `.smi`, `.smiles` | RDKit sanitization, exact isomeric deduplication, retained metadata, explicit hydrogens, seeded conformers only if absent; canonical SDF stays authoritative |
| Pilot | Prepared reference + receptor | Meeko charges/types, Vina poses, chemical identity recovery, receptor-frame symmetry-aware RMSD, physical checks |
| Benchmark | Independent controls or explicit exploratory decision | No controls means unvalidated, never pass; separate held-out evaluation harness |
| Screen | Approved protocol + supported states | Cancellable Vina jobs, checkpointing, per-compound failures, PoseBusters, PLIP when available |
| Prioritize | Scores + validity + chemical identity | One-engine ranking, diversity, inclusion/exclusion rationale; missing ADMET/selectivity/availability remain unknown |
| Handoff | Candidates, evidence, logs, decisions | Markdown/HTML, CSV/JSON, pose/prepared SDF, raw logs, immutable manifest; proposed follow-up assays |

PDBQT alone is not a canonical ligand record. Multi-state enumeration is not silently performed. Input protonation/tautomers are retained; adding hydrogens by valence does not establish physiological protonation. Curated PDBQT must correspond to the supplied PDB and reference coordinate frame. mmCIF originals are preserved and inspected; the executable curated path requires a reviewed matching PDB export.

Upload `identity.json` or use the local identity review form for proprietary/custom targets. It must match project organism, accession, requested isoform/mutations and any supplied sequence. A user-supplied identity stays labeled as such. Source evidence can be imported as a typed JSON array, and measured assay results as CSV with compound ID, endpoint, value, units, relation, assay ID and organism. The interface provides a template.

See [scientific contracts, limitations and verification](docs/science.md) and [retrieval clients and source coverage](docs/retrieval.md).

## Model calls and privacy

In Classic workbench, model advice is optional. Public evidence retrieval and remote inference require disabling **Keep project inputs local**. A second explicit checkbox enables the model. Model identifier, call count and output-token caps are persisted per project. Configure `OPENAI_API_KEY` in the process environment; it is never written into prompts, artifacts or source files. `OPENAI_MODEL` sets the intake default. No key was required for the executed docking example.

The official [Responses function-calling](https://developers.openai.com/api/docs/guides/function-calling) and [structured-output](https://developers.openai.com/api/docs/guides/structured-outputs) interfaces are used through a separate adapter. The model can read existing evidence records using a strict allowlisted tool and return a validated advisory schema. It cannot submit arbitrary shell commands, claim a computation, alter scores, change budgets or approve a gate. API requests use `store=False` and bounded timeouts. Application code validates evidence IDs and requires cited observations to have been read; scientific entailment still needs review.

The timeline shows public requests/responses, tool arguments, status, usage and prompt version/hash. It does not expose private model reasoning. A live GPT-4.1 mini research run made three requests and recorded 14,553 input / 880 output tokens. The API and evidence-reading tool worked, but both returned advice attempts failed citation validation; this is not a successful scientific-advice evaluation. The rejected answer is readable under **Model activity**, separately labeled and excluded from accepted advice. See [live verification](runs/gpt-live-verified/verification.json).

The final permitted request is reserved for an answer. Citation schemas restrict references to opened evidence records; validation feedback can use only the remaining call budget. After a stopped failed attempt, **Correct rejected advice · 1 model call** offers one explicit correction when budget remains; the reservation prevents duplicate submissions. A project with its budget exhausted cannot request a correction. Existing failures and original responses remain in the timeline. The separate [model evaluation harness](evals/README.md) saves actual responses and limited independent rubric flags when explicitly run with credentials.

All files and project metadata stay below `BIOCENTRIC_HOME` (default `.biocentric/`). Local-only mode blocks outbound retrieval and model calls, including silent sequence submission. A local file does not become upload permission. External source content is untrusted evidence, never executable instructions. Use the app as a single-user loopback service; multiuser authentication and network deployment are outside this release.

## CLI

```sh
.venv/bin/python -m biocentric --help
.venv/bin/python -m biocentric adapters
.venv/bin/python -m biocentric example
.venv/bin/python -m biocentric status PROJECT_ID
.venv/bin/python -m biocentric run PROJECT_ID
.venv/bin/python -m biocentric decide PROJECT_ID protocol \
  --value '{"approved":true,"receptor_reviewed":true,"ph_assumption":7.4}' \
  --reason 'Reviewed the recorded receptor, pocket and explicit preparation assumptions'
.venv/bin/python -m biocentric resume PROJECT_ID
.venv/bin/python -m biocentric decide PROJECT_ID exploratory \
  --value '{"approved":true}' --reason 'Exploratory protocol; no held-out screening validation'
.venv/bin/python -m biocentric resume PROJECT_ID
.venv/bin/python -m biocentric export PROJECT_ID
```

`run/resume --foreground` uses the same worker code synchronously for automation. The interface always starts a separate worker. Jobs use atomic SQLite claims and durable checkpoints; a retry reuses the same job rather than silently duplicating work. Original inputs and approved configurations are immutable after preparation; create a new project for a changed scientific run.

`BIOCENTRIC_PRODUCT_NAME` changes the displayed name. `OPENAI_MODEL` changes the intake default, never a running project's recorded model.

## Verification and benchmarking

```sh
.venv/bin/python -m pytest -q
.venv/bin/python examples/run_live.py --output runs/new-live-smoke
```

The core tests exercise target ambiguity, privacy, endpoint/unit fidelity, chemical identity/stereochemistry, unsupported chemistry, timeout/cancellation, idempotent claims/resume, failed gates, hashes/citations, strict model tools and untrusted content. Streamlit AppTests exercise real controls and rendering. Local 3D receptor/pose inspection uses a pinned, vendored viewer; structure files are not sent to a viewer service. See [verification results](docs/verification.md) for actual run details and remaining limits.

The [full application report](runs/integration-verified/report.md) and [scientific smoke report](runs/live-verified/report.md) contain real calculations and linked manifests. No independent active/inactive screening dataset was supplied, so the honest predictive-performance status is **unvalidated/exploratory**.

The independent [benchmark harness](biocentric/benchmark.py) requires dataset version/hash, disjoint tuning/evaluation IDs and scaffold groups, original measured labels, one defined endpoint, frozen criteria, and score provenance. `freeze-benchmark DATASET PARAMETERS OUTPUT` saves the protocol before evaluation. `benchmark DATASET FROZEN_PROTOCOL SCORES` reports precision/top-k enrichment, class balance and uncertainty; it preserves failed criteria. It cannot establish chronology or assay truth merely from a supplied JSON file, and its result does not automatically promote an application run to validated.

## Current boundaries

- Predicted structures can be imported/inspected; prediction, uncertain-pocket validation, automatic loop rebuilding and general target–construct alignment are not implemented.
- Only supported reversible, noncovalent small molecules. Covalent, important-metal, membrane-specific, peptide, macrocycle and other special chemistry need separate validated workflows.
- Local preparation functions for PDBFixer, PDB2PQR/PROPKA and Meeko are inspectable library APIs. The default interface uses reviewed curated receptors; it does not silently repair proteins.
- Live source clients are bounded retrieval, not exhaustive literature review; source ambiguity/errors remain visible. Abstracts are not treated as read full text.
- No validated predictive ADMET, QSAR, molecular dynamics, free-energy, target-generation or retrosynthesis adapter is claimed. No experimental assay is physically executed.
- The initial shortlist is a transparent score/diversity heuristic; quantitative multiobjective optimization, ensemble protocols and automatic cross-docking are future work.

[Optional-module backlog](docs/backlog.md) lists the next extensions and their validation requirements. The full requested specification and supplied scientific catalog are retained in `docs/` alongside the versioned runtime prompt in `biocentric/prompts/`.
