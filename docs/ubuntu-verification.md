# Ubuntu verification — 26 September 2026

The copied application now runs in a fresh Ubuntu environment. Real Linux Vina, GPU GNINA, GPU Uni-Dock and the existing local BioCentric model completed native and application calculations. A new 128-compound screen, BioCentric scoring of those same 128 structures, and a larger 256-compound screen have complete accounting and saved resume evidence. Historical Mac results remain readable and hash verified.

These are software, structure-processing and reference-pose checks. The authorized mouse ABL 1IEP example and unlabeled NCI library do not establish affinity accuracy, enrichment, new binders or a predictively superior method. No independent measured-affinity holdout was evaluated. No paid LLM calls were made.

## Host, installation and preserved history

The [recorded host inventory](../runs/ubuntu-verification/host-inventory.json) identifies Ubuntu 24.04.3 LTS, x86_64, kernel 6.8.0-139-generic, Intel i9-12900K/24 logical CPUs, approximately 62.6 GiB RAM and 48.4 GiB free disk at inventory time. The single RTX 4090 has 24,564 MiB VRAM, driver 535.309.01 and compute capability 8.9. The driver's reported maximum CUDA version is 12.2. Actual CUDA tensor, graph and model computations verified compatibility; the restricted command namespace's initial NVIDIA failure was not a host-driver failure.

The application uses new `.venv-linux` Python 3.12.3 with RDKit 2025.3.6, Meeko 0.7.1, PoseBusters 0.5.1 and Linux Vina 1.2.7. The copied Mac environment, executable and lock remain preserved. BioCentric uses an already existing, separate Linux Python 3.9.7/PyTorch 2.1.0+cu118/PyG 2.3.0 environment; it was inspected and executed, not rebuilt. [Observed model requirements](../requirements-biocentric-observed-linux.txt) record its relevant exact pins. GNINA's isolated CUDA libraries do not modify that environment or the system driver. See [setup](ubuntu-setup.md) and [GPU runtime installation](ubuntu-gpu-engines.md).

The original store had seven projects and 1,508 immutable artifacts. Both requested historical project IDs are present. Backups precede relocation; all original artifact hashes verified. An explicit root resolver preserves original Mac paths and hashed contents, and host/boot/process identity protects against copied worker PIDs. The [relocation report](ubuntu-relocation.md) links verification and backup locations. No Git commit or PR was possible because the supplied `.git` is not a working repository.

## Real configured-engine checks

All configured rows below executed scientific calculations. Installed files or passing synthetic adapter tests alone were not counted as real execution.

| Method | Saved result | Evidence |
|---|---|---|
| Vina 1.2.7, CPU | 144 calculations in the 128 screen; 272 in the 256 screen; zero execution failures | [128 summary](../runs/ubuntu-diverse-docking-verification/summary.json), [256 summary](../runs/ubuntu-diverse-256-verification/summary.json) |
| Uni-Dock 1.1.0 compatibility binary, GPU | Native/application imatinib outputs byte-identical; −9.316 kcal/mol; native 2.115 s | [verification](../runs/ubuntu-gpu-engine-verification/unidock-1iep/verification.json) |
| GNINA 1.3.2 docking, GPU | Native/application outputs byte-identical; empirical −13.33497 kcal/mol, CNN pose 0.9672928452, CNN affinity 8.385184288 pK; native 8.255 s | [verification](../runs/ubuntu-gpu-engine-verification/gnina-docking-1iep-nvrtc/verification.json) |
| GNINA score-only, GPU | Native/application outputs byte-identical; empirical −12.51667 kcal/mol, CNN pose 0.9203819633, CNN affinity 8.3367071152 pK; native 0.751 s | [verification](../runs/ubuntu-gpu-engine-verification/gnina-affinity-final-runtime/verification.json) |
| BioCentric three-model ensemble, GPU | Native 8.1943181356 versus adapter 8.1943178177 pK; largest per-model difference 9.54×10⁻⁷, below predeclared 10⁻⁵ tolerance; identical-request resume reused the pair with zero new predictions | [latest reference equivalence](../runs/linux-biocentric-reference-v2/verification.json) |

GNINA initially failed after docking because `libnvrtc.so.12` was absent; that [failed run](../runs/ubuntu-gpu-engine-verification/gnina-docking-1iep/verification.json) remains preserved. Installing the pinned NVRTC dependency resolved it. The successful docking log still contains a native `ligand outside box` message; score parsing never overrides physical checks. Sampled GNINA process memory reached 530 MiB. Uni-Dock's installed compatibility build is explicitly restricted to supported chemistry and boxes below 30 Å; its smoke check is not a general chemistry guarantee.

BioCentric's [code-supported contract](biocentric-model.md) establishes structure/pocket plus posed-ligand graph inputs, three usable checkpoints, native featurization and inverse scaling. It is distinct from PLANET and from the planning LLM. Its output is a learned PDBbind-scale pK score with mixed training endpoints, not measured Kd/Ki or a probability. Ensemble spread is not calibrated uncertainty. Receptor/ligand bytes, identity, configuration, code and weights are bound to reviewed requests. CPU inference also ran, but differed from GPU by approximately 0.00342 pK; those devices remain separate protocols.

## Frozen comparison and interactive research execution

The [mixed-engine verification](../runs/linux-multiengine-research-v1/local-science-verification.json) completed all nine stages in 168.06 s. Candidate input/runtime/adapter hashes, parameters, criteria and the comparison record were frozen before pilots; the selected protocol was locked before the held-out status and production execution. The receptor-frame, symmetry-aware top-pose cutoff was 2 Å, with physical validity and complete reference accounting required.

| Frozen pose protocol | Reference RMSD | Decision |
|---|---:|---|
| GNINA, exhaustiveness 4 | 4.3906 Å | Ineligible: failed pose recovery |
| Uni-Dock, balance | 0.9523 Å | Passed |
| Vina, exhaustiveness 4 | 0.3355 Å | Selected under the declared rule |
| Vina, exhaustiveness 8 | 0.3355 Å | Passed; tie resolved by protocol ID |

Production retained nine screening poses and then made three genuine BioCentric predictions with zero failures. Completed same-job continuation issued zero new scientific calls. Five durable review answers and 200 registered artifacts are recorded. This runner explicitly mocked advisor, literature and retrieval responses, with no network or paid model calls; its preparation, docking, validation and affinity execution were real. The separate [Vina workflow](../runs/linux-research-comparison-vina/local-science-verification.json) also completed with zero calls on completed-job resume.

The mixed record's protein identifier binds the receptor hash but retains `unknown` species/construct metadata from its mocked proposal. It cannot support target-matched predictive evaluation. The actual authorized ABL source context is recorded separately in the [primary-source review](../runs/ubuntu-verification/primary-source-review.json). Held-out evaluation remains `not_evaluated`, predictive validation `not_established`, and selection scope `software_pose_only`. GNINA's pretrained overlap is unknown. The measured-data comparison API exists, but a scientifically compatible measured cohort and reviewed predictor adapter were not supplied or executed. See [comparison semantics and limits](comparative-evaluation.md).

## Bounded 128- and 256-compound screens

Both fresh Linux workloads select one representative per exact nonchiral Bemis–Murcko scaffold from RDKit's bundled NCI sample, using seeded Morgan-fingerprint MaxMin selection. The source has 4,999 records, 2,536 eligible molecules and 492 eligible scaffolds; selection seed is 42, molecular weight 100–500, heavy atoms 8–35 and rotatable bonds at most eight. The input state, charge and tautomer are retained. No pH/state enumeration, activity labels or inactive controls were invented.

Each baseline used one Vina CPU, exhaustiveness 1, one pose, seed 42 and 60 s per-ligand timeout. Eight-compound seed-43 and exhaustiveness-4 variants used separate protocol fingerprints. These low-depth runs measure bounded software throughput, not exhaustive sampling.

| Accounting/resource | 128 screen | 256 screen |
|---|---:|---:|
| Prepared / requested | 128 / 128 | 256 / 256 |
| Distinct selected scaffolds | 128 | 256 |
| Baseline poses | 128 | 256 |
| Physical passes / retained failures | 121 / 7 | 244 / 12 |
| Execution failures / unaccounted compounds | 0 / 0 | 0 / 0 |
| Additional variant calculations | 16, all physical passes | 16, all physical passes |
| Total actual Vina invocations | 144 | 272 |
| Whole-run wall time | 327.465 s | 618.416 s |
| Python / child CPU time | 11.321 / 408.675 s | 19.239 / 778.621 s |
| Python peak / maximum child RSS | 231.2 / 531.3 MiB | 248.0 / 531.5 MiB |
| Recorded evidence disk size | 10,961,815 bytes | 17,625,352 bytes |
| Resume checkpoints / new Vina calls | 128 / 0 | 256 / 0 |
| Resume wall time | 0.147 s | 0.292 s |

Both stayed within their 1,800 s budgets and used no GPU, network, LLM or PLIP calls. The 256 baseline's 12 physical failures comprise seven non-aromatic ring geometry failures, two internal clashes, one internal energy failure and two protein-distance failures. They remain in all-output exports and outside valid shortlists. The CSV/SDF exports contain exactly 256 all-pose records and 244 valid-pose records. Preparation adversarial controls passed, retaining two supported inputs and explicitly excluding five incompatible cases.

The [independent 256 audit/import](../runs/ubuntu-verification/256-independent-audit-and-import.json) verified 2,221 saved inventory hashes, unique identities, scaffold count, all current pose hashes, 272 start events and 256 resume events with zero new start events. The original runner recorded unchanged pose hashes and modification times on resume. Audit/import itself executed no science. It registered and verified 1,940 artifacts in new completed main-store project `b15deb085ea846b595d31be73f38c38e`; repeat import reused that project. All prior project IDs and both named historical projects remained unchanged. The new 128 main-store project is `e990976e062b4b0bbb9af936dd057eda`.

## BioCentric on the same 128 structures

The [frozen scoring protocol](../runs/ubuntu-biocentric-diverse/frozen_protocol.json) binds the actual Linux 128 Vina poses. Four chunks of 32 run sequentially, with native batch size one and two CPU threads while all three models remain loaded. The [summary](../runs/ubuntu-biocentric-diverse/summary.json) records 128 predictions, zero model failures, 18.545 s of worker time, 32.807 s including resume, peak RSS 925,556 KiB and peak allocated CUDA memory 87,532,032 bytes. Allocated tensor memory is not total process VRAM. Scores ranged from 2.6871 to 7.5104 learned pK. All 128 checkpoints reused with zero new predictions.

All seven physically invalid Vina poses were scored transparently and remained excluded from the valid shortlist; 121 were eligible. [Per-pair comparison JSON](../runs/ubuntu-biocentric-diverse/comparison.json), [CSV](../runs/ubuntu-biocentric-diverse/comparison.csv) and [SDF](../runs/ubuntu-biocentric-diverse/comparison.sdf) retain distinct Vina and BioCentric metrics, exact compound/protein identities, actual structures, pose/checkpoint provenance and failure/eligibility information. No unlike scores were blended into a ranking. These results are imported into the separate affinity-results view of the new 128 project.

The bridge retains per-pair preprocessing/model failures, emits observed progress/resources and supports same-request resume. Current checkpoints additionally verify a digest of each result row. OOM cleanup releases exception-held tensors before collecting objects and clearing the CUDA cache; at batch size one a failed pair remains explicit and later pairs may continue. No physical GPU OOM was induced: that cleanup path has synthetic lifetime tests. Earlier workload evidence retains its original bridge version and hashes; the latest native equivalence check covers the strengthened bridge.

The [actual application interruption/resume test](../runs/linux-biocentric-research-resume-v1/verification.json) used 16 real Linux poses in an isolated store. Cancellation retained two completed predictions. The application resumed the same reviewed request and output directory, reused those two checkpoints with unchanged SHA256 and modification times, and computed only the remaining 14. Final accounting is 16 completed, zero failures or duplicate IDs; 17 original output artifacts were archived and verified. Total test time was 7.886 seconds on the RTX 4090. The final resumed native pass peaked at 79,754,752 bytes of PyTorch CUDA allocation and 921,432 KiB process RSS. This verifies interrupted-work reuse through the application, in addition to the earlier native-only resume checks. Explicit failed-pair retry is covered by deterministic subprocess tests; a physical GPU OOM was not induced.

## Browser review, regression and remaining limits

Recorded browser downloads already match the application's verified bytes: 128-row affinity CSV, ranked-ligand CSV and a ten-molecule actual shortlist SDF. The [download evidence](../runs/ubuntu-verification/browser-downloads.json) includes SHA256s and sizes.

[Browser continuation verification](../runs/linux-browser-review-v1/browser-verification.json) completed all nine stages of project `a6e3ec72712943558a5866a2f1f1299e` on the same job through four worker attempts. Both the exploratory and BioCentric native review were answered through the actual browser. The regular worker produced nine docking poses and three real CUDA BioCentric predictions in 3.318 seconds of native adapter time; all 150 registered artifacts verified. Its protein identity explicitly records P00520, Mus musculus, chain A and the receptor hash. A genuine progress event records pair accounting/resources. The downloaded prediction SDF contains all three actual structures and exactly matches its registered hash. Stage selection, the locked Vina-exhaustiveness-4 review, unavailable holdout status, and return to results were checked in the browser. Advisor stages used explicitly labeled test doubles; no live paid planner verification is claimed.

**Final full regression: 572 passed, zero failed, three existing nonfatal warnings, in 41.21 seconds.** Command: `.venv-linux/bin/python -m pytest -q`; [complete log](../runs/ubuntu-verification/full-regression-final.log). `pip check` reports no broken requirements. Synthetic tests verify software behavior and are not predictive validation.

AutoDock4 remains blocked by missing executable and matching AutoGrid maps. PLANET lacks its separate configured source/environment/weights. MM/PBSA and MM/GBSA lack Amber, matching parameterized topologies and an equilibrated sampled trajectory. PLIP is unavailable in the fresh application environment. No docked SDF was substituted for MD evidence, and none of these methods was declared scientifically validated.

Result tables now paginate at 100 rows; browser verification reached rows 101–128 and returned to page one. The [streamed disk export](../runs/ubuntu-diverse-256-streamed-export/export-manifest.json) writes CSV/JSON incrementally and splits actual poses into four 64-record SDF parts, with a separate ten-pose shortlist. All 256 molecular records and nine output hashes verified. Browser downloads and project snapshots remain bounded in memory.

The exercised scale is 256 docking compounds and 128 model pairs, below the existing 1,000-input cap, 200-row test-set preview, 50 MiB upload/export limits and per-engine budgets. Chunking and durable checkpoints support this bounded workload; this report does not claim arbitrary-size streaming or successful near-limit GPU batches. Complete experimental validation still requires compatible measured target/species/construct/assay/endpoint records, a frozen independent evaluation protocol and resolved training-overlap concerns.

Launch from this checkout with `bash scripts/launch_ubuntu.sh`, then open <http://127.0.0.1:8502>. It uses loopback and does not store an API key. See [setup instructions](ubuntu-setup.md) for installation, overrides and separate environment details.
