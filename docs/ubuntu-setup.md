# Ubuntu continuation and launch

This checkout runs the existing Python/Streamlit application on Ubuntu 24.04 x86_64 with Python 3.12.3. The new application environment is `.venv-linux`; the copied Mac `.venv`, Mach-O `.tools/vina`, historical manifests and Mac dependency lock are preserved. The supplied `.git` directory is not a working Git repository, so this continuation has no Git commit or PR.

```bash
cd /home/sargis/BioGPT/BioGPT
bash scripts/setup_ubuntu.sh
bash scripts/setup_gpu_runtimes.sh
bash scripts/launch_ubuntu.sh
```

Open [BioCentricGPT on loopback](http://127.0.0.1:8502). The launch script sets the local scientific runtimes and never stores an API key. Set `BIOCENTRIC_PORT` to choose another local port or `BIOCENTRIC_HOME` for a separate experiment store. For access from another computer, use an SSH tunnel to port 8502; this single-user app has no remote authentication layer.

The source supports Python 3.11–3.13. The Ubuntu lock records the actual fresh Python 3.12 install and passes `pip check`. This machine lacks Python's `ensurepip`; setup therefore creates the environment with `--without-pip` and uses the system pip's `--python` option to install into it. It does not install system packages or change the GPU driver.

| Runtime | Verified host setup | Evidence |
|---|---|---|
| CPU | Intel i9-12900K, 24 logical CPUs, about 62 GiB RAM | `runs/ubuntu-verification/host-inventory.json` |
| GPU | RTX 4090, 24,564 MiB, driver 535.309.01, compute capability 8.9 | Actual tensor/graph computations and native GPU logs |
| Application science | RDKit 2025.3.6, Meeko 0.7.1, PoseBusters 0.5.1, Vina 1.2.7 Linux | `requirements-lock-ubuntu-x86_64.txt` |
| GNINA / Uni-Dock | GNINA 1.3.2 with isolated CUDA 12.1/cuDNN 9.1 libraries; Uni-Dock 1.1.0 compatibility binary | [GPU installation and actual execution](ubuntu-gpu-engines.md) |
| BioCentric | Existing separate Linux Python 3.9.7, PyTorch 2.1.0+cu118 and graph extensions | [Inspected model contract and real equivalence checks](biocentric-model.md) |

The NVIDIA device is hidden in restricted command execution. `nvidia-smi` failed there but succeeded on the host. Actual model and native GPU calculations were executed with approved host access. Framework imports and the driver version alone are not treated as proof of CUDA compatibility.

BioCentric defaults in the launcher refer to `/home/sargis/dti` and `/home/sargis/anaconda3/bin/python` through the current user's home directory. Override `BIOCENTRIC_MODEL_DIR` and `BIOCENTRIC_MODEL_PYTHON` for a different installation. The original model source is unchanged. Its observed pins document the verified existing environment; a newly rebuilt model environment has not been claimed verified.

History was backed up before migration, all 1,508 supplied artifacts were hash-checked, and explicit root mappings keep completed records/downloads readable. [Relocation details and worker lease safety](ubuntu-relocation.md) describe the evidence. New host calculations use separate run/work directories.

AutoDock4 remains blocked without its executable and matching AutoGrid maps. PLANET has no configured separate checkout/environment/weights. MM/PBSA and MM/GBSA lack the Amber runtime, matching parameterized topologies and sampled trajectory. No docked pose is substituted for MD evidence. PLIP is optional and unavailable in the new application environment; its absence is recorded.

No OpenAI key was copied or configured, and this continuation spent no paid API calls. The existing private terminal key prompt remains available via `.venv-linux/bin/python -m biocentric.launch --with-openai`; it makes its documented connection-check request only when explicitly invoked. Scientific subprocess environments strip credential variables. GPT-4.1-mini remains the configurable planner; the verified local BioCentric model is a separate scientific predictor.

The authorized target is the supplied mouse ABL 1IEP example. The [primary-source review](../runs/ubuntu-verification/primary-source-review.json), frozen comparison records and software controls do not contain an independent measured-affinity dataset. Predictive validation stays unavailable. See [comparative evaluation](comparative-evaluation.md) for the selection/holdout boundary.

## Current scaling and export boundaries

The Linux verification runs are bounded jobs. They do not establish support for an arbitrarily large library. These limits remain enforced:

| Surface | Current bound |
|---|---|
| Project uploads / imported measurements | 50 MiB per upload / 10 MiB measurement import |
| Scientific preparation and Vina library | 1,000 compounds; 100 MiB per science input file |
| Native Uni-Dock batch / other native SDF input | 1,000 ligands / 1,000 valid SDF records |
| Research test-set preview | At most 200 compounds, with explicit exclusions beyond the chosen bound |
| BioCentric | Default 128 pairs, configurable to at most 1,000; 50 MiB input files; 20,000 receptor ATOM records; 256 ligand heavy atoms; one pair per model inference batch |
| Molecular preview and its individual download | 10 MiB per registered artifact |
| Combined pose SDF download | 50 MiB |
| Native output registration | 100 MiB per output file |

Streamlit tables paginate every 100 result rows; downloads still include all recorded rows. Browser CSV/JSON and combined SDF downloads remain bounded, in-memory payloads. For incremental disk output, use:

```bash
.venv-linux/bin/python -m biocentric.cli export-streamed PROJECT_ID runs/new-export --poses-per-part 64
```

This writes CSV/JSON record by record and authenticated SDF parts bounded by both record count and 50 MiB, plus a hash manifest. It publishes a fresh output directory only after all records verify; failure never publishes a partial export or changes original artifacts. The 256-compound run exported four 64-pose SDF parts and a ten-pose shortlist, with all molecular records and nine output hashes verified. Project metadata and derived rankings are still loaded in memory. General library scheduling beyond the enforced input cap and HTTP streaming downloads remain unimplemented; raising constants alone would not solve those limits.

Vina has durable per-compound checkpoints. BioCentric has request- and row-hash-verified per-pair checkpoints, actual progress/resource events, cancellation and timeouts, and explicit failed rows for unfinished pairs after native process failure. After pair exceptions, including OOM, its bridge releases the exception traceback, collects garbage and clears unused CUDA cache. This is failure isolation with a fixed batch size of one; it is not adaptive GPU batch sizing, a guaranteed recovery from every OOM, or an automatic retry of failed checkpoint rows. A changed input, model/runtime hash or request cannot reuse a prior checkpoint as the same calculation. CPU, elapsed-time and per-run input bounds still apply. GPU peak allocation is reported by PyTorch for BioCentric; it is not a measurement of all GPU memory used by other processes.

BioCentric's native bridge and application now retain the original reviewed request and output directory when resuming an interrupted calculation. Before reuse, preflight rechecks the registered plan, input/model/runtime hashes, request bytes and individual checkpoint rows. Completed checkpoints remain unchanged. Original logs/results are archived as immutable artifacts. Existing failed pairs remain failed by default; a separate durable review can authorize retrying failed and unfinished pairs. Tests reject altered inputs, weights, plans or checkpoint evidence before model invocation. The [actual 16-pair GPU test](../runs/linux-biocentric-research-resume-v1/verification.json) cancelled after two completed pairs, reused both unchanged and computed the remaining 14, with zero final failures.
