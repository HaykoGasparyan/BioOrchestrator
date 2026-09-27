# Ubuntu GPU docking and scoring verification

The RTX 4090 / NVIDIA 535.309.01 host completed real Uni-Dock 1.1.0, GNINA 1.3.2 docking and GNINA CNN scoring. Native CLI output and the application's typed `engines.run_engine` path agreed exactly on the recorded one-ligand 1IEP software checks. These are execution checks, not independent affinity benchmarks or evidence of a new binder. Physical pose checks in the comparative workflow remain authoritative.

| Check | Native time | Observed result | Evidence directory |
|---|---:|---|---|
| Uni-Dock 1.1.0 GPU, fast preset | 2.115 s | One imatinib pose, empirical −9.316 kcal/mol; native/application PDBQT hashes identical | `runs/ubuntu-gpu-engine-verification/unidock-1iep/` |
| GNINA 1.3.2 GPU docking, exhaustiveness 1 | 8.255 s | One pose, empirical −13.33497 kcal/mol, CNN pose 0.9672928452, CNN affinity 8.385184288 pK; native/application SDF hashes identical | `runs/ubuntu-gpu-engine-verification/gnina-docking-1iep-nvrtc/` |
| GNINA 1.3.2 GPU score-only, observed reference geometry | 0.751 s | Empirical −12.51667 kcal/mol, CNN pose 0.9203819633, CNN affinity 8.3367071152 pK; native/application SDF hashes identical | `runs/ubuntu-gpu-engine-verification/gnina-affinity-final-runtime/` |

The docked GNINA pose produced a native `ligand outside box` message. A parsed score and a successful process exit do not override that observation, redocking RMSD, or physical checks. Learned pK and CNN pose score remain separate from empirical energy. `CNNvariance` is retained in native output; it is not presented as calibrated predictive uncertainty.

The initial GNINA docking attempt failed after search because `libnvrtc.so.12` was missing. The failed log and empty output remain under `gnina-docking-1iep/`. Adding the pinned NVIDIA CUDA runtime compiler resolved this real execution failure. A successful `--version` or score-only call had not exposed it.

GPU proof includes native Uni-Dock kernel timing and device-memory information; later GNINA verification records process-specific NVIDIA memory samples at approximately 0.1-second intervals. The GNINA launcher executes through the glibc loader, so NVIDIA reports `/lib64/ld-linux-x86-64.so.2` as the process name. Docking runs used up to 530 MiB in those sampled processes. Device-wide usage also includes unrelated desktop applications, so it is not attributed entirely to this calculation.

## Reproducible local installation

Run `bash scripts/setup_ubuntu.sh`, then `bash scripts/setup_gpu_runtimes.sh`. The latter downloads hash-pinned Linux x86_64 native binaries and installs the pinned packages in `requirements-gpu-runtime-ubuntu.txt` into `.tools/linux/cuda12-libs/`; NVRTC 12.1.105 occupies `.tools/linux/cuda12-nvrtc/`. It does not alter system CUDA, the GPU driver, the copied Mac environment, or the BioCentric model environment.

GNINA's [1.3.2 release](https://github.com/gnina/gnina/releases/tag/v1.3.2) requires cuDNN 9. Its [compatibility-build Dockerfile](https://raw.githubusercontent.com/gnina/gnina/v1.3.2/docker/compat/Dockerfile) builds against older CUDA and PyTorch 2.4.1. The [NVIDIA cuDNN 9.1.1 matrix](https://docs.nvidia.com/deeplearning/cudnn/backend/v9.1.1/reference/support-matrix.html) lists CUDA 12.1 and Linux driver ≥525.60.13 for this GPU generation; the [official installation guide](https://docs.nvidia.com/deeplearning/cudnn/backend/v9.1.1/installation/linux.html) documents versioned pip wheels. The verified launcher supplies these isolated libraries and inhibits the binary's old RPATH so `/usr/local/cuda` does not inject incompatible CUDA 11 libraries.

`.tools/linux/gpu-runtime-manifest.json` and the copied `runs/ubuntu-gpu-engine-verification/runtime-manifest.json` record package versions and SHA256s of the launcher, native executables and shared libraries. GNINA's default CNN weights are embedded in the hash-recorded native binary. The application binds these runtime files into its reviewed plan and rechecks them before and after execution. Installing files alone is not marked as a passed GPU test.

The application discovers `.tools/linux/gnina` and `.tools/linux/unidock`; explicit `BIOCENTRIC_GNINA` and `BIOCENTRIC_UNIDOCK` still override discovery. Launch the local app using the project Ubuntu launcher on loopback. The Codex restricted process namespace cannot access the NVIDIA driver, so GPU validation commands here needed approved execution outside that namespace.

## Uni-Dock 1.1.0 compatibility scope

The installed [1.1.0 release binary](https://github.com/dptech-corp/Uni-Dock/releases/tag/1.1.0) is a compatibility runtime, not the newest source release. This application limits that exact bundled binary to box sides below 30 Å, common supported PDBQT atom types and non-macrocyclic ligands with Meeko SMILES metadata. Unsupported atom types, macrocycle glue atoms, rings of eight or more atoms and absent identity metadata fail before execution. The [official Uni-Dock documentation](https://github.com/dptech-corp/Uni-Dock/blob/main/unidock/README.md) explains GPU requirements, search presets and GPU memory constraints. Do not infer reliability on unsupported chemistry or large boxes from this one-ligand smoke result. Actual per-target pose validation is still required before any full experiment.

## Rerunning the bounded checks

`scripts/verify_gpu_engines.py` preserves the bundled 1IEP reference coordinates, adds valence hydrogens without state enumeration, prepares Meeko PDBQT, and records its box and preparation assumptions in the new run's `inputs/` directory. Each run requires a fresh output directory, runs one native CLI calculation and one typed adapter calculation, saves raw logs/output hashes, and records NVIDIA samples. It strips application credentials from subprocesses and caps each scientific invocation at 300 seconds.

```sh
.venv-linux/bin/python scripts/verify_gpu_engines.py --engine unidock --output runs/new-unidock-smoke
.venv-linux/bin/python scripts/verify_gpu_engines.py --engine gnina --output runs/new-gnina-smoke
.venv-linux/bin/python scripts/verify_gpu_engines.py --engine gnina_affinity --output runs/new-gnina-affinity-smoke
```

AutoDock4 still requires compatible maps; the Amber methods still require parameterized matching topologies and sampled MD trajectories; PLANET requires its separate reviewed model installation. These methods were not silently enabled by the GPU setup.

## Fast discovery preset

New starts default to **Fast**, including starts from an older saved workspace.
Existing accepted/running plans retain their original frozen settings. In Docking
settings, Balanced restores the earlier search effort; editing a parameter switches
the selection to Custom. Presets are persisted with the experiment.

| Engine | Fast configuration |
|---|---|
| AutoDock Vina | exhaustiveness 2, 5,000 evaluations per search run, up to 4 CPU threads, 1 output pose |
| GNINA | exhaustiveness 2, 20 Monte Carlo steps, up to 4 CPU threads, GPU CNN scoring, 1 output pose |
| Uni-Dock | GPU `fast` search preset, 1 output pose |

Search parallelism is limited by independent search tasks; adding CPU threads beyond
that does not accelerate a single Vina search. Vina itself remains CPU-only. Fast
search sacrifices sampling accuracy and one output pose reduces downstream scoring
work. Validation gates, pose-identity checks, and measured labels are unchanged.
Native option names were checked against the installed binaries' `--help` and
GNINA `--help_hidden`. The typed adapters reject options intended for another engine.

The v2 Bayesian initialization covers distinct docking engines before repeating one,
then covers pocket methods and scorers. This fixes a max-min-distance tie that left
Uni-Dock out of the first round in the reported six-chain run. V1 plans retain their
original configuration validation and initializer. A pending Uni-Dock chain is not
an execution failure; failed native stages still report their actual error.

`runs/fast-docking-verification/20260927-015400/verification.json` records native
execution on the current 1O41 receptor and control ligand: Vina 0.484 s, Uni-Dock
1.449 s, GNINA 5.171 s. All returned one chemistry-verified pose, could be rescored,
and passed PoseBusters. GPU process samples include Uni-Dock and GNINA's loader.
These are one-ligand preparation/docking timings, not full-benchmark timings or
predictive accuracy claims. Uni-Dock 1.1.0's explicit chemistry/box limitations above
remain in effect.
