# Ligand-first Bayesian pipeline search

New evaluation plans freeze a `search` configuration. Old plans and saved results retain their original execution protocol. The current search domain is the selected, executable pocket × docking × scoring configurations; it does not invent plugins or tune continuous hyperparameters.

Each search round chooses a small candidate set. It executes **crystal control → each measured ligand**, trying all candidates in that round before advancing to the next ligand. Every selected candidate receives the identical frozen measurement set and control. Repeated assay records remain distinct. Immutable docking outputs are reused when pocket, docking settings and chemistry match. BioCentric scores all poses of the current ligand together, rather than waiting for a multi-ligand chunk. Its configured compound batch cap is an upper bound; ligand-first execution uses one ligand group per call (plus the crystal pose when needed).

## Surrogate and acquisition

The initial design uses up to three deterministic, diverse categorical configurations (seed 42). After complete rounds, an exact Gaussian process models a scalar objective over categorical configurations:

- Distance is the fraction of differing pocket/docking/scoring configuration blocks, including frozen parameters.
- Kernel: `k(x,x') = exp(-d(x,x') / 0.7)`; prior mean zero, observation noise variance 0.0025.
- Objective is mean within-assay Spearman rho, with scorer direction handled separately. In v3, failed ligands are skipped. Chains are compared on the intersection of their successful measurements, retaining within-assay grouping. All earlier observations are recomputed on this shared subset after each round. A chain failing checks on its successfully scored poses, lacking finite shared-set rank metrics, or failing crystal recovery gets a **surrogate feasibility penalty −1.25**, not a replacement affinity score. Its original failures and metrics remain visible.
- Expected improvement: `EI = (μ − f_best − ξ) Φ(z) + σ φ(z)`, `z = (μ − f_best − ξ)/σ`, `ξ = 0.01`.
- Later rounds select up to two unused configurations by greedy batch EI with posterior-mean fantasies. Fantasy values are used only to select the second candidate; they are never persisted as observations or displayed as measurements.
- Default trial budget: `min(N, max(6, min(12, ceil(2*sqrt(N)))))`. Thus a pool of nine tries six; a pool of three tries all three. The minimum six allows a winner plus five alternatives when enough candidates exist.

This is a small implementation of GP Bayesian optimization / expected improvement, described in [Brochu, Cora and de Freitas (2010)](https://arxiv.org/abs/1012.2599). It does not establish a global optimum. Surrogate uncertainty is not calibrated binding confidence. Search trace, objective values, acquisition estimates, selected IDs and rounds are stored in `status.json`; runtime fingerprints bind the optimizer, runner and decision policy implementation.

## Completion and final decision

Progress counts complete evaluated tool chains against the frozen search budget. Candidates outside the budget are explicitly `not_selected`, remain gray, and have no calculated scores; interrupted selected work is `not_run`. No score is fabricated for missing evidence. Metrics are finalized only after every requested ligand is accounted for. Failed ligands stay explicit; zero-success chains do not shrink other comparisons. Each chain's coverage, shared-set count and correlation sample size are reported. Separate own-subset metrics are diagnostic and never rank chains against one another.

After search completes, the best available chain is fixed automatically. Shared-subset selection prioritizes finite correlations, passing checks on all successfully scored poses, and comparable crystal RMSD ≤2 Å; it ranks by shared-set Spearman, Kendall, coverage, then candidate ID. Individual ligand failures do not disqualify an otherwise comparable chain. Strict full-completion checks remain separate and visible. Saved pose hashes are verified before locking. If none passes, the fixed best available result retains its failed checks and unavailable metrics. Cancellation does not auto-lock a partial search.

The **Final chain** tab displays only the fixed path with green outlines/connections, its frozen parameters, Spearman/Kendall metrics, crystal RMSD and physical-pass count. **Final results** highlights the fixed winner in green and shows up to five evaluated alternatives, with all candidates and the Bayesian trace available in expanders. **Final visualizations** retains comparisons with the crystal reference.

All statistics describe the same evaluation set used to select the chain. They are **exploratory**, not unbiased held-out predictive validation. Learned-model training overlap also remains unknown. The implementation does not silently split small assay groups or claim a holdout where none was provided.

## Verification

`tests/test_bayesian_search.py` checks the GP posterior against direct linear algebra, acquisition selection, deterministic diversity, ligand-first ordering, six-of-nine budget accounting, cancellation, automatic pose-bound locking, skipped candidates, and rendering of the fixed chain and top alternatives. A small real Vina/BioCentric check is recorded under `runs/bayesian-search-verification`; it tests software execution, not predictive performance.

## GPU execution and reusable artifacts

New GNINA docking/pose-scoring settings default to GPU for its neural scoring;
its docking search retains four CPU threads. Vina remains CPU-only (four threads
by default), Uni-Dock uses its GPU engine, and BioCentric retains CUDA. An explicit
CPU setting is respected, and old frozen plans retain their original settings.
No unsupported GPU fallback substitutes another tool.

Workers share `.biocentric/evaluation-archive/`. Pocket identification is computed
once for the receptor and configured method (including the crystal ligand for a
reference-derived box). Prepared ligands are shared across docking methods.
Docked poses are shared across scoring chains and future runs only when receptor,
prepared receptor, ligand chemistry/coordinates, pocket, docking settings and
runtime signatures match. GNINA CNN pose and affinity use both outputs of one
fixed-pose scoring call. The archive hashes artifacts before reuse, publishes
successful stages atomically under a per-key process lock, and copies artifacts
into each run so the run stays self-contained. Failed/incomplete stages are not
published. It does not import unverified outputs from experiments made before
this feature. Runtime/model file hashes are reused within a worker only while
file identity, size, modification and change timestamps remain unchanged;
recently modified files are always read again.

The displayed Cartesian count includes unavailable selections. For example,
fpocket + DeepPocket, three docking engines and three scorers is **18 selected**,
**9 executable** (DeepPocket has no adapter), and **6 Bayesian trials**. With two
executable pocket methods, all 18 become candidates, with a separate search
budget. Progress counts complete chains on the entire evaluation set, never
individual compounds. Skipped/unavailable chains receive no fabricated scores.

Native verification: `runs/gpu-archive-verification/20260927-014425/verification.json`
records one 1IEP ligand through GNINA GPU docking and both CNN outputs. A fresh
archive took 11.88 s and matching reuse took 0.009 s with identical scores. These
are stage timings from a software check, excluding full-set validation and
initial runtime fingerprinting; they are not an overall benchmark speedup claim.


Final decisions now distinguish passing the strict validation gates from a best-available selection. A terminal run with trustworthy attempted results always retains the best tested candidate, its coverage and available shared-subset metrics. If all chains fail physical or crystal checks, the selected chain remains fixed and green in the UI, while failed checks remain explicit and its evidence color remains yellow. Missing/constant measurements produce unavailable correlations, never fabricated values. A failed crystal control does not hide correlations computed from a complete, identity-checked measured test set. Existing decisions and result files remain immutable.

Search v1/v2 plans retain their frozen objective. New v3 plans use the shared-subset objective and refresh observations each round. Saved decisions remain immutable. Conditioning on scoring success can introduce missingness bias; these metrics are exploratory, and perfect rho on four measurements is not evidence of general predictive performance.


## Correlation-first selection (v4)

At the user's request, new final decisions rank by shared-set within-assay Spearman, then Kendall, measurement coverage and candidate ID. Physical and crystal-recovery checks remain separately visible; they do not override the correlation ranking. A selected chain with failed checks is explicitly an exploratory best-correlation result, not a validated chain.

New search v4 uses the same correlation-first objective (no physical feasibility penalty); v1–v3 retain their frozen objectives and archived traces. Explicitly revising an existing decision preserves its exact original in `decision-history/<sha256>.json`, records the revision reason and prior digest, verifies pose hashes, and leaves plans, measured evidence, scores and search history unchanged. Existing decisions are not automatically replaced.

The requested correction of TYK2 experiment `c05ba69abbb74d39b3e1324ce8dc9b80` selects crystal pocket → GNINA docking → BioCentric affinity: rho 0.236738 versus 0.176481 for the previously selected Vina → GNINA affinity chain. Both correlations use the same 125 rankable measurements across 8 assays within the 133-measurement shared successful set (134 requested). BioCentric has one failed measurement and five scored poses with physical-check failures; crystal RMSD is 1.215 Å. These checks remain reported after correction.
