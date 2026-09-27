# Ubuntu continuation

Run `bash scripts/launch_ubuntu.sh` and open [the loopback app](http://127.0.0.1:8502). The new workflow freezes multiple candidate protocols before their real pilots, locks selection before held-out evaluation, and keeps independent validation explicitly unavailable for the unlabeled ABL example. Affinity predictions appear beside their own docking scores, with separate units, physical status, checkpoint provenance and CSV/JSON/SDF downloads. Historical Mac walkthrough details follow. See [Ubuntu setup](ubuntu-setup.md) and [comparative evaluation](comparative-evaluation.md).

# BioCentricGPT demo

Open the research conversation on port 8502. Existing completed conversations display their saved results without another model or docking call.

The visible sequence is **Research → Test set → Method → Pilot → Screen → Results**:

1. Describe the protein and research question, or prepare the public 1IEP example.
2. Research retrieves literature and available experimental evidence. The test set keeps compatible assays and molecular structures together, with separate tuning and evaluation compounds where the data support that split.
3. The model proposes a tool using the actual test set, documented method requirements and available runtimes. Its main explanation is one or two sentences.
4. Review the inputs in the conversation. The pilot measures pose recovery and physical validity. A reference redock alone does not establish binding-affinity accuracy, screening enrichment or superiority over another engine.
5. Review the pilot and continue the exploratory screen. The ligand list changes from queued/prepared to docking/scored/failed as actual results arrive.
6. Inspect **Ligands & results** and download the structures and score table. The stage cards remain visible alongside results: choose a card to see its short decision, evidence and test-set details, then expand **Exact recorded steps & calls** for the actual requests and outputs.

The workspace follows the current stage while work runs. Clicking a stage keeps that selection during refresh; **Follow current step** resumes automatic following. **Ligands & results** switches the detail panel back to the table. Inline questions stay in the conversation and **Send & continue** resumes the same job. Cards marked **AI** contain recorded model-call activity; **Workflow** cards describe scientific tools and application rules. The interface shows public explanations, not invented internal thoughts.

## Understanding the example

The completed public example contains imatinib, aspirin and caffeine. Imatinib is the cognate bound-ligand control; aspirin and caffeine are untested controls. All three can appear on the exploratory shortlist because the example asks for three diverse physically valid candidates. Being shortlisted does not make them confirmed binders.

For this Vina protocol, lower (more negative) kcal/mol scores rank better. Compare scores only within a compatible engine and protocol. GNINA CNN affinity and other affinity estimates retain their own labels and units; they are not interchangeable with Vina scores.

The current Ubuntu host has verified Vina, GPU GNINA, GPU Uni-Dock and GPU BioCentric execution. AutoDock4, PLANET and Amber methods still need their runtimes and required scientific inputs. An unavailable method is displayed honestly rather than simulated. If no independent assay dataset is available, the test set is a pose/software control and the screen remains exploratory.

The separate **Local validation · 128 diverse ligands** saved run contains the broader software test. It has only its four actual recorded stages, is labeled **no model calls**, and does not imply literature-guided selection or an assay benchmark. Its baseline Vina poses are ranked together; results from different Vina settings remain separate in the validation report. Inspect the [verification record](verification.md) for actual counts and limits.

## Downloads

- The **shortlist SDF** contains the recorded selected docked poses with chemical identity and coordinates.
- The **ranked CSV** contains the displayed ligand scores and selection status.
- **All docked poses** includes registered screening poses; unavailable or integrity-failing files are reported rather than silently substituted.
- An individual ligand's pose download lets you inspect that compound separately.

Before docking finishes, input/prepared files are distinct from docked outputs. Live scores remain provisional until the scientific stage has saved and registered its outputs.
