# BioCentricGPT: discovery workflow, tool catalog, and recommended first release

Research date: 26 September 2026. Scope: target-based, early small-molecule discovery. This is a broad, curated tool landscape, not a claim to list every available program or to have benchmarked them. Official documentation, primary repositories, and research publications are linked beside the relevant statements. Access, licenses, model weights, service terms, and APIs must be checked again when an adapter is implemented.

## What your initial outline covers—and what to add

Your outline captures structure-based virtual screening: investigate a protein, identify chemical starting points, obtain a structure, choose tools, dock, score, and select compounds. That is a useful core for BioCentricGPT. It is one part of discovery, and a computed shortlist remains a set of experimental hypotheses.

Protein-first discovery is one route. Phenotypic discovery can begin with a desired cellular effect and identify the target later. Target-based projects also need a disease hypothesis, desired direction of modulation, modality choice, and evidence that affecting the target will help. The [NCATS discovery overview](https://toolkit.ncats.nih.gov/module/getting-started/understand-the-rd-process/discovery-of-therapeutic-approach/) places target validation, therapeutic approach, screening, and lead selection in this broader context.

| Stage | What BioCentricGPT should help establish | What was missing or implicit in the original outline |
|---|---|---|
| 1. Biological objective | Disease/experimental question, relevant population or model, desired effect, exposure constraints, success criteria | Why this intervention should be useful |
| 2. Target identity and validation | Accession, species, isoform, mutations, causal/functional evidence, tissue/pathway context, liabilities | Whether this is the right target and direction of modulation |
| 3. Modality and tractability | Small molecule versus another approach; inhibitor/agonist/allosteric mechanism; site and assay feasibility | Whether ordinary docking is appropriate |
| 4. Evidence and known chemistry | Literature, patents, measured activity, negative findings, clinical precedents, existing structures | Assay comparability, source provenance, duplicate experiments |
| 5. Assay and control strategy | Primary functional assay, orthogonal confirmation, interference counterscreens, cellular measurements | How computational hypotheses will be tested |
| 6. Structural strategy | Relevant experimental state, template/prediction if needed, site-level uncertainty | Structure quality and biological relevance, beyond simply having coordinates |
| 7. Compound library | Known actives, diverse purchasable compounds, fragments, similarity search, optional designs | Availability, identity, synthesis feasibility, negative controls |
| 8. Preparation | Protonation, tautomers, stereochemistry, cofactors, waters, missing atoms, supported chemistry | A major source of downstream error |
| 9. Protocol validation | Pose recovery, physical validity, known active/inactive enrichment, held-out evaluation | Evidence for choosing a tool and its parameters |
| 10. Screening and prioritization | Docking and/or ligand-based methods; selected rescoring; diversity and uncertainty | Multiple criteria instead of one score |
| 11. Hit confirmation | Compound QC, repeat concentration–response, orthogonal binding/function, selectivity, cellular engagement | What turns a prediction into experimental evidence |
| 12. Hit-to-lead and optimization | Iterative design–make–test–analyze, measured SAR, ADME/PK, efficacy and safety | Feedback loop and developability |
| 13. Development | Broader preclinical evaluation, formulation/manufacturing, clinical development | Beyond the proposed first product |

This table is a planning synthesis, not a mandatory universal sequence. Assay development and biology continue alongside chemistry; a failed experiment can change the target hypothesis. The [Assay Guidance Manual](https://www.ncbi.nlm.nih.gov/books/NBK53196/) provides operational guidance for assay development, screening, and hit characterization.

## Broad tool catalog by purpose

“Core” means recommended for the first release. “Optional” means add when a project requires it. A hosted website, a local program, a data API, and a laboratory method are different integration types. The presence of a website does not imply permission or an interface for automation.

### Evidence, biology, and chemical data

| Purpose | Tools / sources | Suggested use and limits |
|---|---|---|
| Protein identity and sequence | [UniProt](https://www.uniprot.org/help/programmatic_access) | **Core API.** Resolve species/isoform/sequence and retain construct/mutation metadata. |
| Target–disease evidence | [Open Targets](https://platform-docs.opentargets.org/data-access/graphql-api), [GWAS Catalog](https://www.ebi.ac.uk/gwas/) | **Core Open Targets; optional genetics expansion.** Associations support a hypothesis but do not automatically establish causality. |
| Tissue/cell context | [Human Protein Atlas](https://www.proteinatlas.org/about/download), [DepMap](https://depmap.org/portal/) | **Optional datasets.** DepMap is particularly relevant to cancer dependencies; results are context-specific. |
| Pathways and networks | [Reactome](https://reactome.org/dev/content-service), [STRING API](https://string-db.org/help/api/), [STRING MCP](https://string-db.org/help/mcp/) | **Optional.** Functional associations are not necessarily physical interactions. Preserve evidence type. |
| Literature | [Europe PMC REST](https://europepmc.org/RestfulWebService), [PubMed/NCBI E-utilities](https://www.ncbi.nlm.nih.gov/books/NBK25497/) | **Core Europe PMC; optional NCBI adapter.** Retain queries, dates, PMID/DOI, article status and supporting passages. Full text rights vary. |
| Patents and clinical landscape | [SureChEMBL](https://www.ebi.ac.uk/about/news/updates-from-data-resources/surechembl-2-0-is-here/), [ClinicalTrials.gov API information](https://www.nlm.nih.gov/pubs/techbull/ma24/ma24_clinicaltrials_api.html) | **Optional.** Verify current programmatic access. Patent matches are research leads, not a freedom-to-operate determination. |
| Known ligands and bioactivity | [ChEMBL API](https://www.ebi.ac.uk/chembl/api/data/docs), [BindingDB REST](https://ww.bindingdb.org/rwd/bind/BindingDBRESTfulAPI.jsp), [PubChem PUG REST](https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest) | **Core ChEMBL/PubChem; BindingDB as needed.** Keep experimental endpoints, assay context, inequalities and sources. Missing activity is not inactivity. |
| Chemical entity mappings | [ChEBI](https://www.ebi.ac.uk/chebi/tools), PubChem, RDKit | **Optional ontology; core identifier handling.** Identical names do not guarantee identical stereochemistry, salt, or state. |
| Purchasable chemical space | [ZINC22](https://wiki.docking.org/index.php/ZINC22%3AGetting_started), [Enamine REAL](https://enamine.net/compound-collections/real-compounds), uploaded vendor catalogs | **Start with file import.** Distinguish in-stock from make-on-demand and recheck availability before selection. |

ChEMBL's documentation distinguishes assay types, target confidence, activity validity, and potential duplicate records. BioCentricGPT should preserve those fields, and should not count the same experiment imported into several databases as independent corroboration. [ChEMBL data documentation](https://chembl.gitbook.io/chembl-interface-documentation/frequently-asked-questions/chembl-data-questions)

### Structures, preparation, and ligand-based alternatives

| Purpose | Tools / sources | Suggested use and limits |
|---|---|---|
| Experimental structures and quality | [RCSB PDB](https://www.rcsb.org/docs/general-help/assessing-the-quality-of-3d-structures), [RCSB ligand quality](https://www.rcsb.org/docs/general-help/ligand-structure-quality-in-pdb-structures), PDBe/wwPDB services | **Core retrieval and review.** Inspect local pocket quality, assembly, construct and state; global resolution alone is insufficient. |
| Existing predicted structures | [AlphaFold DB](https://alphafold.ebi.ac.uk/faq) | **Optional early adapter.** Assess pocket pLDDT and relative-domain uncertainty using PAE where provided. These are structural confidence measures. |
| Template-based modeling | [SWISS-MODEL](https://swissmodel.expasy.org/), [MODELLER](https://salilab.org/modeller/) | **Optional.** Relevant templates, alignment and state matter. Check commercial permissions separately. |
| Protein prediction and complex prediction | [ColabFold](https://github.com/sokrypton/ColabFold), [Boltz](https://github.com/jwohlwend/boltz), [Chai-1](https://github.com/chaidiscovery/chai-lab), [AlphaFold 3](https://github.com/google-deepmind/alphafold3) | **Later compute adapters.** Separate monomer prediction, complex pose prediction, and affinity models. A plausible predicted complex does not prove binding. |
| Receptor repair/protonation/preparation | [PDBFixer](https://github.com/openmm/pdbfixer), [PDB2PQR/PROPKA](https://pdb2pqr.readthedocs.io/en/latest/), [Meeko](https://github.com/forlilab/Meeko) | **Core supported path.** Keep repairs inspectable and preserve necessary components. These programs perform complementary tasks; no single conversion solves all preparation. |
| Ligand standardization and conformers | [RDKit](https://www.rdkit.org/docs/GettingStartedInPython.html), [Dimorphite-DL](https://github.com/durrantlab/dimorphite_dl), [Gypsum-DL](https://github.com/durrantlab/gypsum_dl), Meeko | **Core RDKit/Meeko; optional enumeration adapters.** Track parent → stereoisomer → tautomer/protonation state → conformer. |
| Binding-site hypotheses | [fpocket](https://github.com/Discngine/fpocket), [P2Rank](https://github.com/rdk/p2rank) | **Optional.** Combine predicted cavities with known ligand or functional evidence. A cavity is not automatically a useful drug-binding site. |
| Similarity and chemical diversity | RDKit fingerprints, substructure search and clustering | **Core baseline.** Useful with known actives and to avoid a shortlist containing near-duplicate compounds. |
| Pharmacophore screening | [Pharmit](https://pharmit.csb.pitt.edu/help.html), [Phase](https://www.schrodinger.com/platform/products/phase/) | **Optional; Phase commercial.** A ligand-based route when structures are limited, or complementary evidence alongside docking. Validate the pharmacophore. |
| Project-specific QSAR | [Chemprop](https://github.com/chemprop/chemprop), [DeepChem/MoleculeNet](https://deepchem.io/tutorials/an-introduction-to-moleculenet/), [TDC benchmarks](https://tdcommons.ai/benchmark/admet_group/overview/) | **Later.** Need sufficient curated data and useful held-out performance versus simple baselines. Framework availability is not model validity. |

### Docking, scoring, simulation, and validation

| Purpose | Tools / sources | Suggested use and limits |
|---|---|---|
| Baseline docking | [AutoDock Vina](https://autodock-vina.readthedocs.io/en/latest/docking_basic.html) + Meeko | **Core.** A practical reproducible baseline; validate for the target and supported chemistry. |
| Additional docking/comparison | [GNINA](https://github.com/gnina/gnina), [UCSF DOCK 6](https://dock.compbio.ucsf.edu/DOCK_6/), [RxDock](https://github.com/rxdock/rxdock/blob/master/docs/user-guide/docking-strategies.rst) | **Optional.** Add an engine when a controlled comparison justifies the maintenance and compute. |
| Commercial docking suites | [Glide](https://www.schrodinger.com/platform/products/glide/), [GOLD](https://www.ccdc.cam.ac.uk/solutions/software/gold/), [OpenEye FRED/HYBRID/POSIT](https://docs.eyesopen.com/applications/oedocking/tutorial/overview.html), [FlexX/HYDE](https://www.biosolveit.de/download/software-downloads/), [MOE](https://www.chemcomp.com/en/Products.htm) | **License-dependent alternatives.** Capabilities for constraints, flexibility, waters, metals and special chemistry differ. Do not assume identical coverage. |
| Learned pose generation | [DiffDock](https://github.com/gcorso/DiffDock), Boltz, Chai-1, AlphaFold 3 complex prediction | **Later comparisons.** DiffDock pose confidence is not a binding-affinity prediction. Evaluate chemical plausibility and novel-target generalization. |
| Pose validity | [PoseBusters](https://github.com/maabuu/posebusters), RDKit chemistry/geometry checks | **Core.** Check geometry, stereochemistry and clashes. Pose recovery and physical plausibility are separate checks. |
| Interaction interpretation | [PLIP](https://plip-tool.biotec.tu-dresden.de/plip-web/plip/index), [ProLIF](https://github.com/chemosim-lab/ProLIF) | **Core PLIP; ProLIF optional.** Contact patterns aid inspection; more hydrogen bonds do not automatically mean stronger binding. |
| Learned rescoring / affinity | GNINA, [Boltz-2](https://github.com/jwohlwend/boltz/blob/main/docs/prediction.md), calibrated target-specific models | **Optional after validation.** Store each metric's exact meaning, scale and direction; keep raw scores separate. |
| MD simulation | [OpenMM](https://docs.openmm.org/latest/userguide/), [GROMACS](https://manual.gromacs.org/current/), [AmberTools](https://ambermd.org/AmberTools.php) | **Later for a stated question.** Need appropriate system setup, force-field coverage, sampling and replication. Short pose stability does not establish affinity. |
| Ligand force fields / parameters | [Open Force Field](https://openforcefield.org/), AmberTools | **Simulation dependencies.** Unsupported chemistry needs specialist parameterization and checks. |
| Trajectory analysis / enhanced sampling | [MDAnalysis](https://www.mdanalysis.org/), [MDTraj](https://mdtraj.org/), ProLIF, [PLUMED](https://www.plumed.org/) | **Later.** Replicates and convergence matter; enhanced sampling requires a justified method. |
| Endpoint energy estimates | AmberTools MMPBSA.py, [gmx_MMPBSA](https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA) | **Optional approximate rescoring.** Setup-sensitive and not automatically more reliable than the original screening method. |
| Alchemical free energies | [OpenFE](https://docs.openfree.energy/en/latest/tutorials/index.html), [FEP+](https://www.schrodinger.com/platform/products/fep/) | **Advanced small-set work.** Distinguish relative and absolute calculations; validate preparation, sampling and applicability to the series. |

The [PoseBusters research paper](https://pubs.rsc.org/se/content/articlelanding/2024/sc/d3sc04185a) demonstrates why low pose RMSD alone can conceal physically implausible predictions. This is a reason to validate geometry, not a claim that its 2024 ranking of methods establishes the best model today.

For screening validation, useful inputs include target-specific curated assay sets and public benchmarks such as [LIT-PCBA](https://pubs.acs.org/doi/10.1021/acs.jcim.0c00155). Benchmark names do not guarantee independence: audit duplication, scaffold overlap, assay semantics and pretrained-model overlap; a [2025 LIT-PCBA audit](https://arxiv.org/abs/2507.21404) specifically raises leakage and redundancy concerns. Label artificial decoys separately from measured inactives.

### Developability, design, experiments, and operations

| Purpose | Tools / sources | Suggested use and limits |
|---|---|---|
| Physicochemical triage | RDKit descriptors, property distributions and structural alerts | **Core.** Molecular weight, lipophilicity estimates and polar surface area are descriptors, not measured ADMET or universal pass/fail rules. |
| ADMET prediction | [ADMET-AI](https://github.com/swansonk14/admet_ai), [SwissADME](https://www.swissadme.ch/), [ADMETlab 3.0](https://admetlab3.scbdd.com/apis/), [ADME@NCATS](https://github.com/ncats/ncats-adme/) | **Optional.** Prefer a supported local model for private chemistry. Record endpoint, applicability, model version and uncertainty where available. |
| Selectivity / off-target hypotheses | ChEMBL/BindingDB cross-target data; [SwissTargetPrediction](https://swisstargetprediction.ch/help.php) | **Evidence first.** Predicted off-targets suggest counterscreens; absent predictions do not prove safety. |
| Analog/SAR analysis | [mmpdb](https://github.com/rdkit/mmpdb), [SwissBioisostere](https://www.swissbioisostere.ch/) | **Later.** Interpret transformations in the context of comparable measured assays. |
| Generative chemistry | [REINVENT4](https://github.com/MolecularAI/REINVENT4) | **Later optional design module.** Validate objectives and synthesis feasibility; optimizing an imperfect score can exploit its errors. |
| Retrosynthesis | [AiZynthFinder](https://github.com/MolecularAI/aizynthfinder), [ASKCOS](https://askcos-docs.mit.edu/) | **Later.** Proposed routes need chemist review, appropriate stock catalogs, and separately verified model/data terms. |
| Assay planning and interference | [NIH Assay Guidance Manual](https://www.ncbi.nlm.nih.gov/books/NBK53196/), [chemical reactivity guidance](https://www.ncbi.nlm.nih.gov/books/NBK326709/), [aggregation guidance](https://ncbi.nlm.nih.gov/books/NBK442297/?report=reader) | **Core planning sources.** Select mechanism-appropriate confirmation and counterscreens; software cannot execute physical assays by itself. |
| Direct binding and functional confirmation | SPR, BLI, ITC, MST, DSF, NMR; biochemical/functional assays | **Experimental methods.** Choice depends on target, affinity range and artifacts. Thermal shifts alone do not establish affinity. [Biacore small-molecule guide](https://cdn.cytivalifesciences.com/api/public/content/digi-33039-original) |
| Cellular target engagement | [CETSA guidance](https://www.ncbi.nlm.nih.gov/books/NBK374282/pdf/Bookshelf_NBK374282.pdf), [NanoBRET](https://ireland.promega.com/resources/technologies/nanoluc-luciferase-enzyme/cellular-target-engagement/) | **Experimental methods.** Separate target engagement from downstream functional efficacy. |
| Laboratory records | [eLabFTW API](https://doc.elabftw.net/docs/usage/api/), [Benchling developer platform](https://help.benchling.com/hc/en-us/articles/9714802977805-Access-Benchling-s-Developer-Platform-and-Documentation-Base) | **Optional ELN/LIMS adapters.** Preserve sample/batch identity, raw measurements, assay metadata, replicates and permissions. |
| Run and dataset tracking | SQLite + immutable manifests initially; [MLflow](https://www.mlflow.org/docs/latest/api_reference/rest-api.html), [DVC](https://doc.dvc.org/) later | **Core provenance; optional platforms.** Persist provenance independently of chat history. |

Early experimental developability work should assess relevant solubility, permeability, metabolic stability, protein binding and safety liabilities, then PK/exposure and efficacy in suitable models as the program advances. These are project-dependent measurement programs. A universal prediction dashboard cannot replace them or guarantee clinical success. The broader progression is covered by [NCATS](https://toolkit.ncats.nih.gov/module/getting-started/understand-the-rd-process/discovery-of-therapeutic-approach/).

## How BioCentricGPT should choose a structure and a protocol

| Situation | Recommended branch |
|---|---|
| Relevant high-quality experimental complex exists | Inspect state and site, retain necessary biology, prepare carefully, then benchmark a docking baseline. |
| Several structures or states exist | Compare biologically plausible states; use an ensemble only when it answers a defined question and improves validation. |
| No suitable experimental structure, good template available | Consider a template-based model with documented alignment and uncertain regions; validate pocket suitability. |
| No suitable template, credible predicted pocket available | Use a predicted structure as an explicit hypothesis; evaluate site confidence, interdomain uncertainty and state dependence. |
| Binding pocket/state is unreliable but known ligands exist | Consider similarity, pharmacophore or validated QSAR approaches; propose evidence that would resolve structural uncertainty. |
| Neither reliable structure nor known chemistry exists | Produce a research/experimental strategy and an exploratory designation; do not invent a validated screening pipeline. |
| Covalent chemistry, important metals, membrane context, unusual ligands or cryptic sites | Require a dedicated supported protocol and specialist assessment. |

There are no universally best docking parameters for a protein. Start from the site's geometry and chemical requirements, define a small comparison, and select settings using controls. Search effort, receptor conformation, protonation, waters and constraints can matter differently. Vina's [FAQ](https://autodock-vina.readthedocs.io/en/latest/faq.html) explains search/protonation issues; its [manual](https://vina.scripps.edu/manual/) recommends target-specific checks with known actives or a bound ligand.

Keep tuning and final evaluation separate. A useful pilot might compare a few defensible receptor/preparation choices and search settings. This is a design suggestion, not a validated numeric recipe. Judge pose recovery, physical validity, enrichment at the experimental testing budget, chemical diversity and uncertainty. If the pilot fails, revise the hypothesis or strategy before scaling.

## What “binding affinity tools” actually mean

| Output category | Examples | Interpretation |
|---|---|---|
| Docking score | Vina, Glide and other scoring functions | Approximate ranking signal within a specific protocol; not a measured Kd. Scores from different functions may not be comparable. |
| Pose confidence | DiffDock confidence, complex structural confidence | Model confidence in a proposed geometry, not proof that binding occurs. |
| Learned bioactivity/affinity estimate | Boltz-2 or a target-specific predictor | A prediction tied to a defined endpoint and training domain; must be validated for intended use. |
| Endpoint energy estimate | MM/GBSA, MM/PBSA | Approximate, setup-sensitive energetic estimate. |
| Alchemical free energy | OpenFE/FEP workflows | A physics-based estimate under an explicit model/protocol, with sampling and parameterization uncertainty. |
| Experimental observation | Kd/Ki or an assay-specific IC50/EC50 | A measured endpoint with assay context; these endpoint types are not interchangeable. |

Vina documents that AutoDock and Vina energy scores cannot be compared directly. [Vina docking documentation](https://autodock-vina.readthedocs.io/en/latest/docking_basic.html)

Boltz-2 has separate outputs for predicted binder probability and `log10(IC50 in μM)`. Its continuous value is intended for comparing active binders during optimization; the probability output serves a different screening task. Do not translate either into measured affinity. If `y = log10(IC50 in μM)`, then `pIC50 = 6 − y` is dimensionless; it is not an energy in kcal/mol, and IC50 is not generally Kd. [Official Boltz repository](https://github.com/jwohlwend/boltz)

Free-energy work is often most useful for selected optimization decisions after establishing credible structures and chemistry. Preparation, sampling and the suitability of the chemical series still determine reliability. [Best practices for alchemical free-energy calculations](https://arxiv.org/abs/2008.03067)

## Recommended compact first release

This is an engineering recommendation, not a claim of universally optimal architecture.

Build one model-driven research assistant around a deterministic scientific workflow. The LLM researches, proposes and explains; typed tool adapters retrieve evidence and run programs; application code validates inputs, enforces gates, computes metrics and records results. Strict tool schemas help constrain calls, but do not establish scientific truth. [OpenAI function calling](https://developers.openai.com/api/docs/guides/function-calling), [Structured Outputs limitations](https://developers.openai.com/api/docs/guides/structured-outputs)

Use the Responses API with one compact agent loop, a configurable model and a versioned runtime prompt. A small Python core, CLI, simple interface, SQLite metadata, artifact files and one worker process are adequate starting choices. OpenAI's [runtime comparison](https://developers.openai.com/api/docs/guides/agents) distinguishes direct Responses integrations from SDK-managed and hosted agent runtimes.

**First milestone:** resolve the target and produce a cited dossier using UniProt, Europe PMC, Open Targets, RCSB PDB, ChEMBL and PubChem; accept private structures and ligand files.

**Second milestone:** implement one supported preparation and docking path with RDKit, Meeko, Vina, PoseBusters and PLIP. Offer curated prepared inputs first, then inspectable PDBFixer/PROPKA/PDB2PQR support. Persist benchmark gates, failures and run manifests.

**Third milestone:** produce a diverse experimental shortlist, explanatory report, downloadable compounds/poses, and laboratory-result import. Start with descriptors and evidence-based flags. Add one predictive ADMET adapter only when its limitations and data handling are explicit.

**Later modules:** GNINA comparison, alternate docking engines, structure prediction, Boltz-2, advanced ADMET/QSAR, MD, free energies, retrosynthesis and generative chemistry. Each needs an actual adapter, documented applicability and tests; a list of links is not an implemented integration.

The key output is an auditable recommendation: “These compounds are worth testing, for these reasons, under these assumptions.” The application should be able to say that the evidence is insufficient or that the benchmark failed.

## Integration details worth checking before commercialization

- **Name:** Use BioCentricGPT consistently throughout the product and its documentation. Keep the displayed name configurable.
- **Separate code, weights and service terms:** [Boltz](https://github.com/jwohlwend/boltz) lists MIT code/model terms; [Chai-1](https://github.com/chaidiscovery/chai-lab) lists Apache 2.0. AlphaFold 3's code and [weights/output terms](https://github.com/google-deepmind/alphafold3/blob/main/WEIGHTS_TERMS_OF_USE.md) are separate; the freely distributed weights have noncommercial restrictions. The official repository links a separate commercial route. Do not infer permissions for data or external services from a code license.
- **Remote-service suitability:** [SwissADME](https://www.swissadme.ch/about.php) gives no confidentiality guarantee. Its public website should not become a silent upload destination for proprietary compounds. [ADMETlab's API page](https://admetlab3.scbdd.com/apis/) currently warns of service instability; keep it optional.
- **Version predictive models:** [ADMET-AI](https://github.com/swansonk14/admet_ai) currently distinguishes its newer repository models from the version used by its original paper/live server. Pin exact model artifacts, not just a familiar tool name.
- **Benchmark and commercial claims:** A provider's benchmark does not establish performance on your target. Audit overlap and validate locally before promoting an optional tool into the default workflow.

## How to use the companion prompts

Use `BioCentricGPT-Codex-build-prompt.md` as the instruction to Codex in your application repository. It requests an executable product with explicit scientific and software acceptance criteria. Use `BioCentricGPT-runtime-system-prompt.md` as the starting behavior specification for the agent inside that product. Keep scientific enforcement in application code as well as the prompt.

The prompts are designed for a focused first release. They do not instruct Codex to implement every tool in this catalog at once.
