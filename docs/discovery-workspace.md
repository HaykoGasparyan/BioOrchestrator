# Discovery workspace

The root Streamlit page is a chat and evidence workspace. Saved project URLs and `?view=workbench` retain the existing execution interface.

Docking and scoring selections are candidate pools. Each selected tool retains independent parameters; changing the pool or parameters invalidates its acceptance. The combinations tab exports a Cartesian candidate plan, including method-specific prerequisites. It does not execute a benchmark or claim compatibility, installed availability, or predictive validation. Custom plugins are session configuration entries, not executable adapters.

After an explicit target and organism, the workspace resolves UniProt identity. Ambiguous matches require selection. It retrieves accession-mapped experimental RCSB entries, organism-matched ChEMBL single-protein binding assays (Kd, Ki, IC50), and a bounded Europe PMC literature search. Structure and assay pagination expose the complete matched source result window; literature coverage is not exhaustive. Activity values, inequality qualifiers, quality flags, assay descriptions and publication links are retained. Exported measurement pages require curation before use as a held-out test set.

Source responses and coordinate hashes are cached under `BIOCENTRIC_HOME/discovery-cache`. Structures use deposited coordinates, without generated ligand poses. The offline 3Dmol viewer supports component focus and model selection. NMR ensembles use PDB MODEL records because the bundled mmCIF parser merges models; unavailable or oversized files show an explicit preview limitation. Displayed models are not automatically distinct biological conformations. Ligands, additives, cofactors and non-target polymer partners require binding-site review.

Implementation: `discovery_studio.py`, `discovery_panels.py`, `discovery_evidence.py` and `discovery_studio.css`.

Validation: `python -m pytest tests/test_discovery_studio.py tests/test_discovery_evidence.py tests/test_research_ui.py tests/test_ui.py -q`.

API references: [RCSB Data](https://data.rcsb.org/), [RCSB Search](https://search.rcsb.org/), [ChEMBL services](https://chembl.gitbook.io/chembl-interface-documentation/web-services/chembl-data-web-services).
