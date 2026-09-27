# Optional module backlog

1. Target/construct alignment and pocket-level confidence review, biological assembly selection, predicted-model PAE/pLDDT inspection, explicit state/protonation ensembles.
2. Human-reviewed PDBFixer/PDB2PQR/Meeko repair diffs in the UI. General receptor preparation must remain distinct from file conversion.
3. Actual curated target-specific active/inactive benchmark datasets, predeclared criteria, scaffold/temporal splits, cross-docking and ligand-similarity baselines. Integrate an audited evaluation into the gate only after dataset provenance and protocol bindings are verified.
4. Locally validated comparison with GNINA or another engine; separate metrics, model weights and terms. Do not average incompatible raw scores.
5. ADMET-AI or target-specific QSAR with endpoint-specific domains, pinned weights and applicability/uncertainty. No default favorable assessment for missing data.
6. AlphaFold DB retrieval; Boltz-2 or other structure/complex predictors only through explicit typed, licensed, hardware-bounded adapters and target-specific evaluation.
7. OpenMM/GROMACS, endpoint energies and OpenFE for explicitly motivated questions, with replicated sampling and force-field coverage checks.
8. ELN/LIMS import/export and richer measured-result comparison, batches, replicates and next-cycle project links.
9. Multiuser authentication, remote worker deployment and resource isolation for a demonstrated need. The current single-user worker is not a container security boundary.
10. Live LLM evaluations for support/contradiction, endpoint interpretation and correct tool choice. Existing tests establish programmatic boundaries, not model scientific accuracy.
