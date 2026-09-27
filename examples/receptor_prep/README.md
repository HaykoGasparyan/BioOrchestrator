# Protein preparation adapter verification

`1crn.pdb` is the unchanged public [PDB 1CRN](https://www.rcsb.org/structure/1CRN) entry downloaded from `https://files.rcsb.org/download/1CRN.pdb` on 2026-09-26. [wwPDB usage policy](https://www.wwpdb.org/about/usage-policy) applies. It is a small standard protein used solely for verifying protein-preparation software; it is not a target-specific screening example.

PDB2PQR 3.7.1/PROPKA 3.5.1, PDBFixer 1.12.0/OpenMM 8.3.1 and Meeko 0.7.1 were executed locally. pH 7.4 is an explicit test assumption. PDBFixer did not fill missing loops or replace residues. PDB2PQR's generated `protonated.pdb` was supplied to Meeko, which generated receptor PDBQT and its JSON record. Outputs/logs are in `runs/adapter-checks/1crn-*`; they are proposed prepared structures requiring inspection, not validated docking receptors.

The official Vina 1IEP tutorial's PDB contains an interrupted SER A 438 record. PDB2PQR and Meeko receptor preparation rejected it; those honest failures are retained in `runs/adapter-checks/pdb2pqr/` and `runs/adapter-checks/meeko/`. The provided curated 1IEP PDBQT remains the verified docking input. The application never repairs/reorders that structure silently.
