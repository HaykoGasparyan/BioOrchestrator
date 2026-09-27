import json
import sys
from pathlib import Path

import pytest

from biocentric.science import ScienceError, derive_box, prepare_library, prepare_receptor, run_process, _symmetry_rmsd
from biocentric.benchmark import dataset_hash, evaluate, freeze_protocol


def test_stereochemistry_duplicates_invalid_and_unsupported_are_explicit(tmp_path):
    pytest.importorskip("rdkit")
    source = tmp_path / "library.smi"
    source.write_text("C[C@H](O)C(=O)O left\nC[C@@H](O)C(=O)O right\nC[C@H](O)C(=O)O duplicate\nCC(O)C(=O)O ambiguous\nCCO.[Na+] salt\n[Fe] metal\nnot_a_smiles invalid\n")
    before = source.read_bytes()
    result = prepare_library(source, tmp_path / "prep")
    assert source.read_bytes() == before
    assert len(result["compounds"]) == 2
    assert len({row["canonical_smiles"] for row in result["compounds"]}) == 2
    assert {row["code"] for row in result["exclusions"]} == {"exact_duplicate", "ambiguous_stereochemistry", "unsupported_salt", "unsupported_element", "invalid_chemistry"}
    assert len(result["artifacts"]) == 3


def test_sdf_metadata_geometry_and_identity_preserved(tmp_path):
    Chem = pytest.importorskip("rdkit.Chem")
    AllChem = pytest.importorskip("rdkit.Chem.AllChem")
    mol = Chem.AddHs(Chem.MolFromSmiles("C[C@H](O)C(=O)O"))
    AllChem.EmbedMolecule(mol, randomSeed=12)
    mol.SetProp("measured_endpoint", "IC50 > 10 uM; private assay")
    source = tmp_path / "source.sdf"
    with Chem.SDWriter(str(source)) as writer:
        writer.write(mol)
    result = prepare_library(source, tmp_path / "out")
    prepared = next(iter(Chem.SDMolSupplier(result["prepared_sdf"], removeHs=False)))
    assert prepared.GetProp("measured_endpoint") == mol.GetProp("measured_endpoint")
    assert result["compounds"][0]["preparation"]["input_3d_preserved"]
    assert result["compounds"][0]["preparation"]["ph"] is None
    assert (abs(prepared.GetConformer().GetPositions() - mol.GetConformer().GetPositions()) < 0.0001).all()


def test_redocking_metric_does_not_align_away_wrong_pocket_placement():
    Chem = pytest.importorskip("rdkit.Chem")
    AllChem = pytest.importorskip("rdkit.Chem.AllChem")
    reference = Chem.AddHs(Chem.MolFromSmiles("CCO"))
    AllChem.EmbedMolecule(reference, randomSeed=7)
    probe = Chem.Mol(reference)
    conformer = probe.GetConformer()
    for atom in range(probe.GetNumAtoms()):
        point = conformer.GetAtomPosition(atom)
        conformer.SetAtomPosition(atom, (point.x + 10, point.y, point.z))
    result = _symmetry_rmsd(probe, reference)
    assert result["rmsd_angstrom"] == pytest.approx(10)
    assert result["aligned_ligand_shape_rmsd_angstrom"] < 0.001


def test_runner_timeout_and_cancellation(tmp_path):
    with pytest.raises(ScienceError, match="timeout") as caught:
        run_process([sys.executable, "-c", "import time; time.sleep(10)"], tmp_path, tmp_path / "timeout.log", 0.1)
    assert caught.value.code == "timeout"
    with pytest.raises(ScienceError) as caught:
        run_process([sys.executable, "-c", "print('never')"], tmp_path, tmp_path / "cancel.log", 1, lambda: True)
    assert caught.value.code == "cancelled"


def test_unknown_receptor_residues_block_before_external_template_lookup(tmp_path):
    source = tmp_path / "unknown.pdb"
    source.write_text("HETATM    1  C1  UNK A   1       1.000   2.000   3.000  1.00 10.00           C  \nEND\n")
    for method in ("pdbfixer", "meeko", "pdb2pqr"):
        with pytest.raises(ScienceError) as caught:
            prepare_receptor(source, tmp_path / method, method, reviewed=True)
        assert caught.value.code == "unsupported_receptor_residue"
        assert not (tmp_path / method / f"{method}.log").exists()


def test_box_uses_bound_ligand_extents(tmp_path):
    Chem = pytest.importorskip("rdkit.Chem")
    AllChem = pytest.importorskip("rdkit.Chem.AllChem")
    mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
    AllChem.EmbedMolecule(mol, randomSeed=1)
    path = tmp_path / "bound.sdf"
    with Chem.SDWriter(str(path)) as writer:
        writer.write(mol)
    box = derive_box(path, 4)
    assert all(8 < side < 15 for side in box["size"])
    assert box["method"] == "bound_ligand_heavy_atom_extents_plus_padding"


def _benchmark_data():
    rows = [{"compound_id": f"c{i}", "partition": "tuning" if i < 2 else "evaluation", "scaffold_group": f"s{i}",
             "label": "active" if i % 2 == 0 else "tested_inactive", "endpoint_definition": "human ABL1 binding IC50 <= 1 uM",
             "source_id": f"private_assay:record_{i}"} for i in range(8)]
    params = {"dataset_version": "test-1", "endpoint_definition": rows[0]["endpoint_definition"], "tuning_ids": ["c0", "c1"],
              "evaluation_ids": [f"c{i}" for i in range(2, 8)], "score_name": "vina_kcal_mol", "top_k": 2,
              "minimum_evaluation_samples": 6, "minimum_active_samples": 2, "minimum_inactive_samples": 2,
              "minimum_precision_at_k": 0.75, "minimum_enrichment_at_k": 1.5, "declared_before_evaluation": True}
    return rows, params


def test_failed_enrichment_report_and_frozen_protocol(tmp_path):
    rows, params = _benchmark_data()
    frozen = freeze_protocol(rows, params, tmp_path / "protocol.json")
    scores = [{"compound_id": f"c{i}", "score": -10 if i % 2 else -5, "score_name": "vina_kcal_mol", "job_id": "job-real", "artifact_sha256": "a" * 64} for i in range(2, 8)]
    result = evaluate(rows, scores, frozen)
    assert result["status"] == "failed"
    assert result["enrichment_factor_at_k"] == 0
    assert result["precision_at_k_95pct_wilson_interval"][1] > 0
    with pytest.raises(FileExistsError):
        freeze_protocol(rows, params, tmp_path / "protocol.json")
    with pytest.raises(ValueError, match="every evaluation"):
        evaluate(rows, scores[:-1], frozen)


def test_scaffold_leakage_rejected(tmp_path):
    rows, params = _benchmark_data()
    rows[2]["scaffold_group"] = rows[0]["scaffold_group"]
    with pytest.raises(ValueError, match="leakage"):
        freeze_protocol(rows, params, tmp_path / "protocol.json")
