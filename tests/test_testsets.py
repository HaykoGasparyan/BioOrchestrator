from copy import deepcopy
import hashlib

import pytest

from biocentric.testsets import build_testset, validate_testset


def assay(ident, smiles, **extra):
    return {"id": "measurement_" + ident, "compound_id": ident, "canonical_smiles": smiles,
            "assay_id": "assay1", "endpoint": "Kd", "value": 10, "units": "nM", "relation": "=",
            "organism": "Mus musculus", "source": "ChEMBL", **extra}


def snapshot(rows):
    return {"project": {"spec": {"organism": "Mus musculus"}}, "data": {"assays": rows}, "artifacts": []}


def test_actual_cohorts_keep_endpoint_species_assay_and_censoring_separate():
    rows = [assay("a", "c1ccccc1", activity_class="tested_active"),
            assay("b", "C1CCCCC1", value=100, activity_class="tested_inactive"),
            assay("c", "CCO", relation=">"),
            assay("d", "CCN", endpoint="IC50"),
            assay("e", "CCO", assay_id="different"),
            assay("f", "CCC", organism="Homo sapiens")]
    result = build_testset(snapshot(rows))
    assert len(result["cohorts"]) == 3
    assert {row["endpoint"] for row in result["records"]} == {"Kd"}
    censored = next(row for row in result["records"] if row["compound_id"] == "c")
    assert censored["relation"] == ">" and censored["partition"] == "excluded" and censored["p_measurement"] is None
    assert all(row["compound_id"] != "f" for row in result["records"])
    assert any("species" in row["reason"] for row in result["exclusions"])
    assert validate_testset(result)


def test_deterministic_scaffold_split_never_shares_generic_scaffolds():
    rows = [assay("a", "c1ccccc1"), assay("b", "Cc1ccccc1"), assay("c", "c1ccncc1"),
            assay("d", "C1CCCC1"), assay("e", "CCO"), assay("f", "CCN")]
    first = build_testset(snapshot(rows), seed=9)
    second = build_testset(snapshot(list(reversed(rows))), seed=9)
    assert first["dataset_sha256"] == second["dataset_sha256"]
    tuning = {row["scaffold_group"] for row in first["records"] if row["partition"] == "tuning"}
    evaluation = {row["scaffold_group"] for row in first["records"] if row["partition"] == "evaluation"}
    assert tuning and evaluation and not tuning & evaluation
    assert all(row["label"] is None for row in first["records"])
    assert first["status"] == "ready_for_review"
    assert first["criteria"]["declared_before_evaluation"]


def test_conflicting_duplicates_are_excluded_not_averaged_or_relabelled():
    rows = [assay("a", "CCO", value=1), assay("b", "CCO", value=2),
            assay("c", "CCN", activity_class="untested"), assay("d", "CCN", activity_class="untested")]
    result = build_testset(snapshot(rows))
    assert len(result["records"]) == 1
    assert result["records"][0]["duplicate_compound_ids"] == ["c", "d"]
    assert result["records"][0]["label"] is None
    assert len(result["exclusions"]) == 2
    assert all("Conflicting" in row["reason"] for row in result["exclusions"])


def test_reference_only_is_hash_bound_software_pilot_with_no_fake_affinity(tmp_path):
    path = tmp_path / "reference.sdf"
    path.write_text("registered reference content")
    snap = snapshot([])
    snap["data"]["example"] = True
    snap["artifacts"] = [{"id": "ref1", "kind": "reference", "path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}]
    result = build_testset(snap)
    assert result["status"] == "software_only"
    assert result["pose_records"][0]["activity_label"] is None
    assert not result["evaluation_ids"] and not result["records"]
    assert result["criteria"]["pose_recovery"]["maximum"] == 2.0
    assert result["criteria"]["pose_recovery"]["require_physical_validity"]
    assert result["questions"]
    path.write_text("changed reference")
    with pytest.raises(ValueError, match="hash mismatch"):
        build_testset(snap)


def test_changed_criteria_data_or_manifest_cannot_pass_freeze_validation():
    result = build_testset(snapshot([assay("a", "CCO")]))
    for key in ("criteria", "records", "limitations"):
        changed = deepcopy(result)
        if key == "criteria":
            changed[key]["pose_recovery"]["maximum"] = 20
        elif key == "records":
            changed[key][0]["value"] = 900
        else:
            changed[key] = []
        with pytest.raises(ValueError, match="hash mismatch"):
            validate_testset(changed)


def test_missing_smiles_nonfinite_measurements_and_unknown_species_stay_excluded():
    result = build_testset(snapshot([assay("a", ""), assay("b", "CCC", value="nan"), assay("c", "CCN", organism=None)]))
    assert not result["records"] and len(result["exclusions"]) == 3


def test_conflicting_compound_identity_cannot_leak_between_partitions():
    result = build_testset(snapshot([assay("a", "CCO"), assay("a", "c1ccccc1")]))
    assert not result["records"] and len(result["exclusions"]) == 2
    assert all("Conflicting structures" in row["reason"] for row in result["exclusions"])
    assert result["criteria"]["affinity"]["maximum_mae_p_units"] is None
    assert result["criteria"]["screening"]["status"] == "requires_separate_scientist_approved_protocol"


def test_imported_csv_original_smiles_are_valid_structure_provenance():
    row = assay("a", "", source="user_import", original={"smiles": "CCO"})
    row.pop("id")
    result = build_testset(snapshot([row]))
    assert result["records"][0]["canonical_smiles"] == "CCO"
    assert result["records"][0]["source_ids"][0].startswith("user_measurement_")


def test_latest_registered_reference_is_the_only_frozen_reference(tmp_path):
    snap = snapshot([])
    for number in range(2):
        path = tmp_path / f"reference-{number}.sdf"
        path.write_text(f"registered reference {number}")
        snap["artifacts"].append({"id": f"reference{number}", "kind": "reference", "path": str(path),
                                  "sha256": hashlib.sha256(path.read_bytes()).hexdigest()})
    frozen = build_testset(snap)
    assert len(frozen["pose_records"]) == 1 and len(frozen["input_artifacts"]) == 1
    assert frozen["pose_records"][0]["artifact_id"] == "reference1"
    assert frozen["input_artifacts"][0]["sha256"] == snap["artifacts"][-1]["sha256"]


def test_model_selection_view_withholds_heldout_values_and_raw_nested_provider_payloads():
    from biocentric.testsets import model_selection_view
    rows=[assay('a','c1ccccc1',value=12.34567),assay('b','CCO',value=765.4321),assay('c','C1CCC1',value=98.7654)]
    for row in rows:
        row['raw_provider_document']={'DO_NOT_EXPOSE_RAW_TABLE':[12.34567,765.4321,98.7654]}
    frozen=build_testset(snapshot(rows))
    original=deepcopy(frozen)
    view=model_selection_view(frozen)
    assert view['manifest_sha256']==frozen['manifest_sha256'] and frozen==original
    for row in view['records']:
        assert 'original' not in row
        if row['partition']=='evaluation':
            assert set(row)=={'compound_id','canonical_smiles','scaffold_group','partition'}
        else:
            assert row['p_measurement'] is not None and 'value' in row
    import json
    assert 'DO_NOT_EXPOSE_RAW_TABLE' not in json.dumps(view)
    assert view['view_type']=='method_selection_with_final_holdout_measurements_withheld'
    assert validate_testset(frozen)
