"""Synthetic invariants only; these examples are not measured affinity evidence."""
from copy import deepcopy
import json

import pytest

from biocentric.comparative import (CandidateProtocol, ComparisonCriteria, content_hash, evaluate_locked,
                                    freeze_comparison, run_comparison, select_and_lock)
from biocentric.testsets import build_testset


def data():
    rows = [{"id": "record" + str(i), "compound_id": "c" + str(i), "smiles": "C1" + "C" * i + "C1",
             "assay_id": "assay", "target_id": "P00520", "organism": "Mus musculus", "construct": "synthetic_test_construct",
             "endpoint": "Kd", "units": "nM", "value": float(10 ** i), "relation": "="} for i in range(1, 9)]
    return build_testset({"data": {"measurements": rows}, "artifacts": []}, organism="Mus musculus")


def candidate(ident="a", **extra):
    return {"protocol_id": ident, "engine_id": "synthetic_fixture_only", "metric": "pKd", "units": "-log10(mol/L)",
            "direction": "higher_is_better", "parameters": {"fixture": ident}, "input_hashes": {"protein": "a" * 64},
            "provenance": {"kind": "test_double_not_a_runtime"}, **extra}


def criteria(**extra):
    return {"mode": "affinity_regression", "minimum_tuning_samples": 3, "minimum_evaluation_samples": 3,
            "minimum_evaluation_scaffolds": 2, "maximum_mae_p_units": .5, "minimum_spearman_rho": .8,
            "scientist_reviewed": True, "bootstrap_samples": 100, **extra}


def freeze(tmp_path, **extra):
    dataset = data()
    frozen = freeze_comparison(dataset, [candidate("a"), candidate("b")], criteria(), tmp_path / "frozen.json")
    return dataset, frozen


def predictions(dataset, frozen, partition="tuning", lock=None):
    candidates = [lock["selected_protocol"]] if lock else frozen["candidates"]
    return [{"protocol_id": protocol["protocol_id"], "compound_id": row["compound_id"], "partition": partition,
             "status": "completed", "value": row["p_measurement"] + (.1 if protocol["protocol_id"] == "a" else .7),
             "metric": protocol["metric"], "units": protocol["units"], "direction": protocol["direction"],
             "protein_id": row["cohort"]["target_id"], "canonical_smiles": row["canonical_smiles"], "job_id": "synthetic_test_job",
             "artifact_sha256": "b" * 64, "comparison_sha256": frozen["comparison_sha256"],
             "selection_sha256": lock["selection_sha256"] if lock else None}
            for protocol in candidates for row in dataset["records"] if row["partition"] == partition]


def test_tuning_only_selection_locks_before_one_time_heldout_evaluation(tmp_path):
    dataset, frozen = freeze(tmp_path)
    lock = select_and_lock(dataset, predictions(dataset, frozen), frozen, tmp_path / "selection.json")
    assert lock["selected_protocol_id"] == "a" and lock["heldout_access"] == "not_consumed"
    heldout = predictions(dataset, frozen, "evaluation", lock)
    result = evaluate_locked(dataset, heldout, frozen, lock, tmp_path / "evaluation.json")
    assert result["status"] == "passed" and result["endpoint"] == "Kd"
    assert result["metrics"]["mae_p_units"] == pytest.approx(.1)
    assert result["metrics"]["spearman_rho"] == pytest.approx(1)
    assert result["uncertainty"]["method"] == "scaffold_cluster_bootstrap"
    assert result["training_overlap"] == "unknown" and not result["reselection_allowed"]
    with pytest.raises(FileExistsError):
        evaluate_locked(dataset, heldout, frozen, lock, tmp_path / "evaluation.json")
    with pytest.raises(FileExistsError):
        select_and_lock(dataset, predictions(dataset, frozen), frozen, tmp_path / "selection.json")


@pytest.mark.parametrize("mutation", ["holdout", "duplicate", "missing", "nan", "boolean", "identity", "protein", "metric", "units", "hash"])
def test_selection_rejects_contamination_substitution_and_incomplete_pairs(tmp_path, mutation):
    dataset, frozen = freeze(tmp_path)
    rows = predictions(dataset, frozen)
    if mutation == "holdout": rows = predictions(dataset, frozen, "evaluation")
    elif mutation == "duplicate": rows.append(deepcopy(rows[0]))
    elif mutation == "missing": rows.pop()
    elif mutation == "nan": rows[0]["value"] = float("nan")
    elif mutation == "boolean": rows[0]["value"] = True
    elif mutation == "identity": rows[0]["canonical_smiles"] = "CCO"
    elif mutation == "protein": rows[0]["protein_id"] = "HUMAN_SUBSTITUTION"
    elif mutation == "metric": rows[0]["metric"] = "binding_probability"
    elif mutation == "units": rows[0]["units"] = "kcal/mol"
    elif mutation == "hash": rows[0]["comparison_sha256"] = "c" * 64
    with pytest.raises(ValueError):
        select_and_lock(dataset, rows, frozen, tmp_path / "lock.json")
    assert not (tmp_path / "lock.json").exists()


def test_failure_is_accounted_and_cannot_improve_candidate_by_selective_exclusion(tmp_path):
    dataset, frozen = freeze(tmp_path)
    rows = predictions(dataset, frozen)
    rows[0].update(status="failed", value=None, reason="Synthetic OOM fixture")
    lock = select_and_lock(dataset, rows, frozen, tmp_path / "lock.json")
    assert lock["selected_protocol_id"] == "b"
    assert not lock["candidate_results"][0]["eligible"]
    heldout = predictions(dataset, frozen, "evaluation", lock)
    heldout[0].update(status="excluded", value=None, reason="Unsupported chemistry fixture")
    result = evaluate_locked(dataset, heldout, frozen, lock, tmp_path / "report.json")
    assert result["status"] == "failed" and result["metrics"] is None
    assert len(result["failed_pairs"]) == 1


def test_holdout_requires_prior_lock_and_never_evaluates_unselected_candidate(tmp_path):
    dataset, frozen = freeze(tmp_path)
    lock = select_and_lock(dataset, predictions(dataset, frozen), frozen, tmp_path / "lock.json")
    for alteration in ("no_lock", "wrong_lock", "unselected", "tuning"):
        heldout = predictions(dataset, frozen, "evaluation", lock)
        if alteration == "no_lock": heldout[0]["selection_sha256"] = None
        elif alteration == "wrong_lock": heldout[0]["selection_sha256"] = "c" * 64
        elif alteration == "unselected": heldout[0]["protocol_id"] = "b"
        else: heldout[0]["partition"] = "tuning"
        with pytest.raises(ValueError):
            evaluate_locked(dataset, heldout, frozen, lock, tmp_path / "report.json")


@pytest.mark.parametrize("endpoint", ["IC50", "EC50"])
def test_assay_response_endpoint_is_not_relabelled_as_affinity(tmp_path, endpoint):
    dataset = data()
    # Rebuild from preserved source records, preserving the intentionally changed endpoint.
    originals = [{**row["original"], "endpoint": endpoint} for row in dataset["records"]]
    dataset = build_testset({"data": {"measurements": originals}}, organism="Mus musculus")
    with pytest.raises(ValueError, match="Kd or Ki"):
        freeze_comparison(dataset, [candidate()], criteria(), tmp_path / "frozen.json")


def test_unknown_target_construct_and_known_pretraining_overlap_block_affinity(tmp_path):
    dataset = data()
    member = next(row for row in dataset["records"] if row["partition"] == "evaluation")
    with pytest.raises(ValueError, match="pretrained"):
        freeze_comparison(dataset, [candidate(known_training_compound_ids=[member["compound_id"]])], criteria(), tmp_path / "frozen.json")
    originals = [{**row["original"], "construct": None} for row in dataset["records"]]
    dataset = build_testset({"data": {"measurements": originals}}, organism="Mus musculus")
    with pytest.raises(ValueError, match="construct"):
        freeze_comparison(dataset, [candidate()], criteria(), tmp_path / "frozen.json")


def test_unlike_scores_and_missing_reviewed_cutoffs_cannot_be_frozen(tmp_path):
    for change in ({"metric": "pKi"}, {"metric": "DTI_probability"}, {"units": "kcal/mol"}, {"direction": "lower_is_better"}):
        with pytest.raises(ValueError, match="endpoint"):
            freeze_comparison(data(), [candidate(**change)], criteria(), tmp_path / "frozen.json")
    with pytest.raises(ValueError, match="scientist-approved"):
        freeze_comparison(data(), [candidate()], criteria(maximum_mae_p_units=None), tmp_path / "frozen.json")


def test_software_only_reference_comparison_never_creates_predictive_evidence(tmp_path):
    reference = tmp_path / "reference.sdf"
    reference.write_text("Registered test fixture, not a real molecular result")
    import hashlib
    dataset = build_testset({"data": {}, "artifacts": [{"id": "ref", "kind": "reference", "path": str(reference),
                                                       "sha256": hashlib.sha256(reference.read_bytes()).hexdigest()}]})
    candidates = [candidate(key, metric="top_1_symmetry_aware_rmsd_angstrom", units="angstrom", direction="lower_is_better") for key in ("a", "b")]
    frozen = freeze_comparison(dataset, candidates, {"mode": "pose_software", "scientist_reviewed": True}, tmp_path / "freeze.json")
    rows = [{"protocol_id": key, "compound_id": dataset["pose_records"][0]["id"], "partition": "tuning", "status": "completed",
             "value": .5, "metric": candidates[0]["metric"], "units": "angstrom", "direction": "lower_is_better", "physical_validity": True,
             "job_id": "test", "artifact_sha256": "d" * 64, "comparison_sha256": frozen["comparison_sha256"]} for key in ("b", "a")]
    lock = select_and_lock(dataset, rows, frozen, tmp_path / "selection.json")
    assert lock["selected_protocol_id"] == "a"  # predeclared stable tie break, independent of result order
    result = evaluate_locked(dataset, [], frozen, lock, tmp_path / "evaluation.json")
    assert result["status"] == "not_evaluated" and result["metrics"] is None
    assert result["predictive_validation"] == "not_established" and lock["selection_scope"] == "software_pose_only"


def test_no_passing_candidate_cannot_proceed_to_final_evaluation(tmp_path):
    dataset, frozen = freeze(tmp_path)
    rows = predictions(dataset, frozen)
    for row in rows:
        row.update(status="failed", value=None, reason="Synthetic calculation failure")
    lock = select_and_lock(dataset, rows, frozen, tmp_path / "lock.json")
    assert lock["status"] == "no_passing_protocol"
    with pytest.raises(ValueError, match="passing locked"):
        evaluate_locked(dataset, [], frozen, lock, tmp_path / "report.json")


def test_callback_orchestration_hides_labels_locks_before_holdout_and_resumes(tmp_path):
    dataset = data()
    truth = {row["compound_id"]: row for row in dataset["records"]}
    calls = []
    def runner(protocol, request):
        calls.append((protocol["protocol_id"], request["partition"]))
        assert "p_measurement" not in json.dumps(request) and "value_molar" not in json.dumps(request)
        if request["partition"] == "evaluation":
            assert (tmp_path / "selection-lock.json").is_file()
            assert request["selection_sha256"] and protocol["protocol_id"] == "a"
        return [{"protocol_id": protocol["protocol_id"], "compound_id": pair["compound_id"], "partition": request["partition"],
                 "status": "completed", "value": truth[pair["compound_id"]]["p_measurement"] + (.1 if protocol["protocol_id"] == "a" else .7),
                 "metric": protocol["metric"], "units": protocol["units"], "direction": protocol["direction"],
                 "canonical_smiles": pair["canonical_smiles"], "protein_id": pair["protein_id"], "job_id": "fixture",
                 "artifact_sha256": "a" * 64, "comparison_sha256": request["comparison_sha256"], "selection_sha256": request["selection_sha256"]}
                for pair in request["pairs"]]
    result = run_comparison(dataset, [candidate("b"), candidate("a")], criteria(), tmp_path, runner)
    assert calls == [("a", "tuning"), ("b", "tuning"), ("a", "evaluation")]
    repeated = run_comparison(dataset, [candidate("a"), candidate("b")], criteria(), tmp_path, runner)
    assert repeated == result and len(calls) == 3
    with pytest.raises(ValueError, match="cannot change"):
        run_comparison(dataset, [candidate("a")], criteria(), tmp_path, runner)
    checkpoint = tmp_path / "predictions" / "a-evaluation.json"
    changed = json.loads(checkpoint.read_text()); changed["results"][0]["value"] = 999
    checkpoint.write_text(json.dumps(changed))
    with pytest.raises(ValueError, match="integrity"):
        run_comparison(dataset, [candidate("a"), candidate("b")], criteria(), tmp_path, runner)
