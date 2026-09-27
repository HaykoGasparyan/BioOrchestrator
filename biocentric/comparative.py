"""Frozen multi-protocol selection and a single independent held-out evaluation.

Selection consumes tuning rows only. Predictions must carry the freeze hash;
held-out predictions additionally carry the selection-lock hash. These are
application provenance barriers, not a claim that a caller cannot have viewed
external labels. A reference redock can select an operating protocol only.
"""
from __future__ import annotations

from collections import defaultdict
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import random
from typing import Literal

from pydantic import BaseModel, ConfigDict, Field, model_validator

from .testsets import validate_testset


HASH_PATTERN = r"^[a-f0-9]{64}$"


def content_hash(value) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()).hexdigest()


class CandidateProtocol(BaseModel):
    model_config = ConfigDict(extra="forbid")
    protocol_id: str = Field(min_length=1, max_length=100, pattern=r"^[A-Za-z0-9_.-]+$")
    engine_id: str = Field(min_length=1)
    metric: str = Field(min_length=1)
    units: str = Field(min_length=1)
    direction: Literal["lower_is_better", "higher_is_better"]
    parameters: dict
    input_hashes: dict[str, str]
    provenance: dict
    training_overlap: Literal["unknown", "not_applicable", "checked"] = "unknown"
    known_training_compound_ids: list[str] = Field(default_factory=list)
    known_training_scaffold_groups: list[str] = Field(default_factory=list)

    @model_validator(mode="after")
    def bindings(self):
        import re
        if not self.input_hashes or any(not re.fullmatch(HASH_PATTERN, value) for value in self.input_hashes.values()):
            raise ValueError("Every candidate needs reviewed input hashes")
        if not self.provenance:
            raise ValueError("Engine/code/weight provenance must be recorded, including unavailable information")
        # Ensure protocol values cannot contain NaN or unhashable objects.
        content_hash(self.model_dump())
        return self


class ComparisonCriteria(BaseModel):
    model_config = ConfigDict(extra="forbid")
    mode: Literal["pose_software", "affinity_regression"]
    maximum_pose_rmsd_angstrom: float = Field(default=2.0, gt=0, allow_inf_nan=False)
    require_physical_validity: Literal[True] = True
    require_complete_accounting: Literal[True] = True
    minimum_tuning_samples: int = Field(default=1, ge=1)
    minimum_evaluation_samples: int | None = Field(default=None, ge=3)
    minimum_evaluation_scaffolds: int | None = Field(default=None, ge=2)
    maximum_mae_p_units: float | None = Field(default=None, ge=0, allow_inf_nan=False)
    minimum_spearman_rho: float | None = Field(default=None, ge=-1, le=1, allow_inf_nan=False)
    bootstrap_samples: int = Field(default=1000, ge=100, le=10000)
    seed: int = 42
    scientist_reviewed: Literal[True]

    @model_validator(mode="after")
    def predictive_cutoffs(self):
        if self.mode == "affinity_regression" and any(getattr(self, key) is None for key in (
            "minimum_evaluation_samples", "minimum_evaluation_scaffolds", "maximum_mae_p_units", "minimum_spearman_rho"
        )):
            raise ValueError("Affinity evaluation requires explicit scientist-approved sample-size and performance criteria")
        return self


class PairResult(BaseModel):
    model_config = ConfigDict(extra="forbid", strict=True)
    protocol_id: str
    compound_id: str
    partition: Literal["tuning", "evaluation"]
    status: Literal["completed", "failed", "excluded"]
    value: float | None = Field(default=None, allow_inf_nan=False)
    metric: str
    units: str
    direction: Literal["lower_is_better", "higher_is_better"]
    canonical_smiles: str | None = None
    protein_id: str | None = None
    physical_validity: bool | None = None
    reason: str | None = None
    job_id: str = Field(min_length=1)
    artifact_sha256: str = Field(pattern=HASH_PATTERN)
    comparison_sha256: str = Field(pattern=HASH_PATTERN)
    selection_sha256: str | None = Field(default=None, pattern=HASH_PATTERN)

    @model_validator(mode="after")
    def outcome(self):
        if self.status == "completed" and self.value is None:
            raise ValueError("Completed result requires a finite measured calculation")
        if self.status != "completed" and (not self.reason or self.value is not None):
            raise ValueError("Failure/exclusion requires a reason and must not carry a fabricated score")
        return self


def _envelope(record, hash_key):
    return {**record, hash_key: content_hash(record)}


def _verify(record, hash_key):
    if not record or record.get(hash_key) != content_hash({k: v for k, v in record.items() if k != hash_key}):
        raise ValueError(f"{hash_key} integrity mismatch")


def _exclusive(destination, record):
    path = Path(destination)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8") as handle:
        json.dump(record, handle, indent=2, allow_nan=False)
    return record


def _dataset(testset, mode):
    validate_testset(testset)
    if mode == "pose_software":
        rows = [{"compound_id": row["id"], "partition": "tuning", "scaffold_group": None} for row in testset["pose_records"]]
        if not rows:
            raise ValueError("A software pose comparison requires an actual registered reference")
        return rows
    rows = [row for row in testset["records"] if row["partition"] in {"tuning", "evaluation"}]
    if not rows or not testset["tuning_ids"] or not testset["evaluation_ids"]:
        raise ValueError("Affinity comparison requires independent tuning and evaluation measurements")
    cohorts = {content_hash(row["cohort"]) for row in rows}
    if len(cohorts) != 1:
        raise ValueError("Do not pool target, species, construct, assay, endpoint or conditions")
    for row in rows:
        if row["endpoint"] not in {"Kd", "Ki"} or row["relation"] != "=":
            raise ValueError("Affinity regression requires exact endpoint-matched Kd or Ki; IC50/EC50 and censoring are unsupported")
        if any(not row["cohort"].get(key) for key in ("target_id", "organism", "construct", "assay_id")):
            raise ValueError("Affinity comparison requires target, species, construct and assay identity")
        if not row.get("source_ids") or not row.get("canonical_smiles"):
            raise ValueError("Measured records require structure and source provenance")
        if row.get("p_measurement") is None or not math.isfinite(row["p_measurement"]):
            raise ValueError("A finite endpoint-matched p measurement is required")
    tuning = [row for row in rows if row["partition"] == "tuning"]
    evaluation = [row for row in rows if row["partition"] == "evaluation"]
    for key in ("compound_id", "canonical_smiles", "scaffold_group"):
        if {row[key] for row in tuning} & {row[key] for row in evaluation}:
            raise ValueError(f"Chemical/scaffold leakage across partitions: {key}")
    return rows


def freeze_comparison(testset: dict, candidates: list[dict], criteria: dict, destination: str | Path) -> dict:
    """Exclusive, pre-prediction freeze of every candidate and selection rule."""
    criteria = ComparisonCriteria.model_validate(criteria)
    candidates = [CandidateProtocol.model_validate(row).model_dump() for row in candidates]
    ids = [row["protocol_id"] for row in candidates]
    if not 1 <= len(ids) <= 16 or len(ids) != len(set(ids)):
        raise ValueError("Declare 1–16 distinct candidate protocol IDs")
    rows = _dataset(testset, criteria.mode)
    expected_metric = "top_1_symmetry_aware_rmsd_angstrom" if criteria.mode == "pose_software" else "p" + rows[0]["endpoint"]
    expected_units = "angstrom" if criteria.mode == "pose_software" else "-log10(mol/L)"
    expected_direction = "lower_is_better" if criteria.mode == "pose_software" else "higher_is_better"
    evaluation = [row for row in rows if row["partition"] == "evaluation"]
    for candidate in candidates:
        if (candidate["metric"], candidate["units"], candidate["direction"]) != (expected_metric, expected_units, expected_direction):
            raise ValueError("Candidate output metric, units and direction must match the frozen endpoint; unlike raw scores cannot be pooled")
        if ({row["compound_id"] for row in evaluation} & set(candidate["known_training_compound_ids"]) or
                {row["scaffold_group"] for row in evaluation} & set(candidate["known_training_scaffold_groups"])):
            raise ValueError("Known pretrained training overlap contaminates the held-out evaluation; revise inputs before freezing")
    result = {"schema_version": 1, "created_at": datetime.now(timezone.utc).isoformat(),
              "testset_sha256": testset["manifest_sha256"], "dataset_sha256": testset["dataset_sha256"],
              "candidates": sorted(candidates, key=lambda row: row["protocol_id"]), "criteria": criteria.model_dump(),
              "tuning_ids": sorted(row["compound_id"] for row in rows if row["partition"] == "tuning"),
              "evaluation_ids": sorted(row["compound_id"] for row in evaluation),
              "selection_rule": "Complete passing candidates only; lowest tuning mean RMSD or MAE, then protocol_id. No held-out scores are accepted for selection.",
              "predictive_scope": "not_established" if criteria.mode == "pose_software" else "endpoint_matched_independent_evaluation_pending",
              "limitations": ["Model training overlap remains unknown unless checked against exact training data and checkpoint.",
                              "Scaffold separation does not eliminate all analog leakage.",
                              "A reference redock and unlabeled compounds cannot establish affinity accuracy or screening enrichment."]}
    return _exclusive(destination, _envelope(result, "comparison_sha256"))


def validate_comparison(testset, frozen):
    _verify(frozen, "comparison_sha256")
    if frozen["testset_sha256"] != testset["manifest_sha256"] or frozen["dataset_sha256"] != testset["dataset_sha256"]:
        raise ValueError("Comparison does not bind this frozen test set")
    return _dataset(testset, frozen["criteria"]["mode"])


def _results(results, frozen, records, partition, protocols, selection_sha256=None):
    parsed = [PairResult.model_validate(row).model_dump() for row in results]
    expected = {(protocol["protocol_id"], row["compound_id"]) for protocol in protocols for row in records if row["partition"] == partition}
    actual = [(row["protocol_id"], row["compound_id"]) for row in parsed]
    if len(actual) != len(set(actual)) or set(actual) != expected:
        raise ValueError("Every requested protocol/compound pair must appear exactly once, including explicit failures")
    candidates = {row["protocol_id"]: row for row in protocols}
    identities = {row["compound_id"]: row for row in records}
    for result in parsed:
        candidate = candidates[result["protocol_id"]]
        if result["partition"] != partition or result["comparison_sha256"] != frozen["comparison_sha256"]:
            raise ValueError("Wrong partition or comparison hash; held-out results cannot enter tuning selection")
        if result["selection_sha256"] != selection_sha256:
            raise ValueError("Held-out predictions must bind the prior selection lock; tuning results must not bind a later lock")
        if any(result[key] != candidate[key] for key in ("metric", "units", "direction")):
            raise ValueError("Prediction metric/units/direction differ from the candidate protocol")
        if frozen["criteria"]["mode"] == "affinity_regression":
            identity = identities[result["compound_id"]]
            if result["canonical_smiles"] != identity["canonical_smiles"] or result["protein_id"] != identity["cohort"]["target_id"]:
                raise ValueError("Prediction protein or compound identity differs from the exact measured pair")
        elif result["status"] == "completed" and result["value"] < 0:
            raise ValueError("Receptor-frame RMSD cannot be negative")
    return sorted(parsed, key=lambda row: (row["protocol_id"], row["compound_id"]))


def select_and_lock(testset: dict, tuning_results: list[dict], frozen: dict, destination: str | Path) -> dict:
    """Select on tuning alone, persist the decision before requesting a holdout."""
    records = validate_comparison(testset, frozen)
    results = _results(tuning_results, frozen, records, "tuning", frozen["candidates"])
    criteria = frozen["criteria"]
    targets = {row["compound_id"]: row for row in records}
    summaries = []
    for candidate in frozen["candidates"]:
        rows = [row for row in results if row["protocol_id"] == candidate["protocol_id"]]
        complete = bool(rows) and all(row["status"] == "completed" for row in rows)
        checks = {"all_pairs_completed": complete, "minimum_tuning_samples": len(rows) >= criteria["minimum_tuning_samples"]}
        if criteria["mode"] == "pose_software":
            checks["physical_validity"] = complete and all(row["physical_validity"] is True for row in rows)
            checks["pose_recovery"] = complete and all(row["value"] <= criteria["maximum_pose_rmsd_angstrom"] for row in rows)
            objective = sum(row["value"] for row in rows) / len(rows) if complete else None
        else:
            objective = sum(abs(row["value"] - targets[row["compound_id"]]["p_measurement"]) for row in rows) / len(rows) if complete else None
        summaries.append({"protocol_id": candidate["protocol_id"], "engine_id": candidate["engine_id"],
                          "checks": checks, "eligible": all(checks.values()), "tuning_objective": objective,
                          "failed_pairs": [row for row in rows if row["status"] != "completed"],
                          "training_overlap": candidate["training_overlap"]})
    eligible = sorted((row for row in summaries if row["eligible"]), key=lambda row: (row["tuning_objective"], row["protocol_id"]))
    selected = eligible[0]["protocol_id"] if eligible else None
    lock = {"schema_version": 1, "created_at": datetime.now(timezone.utc).isoformat(),
            "comparison_sha256": frozen["comparison_sha256"], "testset_sha256": testset["manifest_sha256"],
            "tuning_results_sha256": content_hash(results), "selected_protocol_id": selected,
            "selected_protocol": next((row for row in frozen["candidates"] if row["protocol_id"] == selected), None),
            "status": "locked" if selected else "no_passing_protocol", "candidate_results": summaries,
            "selection_scope": "software_pose_only" if criteria["mode"] == "pose_software" else "tuning_affinity_error_only",
            "heldout_access": "not_consumed", "predictive_validation": "not_established",
            "reason": (f"Locked {selected} using the predeclared tuning rule." if selected else "No candidate passed the frozen tuning checks; execution remains blocked."),
            "comparison_scope": ("One protocol tested; no comparative winner established." if len(summaries) == 1 else
                                 "Compared operational pose protocols; no affinity or enrichment superiority established." if criteria["mode"] == "pose_software" else
                                 "Compared endpoint-matched predictions on tuning data; independent evaluation remains pending.")}
    return _exclusive(destination, _envelope(lock, "selection_sha256"))


def _rank(values):
    ordered = sorted(range(len(values)), key=lambda index: values[index])
    ranks = [0.0] * len(values)
    start = 0
    while start < len(ordered):
        end = start + 1
        while end < len(ordered) and values[ordered[end]] == values[ordered[start]]:
            end += 1
        for index in ordered[start:end]:
            ranks[index] = (start + end - 1) / 2
        start = end
    return ranks


def _correlation(left, right):
    if len(left) < 2:
        return None
    a, b = sum(left) / len(left), sum(right) / len(right)
    numerator = sum((x-a)*(y-b) for x, y in zip(left, right))
    denominator = math.sqrt(sum((x-a)**2 for x in left) * sum((y-b)**2 for y in right))
    return numerator / denominator if denominator else None


def _metrics(rows, truth):
    actual = [truth[row["compound_id"]]["p_measurement"] for row in rows]
    predicted = [row["value"] for row in rows]
    errors = [pred-y for pred, y in zip(predicted, actual)]
    return {"mae_p_units": sum(abs(error) for error in errors) / len(errors),
            "rmse_p_units": math.sqrt(sum(error*error for error in errors) / len(errors)),
            "mean_error_p_units": sum(errors) / len(errors),
            "pearson_r": _correlation(actual, predicted), "spearman_rho": _correlation(_rank(actual), _rank(predicted))}


def _uncertainty(rows, truth, criteria):
    groups = defaultdict(list)
    for row in rows:
        groups[truth[row["compound_id"]]["scaffold_group"]].append(row)
    if len(groups) < 2:
        return {"method": "scaffold_cluster_bootstrap", "intervals_95pct": None, "reason": "At least two distinct scaffold groups are required"}
    rng = random.Random(criteria["seed"])
    keys = sorted(groups)
    samples = defaultdict(list)
    for _ in range(criteria["bootstrap_samples"]):
        sample = [row for key in rng.choices(keys, k=len(keys)) for row in groups[key]]
        for name, value in _metrics(sample, truth).items():
            if value is not None and math.isfinite(value):
                samples[name].append(value)
    intervals = {}
    for name in _metrics(rows, truth):
        values = sorted(samples[name])
        intervals[name] = [values[int((len(values)-1)*.025)], values[int((len(values)-1)*.975)]] if values else None
    return {"method": "scaffold_cluster_bootstrap", "seed": criteria["seed"], "samples": criteria["bootstrap_samples"],
            "intervals_95pct": intervals, "finite_replicates": {key: len(values) for key, values in samples.items()},
            "limitations": "Descriptive percentile intervals; few scaffolds and remaining analog dependence limit coverage. No model epistemic uncertainty is implied."}


def evaluate_locked(testset: dict, heldout_results: list[dict], frozen: dict, lock: dict, destination: str | Path) -> dict:
    """One immutable final evaluation; changing results requires a new study."""
    records = validate_comparison(testset, frozen)
    _verify(lock, "selection_sha256")
    if lock["comparison_sha256"] != frozen["comparison_sha256"] or lock["testset_sha256"] != testset["manifest_sha256"]:
        raise ValueError("Selection lock does not bind this comparison")
    selected = next((row for row in frozen["candidates"] if row["protocol_id"] == lock["selected_protocol_id"]), None)
    if selected is None or lock["status"] != "locked" or selected != lock["selected_protocol"]:
        raise ValueError("A passing locked protocol is required before held-out evaluation")
    criteria = frozen["criteria"]
    base = {"schema_version": 1, "created_at": datetime.now(timezone.utc).isoformat(),
            "comparison_sha256": frozen["comparison_sha256"], "selection_sha256": lock["selection_sha256"],
            "selected_protocol_id": selected["protocol_id"], "training_overlap": selected["training_overlap"],
            "reselection_allowed": False}
    if criteria["mode"] == "pose_software":
        if heldout_results:
            raise ValueError("An unlabeled/reference-only study cannot consume supposed predictive held-out results")
        report = {**base, "status": "not_evaluated", "predictive_validation": "not_established", "metrics": None,
                  "reason": "No approved independent measured-affinity protocol was frozen. Reference-pose selection is software-only; no holdout was evaluated."}
    else:
        results = _results(heldout_results, frozen, records, "evaluation", [selected], lock["selection_sha256"])
        truth = {row["compound_id"]: row for row in records}
        complete = bool(results) and all(row["status"] == "completed" for row in results)
        metrics = _metrics(results, truth) if complete else None
        checks = {"complete_predictions": complete, "minimum_evaluation_samples": len(results) >= criteria["minimum_evaluation_samples"],
                  "minimum_evaluation_scaffolds": len({truth[row["compound_id"]]["scaffold_group"] for row in results}) >= criteria["minimum_evaluation_scaffolds"],
                  "maximum_mae_p_units": bool(metrics) and metrics["mae_p_units"] <= criteria["maximum_mae_p_units"],
                  "minimum_spearman_rho": bool(metrics) and metrics["spearman_rho"] is not None and metrics["spearman_rho"] >= criteria["minimum_spearman_rho"]}
        report = {**base, "status": "passed" if all(checks.values()) else "failed",
                  "predictive_validation": "passed_frozen_endpoint_criteria" if all(checks.values()) else "not_established",
                  "endpoint": records[0]["endpoint"], "cohort": records[0]["cohort"], "n_evaluation": len(results),
                  "metrics": metrics, "checks": checks, "failed_criteria": [key for key, passed in checks.items() if not passed],
                  "failed_pairs": [row for row in results if row["status"] != "completed"], "results_sha256": content_hash(results),
                  "uncertainty": _uncertainty(results, truth, criteria) if complete else None,
                  "limitations": frozen["limitations"]}
    return _exclusive(destination, _envelope(report, "evaluation_sha256"))


def run_comparison(testset, candidates, criteria, output_dir, runner, *, emit=None, cancelled=None):
    """Execute bounded protocols through a supplied real typed-adapter runner.

    ``runner(candidate, request)`` returns PairResult dictionaries. Requests
    contain exact identities, partition/freeze/selection hashes, and no measured
    labels. Tuning is requested for every frozen candidate; the final partition
    is requested once, for the locked candidate only. Engine-specific chemistry,
    checkpoint equivalence, cancellation and subprocess bounds remain adapter
    responsibilities. Completed pair tables are immutable and reused on resume.
    """
    out = Path(output_dir)
    frozen_path = out / "frozen-comparison.json"
    if frozen_path.exists():
        frozen = json.loads(frozen_path.read_text())
        validate_comparison(testset, frozen)
        expected = sorted([CandidateProtocol.model_validate(row).model_dump() for row in candidates], key=lambda row: row["protocol_id"])
        if frozen["candidates"] != expected or frozen["criteria"] != ComparisonCriteria.model_validate(criteria).model_dump():
            raise ValueError("A resumed comparison cannot change candidate protocols or criteria")
    else:
        frozen = freeze_comparison(testset, candidates, criteria, frozen_path)
    records = validate_comparison(testset, frozen)

    def predict(candidate, partition, lock=None):
        if cancelled and cancelled():
            raise InterruptedError("Comparison cancelled; completed protocol results remain reusable")
        path = out / "predictions" / (candidate["protocol_id"] + "-" + partition + ".json")
        if path.exists():
            checkpoint = json.loads(path.read_text())
            _verify(checkpoint, "checkpoint_sha256")
            results = _results(checkpoint["results"], frozen, records, partition, [candidate], lock["selection_sha256"] if lock else None)
            if emit:
                emit({"kind": "checkpoint", "protocol_id": candidate["protocol_id"], "partition": partition, "status": "reused"})
            return results
        request = {"partition": partition, "comparison_sha256": frozen["comparison_sha256"],
                   "selection_sha256": lock["selection_sha256"] if lock else None,
                   "pairs": [{"compound_id": row["compound_id"], "canonical_smiles": row.get("canonical_smiles"),
                              "protein_id": row.get("cohort", {}).get("target_id"),
                              "cohort_identity": {key: row.get("cohort", {}).get(key) for key in ("target_id", "organism", "construct")}}
                             for row in records if row["partition"] == partition]}
        if emit:
            emit({"kind": "comparative_prediction", "protocol_id": candidate["protocol_id"], "partition": partition,
                  "status": "running", "pair_count": len(request["pairs"])})
        results = _results(runner(candidate, request), frozen, records, partition, [candidate], request["selection_sha256"])
        _exclusive(path, _envelope({"results": results}, "checkpoint_sha256"))
        if emit:
            emit({"kind": "comparative_prediction", "protocol_id": candidate["protocol_id"], "partition": partition,
                  "status": "completed", "pair_count": len(results), "failed_pairs": sum(row["status"] != "completed" for row in results)})
        return results

    tuning = [row for candidate in frozen["candidates"] for row in predict(candidate, "tuning")]
    lock_path = out / "selection-lock.json"
    if lock_path.exists():
        lock = json.loads(lock_path.read_text())
        _verify(lock, "selection_sha256")
        if lock["comparison_sha256"] != frozen["comparison_sha256"] or lock["tuning_results_sha256"] != content_hash(sorted(tuning, key=lambda row: (row["protocol_id"], row["compound_id"]))):
            raise ValueError("Selection lock no longer matches the completed tuning predictions")
    else:
        lock = select_and_lock(testset, tuning, frozen, lock_path)
    if lock["status"] != "locked":
        return {"comparison": frozen, "selection_lock": lock, "heldout_evaluation": None}
    if emit:
        emit({"kind": "protocol_selection", "status": "locked", "selection_sha256": lock["selection_sha256"], "protocol_id": lock["selected_protocol_id"]})
    heldout = [] if frozen["criteria"]["mode"] == "pose_software" else predict(lock["selected_protocol"], "evaluation", lock)
    evaluation_path = out / "heldout-evaluation.json"
    if evaluation_path.exists():
        report = json.loads(evaluation_path.read_text())
        _verify(report, "evaluation_sha256")
        if report["selection_sha256"] != lock["selection_sha256"] or (heldout and report["results_sha256"] != content_hash(heldout)):
            raise ValueError("Final evaluation is immutable; changed predictions cannot trigger reselection")
    else:
        report = evaluate_locked(testset, heldout, frozen, lock, evaluation_path)
    return {"comparison": frozen, "selection_lock": lock, "heldout_evaluation": report}
