"""Deterministic scientific evaluation independent of LLM explanations.

This module never tunes scores. The evaluation partition and pass criteria must
be committed before a score file is evaluated. A tiny redock is not enrichment.
"""
from __future__ import annotations

import hashlib
import json
import math
import re
from pathlib import Path
from typing import Literal

from pydantic import BaseModel, ConfigDict, Field, model_validator


class BenchmarkRecord(BaseModel):
    model_config = ConfigDict(extra="forbid")
    compound_id: str
    partition: Literal["tuning", "evaluation"]
    scaffold_group: str
    label: Literal["active", "tested_inactive"]
    endpoint_definition: str
    source_id: str
    measured_at: str | None = None
    quantitative_value: float | None = None
    relation: str = "="
    unit: str | None = None


class BenchmarkProtocol(BaseModel):
    model_config = ConfigDict(extra="forbid")
    dataset_version: str
    dataset_sha256: str = Field(pattern=r"^[0-9a-f]{64}$")
    endpoint_definition: str
    tuning_ids: list[str]
    evaluation_ids: list[str]
    score_name: str
    lower_is_better: bool = True
    top_k: int = Field(ge=1)
    minimum_evaluation_samples: int = Field(default=20, ge=2)
    minimum_active_samples: int = Field(default=5, ge=1)
    minimum_inactive_samples: int = Field(default=5, ge=1)
    minimum_precision_at_k: float = Field(ge=0, le=1)
    minimum_enrichment_at_k: float = Field(ge=0)
    confidence_level: float = Field(default=0.95, ge=0.95, le=0.95)
    pretrained_training_overlap: str = "unknown"
    temporal_separation: str = "not_available"
    declared_before_evaluation: bool

    @model_validator(mode="after")
    def check_disjoint(self):
        if set(self.tuning_ids) & set(self.evaluation_ids):
            raise ValueError("Tuning and evaluation compounds overlap.")
        if len(set(self.tuning_ids)) != len(self.tuning_ids) or len(set(self.evaluation_ids)) != len(self.evaluation_ids):
            raise ValueError("Repeated IDs in partitions.")
        if not self.declared_before_evaluation:
            raise ValueError("Pass criteria must be declared before inspecting held-out scores.")
        return self


def _hash(data) -> str:
    return hashlib.sha256(json.dumps(data, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()).hexdigest()


def dataset_hash(records: list[dict] | list[BenchmarkRecord]) -> str:
    rows = [row.model_dump() if isinstance(row, BenchmarkRecord) else BenchmarkRecord.model_validate(row).model_dump() for row in records]
    return _hash(sorted(rows, key=lambda row: row["compound_id"]))


def _validate(records: list[dict], protocol: BenchmarkProtocol) -> list[BenchmarkRecord]:
    rows = [BenchmarkRecord.model_validate(row) for row in records]
    ids = [row.compound_id for row in rows]
    if len(set(ids)) != len(ids):
        raise ValueError("Duplicate compound identifiers in benchmark dataset.")
    if dataset_hash(rows) != protocol.dataset_sha256:
        raise ValueError("Benchmark dataset hash does not match the frozen protocol.")
    tuning = [row for row in rows if row.partition == "tuning"]
    evaluation = [row for row in rows if row.partition == "evaluation"]
    if set(row.compound_id for row in tuning) != set(protocol.tuning_ids) or set(row.compound_id for row in evaluation) != set(protocol.evaluation_ids):
        raise ValueError("Dataset partitions differ from frozen protocol.")
    if {row.scaffold_group for row in tuning} & {row.scaffold_group for row in evaluation}:
        raise ValueError("Scaffold/analog leakage across tuning and evaluation partitions.")
    if any(not row.scaffold_group.strip() or not row.source_id.strip() for row in rows):
        raise ValueError("Scaffold groups and measurement source identifiers are required.")
    if any(row.endpoint_definition != protocol.endpoint_definition for row in rows):
        raise ValueError("Incompatible endpoints cannot share this classification benchmark.")
    return evaluation


def freeze_protocol(records: list[dict], parameters: dict, destination: Path) -> dict:
    """Exclusive creation freezes a reviewable protocol before score evaluation."""
    protocol = BenchmarkProtocol.model_validate({**parameters, "dataset_sha256": dataset_hash(records)})
    _validate(records, protocol)
    data = protocol.model_dump()
    envelope = {"protocol": data, "protocol_sha256": _hash(data)}
    Path(destination).parent.mkdir(parents=True, exist_ok=True)
    with Path(destination).open("x") as handle:
        json.dump(envelope, handle, indent=2)
    return envelope


def _wilson(successes: int, n: int) -> list[float] | None:
    if n == 0:
        return None
    z = 1.959963984540054
    p = successes / n
    denominator = 1 + z*z/n
    center = (p + z*z/(2*n)) / denominator
    half = z * math.sqrt(p*(1-p)/n + z*z/(4*n*n)) / denominator
    return [max(0.0, center - half), min(1.0, center + half)]


def evaluate(records: list[dict], scores: list[dict], frozen_protocol: dict) -> dict:
    """Evaluate a frozen held-out classifier ranking; preserve all failed gates."""
    data = frozen_protocol["protocol"]
    if frozen_protocol["protocol_sha256"] != _hash(data):
        raise ValueError("Frozen protocol hash mismatch.")
    protocol = BenchmarkProtocol.model_validate(data)
    rows = _validate(records, protocol)
    expected = set(protocol.evaluation_ids)
    if {row["compound_id"] for row in scores} != expected or len(scores) != len(expected):
        raise ValueError("Scores must contain every evaluation compound exactly once; failures cannot be silently dropped.")
    if any(isinstance(row.get("score"), bool) or not isinstance(row.get("score"), (int, float)) or not math.isfinite(row["score"]) for row in scores):
        raise ValueError("Benchmark scores must all be finite computed numbers.")
    if any(row.get("score_name") != protocol.score_name or not row.get("job_id") or not re.fullmatch(r"[a-f0-9]{64}", str(row.get("artifact_sha256", ""))) for row in scores):
        raise ValueError("Every score needs the declared metric, completed-job reference and artifact hash.")
    labels = {row.compound_id: row.label == "active" for row in rows}
    ordered = sorted(scores, key=lambda row: (row["score"] if protocol.lower_is_better else -row["score"], row["compound_id"]))
    n, active = len(rows), sum(labels.values())
    k = min(protocol.top_k, n)
    hits = sum(labels[row["compound_id"]] for row in ordered[:k])
    precision = hits / k if k else None
    prevalence = active / n if n else None
    enrichment = precision / prevalence if prevalence and precision is not None else None
    uncertainty = _wilson(hits, k)
    checks = {"evaluation_sample_size": n >= protocol.minimum_evaluation_samples,
              "active_sample_size": active >= protocol.minimum_active_samples,
              "inactive_sample_size": n-active >= protocol.minimum_inactive_samples,
              "full_requested_top_k": n >= protocol.top_k,
              "precision_at_k": precision is not None and precision >= protocol.minimum_precision_at_k,
              "enrichment_at_k": enrichment is not None and enrichment >= protocol.minimum_enrichment_at_k}
    return {"status": "passed" if all(checks.values()) else "failed", "scope": "independent held-out ranking under declared endpoint",
            "dataset_version": protocol.dataset_version, "dataset_sha256": protocol.dataset_sha256,
            "protocol_sha256": frozen_protocol["protocol_sha256"], "scores_sha256": _hash(scores),
            "n_evaluation": n, "n_active": active, "n_tested_inactive": n-active, "class_balance_active_fraction": prevalence,
            "top_k": k, "hits_at_k": hits, "precision_at_k": precision, "enrichment_factor_at_k": enrichment,
            "precision_at_k_95pct_wilson_interval": uncertainty,
            "enrichment_95pct_interval_conditional_on_observed_prevalence": [v/prevalence for v in uncertainty] if uncertainty and prevalence else None,
            "uncertainty_limitations": "Wilson binomial approximation; correlated analogs and finite-library sampling can violate independence. Enrichment interval conditions on observed prevalence.",
            "checks": checks, "failed_criteria": [key for key, passed in checks.items() if not passed],
            "tie_policy": "score then compound_id; tied boundary counts reported", "ties_at_boundary": sum(row["score"] == ordered[k-1]["score"] for row in ordered) if k else 0,
            "training_overlap": protocol.pretrained_training_overlap, "temporal_separation": protocol.temporal_separation,
            "affinity_error": None, "affinity_error_reason": "Classification endpoint only; no mixing IC50, Ki, Kd or EC50.",
            "ranked_compound_ids": [row["compound_id"] for row in ordered]}


def smoke_report(docking_result: dict) -> dict:
    controls = [row for row in docking_result.get("poses", []) if row.get("redocking", {}).get("rmsd_angstrom") is not None]
    top = [row for row in controls if row["pose_number"] == 1]
    return {"status": "unvalidated_exploratory", "software_smoke_passed": bool(docking_result.get("poses")),
            "sample_size_bound_controls": len(top), "redocking_rmsd_angstrom": [row["redocking"]["rmsd_angstrom"] for row in top],
            "physical_validity_passed_top_pose": [row["validation"].get("passed") for row in top],
            "independent_tested_actives": 0, "independent_tested_inactives": 0,
            "enrichment": None, "precision_at_k": None, "uncertainty": None,
            "reason": "A tutorial redock and untested chemical controls verify operation only. No independent held-out activity dataset was supplied.",
            "screening_predictivity_gate": "not_passed_missing_controls", "training_overlap": "not_applicable_to_classical_Vina; unknown_for_other_models"}
