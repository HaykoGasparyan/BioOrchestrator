"""Versioned, JSON-serializable contracts at application trust boundaries."""
from __future__ import annotations
import math
import re
from datetime import datetime, timezone
from typing import Any, Literal
from pydantic import BaseModel, ConfigDict, Field, field_validator


def utc() -> str:
    return datetime.now(timezone.utc).isoformat()


class Record(BaseModel):
    model_config = ConfigDict(extra="forbid", allow_inf_nan=False)


class ProjectSpec(Record):
    name: str = Field(min_length=1, max_length=160)
    target: str = Field(min_length=1, max_length=20000)
    organism: str = Field(min_length=1, max_length=160)
    context: str = Field(default="", max_length=4000)
    intended_effect: str | None = None
    isoform: str | None = None
    mutations: list[str] = Field(default_factory=list)
    shortlist_size: int = Field(default=10, ge=1, le=1000)
    cpu: int = Field(default=1, ge=1, le=16)
    has_gpu: bool = False
    commercial_tools: list[str] = Field(default_factory=list)
    compute_budget_seconds: int = Field(default=600, ge=10, le=86400)
    local_only: bool = True
    allow_model: bool = False
    model: str = "gpt-4.1-mini"
    max_model_calls: int = Field(default=3, ge=0, le=10)
    max_model_output_tokens: int = Field(default=1500, ge=200, le=4000)
    workflow_kind: str = "noncovalent"
    seed: int = Field(default=42, ge=1, le=2147483647)
    mode: Literal["live", "offline_fixture"] = "live"
    workflow_mode: Literal["classic", "research"] = "classic"

    @field_validator("name", "target", "organism")
    @classmethod
    def not_blank(cls, value: str) -> str:
        if not value.strip():
            raise ValueError("Must not be blank")
        return value.strip()


class TargetIdentity(Record):
    accession: str | None = None
    name: str
    organism: str
    sequence: str | None = None
    isoform: str | None = None
    mutations: list[str] = Field(default_factory=list)
    construct: str | None = None
    domain: str | None = None


class EvidenceRecord(Record):
    id: str
    source: str
    source_id: str
    url: str
    retrieved_at: str
    evidence_type: str
    text: str
    location: str
    claim_kind: str = "retrieved_observation"
    instruction_trust: Literal["untrusted_data"] = "untrusted_data"


class AssayMeasurement(Record):
    compound_id: str
    endpoint: Literal["Kd", "Ki", "IC50", "EC50", "percent_inhibition", "other"]
    value: float
    units: str
    relation: Literal["=", "<", "<=", ">", ">=", "~"] = "="
    assay_id: str
    organism: str
    construct: str | None = None
    assay_format: str | None = None
    conditions: str | None = None
    publication: str | None = None
    original: dict[str, Any] = Field(default_factory=dict)
    value_molar: float | None = None
    source: str = "user_import"
    status: Literal["experimentally_confirmed_result_imported_by_user"] = "experimentally_confirmed_result_imported_by_user"

    @field_validator("value")
    @classmethod
    def valid_measurement(cls, value: float) -> float:
        if not math.isfinite(value) or value < 0:
            raise ValueError("Measurement must be finite and nonnegative")
        return value


class StructureCandidate(Record):
    id: str
    source: str
    format: Literal["pdb", "mmcif", "pdbqt"]
    path: str | None = None
    accession: str | None = None
    quality: dict[str, Any] = Field(default_factory=dict)
    uncertain_regions: list[str] = Field(default_factory=list)
    supported: bool | None = None


class PreparationRecord(Record):
    input_hash: str
    output_hash: str
    method: str
    version: str
    ph_assumption: float | None = None
    changes: list[str] = Field(default_factory=list)
    warnings: list[str] = Field(default_factory=list)
    parent_id: str | None = None
    state_id: str | None = None
    job_id: str


class CompoundParent(Record):
    id: str
    smiles: str
    source_id: str | None = None
    original_artifact: str


class CompoundState(Record):
    id: str
    parent_id: str
    smiles: str
    sdf_artifact: str
    assumptions: list[str] = Field(default_factory=list)


class ToolSpec(Record):
    name: str
    role: str
    documentation: str
    pinned_version: str
    inputs: list[str]
    outputs: list[str]
    dependencies: list[str] = Field(default_factory=list)
    hardware: str = "CPU"
    code_terms: str = "See official documentation"
    data_terms: str = "Source-specific terms; retain provenance"
    model_terms: str | None = None
    external_data: str = "None"
    cost_class: str = "local"
    implementation_status: str
    validation_status: str = "unverified"
    health_check: str


class PipelinePlan(Record):
    strategy: Literal["structure_based", "ligand_based", "further_experiments"]
    rationale: str
    assumptions: list[str]
    prerequisites: list[str]
    stages: list[str]
    status: str = "proposed"


class BenchmarkResult(Record):
    status: Literal["passed", "failed", "exploratory", "not_evaluated"]
    criteria: dict[str, Any]
    metrics: dict[str, Any]
    dataset_sha256: str | None = None
    tuning_ids: list[str] = Field(default_factory=list)
    evaluation_ids: list[str] = Field(default_factory=list)
    limitations: list[str] = Field(default_factory=list)


class JobRecord(Record):
    id: str
    project_id: str
    status: Literal["queued", "running", "blocked", "completed", "failed", "cancelled"]
    created_at: str
    updated_at: str
    attempts: int = 0
    cancel_requested: bool = False
    lease_until: float | None = None
    error: str | None = None


class PoseResult(Record):
    compound_id: str
    state_id: str
    job_id: str
    artifact: str
    engine: str
    version: str
    score: float
    units: str
    direction: Literal["lower_is_better", "higher_is_better"]
    validity: dict[str, Any]


class CandidateAssessment(Record):
    compound_id: str
    rank: int
    selected: bool
    rationale: str
    metrics: dict[str, Any]
    missing_evidence: list[str]
    job_id: str
    artifact_ids: list[str]
    status: str = "exploratory"


class RunManifest(Record):
    schema_version: str = "1.0"
    run_id: str
    project_id: str
    created_at: str
    status: str
    spec: dict[str, Any]
    input_hashes: dict[str, str]
    artifacts: list[dict[str, Any]]
    jobs: list[dict[str, Any]]
    decisions: list[dict[str, Any]]
    environment: dict[str, Any]
    limitations: list[str]


def normalize_measurement(row: dict) -> AssayMeasurement:
    clean = {k: (None if v == "" else v) for k, v in row.items()}
    allowed = set(AssayMeasurement.model_fields) - {"original", "value_molar", "status", "source"}
    rec = AssayMeasurement(**{k: v for k, v in clean.items() if k in allowed}, original=row)
    factors = {"M": 1.0, "mM": 1e-3, "uM": 1e-6, "µM": 1e-6, "μM": 1e-6, "nM": 1e-9, "pM": 1e-12}
    if rec.endpoint in {"Kd", "Ki", "IC50", "EC50"}:
        if rec.units not in factors:
            raise ValueError(f"Unsupported concentration unit {rec.units!r}; original endpoint retained")
        rec.value_molar = rec.value * factors[rec.units]
    return rec
