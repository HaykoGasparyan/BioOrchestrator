"""Build reviewable, immutable pilot manifests from real measurements and files.

No activity labels or assay thresholds are invented. The manifest is a preview,
not an evaluation result, and its hashes bind criteria before scores are read.
"""
from __future__ import annotations

from collections import defaultdict
from copy import deepcopy
import hashlib
import json
import math
from pathlib import Path


def _hash(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()).hexdigest()


def _scaffold(smiles):
    from rdkit import Chem
    from rdkit.Chem.Scaffolds import MurckoScaffold
    mol = Chem.MolFromSmiles(smiles)
    if mol is None or "." in smiles:
        raise ValueError("Invalid or multicomponent structure; requires curated parent structure")
    canonical = Chem.MolToSmiles(mol, isomericSmiles=True)
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    # Acyclic compounds stay together rather than receiving separate empty scaffolds.
    generic = Chem.MolToSmiles(MurckoScaffold.MakeScaffoldGeneric(scaffold)) if scaffold.GetNumAtoms() else "ACYCLIC"
    return canonical, "scaffold_" + hashlib.sha256(generic.encode()).hexdigest()[:16]


def _measurement(row, organism, compound_smiles):
    endpoint = row.get("endpoint")
    if isinstance(row.get("value"), bool):
        raise ValueError("A boolean is not a quantitative measurement")
    value = float(row.get("value"))
    relation = row.get("relation") or "="
    units = row.get("units")
    factors = {"M": 1, "mM": 1e-3, "uM": 1e-6, "µM": 1e-6, "μM": 1e-6, "nM": 1e-9, "pM": 1e-12}
    if endpoint not in {"Kd", "Ki", "IC50", "EC50"} or units not in factors:
        raise ValueError("Unsupported endpoint or units; no conversion or endpoint pooling applied")
    if not math.isfinite(value) or value <= 0 or relation not in {"=", "<", "<=", ">", ">=", "~"}:
        raise ValueError("A finite positive measurement and supported relation are required")
    if not row.get("organism") or not organism or row["organism"].strip().casefold() != organism.strip().casefold():
        raise ValueError("Measurement species is missing or differs from the project species")
    compound_id, assay = row.get("compound_id"), row.get("assay_id")
    source_id = row.get("original_record_id") or row.get("evidence_id") or row.get("id") or row.get("source_id")
    if not source_id and row.get("source") == "user_import":
        source_id = "user_measurement_" + _hash(row)[:16]
    if not compound_id or not assay or not source_id:
        raise ValueError("Compound, assay and measurement provenance identifiers are required")
    original = row.get("original") or {}
    smiles = row.get("canonical_smiles") or row.get("smiles") or original.get("canonical_smiles") or original.get("smiles") or compound_smiles.get(compound_id)
    if not smiles:
        raise ValueError("No matched molecular structure; scaffold separation cannot be verified")
    canonical, scaffold = _scaffold(smiles)
    value_molar = value * factors[units]
    if not math.isfinite(value_molar) or value_molar <= 0:
        raise ValueError("Converted measurement is outside the finite positive molar range")
    label = {"tested_active": "active", "tested_inactive": "tested_inactive"}.get(row.get("activity_class"))
    cohort = {key: row.get(key) for key in ("target_id", "assay_id", "organism", "endpoint", "construct", "assay_format", "conditions")}
    return {"compound_id": str(compound_id), "canonical_smiles": canonical, "scaffold_group": scaffold,
            "source_ids": [str(source_id)], "source_url": row.get("url"), "assay_id": assay,
            "organism": row["organism"], "endpoint": endpoint, "value": value, "units": units,
            "relation": relation, "value_molar": value_molar,
            "p_measurement": -math.log10(value_molar) if relation == "=" else None,
            "label": label, "label_basis": "explicit provider tested activity annotation" if label else "unclassified measured activity",
            "cohort": cohort, "original": deepcopy(row),
            "quantitative_eligible": relation == "=", "partition": "excluded" if relation != "=" else "unassigned"}


def _references(snapshot):
    result, limitations = [], []
    latest = {}
    for artifact in snapshot.get("artifacts", []):
        if artifact.get("kind") in {"reference", "receptor", "structure", "library"}:
            latest[artifact["kind"]] = artifact
    for kind, artifact in latest.items():
        path = Path(artifact.get("path", ""))
        if not path.is_file():
            limitations.append(f"Registered {kind} artifact is missing")
            continue
        actual = hashlib.sha256(path.read_bytes()).hexdigest()
        if not artifact.get("sha256") or actual != artifact["sha256"]:
            raise ValueError(f"Registered {kind} artifact hash mismatch")
        result.append({"kind": kind, "artifact_id": artifact.get("id"), "name": artifact.get("name", path.name),
                       "path": str(path), "sha256": actual})
    return result, limitations


def build_testset(snapshot, *, seed=42, max_compounds=24, organism=None):
    """Freeze one compatible quantitative cohort plus a reference-pose pilot.

    Missing/censored/conflicting labels remain visible. At least two generic
    scaffolds are needed to split; imported thresholds are never inferred.
    """
    if not isinstance(max_compounds, int) or not 2 <= max_compounds <= 200:
        raise ValueError("Test-set preview is bounded to 2–200 compounds")
    data = snapshot.get("data", snapshot)
    spec = snapshot.get("project", {}).get("spec", data.get("spec", {}))
    organism = organism or spec.get("organism") or data.get("target_identity", {}).get("organism")
    compound_smiles = {}
    for row in data.get("compounds", []):
        ident = row.get("compound_id") or row.get("chembl_id") or row.get("id")
        smiles = row.get("canonical_smiles") or row.get("smiles") or row.get("properties", {}).get("ConnectivitySMILES") or row.get("properties", {}).get("CanonicalSMILES")
        if ident and smiles:
            compound_smiles[ident] = smiles
    groups, exclusions = defaultdict(list), []
    for row in [*data.get("measurements", []), *data.get("assays", [])]:
        try:
            record = _measurement(row, organism, compound_smiles)
            groups[_hash(record["cohort"])].append(record)
        except (ValueError, TypeError, KeyError) as exc:
            exclusions.append({"compound_id": row.get("compound_id"), "source_id": row.get("id") or row.get("source_id"),
                               "reason": str(exc), "original": deepcopy(row)})
    cleaned = {}
    for cohort_id, rows in groups.items():
        identity_structures = defaultdict(set)
        for row in rows:
            identity_structures[row["compound_id"]].add(row["canonical_smiles"])
        conflicts = {ident for ident, structures in identity_structures.items() if len(structures) > 1}
        by_structure = defaultdict(list)
        for row in rows:
            if row["compound_id"] in conflicts:
                exclusions.append({"compound_id": row["compound_id"],
                                   "reason": "Conflicting structures for one compound identifier; identity review required",
                                   "original": row["original"]})
                continue
            by_structure[row["canonical_smiles"]].append(row)
        unique = []
        for duplicates in by_structure.values():
            values = {(row["value_molar"], row["relation"], row["label"]) for row in duplicates}
            if len(values) != 1:
                exclusions.extend({"compound_id": row["compound_id"], "reason": "Conflicting replicate measurements in one assay; no averaging or preferred-label selection", "original": row["original"]} for row in duplicates)
                continue
            first = sorted(duplicates, key=lambda row: row["compound_id"])[0]
            first["source_ids"] = sorted({source for row in duplicates for source in row["source_ids"]})
            first["duplicate_compound_ids"] = sorted({row["compound_id"] for row in duplicates})
            unique.append(first)
        cleaned[cohort_id] = unique
    ranked = sorted(cleaned, key=lambda key: (-sum(row["quantitative_eligible"] for row in cleaned[key]), -len(cleaned[key]), key))
    selected_id = ranked[0] if ranked else None
    selected = cleaned.get(selected_id, [])
    ordered = sorted(selected, key=lambda row: (_hash([seed, row["scaffold_group"]]), row["compound_id"]))
    selected = ordered[:max_compounds]
    for cohort_id, rows in cleaned.items():
        omitted = rows if cohort_id != selected_id else ordered[max_compounds:]
        reason = "Different assay/target/species/construct/endpoint cohort; not pooled" if cohort_id != selected_id else "Outside the bounded test-set preview; retained as an explicit exclusion"
        exclusions.extend({"compound_id": row["compound_id"], "reason": reason, "original": row["original"]} for row in omitted)
    scaffolds = sorted({row["scaffold_group"] for row in selected if row["quantitative_eligible"]}, key=lambda group: _hash([seed, group]))
    heldout_scaffolds = set(scaffolds[len(scaffolds)//2:]) if len(scaffolds) >= 2 else set()
    for row in selected:
        if row["quantitative_eligible"]:
            row["partition"] = "evaluation" if row["scaffold_group"] in heldout_scaffolds else "tuning"
    tuning = [row["compound_id"] for row in selected if row["partition"] == "tuning"]
    evaluation = [row["compound_id"] for row in selected if row["partition"] == "evaluation"]
    references, limitations = _references(snapshot)
    reference = next((row for row in references if row["kind"] == "reference"), None)
    pose_records = [{"id": "pose_reference_" + reference["sha256"][:16], "artifact_id": reference["artifact_id"],
                     "sha256": reference["sha256"], "kind": "observed_pose_reference",
                     "partition": "pose_recovery", "activity_label": None,
                     "scope": "Pose recovery only; no independent affinity or enrichment label"}] if reference else []
    limitations += ["Preview only; no pilot scores have been computed by test-set construction.",
                    "One assay/endpoint/species/construct/conditions cohort is retained; no endpoint pooling or IC50-to-Kd conversion.",
                    "Generic Murcko scaffold separation limits scaffold leakage; it does not prove absence of all analog leakage.",
                    "Training-set overlap for learned affinity methods is unknown until checked against exact weights/data.",
                    "Numeric activity cutoffs were not invented; untested molecules are never labeled inactive."]
    if not evaluation:
        limitations.append("No independent scaffold-disjoint quantitative evaluation set is available.")
    if reference:
        limitations.append("Reference redocking tests software and pose recovery, not affinity prediction or screening enrichment.")
    if data.get("example"):
        limitations.append("The bundled 1IEP reference and chemical controls are a software demonstration, not an independent labeled benchmark.")
    criteria = {"declared_before_evaluation": True, "provenance": "BioCentricGPT predeclared proof-of-concept criteria v1; scientist approval required",
                "pose_recovery": {"metric": "top_1_symmetry_aware_rmsd_angstrom", "maximum": 2.0,
                                  "require_physical_validity": True, "require_all_reference_poses_scored": True},
                "affinity": {"endpoint": selected[0]["endpoint"] if selected else None,
                             "status": "requires_separate_scientist_approved_protocol",
                             "minimum_evaluation_compounds": None, "minimum_distinct_scaffold_groups": None,
                             "maximum_mae_p_units": None, "minimum_spearman_rho": None,
                             "eligible_relations": ["="], "same_assay_required": True,
                             "scope": "Affinity validation requires exact matched Kd or Ki; IC50/EC50 are assay-response ranking only"},
                "screening": {"status": "requires_separate_scientist_approved_protocol",
                              "minimum_tested_actives": None, "minimum_tested_inactives": None,
                              "minimum_evaluation_compounds": None, "minimum_enrichment_at_k": None,
                              "top_k": None, "requires_independent_measured_labels": True},
                "on_missing_or_failed_predictive_benchmark": "Require explicit exploratory approval; do not claim validated affinity or screening"}
    cohorts = [{"id": key, "cohort": rows[0]["cohort"] if rows else {}, "count": len(rows),
                "quantitative_count": sum(row["quantitative_eligible"] for row in rows), "selected": key == selected_id}
               for key, rows in sorted(cleaned.items())]
    dataset = {"records": selected, "pose_records": pose_records, "input_artifacts": references, "seed": seed,
               "cohort_id": selected_id, "tuning_ids": tuning, "evaluation_ids": evaluation}
    result = {"schema_version": "1.0", "status": "ready_for_review" if tuning and evaluation else "software_only",
              "scope": "quantitative_assay_preview" if tuning and evaluation else "pose_recovery_or_software_only",
              **dataset, "dataset_sha256": _hash(dataset), "criteria": criteria, "criteria_sha256": _hash(criteria),
              "cohorts": cohorts, "exclusions": exclusions, "limitations": limitations,
              "questions": [] if evaluation else [{"id": "exploratory_scope", "question": "No independent labeled benchmark is available. Continue with pose recovery and explicitly exploratory screening?",
                    "options": ["Continue as exploratory", "Provide matched assay data"]}]}
    result["manifest_sha256"] = _hash(result)
    return result


def validate_testset(manifest):
    """Detect post-freeze edits before a planner or evaluator consumes a manifest."""
    data = {key: manifest[key] for key in ("records", "pose_records", "input_artifacts", "seed", "cohort_id", "tuning_ids", "evaluation_ids")}
    if _hash(data) != manifest.get("dataset_sha256") or _hash(manifest["criteria"]) != manifest.get("criteria_sha256"):
        raise ValueError("Frozen test-set data or criteria hash mismatch")
    if _hash({key: value for key, value in manifest.items() if key != "manifest_sha256"}) != manifest.get("manifest_sha256"):
        raise ValueError("Frozen test-set manifest hash mismatch")
    rows = manifest["records"]
    ids = [row["compound_id"] for row in rows]
    if len(ids) != len(set(ids)):
        raise ValueError("Frozen test set contains duplicate compound identities")
    partitions = {name: [row for row in rows if row["partition"] == name] for name in ("tuning", "evaluation")}
    for name, members in partitions.items():
        declared = manifest[name + "_ids"]
        if len(declared) != len(set(declared)) or set(declared) != {row["compound_id"] for row in members}:
            raise ValueError("Frozen test-set partition IDs disagree with measured records")
        if any(row["relation"] != "=" or row.get("p_measurement") is None for row in members):
            raise ValueError("Censored measurements cannot enter exact quantitative partitions")
    for key in ("compound_id", "canonical_smiles", "scaffold_group"):
        if {row[key] for row in partitions["tuning"]} & {row[key] for row in partitions["evaluation"]}:
            raise ValueError("Chemical/scaffold leakage across test-set partitions")
    return True


def model_selection_view(manifest):
    """Give advisors scientific scope and tuning data without final-holdout labels.

    Hashes identify the intact stored manifest, not this deliberately redacted
    view. Raw original records are withheld from every partition because provider
    payloads may embed larger tables. Exclusion originals are also withheld.
    Final metrics may be supplied separately only after deterministic locking.
    """
    validate_testset(manifest)
    view=deepcopy(manifest)
    identity_fields=('compound_id','canonical_smiles','scaffold_group','partition')
    tuning_fields=identity_fields+('source_ids','assay_id','organism','endpoint','value','units','relation',
                                  'value_molar','p_measurement','label','label_basis','cohort','quantitative_eligible')
    view['records']=[{key:deepcopy(row[key]) for key in (tuning_fields if row['partition']=='tuning' else identity_fields) if key in row}
                     for row in manifest['records']]
    view['exclusions']=[{key:deepcopy(row[key]) for key in ('compound_id','reason') if key in row}
                        for row in manifest.get('exclusions',[])]
    view['view_type']='method_selection_with_final_holdout_measurements_withheld'
    view['hash_scope']='Hashes refer to the intact immutable manifest; this is a redacted advisor view, not a replacement manifest.'
    return view
