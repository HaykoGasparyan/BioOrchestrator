"""Bounded, read-only scientific retrieval with explicit provenance and privacy.

Public records are JSON dictionaries for persistence by the workflow service.
Retrieved content is always data; no returned field is evaluated or dispatched.
"""
from __future__ import annotations

import hashlib
import json
import re
import threading
import time
from copy import deepcopy
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable
from urllib.parse import urlparse

import httpx

ACCESSION = re.compile(r"(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(?:-[1-9][0-9]*)?\Z")
AA = frozenset("ACDEFGHIKLMNPQRSTVWYBXZJUO")
ALLOWED_HOSTS = frozenset({"rest.uniprot.org", "www.ebi.ac.uk", "api.platform.opentargets.org", "search.rcsb.org", "data.rcsb.org", "pubchem.ncbi.nlm.nih.gov"})
_last_request: dict[str, float] = {}
_rate_lock = threading.Lock()
_TAXA = {"human": "9606", "homo sapiens": "9606", "mouse": "10090", "mus musculus": "10090", "rat": "10116", "rattus norvegicus": "10116"}


class RetrievalError(RuntimeError):
    """An explicit source error, including a disabled external request."""


class TargetMappingError(RetrievalError):
    """An identity contradiction that blocks all target-specific retrieval."""


def _utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def _digest(value: Any) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, ensure_ascii=True).encode()).hexdigest()


def _literal(value: str, limit: int = 300) -> str:
    """Quote user text as a search literal, preventing query-language injection."""
    value = " ".join(value.split())[:limit]
    return '"' + value.replace("\\", " ").replace('"', " ") + '"'


def _organism_matches(organism: str, record: dict) -> bool:
    needle = organism.strip().lower()
    if not needle:
        return False
    taxon = _TAXA.get(needle, needle if needle.isdigit() else None)
    if taxon:
        return str(record.get("taxonId")) == taxon
    return needle in {str(record.get("scientificName", "")).lower(), str(record.get("commonName", "")).lower()}


def _evidence(source: str, source_id: str, url: str, retrieved_at: str, evidence_type: str, text: str, location: str) -> dict:
    return {"id": "ev_" + _digest([source, source_id, location, retrieved_at, text])[:20], "source": source,
            "source_id": source_id, "url": url, "retrieved_at": retrieved_at,
            "evidence_type": evidence_type, "text": text, "location": location,
            "instruction_trust": "untrusted_data", "claim_kind": "retrieved_observation"}


def deduplicate_assays(records: list[dict]) -> list[dict]:
    """Merge only the same source measurement; equal values alone are not duplicates.

    Cross-database imports can set original_record_id to e.g. ChEMBL:12345.
    Conflicting versions of that record remain separate and are marked for review.
    """
    grouped: dict[str, dict] = {}
    seen_origins: dict[str, list[dict]] = {}
    for index, raw in enumerate(records):
        record = deepcopy(raw)
        origin = record.get("original_record_id")
        if not origin:
            origin = f"{record['source']}:{record['source_id']}" if record.get("source") and record.get("source_id") else f"unidentified-record:{index}"
        # Assay, endpoint, units and original relation must match. Never average.
        fields = ("compound_id", "target_id", "assay_id", "endpoint", "value", "units", "relation", "organism", "construct", "assay_format", "conditions", "publication")
        signature = _digest([origin, {key: record.get(key) for key in fields}])
        lineage = record.get("lineage") or [{"source": record.get("source"), "source_id": record.get("source_id"), "url": record.get("url")}]
        if signature in grouped:
            current = grouped[signature]
            for item in lineage:
                if item not in current["lineage"]:
                    current["lineage"].append(item)
            current.setdefault("imported_records", []).append(record.get("raw", record))
        else:
            record["lineage"] = lineage
            grouped[signature] = record
            if origin in seen_origins:
                record["duplicate_conflict"] = True
                for previous in seen_origins[origin]:
                    previous["duplicate_conflict"] = True
            seen_origins.setdefault(origin, []).append(record)
    return list(grouped.values())


class RetrievalClient:
    """One target at a time, finite result windows, immutable response cache.

    ``emit(name, inputs, outputs)`` is called before and after every request.
    ``outputs['status']`` is started/completed/cached/error/blocked. The callback
    receives request fields, response metadata and counts, never authorization
    headers. Set ``allow_external=False`` to allow cache reads only.
    """

    def __init__(self, cache_dir: Path, allow_external: bool = True, emit: Callable | None = None):
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.allow_external = allow_external
        self.emit = emit
        self.search_log: list[dict] = []
        self.timeout = 20.0
        self.max_attempts = 3
        self.min_interval = 1.0  # Conservative; PubChem permits <=5 requests/s.
        self.max_bytes = 8_000_000
        self._transport = None  # Injectable httpx transport for network-free tests.

    def _event(self, source: str, inputs: dict, outputs: dict):
        event = {"source": source, "inputs": inputs, **outputs, "timestamp": _utc()}
        self.search_log.append(event)
        if self.emit:
            self.emit(f"{source}.retrieve", inputs, outputs)

    def _request(self, source: str, url: str, params: dict | None = None, payload: dict | None = None) -> dict:
        parsed = urlparse(url)
        if parsed.scheme != "https" or parsed.hostname not in ALLOWED_HOSTS or parsed.username or parsed.password or parsed.port not in (None, 443):
            raise RetrievalError("Request destination is not allowlisted HTTPS")
        method = "POST" if payload is not None else "GET"
        inputs = {"method": method, "url": url, "params": params, "json": payload}
        cache_key = _digest(inputs)
        cache_file = self.cache_dir / (cache_key + ".json")
        if cache_file.exists():
            try:
                envelope = json.loads(cache_file.read_text())
                if (envelope.get("request") == inputs and isinstance(envelope.get("data"), dict)
                        and envelope.get("retrieved_at") and envelope.get("data_sha256") == _digest(envelope["data"])):
                    self._event(source, inputs, {"status": "cached", "cache_key": cache_key, "retrieved_at": envelope["retrieved_at"]})
                    return envelope
            except (OSError, ValueError):
                pass  # Corrupt cache never counts as retrieved evidence.
        if not self.allow_external:
            self._event(source, inputs, {"status": "blocked", "reason": "local_only_cache_miss"})
            raise RetrievalError(f"{source}: external retrieval disabled and no cached response exists")
        self._event(source, inputs, {"status": "started", "external_submission": True})
        last_error = "Request did not complete"
        with httpx.Client(timeout=self.timeout, follow_redirects=False, transport=self._transport,
                          headers={"User-Agent": "BioCentricGPT/0.1 (read-only research client)", "Accept": "application/json"}) as client:
            for attempt in range(self.max_attempts):
                with _rate_lock:
                    delay = self.min_interval - (time.monotonic() - _last_request.get(parsed.hostname, 0))
                    if delay > 0:
                        time.sleep(delay)
                    _last_request[parsed.hostname] = time.monotonic()
                try:
                    with client.stream(method, url, params=params, json=payload) as response:
                        if response.status_code in (429, 500, 502, 503, 504):
                            last_error = f"HTTP {response.status_code}"
                            retry_after = response.headers.get("Retry-After", "")
                            # Do not retry early when the source asks for a long wait.
                            if retry_after and (not retry_after.isdigit() or int(retry_after) > 5):
                                last_error += "; Retry-After requires a later run"
                                break
                            if attempt + 1 < self.max_attempts:
                                time.sleep(max(float(retry_after or 0), 0.5 * 2 ** attempt))
                                continue
                        if response.status_code not in (200, 204):
                            last_error = f"HTTP {response.status_code}"
                            break
                        content = bytearray()
                        for chunk in response.iter_bytes():
                            content.extend(chunk)
                            if len(content) > self.max_bytes:
                                raise RetrievalError("Response exceeds configured size limit")
                        data = json.loads(content) if content else {}
                        if not isinstance(data, dict):
                            raise RetrievalError("Source returned an unexpected top-level JSON type")
                        if isinstance(data, dict) and data.get("errors"):
                            raise RetrievalError("GraphQL returned errors: " + json.dumps(data["errors"])[:400])
                        envelope = {"request": inputs, "retrieved_at": _utc(), "source": source,
                                    "data": data, "body_sha256": hashlib.sha256(content).hexdigest(),
                                    "data_sha256": _digest(data), "cache_schema": 1,
                                    "request_url": str(response.url), "status_code": response.status_code,
                                    "service_version": response.headers.get("X-UniProt-Release")}
                        # Same key always retains its original retrieval time. Explicit cache
                        # removal is required to refresh; project evidence stays reproducible.
                        temp = cache_file.with_suffix(".tmp." + str(threading.get_ident()))
                        temp.write_text(json.dumps(envelope, ensure_ascii=False))
                        temp.replace(cache_file)
                        self._event(source, inputs, {"status": "completed", "cache_key": cache_key,
                                    "retrieved_at": envelope["retrieved_at"], "bytes": len(content), "attempt": attempt + 1})
                        return envelope
                except (httpx.HTTPError, ValueError, RetrievalError) as exc:
                    last_error = f"{type(exc).__name__}: {str(exc)[:250]}"
                    if isinstance(exc, (ValueError, RetrievalError)):
                        break
                    if attempt + 1 < self.max_attempts:
                        time.sleep(0.5 * 2 ** attempt)
        self._event(source, inputs, {"status": "error", "error": last_error})
        raise RetrievalError(f"{source}: {last_error}")

    @staticmethod
    def _identity(record: dict, timestamp: str) -> dict:
        description = record.get("proteinDescription", {})
        name = description.get("recommendedName", {}).get("fullName", {}).get("value")
        if not name:
            names = description.get("submissionNames", [])
            name = names[0].get("fullName", {}).get("value") if names else record.get("uniProtkbId")
        accession = record.get("primaryAccession")
        domains = [feature for feature in record.get("features", []) if feature.get("type") in {"Domain", "Region"}]
        return {"accession": accession, "name": name, "organism": record.get("organism", {}).get("scientificName"),
                "taxon_id": record.get("organism", {}).get("taxonId"), "sequence": record.get("sequence", {}).get("value"),
                "isoform": "canonical", "mutations": [], "construct": None, "domain": None,
                "annotated_domains": domains, "reviewed": record.get("entryType") == "UniProtKB reviewed (Swiss-Prot)",
                "gene_names": [gene.get("geneName", {}).get("value") for gene in record.get("genes", []) if gene.get("geneName")],
                "url": f"https://www.uniprot.org/uniprotkb/{accession}/entry", "retrieved_at": timestamp,
                "identity_notes": ["Canonical reference sequence; experimental construct, mutation state and domain require confirmation."],
                "cross_references": record.get("uniProtKBCrossReferences", [])}

    def resolve_target(self, query: str, organism: str) -> dict:
        start = len(self.search_log)
        query, organism = query.strip(), organism.strip()
        limitations: list[str] = []
        base = {"status": "not_found", "candidates": [], "limitations": limitations, "search_log": []}
        if not query or not organism:
            limitations.append("Both a target and an organism are required for species-safe resolution.")
            return base
        lines = query.splitlines()
        if sum(line.lstrip().startswith(">") for line in lines) > 1:
            limitations.append("Multiple FASTA records were supplied. Select one target sequence; records are never concatenated.")
            return base
        sequence = "".join(line.strip() for line in lines if not line.lstrip().startswith(">"))
        if query.startswith(">") or (len(sequence) >= 30 and set(sequence.upper()) <= AA):
            if not sequence or not set(sequence.upper()) <= AA:
                limitations.append("The supplied sequence contains unsupported characters.")
            else:
                base["candidates"] = [{"accession": None, "name": "User-supplied sequence", "organism": organism,
                    "sequence": sequence.upper(), "sequence_sha256": hashlib.sha256(sequence.upper().encode()).hexdigest(),
                    "isoform": None, "mutations": [], "construct": None, "domain": None}]
                base["status"] = "ambiguous"
                limitations.append("Sequence retained locally. Identity/isoform mapping has not been established; provide an accession or confirm a custom target. No sequence was submitted externally.")
            return base
        if len(query) > 300:
            limitations.append("Target names are limited to 300 characters; use a FASTA sequence or accession.")
            return base
        accession = query.upper() if ACCESSION.fullmatch(query.upper()) else None
        if accession and "-" in accession:
            limitations.append("Explicit isoform accession supplied. Canonical substitution is prohibited; import or verify the specific isoform sequence before proceeding.")
            return base
        try:
            if accession:
                result = self._request("UniProt", f"https://rest.uniprot.org/uniprotkb/{accession}.json")
                records = [result["data"]]
            else:
                taxon = _TAXA.get(organism.lower(), organism if organism.isdigit() else None)
                species = f"organism_id:{taxon}" if taxon else f"organism_name:{_literal(organism)}"
                result = self._request("UniProt", "https://rest.uniprot.org/uniprotkb/search",
                                       {"query": f"({_literal(query)}) AND ({species})", "format": "json", "size": 20})
                records = result["data"].get("results", [])
                if len(records) == 20:
                    limitations.append("Only the first 20 identity matches were retrieved; refine the name or provide an accession.")
            candidates = [self._identity(record, result["retrieved_at"]) for record in records if _organism_matches(organism, record.get("organism", {}))]
            base["candidates"] = candidates
            base["status"] = "resolved" if len(candidates) == 1 else "ambiguous" if candidates else "not_found"
            if not candidates:
                limitations.append("No record exactly matched the supplied organism; species were not substituted.")
            elif len(candidates) > 1:
                limitations.append("Multiple accessions match. Choose the correct target before target-specific calculations.")
        except RetrievalError as exc:
            limitations.append(str(exc))
        base["search_log"] = self.search_log[start:]
        return base

    def dossier(self, identity: dict, context: str = "") -> dict:
        identity = deepcopy(identity)  # Enrichment must not mutate persisted run inputs.
        result: dict[str, Any] = {"evidence": [], "assays": [], "structures": [], "compounds": [], "search_log": [], "limitations": [
            "This is a bounded evidence sample, not a systematic review. Association scores are not proof that modulation is beneficial.",
            "Abstracts and database annotations are retrieved; inaccessible full text is not represented as read.",
            "Structure candidates require assembly, construct, ligand, pocket and local-quality inspection before use."]}
        if len(" ".join(context.split())) > 300:
            result["limitations"].append("Literature context query was truncated to its first 300 characters.")
        accession = identity.get("accession") or ""
        if not ACCESSION.fullmatch(accession) or "-" in accession:
            result["limitations"].append("A resolved canonical accession is required for this public-source dossier; custom targets and isoforms need explicit mapping.")
            return result
        start = len(self.search_log)
        for name, operation in [("UniProt", lambda: self._uniprot(identity, result)),
                                ("Europe PMC", lambda: self._literature(identity, context, result)),
                                ("Open Targets", lambda: self._open_targets(identity, result)),
                                ("RCSB PDB", lambda: self._structures(identity, result)),
                                ("ChEMBL", lambda: self._chembl(identity, result)),
                                ("PubChem", lambda: self._pubchem(result))]:
            try:
                operation()
            except TargetMappingError as exc:
                result["limitations"].append(f"Target mapping failed: {exc}; further target-specific retrieval blocked.")
                break
            except (RetrievalError, KeyError, TypeError, ValueError) as exc:
                result["limitations"].append(f"{name} incomplete: {exc}")
        result["assays"] = deduplicate_assays(result["assays"])
        result["search_log"] = self.search_log[start:]
        result["coverage"] = {name: sum(record["source"] == name for record in result["evidence"]) for name in ("UniProt", "Europe PMC", "Open Targets", "RCSB PDB", "ChEMBL", "PubChem")}
        return result

    def _uniprot(self, identity: dict, result: dict):
        accession = identity["accession"]
        envelope = self._request("UniProt", f"https://rest.uniprot.org/uniprotkb/{accession}.json")
        record = envelope["data"]
        if not _organism_matches(identity["organism"], record.get("organism", {})):
            raise TargetMappingError("Dossier accession does not match selected organism")
        for index, comment in enumerate(record.get("comments", [])):
            kind = comment.get("commentType", "annotation")
            result["evidence"].append(_evidence("UniProt", accession, f"https://www.uniprot.org/uniprotkb/{accession}/entry", envelope["retrieved_at"], f"database_annotation:{kind}", json.dumps(comment, ensure_ascii=False), f"comments[{index}]"))
        identity["cross_references"] = record.get("uniProtKBCrossReferences", [])
        identity["taxon_id"] = record.get("organism", {}).get("taxonId")

    def _literature(self, identity: dict, context: str, result: dict):
        gene = next(iter(identity.get("gene_names", [])), None)
        search = f"({_literal(identity['accession'])} OR TITLE_ABS:{_literal(identity.get('name') or gene or identity['accession'])})"
        if context.strip():
            search += f" AND TITLE_ABS:{_literal(context)}"
        envelope = self._request("Europe PMC", "https://www.ebi.ac.uk/europepmc/webservices/rest/search", {"query": search, "format": "json", "resultType": "core", "pageSize": 10})
        records = envelope["data"].get("resultList", {}).get("result", [])
        result["limitations"].append(f"Europe PMC: retrieved {len(records)} of {envelope['data'].get('hitCount', 'unknown')} search hits; name matches can refer to other species and require review.")
        for row in records:
            source_id = f"{row.get('source', 'MED')}:{row.get('id')}"
            text = (row.get("title") or "") + "\n" + (row.get("abstractText") or "Abstract unavailable.")
            result["evidence"].append(_evidence("Europe PMC", source_id, f"https://europepmc.org/article/{row.get('source', 'MED')}/{row.get('id')}", envelope["retrieved_at"], "literature_abstract", text, "title; abstractText"))

    def _open_targets(self, identity: dict, result: dict):
        if str(identity.get("taxon_id", "")) != "9606" and identity.get("organism", "").lower() not in {"human", "homo sapiens"}:
            result["limitations"].append("Open Targets human disease associations skipped for nonhuman target; no ortholog was silently substituted.")
            return
        genes: set[str] = set()
        for ref in identity.get("cross_references", []):
            if ref.get("database") == "Ensembl":
                for prop in ref.get("properties", []):
                    value = prop.get("value", "").split(".")[0]
                    if prop.get("key") == "GeneId" and re.fullmatch(r"ENSG[0-9]+", value):
                        genes.add(value)
        if len(genes) != 1:
            result["limitations"].append("Open Targets requires exactly one verified human Ensembl gene mapping; no target was inferred by name.")
            return
        gene = next(iter(genes))
        query = "query Target($id: String!) { target(ensemblId: $id) { id approvedSymbol approvedName functionDescriptions tractability { label modality value } associatedDiseases(page: {index: 0, size: 10}) { count rows { disease { id name } score } } } }"
        envelope = self._request("Open Targets", "https://api.platform.opentargets.org/api/v4/graphql", payload={"query": query, "variables": {"id": gene}})
        target = envelope["data"].get("data", {}).get("target")
        if not target:
            result["limitations"].append("Open Targets returned no target annotation.")
            return
        url = f"https://platform.opentargets.org/target/{gene}"
        for field in ("functionDescriptions", "tractability"):
            result["evidence"].append(_evidence("Open Targets", gene, url, envelope["retrieved_at"], "database_annotation", json.dumps(target.get(field)), f"data.target.{field}"))
        for index, row in enumerate((target.get("associatedDiseases") or {}).get("rows", [])):
            result["evidence"].append(_evidence("Open Targets", gene, url, envelope["retrieved_at"], "disease_association_not_causation", json.dumps(row), f"data.target.associatedDiseases.rows[{index}]"))
        result["limitations"].append("Open Targets: first 10 associated diseases and tractability annotations; individual genetics studies, direction of effect and safety evidence need separate review.")

    def _structures(self, identity: dict, result: dict):
        query = {"query": {"type": "terminal", "service": "text", "parameters": {"attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession", "operator": "exact_match", "value": identity["accession"]}}, "return_type": "entry", "request_options": {"paginate": {"start": 0, "rows": 4}, "results_content_type": ["experimental"]}}
        envelope = self._request("RCSB PDB", "https://search.rcsb.org/rcsbsearch/v2/query", payload=query)
        matches = envelope["data"].get("result_set", [])
        result["limitations"].append(f"RCSB PDB: first {len(matches)} of {envelope['data'].get('total_count', 0)} experimental entries; no automatic receptor selection.")
        for match in matches:
            pdb = match["identifier"]
            if not re.fullmatch(r"[A-Za-z0-9]{4}|pdb_[A-Za-z0-9]{8}", pdb):
                continue
            data = self._request("RCSB PDB", f"https://data.rcsb.org/rest/v1/core/entry/{pdb}")
            entry = data["data"]
            candidate = {"id": pdb, "pdb_id": pdb, "source": "RCSB PDB", "source_id": pdb,
                         "url": f"https://www.rcsb.org/structure/{pdb}", "accession": identity["accession"],
                         "title": entry.get("struct", {}).get("title"), "methods": [row.get("method") for row in entry.get("exptl", [])],
                         "resolution_angstrom": entry.get("rcsb_entry_info", {}).get("resolution_combined"),
                         "retrieved_at": data["retrieved_at"], "review_status": "requires_structure_review",
                         "assembly_ids": entry.get("rcsb_entry_container_identifiers", {}).get("assembly_ids", []),
                         "coordinate_download_url": f"https://files.rcsb.org/download/{pdb}.cif",
                         "limitations": ["Sequence mapping, mutations, biological assembly, missing pocket residues, cofactors and ligand state are not validated by this metadata retrieval."]}
            result["structures"].append(candidate)
            result["evidence"].append(_evidence("RCSB PDB", pdb, candidate["url"], data["retrieved_at"], "experimental_structure_metadata", json.dumps(candidate), "struct; exptl; rcsb_entry_info; rcsb_entry_container_identifiers"))

    def _chembl(self, identity: dict, result: dict):
        envelope = self._request("ChEMBL", "https://www.ebi.ac.uk/chembl/api/data/target.json", {"target_components__accession": identity["accession"], "target_type": "SINGLE PROTEIN", "limit": 3})
        targets = envelope["data"].get("targets", [])
        if not targets:
            result["limitations"].append("ChEMBL: no single-protein target mapping returned; absence is not evidence of lack of ligands.")
        for target in targets:
            if identity.get("taxon_id") and str(target.get("tax_id")) != str(identity["taxon_id"]):
                result["limitations"].append(f"ChEMBL target {target.get('target_chembl_id')} skipped because taxonomy did not match.")
                continue
            target_id = target["target_chembl_id"]
            activities = self._request("ChEMBL", "https://www.ebi.ac.uk/chembl/api/data/activity.json", {"target_chembl_id": target_id, "limit": 40, "offset": 0})
            rows = activities["data"].get("activities", [])
            assay_ids = sorted({row["assay_chembl_id"] for row in rows if row.get("assay_chembl_id")})
            assay_metadata = {}
            if assay_ids:
                meta = self._request("ChEMBL", "https://www.ebi.ac.uk/chembl/api/data/assay.json", {"assay_chembl_id__in": ",".join(assay_ids), "limit": 40})
                assay_metadata = {row["assay_chembl_id"]: row for row in meta["data"].get("assays", [])}
            for row in rows:
                measurement = self.parse_chembl_activity(row, assay_metadata.get(row.get("assay_chembl_id"), {}), activities["retrieved_at"])
                result["assays"].append(measurement)
                evidence = _evidence("ChEMBL", str(row.get("activity_id")), measurement["url"], activities["retrieved_at"], "assay_measurement", json.dumps({key: measurement[key] for key in ("compound_id", "target_id", "endpoint", "relation", "value", "units", "organism", "assay_id")}), f"activities[activity_id={row.get('activity_id')}] (original measurement fields)")
                measurement["evidence_id"] = evidence["id"]
                result["evidence"].append(evidence)
            total = activities["data"].get("page_meta", {}).get("total_count", "unknown")
            result["limitations"].append(f"ChEMBL {target_id}: {len(rows)} of {total} activity records. No endpoint pooling, IC50-to-Kd conversion, or activity threshold is applied.")

    @staticmethod
    def parse_chembl_activity(row: dict, assay: dict, retrieved_at: str) -> dict:
        source_id = str(row.get("activity_id"))
        url = f"https://www.ebi.ac.uk/chembl/explore/activities/ACTIVITY_ID:{source_id}"
        record = {"id": "chembl_activity_" + source_id, "source": "ChEMBL", "source_id": source_id,
            "original_record_id": "ChEMBL:" + source_id, "url": url, "retrieved_at": retrieved_at,
            "compound_id": row.get("molecule_chembl_id"), "canonical_smiles": row.get("canonical_smiles"),
            "target_id": row.get("target_chembl_id"), "assay_id": row.get("assay_chembl_id"),
            "endpoint": row.get("type"), "value": row.get("value"), "units": row.get("units"), "relation": row.get("relation"),
            "standardized": {"endpoint": row.get("standard_type"), "value": row.get("standard_value"), "units": row.get("standard_units"), "relation": row.get("standard_relation"), "provider": "ChEMBL"},
            "organism": assay.get("assay_organism") or row.get("assay_organism"),
            "construct": assay.get("variant_sequence") or None, "assay_format": assay.get("assay_type") or row.get("assay_type"),
            "conditions": {"description": assay.get("description") or row.get("assay_description"), "parameters": assay.get("assay_parameters"), "cell_type": assay.get("assay_cell_type"), "tissue": assay.get("assay_tissue"), "subcellular_fraction": assay.get("assay_subcellular_fraction")},
            "publication": row.get("document_chembl_id"), "document_year": row.get("document_year"),
            "activity_comment": row.get("activity_comment"), "data_validity_comment": row.get("data_validity_comment"),
            "activity_class": "measured_activity_unclassified", "raw": deepcopy(row), "raw_assay": deepcopy(assay),
            "lineage": [{"source": "ChEMBL", "source_id": source_id, "url": url}]}
        # An absent measurement is not labelled inactive, and numerical cutoffs are
        # project-specific. Only preserve explicit provider activity annotations.
        annotation = str(row.get("activity_comment") or "").strip().lower()
        if annotation in {"active", "inactive"}:
            record["activity_class"] = "tested_" + annotation
        return record

    def _pubchem(self, result: dict):
        initial_count = len(result["compounds"])
        ids = sorted({row["compound_id"] for row in result["assays"] if row.get("compound_id") and re.fullmatch(r"CHEMBL[0-9]+", row["compound_id"])})[:3]
        if not ids:
            result["limitations"].append("PubChem: no public ligand identifiers were available for compound cross-reference; private structures/SMILES were not submitted.")
            return
        molecules = self._request("ChEMBL", "https://www.ebi.ac.uk/chembl/api/data/molecule.json", {"molecule_chembl_id__in": ",".join(ids), "limit": 3})
        for molecule in molecules["data"].get("molecules", []):
            key = (molecule.get("molecule_structures") or {}).get("standard_inchi_key")
            if not key or not re.fullmatch(r"[A-Z]{14}-[A-Z]{10}-[A-Z]", key):
                continue
            try:
                envelope = self._request("PubChem", f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{key}/property/MolecularFormula,MolecularWeight,InChIKey,CanonicalSMILES/JSON")
            except RetrievalError as exc:
                result["limitations"].append(f"PubChem compound {molecule.get('molecule_chembl_id')} unavailable: {exc}; other sampled ligands are still checked.")
                continue
            for row in envelope["data"].get("PropertyTable", {}).get("Properties", []):
                if row.get("InChIKey") != key:
                    result["limitations"].append("PubChem returned a nonmatching InChIKey; compound identity was not merged.")
                    continue
                cid = str(row["CID"])
                compound = {"chembl_id": molecule["molecule_chembl_id"], "pubchem_cid": cid, "inchikey": key,
                            "properties": row, "retrieved_at": envelope["retrieved_at"], "url": f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}"}
                result["compounds"].append(compound)
                result["evidence"].append(_evidence("PubChem", cid, compound["url"], envelope["retrieved_at"], "chemical_identity_properties", json.dumps(compound), f"PropertyTable.Properties[CID={cid}]"))
        result["limitations"].append("PubChem cross-references at most three known public ChEMBL compounds by exact InChIKey. Chemical properties are not ADMET measurements or affinity evidence.")
        if len(result["compounds"]) == initial_count:
            result["limitations"].append("PubChem: no exact compound cross-reference was returned for the sampled public ligands.")
