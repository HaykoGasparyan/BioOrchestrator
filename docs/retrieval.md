# Scientific retrieval

`biocentric.retrieval.RetrievalClient(cache_dir, allow_external=True, emit=None)` implements six read-only public interfaces. `resolve_target(query, organism)` returns `status`, `candidates`, `limitations`, and `search_log`. `dossier(identity, context="")` returns evidence, original assay measurements, structure candidates, compound cross-references, coverage counts, a search log, and limitations. These are bounded evidence samples, not systematic reviews.

## Interfaces and coverage

| Source | Implemented request | Coverage and scientific limits |
| --- | --- | --- |
| UniProt | UniProtKB accession JSON or species-constrained search | Reference sequence, identity and annotation comments; up to 20 name matches. Multiple accessions require selection. Isoform accessions are explicitly unresolved until the specific sequence is verified. Canonical sequence is never substituted silently. |
| Europe PMC | REST `search`, `resultType=core` | First 10 title/abstract records for accession/protein name and optional context. No inaccessible full text is claimed as read; name matches and species relevance need review. |
| Open Targets | GraphQL target query | One human Ensembl gene verified from UniProt cross-references; function, tractability and first 10 disease associations. Nonhuman targets and ambiguous gene mappings are skipped explicitly. Association is not causation or beneficial direction of effect. |
| RCSB PDB | Search API accession query and Data API entry metadata | First four experimental structures, methods, resolution and assembly identifiers. No automatic receptor selection, coordinate download, pocket validation or construct validation is implied. |
| ChEMBL | Target, activity, assay, molecule JSON endpoints | Up to three single-protein target mappings and 40 activities per target, with associated assay metadata. Original endpoint, value, relation, units, assay, organism, publication, conditions and full raw records remain available. |
| PubChem | PUG REST exact InChIKey properties | Cross-reference up to three public ChEMBL compounds; returned InChIKey must match. Does not submit private chemistry, fetch all PubChem bioassays, or infer activity from compound presence. |

Official interface references checked during implementation:

- [UniProt API queries](https://www.uniprot.org/help/api_queries) and [accession syntax](https://www.uniprot.org/help/accession_numbers).
- [Europe PMC REST web service](https://europepmc.org/RestfulWebService).
- [Open Targets GraphQL API](https://platform-docs.opentargets.org/data-access/graphql-api).
- [RCSB API overview](https://www.rcsb.org/docs/programmatic-access/web-apis-overview), [Search API](https://search.rcsb.org/), and [Data API](https://data.rcsb.org/).
- [ChEMBL API](https://www.ebi.ac.uk/chembl/api/data/docs), [target mapping example](https://chembl.github.io/using-new-chembl-web-services/), and [web-service resources](https://chembl.gitbook.io/chembl-interface-documentation/web-services/chembl-data-web-services).
- [PubChem PUG REST specification](https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest) and [usage policy/tutorial](https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest-tutorial).

The application code license does not grant rights to third-party datasets, abstracts or service outputs. Preserve source attribution and follow each provider's terms when redistributing a dossier; this client does not bypass restrictions or retrieve paywalled content.

## Preprocessing and postprocessing

Target names are bounded to 300 characters and quoted as search literals. Organism matching is checked after retrieval as well as constrained in the query. Common human/mouse/rat names map to explicit taxonomy identifiers; other species require an exact returned scientific/common name or taxon ID. A species contradiction blocks downstream dossier requests.

FASTA and sufficiently long amino-acid sequences are recognized locally, validated against the supported residue alphabet, normalized to uppercase, and hashed. They are not submitted to any sequence-search service. A sequence-only target remains unresolved until an accession is supplied or a custom-target route establishes identity. Mutation, construct and domain information are not invented from the reference sequence.

ChEMBL's original `type`, `value`, `units`, and `relation` populate the primary measurement fields without conversion. Provider-standardized fields are retained in a separate `standardized` object. Kd, Ki, IC50 and EC50 are never pooled; IC50 is never converted to Kd. An absent activity or unclassified measurement is never labelled inactive. Identical stable source-record IDs can be merged while retaining lineage; conflicting endpoint/context versions remain separate and flagged. Equal numerical values with different record IDs are not sufficient evidence of duplication.

Each evidence record has a stable ID, provider, source identifier, provider URL, original retrieval time, evidence type, supporting text, and JSON field location. `instruction_trust="untrusted_data"` identifies imported content. Retrieved text is never executed, used as a tool name, or followed as an instruction. The application/model adapter must preserve that trust boundary when summarizing evidence.

## Transparency, caching and privacy

Every request calls `emit(name, inputs, outputs)` before and after the operation. Names are e.g. `UniProt.retrieve`. Output statuses are `started`, `completed`, `cached`, `error`, and `blocked`; request method, endpoint, query, response size, attempt number, cache identifier, original retrieval time and errors are available. These events can be rendered as tool-call cards while a worker runs. HTTP requests are retrieval calls, not LLM model calls.

Public identifier, protein-name, organism and optional disease-context searches are externally disclosed when external access is enabled. Retrieved public compound identifiers may be sent to PubChem. No credentials or authorization headers are logged. The client itself never submits user sequences, uploaded files or private SMILES. `allow_external=False` permits cached reads only; cache misses become explicit limitations without any network request.

HTTPS hosts are allowlisted, redirects are refused, and response bodies are capped at 8 MB. A process-wide per-host lock spaces requests at least one second apart, more conservatively than PubChem's limit of five per second. Requests have a 20-second timeout and at most three attempts for transient network/429/5xx failures. Long or nonnumeric `Retry-After` values stop the operation for a later run instead of retrying early. Multiple separate worker processes must be coordinated by the application to maintain aggregate service limits.

The cache key hashes the complete method/endpoint/parameters/body. Cached JSON envelopes preserve the original response data, request, retrieval timestamp, response-body hash, a verified canonical-data hash, status and available service release. Cache hits retain their original timestamp; altered or corrupt cached data is rejected. There is no hidden freshness claim: an operator can remove a cache entry to refresh it, while already persisted project evidence keeps its original provenance. Failed requests and GraphQL error responses are not cached as successful evidence.

## Validation

`tests/test_retrieval.py` uses `httpx.MockTransport` for deterministic species and target-ambiguity gates, isoform handling, local-only and sequence privacy, original endpoint/unit/relation preservation, lineage deduplication, conflicting measurements, timeout/retry bounds, cache provenance, allowlisted destinations, redirect refusal, source errors and malicious retrieved instructions. Live API availability and coverage are recorded separately from these network-free tests.

Executed on 26 September 2026: `.venv/bin/python -m pytest tests/test_retrieval.py -q` → **26 passed**.

The actual `RetrievalClient` was also exercised against the live public services with human SRC, UniProt **P12931**:

| Live path | Observed result |
| --- | --- |
| UniProt accession resolution and annotations | Resolved human SRC and returned 35 annotation evidence records. |
| Europe PMC core search | Retrieved 10 of 475 title/abstract search hits. |
| RCSB search plus entry data | Retrieved metadata for 1A07, 1A08, 1A09 and 1A1A from 78 accession-mapped entries. These were not selected or validated as docking receptors. |
| ChEMBL target, activity, assay and molecule APIs | Mapped to CHEMBL267; preserved 40 of 15,472 activity records and their source assay metadata. |
| PubChem exact public-ligand cross-reference | Two compounds matched. CHEMBL100174 returned HTTP 404; the client recorded the missing result and continued with the third sampled compound. |
| Open Targets automatic mapping gate | Correctly stopped because the current UniProt SRC entry contains two distinct Ensembl gene IDs, ENSG00000197122 and ENSG00000291971. No gene was selected silently. |
| Open Targets actual GraphQL client and parser | A separate explicit ENSG00000197122 query returned 12 evidence records: function, tractability and 10 disease associations. This test does not resolve the ambiguous automatic mapping. |

The raw live validation output and cached source responses were retained in `/tmp/biocentric-retrieval-live-result.json` and `/tmp/biocentric-retrieval-live/` during the build. They are live retrieval artifacts, not a redistributed fixture or a scientific screening benchmark. All six adapters therefore have executed live paths; unresolved identity mappings and absent database records remain explicit limitations.
