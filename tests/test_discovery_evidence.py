from copy import deepcopy
import hashlib
from types import SimpleNamespace

import pytest

from biocentric import discovery_evidence as evidence
from biocentric.retrieval import RetrievalError

IDENTITY = {'accession': 'P00519', 'taxon_id': 9606, 'organism': 'Homo sapiens'}


def envelope(data):
    return {'data': data, 'retrieved_at': '2026-09-26T00:00:00Z'}


def test_structure_pagination_is_experimental_and_accession_verified(monkeypatch):
    calls = []
    def request(source, url, params=None, payload=None):
        calls.append((url, params, payload))
        if 'search.rcsb' in url:
            return envelope({'result_set': [{'identifier': '1OPL'}, {'identifier': '2BAD'}], 'total_count': 85})
        row = {'rcsb_id': '1OPL', 'struct': {'title': 'Abl kinase'}, 'exptl': [{'method': 'X-RAY DIFFRACTION'}],
               'rcsb_entry_info': {'resolution_combined': [3.42], 'deposited_model_count': 1},
               'polymer_entities': [{'rcsb_polymer_entity': {'pdbx_mutation': 'D382N'},
                 'rcsb_polymer_entity_container_identifiers': {'uniprot_ids': ['P00519'], 'auth_asym_ids': ['A']}}],
               'nonpolymer_entities': [{'nonpolymer_comp': {'chem_comp': {'id': 'P16', 'name': 'Inhibitor', 'formula_weight': 427.3}}}]}
        wrong = deepcopy(row)
        wrong['rcsb_id'] = '2BAD'
        wrong['polymer_entities'][0]['rcsb_polymer_entity_container_identifiers']['uniprot_ids'] = ['P00520']
        return envelope({'data': {'entries': [row, wrong]}})
    monkeypatch.setattr(evidence, 'client', lambda root: SimpleNamespace(_request=request))
    result = evidence.structure_page(IDENTITY, '.', page=2)
    query = calls[0][2]
    assert query['request_options']['paginate'] == {'start': 48, 'rows': 24}
    assert query['request_options']['results_content_type'] == ['experimental']
    assert query['query']['parameters']['value'] == 'P00519'
    assert [r['id'] for r in result['rows']] == ['1OPL']
    assert result['rows'][0]['compounds'][0]['id'] == 'P16'
    assert result['rows'][0]['mutations'] == ['D382N']
    assert result['total'] == 85


def test_binding_data_keeps_relations_quality_and_publication_and_filters_species(monkeypatch):
    calls = []
    def request(source, url, params=None, payload=None):
        calls.append((url, params))
        if url.endswith('/target.json'):
            return envelope({'targets': [{'target_chembl_id': 'CHEMBL1', 'tax_id': 9606, 'target_type': 'SINGLE PROTEIN', 'target_components': [{'accession': 'P00519'}]},
                                         {'target_chembl_id': 'CHEMBL2', 'tax_id': 10090, 'target_type': 'SINGLE PROTEIN', 'target_components': [{'accession': 'P00519'}]}]})
        if url.endswith('/activity.json'):
            return envelope({'activities': [{'target_chembl_id': 'CHEMBL1', 'assay_type': 'B', 'activity_id': 10, 'molecule_chembl_id': 'CHEMBL3',
                  'standard_type': 'IC50', 'standard_relation': '>', 'standard_value': '10000', 'standard_units': 'nM',
                  'document_chembl_id': 'CHEMBL4', 'data_validity_comment': 'Outside typical range', 'assay_variant_mutation': 'T315I'},
                 {'target_chembl_id': 'CHEMBL2', 'assay_type': 'B'}], 'page_meta': {'total_count': 500}})
        return envelope({'documents': [{'document_chembl_id': 'CHEMBL4', 'pubmed_id': 123456, 'title': 'Measured inhibitors'}]})
    monkeypatch.setattr(evidence, 'client', lambda root: SimpleNamespace(_request=request))
    result = evidence.binding_page(IDENTITY, '.', page=3)
    assert calls[1][1]['offset'] == 150
    assert calls[1][1]['target_chembl_id__in'] == 'CHEMBL1'
    assert calls[1][1]['standard_type__in'] == 'Kd,Ki,IC50'
    assert len(result['rows']) == 1
    r = result['rows'][0]
    assert (r['endpoint'], r['relation'], r['value'], r['units']) == ('IC50', '>', '10000', 'nM')
    assert r['validity'] == 'Outside typical range'
    assert r['raw']['assay_variant_mutation'] == 'T315I'
    assert r['publication_url'] == 'https://pubmed.ncbi.nlm.nih.gov/123456/'


def test_partial_outage_preserves_other_sources(monkeypatch):
    monkeypatch.setattr(evidence, 'structure_page', lambda *args: {'rows': [{'id': '1OPL'}], 'total': 1})
    def unavailable(*args):
        raise RetrievalError('ChEMBL timed out')
    monkeypatch.setattr(evidence, 'binding_page', unavailable)
    monkeypatch.setattr(evidence, 'literature', lambda *args: {'evidence': []})
    result = evidence.overview(IDENTITY, '.')
    assert result['structures']['rows'][0]['id'] == '1OPL'
    assert result['binding']['total'] is None
    assert 'timed out' in result['binding']['error']


def test_ambiguous_target_does_not_silently_select_an_accession(monkeypatch):
    candidates = [{'accession': 'P00519', 'gene_names': ['ABL1'], 'reviewed': True},
                  {'accession': 'P42684', 'gene_names': ['ABL2'], 'reviewed': True}]
    monkeypatch.setattr(evidence, 'client', lambda root: SimpleNamespace(resolve_target=lambda *a: {'status': 'ambiguous', 'candidates': candidates}))
    assert evidence.resolve_identity('ABL', 'Human', '.')['status'] == 'ambiguous'
    assert evidence.resolve_identity('ABL1', 'Human', '.')['candidates'] == [candidates[0]]


def test_experimental_viewer_is_offline_and_binds_coordinates_to_entry():
    text = 'data_1OPL\n_atom_site.id\n</script><script>bad()</script>'
    geometry = {'id': '1OPL', 'text': text, 'sha256': hashlib.sha256(text.encode()).hexdigest()}
    row = {'id': '1OPL', 'target_chains': ['A']}
    markup = evidence.experimental_viewer(row, geometry, 'P16')
    assert "connect-src 'none'" in markup
    assert 'addModelsAsFrames' in markup
    assert 'setFrame(Number' in markup
    assert '</script><script>bad()' not in markup
    with pytest.raises(ValueError, match='do not match'):
        evidence.experimental_viewer({'id': '1IEP', 'target_chains': ['A']}, geometry)
    geometry['text'] += 'modified'
    with pytest.raises(ValueError, match='do not match'):
        evidence.experimental_viewer(row, geometry)


@pytest.mark.parametrize('pdb', ['../../secret', 'https://other.example/a', '1OPL?x=y'])
def test_coordinate_fetch_rejects_non_identifier_inputs(tmp_path, pdb):
    with pytest.raises(ValueError, match='Invalid PDB'):
        evidence.coordinates(pdb, tmp_path)


def test_nmr_coordinate_download_and_hash_checked_cache(tmp_path, monkeypatch):
    import httpx
    import json
    calls = []
    pdb = 'HEADER    NMR EXAMPLE                             1AB2\nMODEL        1\nATOM      1  CA  ALA A   1       1.000   2.000   3.000\nENDMDL\nMODEL        2\nATOM      1  CA  ALA A   1       2.000   3.000   4.000\nENDMDL\n'
    def handler(request):
        calls.append(str(request.url))
        return httpx.Response(200, text=pdb)
    original = httpx.Client
    monkeypatch.setattr(evidence.httpx, 'Client', lambda **kw: original(transport=httpx.MockTransport(handler), **kw))
    first = evidence.coordinates('1AB2', tmp_path, 'pdb')
    assert first['format'] == 'pdb'
    assert first['text'].count('MODEL ') == 2
    assert calls == ['https://files.rcsb.org/download/1AB2.pdb']
    assert evidence.coordinates('1AB2', tmp_path, 'pdb') == first
    assert len(calls) == 1
    path = tmp_path / 'discovery-cache/coordinates/1AB2.pdb.json'
    altered = json.loads(path.read_text())
    altered['text'] = 'corrupted coordinates'
    path.write_text(json.dumps(altered))
    assert evidence.coordinates('1AB2', tmp_path, 'pdb')['text'] == pdb
    assert len(calls) == 2


def test_pdb_resolution_retains_requested_structure(monkeypatch):
    fake = SimpleNamespace(_request=lambda *a, **k: envelope({'data': {'entries': [{
        'polymer_entities': [{'rcsb_polymer_entity_container_identifiers': {'uniprot_ids': ['P12931']}}]}]}}),
        resolve_target=lambda *a: {'candidates': [{'accession': 'P12931'}]})
    monkeypatch.setattr(evidence, 'client', lambda root: fake)
    result = evidence.resolve_identity('1a1b', 'Human', '.')
    assert result['status'] == 'resolved'
    assert result['candidates'][0]['requested_pdb'] == '1A1B'
