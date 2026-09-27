import csv
import io
import math

from biocentric.affinity_results import affinity_rows, affinity_csv


def snapshot(results):
    return {'data': {'research': {'production': {'affinity': {
        'engine': 'biocentric_affinity', 'results': results,
        'native_result': {'provenance': {'weights_sha256': 'a'*64}, 'parameters': {'device': 'cuda'}}}}}}}


def test_affinity_retains_semantics_and_missing_values_without_docking_ranks():
    data = snapshot([
        {'compound_id': 'one', 'protein_id': 'mouse-ABL-A', 'status': 'completed',
         'metrics': [{'name': 'BioCentric_affinity', 'value': 6.2, 'unit': 'pK', 'direction': 'higher_is_better'}]},
        {'compound_id': 'two', 'status': 'failed', 'failure_reason': 'unsupported chemistry', 'metrics': []},
        {'compound_id': 'three', 'status': 'completed', 'metrics': [{'name': 'bad', 'value': math.nan}]},
    ])
    rows = affinity_rows(data)
    assert rows[0]['score'] == 6.2 and rows[0]['checkpoint_sha256'] == 'a'*64
    assert rows[0]['protein_id'] == 'mouse-ABL-A' and rows[0]['direction'] == 'higher_is_better'
    assert rows[1]['score'] is None and rows[1]['failure_reason'] == 'unsupported chemistry'
    assert rows[2]['score'] is None and rows[2]['status'] == 'invalid_score'
    assert all('rank' not in row for row in rows)
    csv_rows = list(csv.DictReader(io.StringIO(affinity_csv(rows).decode())))
    assert csv_rows[1]['score'] == ''
    assert csv_rows[0]['predictive_validation'] == 'not_established'


def test_excluded_inputs_remain_visible_and_formula_cells_are_escaped():
    data = snapshot([])
    run = data['data']['research']['production']['affinity']
    run['metadata'] = {'frozen_input_binding': {'exclusions': [{'record_index': 4, 'canonical_smiles': 'CC', 'reason': '=unsafe'}]}}
    rows = affinity_rows(data)
    assert rows[0]['status'] == 'excluded' and rows[0]['score'] is None
    assert "'=unsafe" in affinity_csv(rows).decode()
