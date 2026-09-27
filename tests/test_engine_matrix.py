"""Deterministic native-format fixtures, never native-engine benchmarks.

All subprocess execution is disabled unless a test supplies a fake runner.
Molecules/topologies here exercise input contracts, not scientific validity.
"""
from pathlib import Path

import pytest
from pydantic import ValidationError

from biocentric import engines, science


@pytest.fixture(autouse=True)
def no_native_execution(monkeypatch):
    monkeypatch.setattr(engines, '_executable', lambda name: None)
    monkeypatch.setattr(engines.platform, 'system', lambda: 'Darwin')
    monkeypatch.delenv('BIOCENTRIC_PLANET_DIR', raising=False)
    monkeypatch.delenv('BIOCENTRIC_MODEL_DIR', raising=False)
    monkeypatch.delenv('BIOCENTRIC_MODEL_PYTHON', raising=False)
    monkeypatch.setattr(science, 'run_process', lambda *args, **kwargs: pytest.fail('A fixture must explicitly mock execution'))


def _pdb(path, *, ligand=False, alpha_carbon=False):
    atom = ' CA  ALA' if alpha_carbon else ' C   LIG'
    path.write_text(f'ATOM      1 {atom} A   1       1.000   2.000   3.000  1.00  0.00     0.000 C\n' +
                    ('ROOT\nENDROOT\nTORSDOF 0\n' if ligand else ''))
    return str(path)


def _sdf(path, records=None, *, is3d=True):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    with Chem.SDWriter(str(path)) as writer:
        for name, smiles, properties in records or [('ethanol', 'CCO', {})]:
            mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
            if is3d:
                assert AllChem.EmbedMolecule(mol, randomSeed=42) == 0
            else:
                AllChem.Compute2DCoords(mol)
            mol.SetProp('_Name', name)
            for key, value in properties.items():
                mol.SetProp(key, str(value))
            writer.write(mol)
    return str(path)


def _topology(path, names):
    path.write_text('%VERSION VERSION_STAMP = V0001.000\n%FLAG POINTERS\n%FORMAT(10I8)\n' +
                    f'{len(names):8}\n%FLAG ATOM_NAME\n%FORMAT(20a4)\n' +
                    ''.join(f'{name:4}' for name in names) + '\n%FLAG CHARGE\n%FORMAT(5E16.8)\n' +
                    ''.join(f'{0.:16.8E}' for _ in names) + '\n')
    return str(path)


def _request(engine, folder):
    folder.mkdir(parents=True, exist_ok=True)
    box = {'center': [1, 2, 3], 'size': [12, 14, 16]}
    if engine in {'mmpbsa', 'mmgbsa'}:
        inp = {role: _topology(folder / (role + '.prmtop'), names) for role, names in
               [('solvated_topology', ['CA', 'N', 'C1', 'O1', 'OW']),
                ('complex_topology', ['CA', 'N', 'C1', 'O1']),
                ('receptor_topology', ['CA', 'N']), ('ligand_topology', ['C1', 'O1'])]}
        trajectory = folder / 'trajectory.mdcrd'
        trajectory.write_text('Synthetic contract fixture, never a real MD ensemble\n1 2 3\n')
        inp['trajectory'] = str(trajectory)
        inp['review'] = {'equilibrated_trajectory': True, 'matching_topologies_and_atom_order': True,
                         'force_field': 'ff14SB', 'ligand_parameterization': 'GAFF2 fixture',
                         'sampling_notes': 'Synthetic review attestation for contract tests only.',
                         'reviewed_input_sha256': {role: science.sha256(Path(path)) for role, path in inp.items()}}
        return inp, {'end_frame': 6, 'interval': 2}
    if engine == 'planet':
        return {'receptor': _pdb(folder / 'protein.pdb', alpha_carbon=True),
                'pocket_reference': _sdf(folder / 'reference.sdf'),
                'ligands': _sdf(folder / 'query.sdf', [('ethanol', 'CCO', {}), ('ethylamine', 'CCN', {})], is3d=False),
                'preparation_reviewed': True}, {}
    if engine == 'biocentric_affinity':
        return {'receptor': _pdb(folder / 'protein.pdb', alpha_carbon=True),
                'ligands': _sdf(folder / 'query.sdf', [('ethanol', 'CCO', {}), ('ethylamine', 'CCN', {})]),
                'protein_id': 'fixture_protein_chain_A', 'preparation_reviewed': True}, {'device': 'cpu'}
    inp = {'receptor': _pdb(folder / ('receptor.pdb' if engine.startswith('gnina') else 'receptor.pdbqt')),
           'preparation_reviewed': True}
    if engine == 'unidock':
        inp['ligands'] = [_pdb(folder / (name + '.pdbqt'), ligand=True) for name in ('ligand_a', 'ligand_b')]
        return inp, {'box': box, 'num_modes': 2}
    inp['ligand'] = (_sdf(folder / 'query.sdf') if engine.startswith('gnina') else
                     _pdb(folder / 'query.pdbqt', ligand=True))
    if engine == 'autodock4':
        header = ('GRID_PARAMETER_FILE protein.gpf\nGRID_DATA_FILE protein.maps.fld\nMACROMOLECULE receptor.pdbqt\n'
                  'SPACING 0.375\nNELEMENTS 40 40 40\nCENTER 1 2 3\n0.0\n')
        for name in ('C', 'e', 'd'):
            (folder / f'protein.{name}.map').write_text(header)
        (folder / 'protein.maps.fld').write_text('# AVS field fixture\n')
        inp.update(grid_field=str(folder / 'protein.maps.fld'), atom_maps={'C': str(folder / 'protein.C.map')},
                   electrostatic_map=str(folder / 'protein.e.map'), desolvation_map=str(folder / 'protein.d.map'),
                   maps_match_receptor_reviewed=True)
        return inp, {'ga_runs': 2, 'ga_num_evals': 1000}
    return inp, {} if engine == 'gnina_affinity' else {'box': box, 'exhaustiveness': 1, 'num_modes': 2}


def _ready_catalog(monkeypatch):
    original = engines.engine_catalog
    def catalog():
        return [{**row, 'available': True, 'blockers': [], 'executable': '/fixture/bin/' + row['id'],
                 'hardware': {**row['hardware'], 'gpu_ready': True}} for row in original()]
    monkeypatch.setattr(engines, 'engine_catalog', catalog)


@pytest.mark.parametrize('engine', engines.ENGINE_IDS)
def test_every_engine_has_bounded_reviewable_plan_without_execution(engine, tmp_path):
    inp, params = _request(engine, tmp_path / 'input')
    out = tmp_path / 'output'
    plan = engines.plan_engine_run(engine, inp, params, out)
    assert plan['status'] == 'blocked' and plan['blockers']
    assert not out.exists()
    assert len(plan['commands']) == 1 and isinstance(plan['commands'][0], list)
    assert plan['validated_inputs'] and all(len(record['sha256']) == 64 for record in plan['validated_inputs'].values())
    assert 0 < plan['parameters']['timeout_seconds'] <= 3600
    for metric in plan['output_metrics']:
        assert metric['experimental_measurement'] is False
        assert metric['unit'] in {'kcal/mol', 'pK', 'dimensionless'}
        assert metric['direction'] in {'lower_is_better', 'higher_is_better'}


@pytest.mark.parametrize('engine', engines.ENGINE_IDS)
def test_every_engine_rejects_unbounded_or_arbitrary_options(engine, tmp_path):
    inp, params = _request(engine, tmp_path / 'input')
    for extra in ({'extra_args': ['--anything']}, {'timeout_seconds': float('inf')}, {'timeout_seconds': 3601}):
        with pytest.raises(ValidationError):
            engines.plan_engine_run(engine, inp, {**params, **extra}, tmp_path / 'output')


@pytest.mark.parametrize('engine', engines.ENGINE_IDS)
def test_every_engine_rejects_missing_scientific_review(engine, tmp_path):
    inp, params = _request(engine, tmp_path / 'input')
    if engine in {'mmpbsa', 'mmgbsa'}:
        inp['review']['equilibrated_trajectory'] = False
    else:
        inp['preparation_reviewed'] = False
    with pytest.raises(ValidationError):
        engines.plan_engine_run(engine, inp, params, tmp_path / 'output')


@pytest.mark.parametrize('engine', engines.ENGINE_IDS)
@pytest.mark.parametrize('damage', ['empty', 'wrong_extension'])
def test_native_file_contract_rejects_wrong_or_empty_input(engine, damage, tmp_path):
    inp, params = _request(engine, tmp_path / 'input')
    role = 'trajectory' if engine in {'mmpbsa', 'mmgbsa'} else 'receptor'
    original = Path(inp[role])
    if damage == 'empty':
        original.write_bytes(b'')
    else:
        replacement = original.with_suffix('.txt')
        replacement.write_bytes(original.read_bytes())
        inp[role] = str(replacement)
    with pytest.raises(ValueError, match='requires a nonempty'):
        engines.plan_engine_run(engine, inp, params, tmp_path / 'output')


@pytest.mark.parametrize('engine', ['vina', 'unidock'])
def test_vina_family_preserves_multiple_poses_and_scientific_notation(engine, tmp_path):
    path = tmp_path / 'poses.pdbqt'
    path.write_text('MODEL 1\nREMARK VINA RESULT: -8.125e+00 0.000 0.000\nENDMDL\n'
                    'MODEL 2\nREMARK VINA RESULT: +.25E1 1.000 2.000\nENDMDL\n')
    rows = engines.parse_engine_output(engine, path)
    assert [row['pose_index'] for row in rows] == [1, 2]
    assert [row['metrics'][0]['value'] for row in rows] == [-8.125, 2.5]
    assert all(row['source']['sha256'] == science.sha256(path) for row in rows)
    assert {row['metrics'][0]['name'] for row in rows} == {'vina_score' if engine == 'vina' else 'unidock_vina_score'}


@pytest.mark.parametrize('engine', ['vina', 'unidock'])
@pytest.mark.parametrize('bad_line', ['-7bad 0 0', 'nan 0 0', '1e999 0 0', '-7 0', '-7 0 nan'])
def test_vina_family_rejects_partial_or_malformed_pose_in_otherwise_valid_output(engine, bad_line, tmp_path):
    path = tmp_path / 'poses.pdbqt'
    path.write_text('REMARK VINA RESULT: -8 0 0\nREMARK VINA RESULT: ' + bad_line + '\n')
    with pytest.raises(ValueError):
        engines.parse_engine_output(engine, path)


def test_autodock_repeated_runs_preserve_order_ignore_nonpose_summary(tmp_path):
    path = tmp_path / 'docking.dlg'
    path.write_text('USER Estimated Free Energy of Binding = -999 kcal/mol\n'
                    'DOCKED: MODEL 1\nDOCKED: USER Estimated Free Energy of Binding = -8.25E+00 kcal/mol [=(1)+(2)+(3)-(4)]\n'
                    'DOCKED: ENDMDL\nDOCKED: MODEL 2\nDOCKED: USER Estimated Free Energy of Binding = +.50 kcal/mol\nDOCKED: ENDMDL\n')
    rows = engines.parse_engine_output('autodock4', path)
    assert [row['pose_index'] for row in rows] == [1, 2]
    assert [row['metrics'][0]['value'] for row in rows] == [-8.25, .5]
    assert all(row['metrics'][0]['name'] == 'autodock4_score' for row in rows)


@pytest.mark.parametrize('value', ['nan kcal/mol', '-7 kJ/mol', '-7bad kcal/mol', '1e999 kcal/mol'])
def test_autodock_does_not_silently_drop_invalid_run(value, tmp_path):
    path = tmp_path / 'docking.dlg'
    path.write_text('DOCKED: USER Estimated Free Energy of Binding = -8 kcal/mol\n'
                    'DOCKED: USER Estimated Free Energy of Binding = ' + value + '\n')
    with pytest.raises(ValueError):
        engines.parse_engine_output('autodock4', path)


@pytest.mark.parametrize('engine', ['gnina', 'gnina_affinity'])
def test_gnina_multiple_poses_keep_three_distinct_score_meanings(engine, tmp_path):
    path = tmp_path / 'poses.sdf'
    _sdf(path, [('first', 'CCO', {'minimizedAffinity': -8, 'CNNscore': 0, 'CNNaffinity': 5}),
                ('second', 'CCO', {'minimizedAffinity': -6, 'CNNscore': 1, 'CNNaffinity': 7})])
    rows = engines.parse_engine_output(engine, path)
    assert [row['name'] for row in rows] == ['first', 'second']
    assert [row['pose_index'] for row in rows] == [1, 2]
    values = [{metric['name']: metric for metric in row['metrics']} for row in rows]
    assert values[0]['minimizedAffinity']['value'] < values[1]['minimizedAffinity']['value']
    assert values[0]['CNNaffinity']['value'] < values[1]['CNNaffinity']['value']
    assert values[0]['minimizedAffinity']['direction'] == 'lower_is_better'
    assert values[0]['CNNaffinity']['direction'] == values[0]['CNNscore']['direction'] == 'higher_is_better'


@pytest.mark.parametrize('engine', ['gnina', 'gnina_affinity'])
@pytest.mark.parametrize('field,bad', [('CNNscore', None), ('CNNaffinity', None), ('minimizedAffinity', None),
                                    ('CNNscore', '-.001'), ('CNNscore', '1.001'), ('CNNscore', 'nan'),
                                    ('CNNaffinity', 'inf'), ('minimizedAffinity', '-inf'), ('CNNaffinity', 'not-a-number')])
def test_gnina_rejects_bad_second_pose_instead_of_partial_success(engine, field, bad, tmp_path):
    path = tmp_path / 'poses.sdf'
    properties = {'minimizedAffinity': -7, 'CNNscore': .8, 'CNNaffinity': 6}
    broken = {**properties, field: bad}
    if bad is None:
        del broken[field]
    _sdf(path, [('valid', 'CCO', properties), ('broken', 'CCN', broken)])
    with pytest.raises(ValueError):
        engines.parse_engine_output(engine, path)


@pytest.mark.parametrize('engine', ['gnina', 'gnina_affinity'])
def test_gnina_invalid_molecular_record_is_not_skipped(engine, tmp_path):
    path = tmp_path / 'poses.sdf'
    _sdf(path, [('valid', 'CCO', {'minimizedAffinity': -7, 'CNNscore': .8, 'CNNaffinity': 6})])
    with path.open('a') as stream:
        stream.write('Broken record\nnot a molblock\n\ninvalid\n$$$$\n')
    with pytest.raises(ValueError, match='invalid SDF'):
        engines.parse_engine_output(engine, path)


def _amber(engine, values, units='kcal/mole', *, alternate=False):
    heading = 'POISSON BOLTZMANN' if engine == 'mmpbsa' else 'GENERALIZED BORN'
    label = 'DELTA G binding =' if alternate else 'DELTA TOTAL'
    return (f'All units are reported in {units}.\n{heading}:\nComplex:\nTOTAL -9999 1 1\n'
            f'Differences (Complex - Receptor - Ligand):\n{label} {values}\n')


@pytest.mark.parametrize('engine', ['mmpbsa', 'mmgbsa'])
@pytest.mark.parametrize('alternate', [False, True])
def test_amber_scientific_notation_and_uncertainty_remain_endpoint_estimates(engine, alternate, tmp_path):
    path = tmp_path / 'FINAL_RESULTS.dat'
    path.write_text(_amber(engine, '-2.125E1 ' + ('+/- ' if alternate else '') + '3.0e0 4.0e-1', alternate=alternate))
    row = engines.parse_engine_output(engine, path)[0]
    metric = row['metrics'][0]
    assert metric['value'] == -21.25 and metric['standard_deviation'] == 3 and metric['standard_error'] == .4
    assert metric['experimental_measurement'] is False and row['entropy_included'] is False


@pytest.mark.parametrize('engine', ['mmpbsa', 'mmgbsa'])
@pytest.mark.parametrize('values', ['-12.5 3.2', '-12.5 -3 .2', '-12.5 3 -.2', '-12.5 1e999 .2',
                                  '-12.5 3 1e999', 'nan 3 .2', '-12.5 3 .2extra'])
def test_amber_rejects_truncated_nonfinite_or_invalid_uncertainties(engine, values, tmp_path):
    path = tmp_path / 'FINAL_RESULTS.dat'
    path.write_text(_amber(engine, values))
    with pytest.raises(ValueError):
        engines.parse_engine_output(engine, path)


@pytest.mark.parametrize('engine', ['mmpbsa', 'mmgbsa'])
@pytest.mark.parametrize('units', ['kJ/mol', 'kcal', 'unknown'])
def test_amber_never_relabels_wrong_energy_units_as_kcal(engine, units, tmp_path):
    path = tmp_path / 'FINAL_RESULTS.dat'
    path.write_text(_amber(engine, '-20 3 .2', units))
    with pytest.raises(ValueError, match='units'):
        engines.parse_engine_output(engine, path)


@pytest.mark.parametrize('engine', ['mmpbsa', 'mmgbsa'])
def test_amber_requires_requested_method_and_unambiguous_difference(engine, tmp_path):
    path = tmp_path / 'FINAL_RESULTS.dat'
    other = 'mmgbsa' if engine == 'mmpbsa' else 'mmpbsa'
    path.write_text(_amber(other, '-20 3 .2'))
    with pytest.raises(ValueError, match='requested'):
        engines.parse_engine_output(engine, path)
    path.write_text(_amber(engine, '-20 3 .2') + _amber(engine, '-10 4 .3'))
    with pytest.raises(ValueError, match='ambiguous|multiple'):
        engines.parse_engine_output(engine, path)


@pytest.mark.parametrize('engine', ['mmpbsa', 'mmgbsa'])
def test_amber_missing_units_and_duplicate_totals_require_review(engine, tmp_path):
    path = tmp_path / 'FINAL_RESULTS.dat'
    path.write_text(_amber(engine, '-20 3 .2').split('\n', 1)[1])
    with pytest.raises(ValueError, match='units'):
        engines.parse_engine_output(engine, path)
    path.write_text(_amber(engine, '-20 3 .2', units='kcal/mol'))
    assert engines.parse_engine_output(engine, path)[0]['metrics'][0]['value'] == -20
    with path.open('a') as stream:
        stream.write('DELTA TOTAL -19 2 .1\n')
    with pytest.raises(ValueError, match='ambiguous|multiple'):
        engines.parse_engine_output(engine, path)


def test_planet_native_csv_preserves_quoted_names_and_molecular_rows(tmp_path):
    path = tmp_path / 'result.csv'
    path.write_text(',mol_name,SMILES,PLANET_affinity\n0,"Ligand, one",CCO,4.2\n1,ligand_two,NCC,6.7e0\n')
    rows = engines.parse_engine_output('planet', path)
    assert [row['name'] for row in rows] == ['Ligand, one', 'ligand_two']
    assert [row['smiles'] for row in rows] == ['CCO', 'NCC']
    assert all(row['generated_docking_pose'] is False for row in rows)
    assert all(row['metrics'][0]['unit'] == 'pK' and row['metrics'][0]['direction'] == 'higher_is_better' for row in rows)


@pytest.mark.parametrize('row', ['bad,not-smiles,6.0', 'blank,,6.0', 'broken,CCO', 'overflow,CCO,1e999',
                               'wrong,CCO,5,unexpected_column'])
def test_planet_rejects_malformed_or_invalid_molecular_row_after_valid_row(row, tmp_path):
    path = tmp_path / 'result.csv'
    path.write_text('mol_name,SMILES,PLANET_affinity\nvalid,CCO,5.1\n' + row + '\n')
    with pytest.raises(ValueError):
        engines.parse_engine_output('planet', path)


@pytest.mark.parametrize('problem', ['wrong_identity', 'duplicate_identity', 'missing_sdf', 'wrong_sdf_identity'])
def test_planet_fake_execution_requires_complete_matching_molecules(problem, tmp_path, monkeypatch):
    inp, params = _request('planet', tmp_path / 'input')
    _ready_catalog(monkeypatch)
    def runner(argv, cwd, log, timeout, cancelled):
        names = 'ethanol,CCO,5\nethylamine,CCN,6\n'
        if problem == 'wrong_identity':
            names = 'ethanol,CCO,5\nsubstituted,CCC,6\n'
        elif problem == 'duplicate_identity':
            names = 'ethanol,CCO,5\nethanol_again,CCO,6\n'
        (cwd / 'result.csv').write_text('mol_name,SMILES,PLANET_affinity\n' + names)
        if problem != 'missing_sdf':
            (cwd / 'result.sdf').write_bytes(Path(inp['ligands']).read_bytes())
        if problem == 'wrong_sdf_identity':
            _sdf(cwd / 'result.sdf', [('ethanol', 'CCO', {}), ('substitute', 'CCC', {})], is3d=False)
        log.write_text('Fake runner; no PLANET model was loaded or executed.\n')
        return {'returncode': 0}
    monkeypatch.setattr(science, 'run_process', runner)
    result = engines.run_engine('planet', inp, params, tmp_path / 'output')
    assert result['status'] == 'failed' and result['results'] == []
    assert result['executed'] is True  # fake process boundary only, not a scientific runtime test


def test_planet_fake_execution_accepts_equivalent_smiles_and_complete_2d_sdf(tmp_path, monkeypatch):
    inp, params = _request('planet', tmp_path / 'input')
    _ready_catalog(monkeypatch)
    def runner(argv, cwd, log, timeout, cancelled):
        (cwd / 'result.csv').write_text('mol_name,SMILES,PLANET_affinity\nethylamine,NCC,6\nethanol,OCC,5\n')
        (cwd / 'result.sdf').write_bytes(Path(inp['ligands']).read_bytes())
        log.write_text('Fake runner; fixtures only, no model was executed.\n')
        return {'returncode': 0}
    monkeypatch.setattr(science, 'run_process', runner)
    result = engines.run_engine('planet', inp, params, tmp_path / 'output')
    assert result['status'] == 'completed' and len(result['results']) == 2
    assert all(row['generated_docking_pose'] is False for row in result['results'])
    assert 'Not established' in result['scientific_validation']


@pytest.mark.parametrize('input_smiles,output_smiles', [('F[C@H](Cl)Br', 'F[C@@H](Cl)Br'), ('CC[NH3+]', 'CCN')])
def test_planet_identity_validation_preserves_stereochemistry_and_charge(input_smiles, output_smiles, tmp_path, monkeypatch):
    inp, params = _request('planet', tmp_path / 'input')
    _sdf(Path(inp['ligands']), [('reviewed-state', input_smiles, {})], is3d=False)
    _ready_catalog(monkeypatch)
    def runner(argv, cwd, log, timeout, cancelled):
        (cwd / 'result.csv').write_text('mol_name,SMILES,PLANET_affinity\nchanged-state,' + output_smiles + ',5\n')
        (cwd / 'result.sdf').write_bytes(Path(inp['ligands']).read_bytes())
        log.write_text('Fake output with changed stereochemistry or charge.\n')
        return {'returncode': 0}
    monkeypatch.setattr(science, 'run_process', runner)
    result = engines.run_engine('planet', inp, params, tmp_path / 'output')
    assert result['status'] == 'failed' and result['results'] == []
    assert 'molecular identities' in result['error']['message']


def test_unidock_fake_batch_preserves_per_file_sources_and_multiple_poses(tmp_path, monkeypatch):
    inp, params = _request('unidock', tmp_path / 'input')
    _ready_catalog(monkeypatch)
    def runner(argv, cwd, log, timeout, cancelled):
        for ligand in inp['ligands']:
            (cwd / 'poses' / (Path(ligand).stem + '_out.pdbqt')).write_text(
                'REMARK VINA RESULT: -8 0 0\nREMARK VINA RESULT: -7 1 2\n')
        log.write_text('Fake runner; no GPU computation.\n')
        return {'returncode': 0}
    monkeypatch.setattr(science, 'run_process', runner)
    result = engines.run_engine('unidock', inp, params, tmp_path / 'output')
    assert result['status'] == 'completed' and len(result['results']) == 4
    assert len({row['source']['path'] for row in result['results']}) == 2
    assert [row['pose_index'] for row in result['results']] == [1, 2, 1, 2]


def test_unidock_fake_batch_cannot_report_partial_library_as_completed(tmp_path, monkeypatch):
    inp, params = _request('unidock', tmp_path / 'input')
    _ready_catalog(monkeypatch)
    def runner(argv, cwd, log, timeout, cancelled):
        (cwd / 'poses' / 'ligand_a_out.pdbqt').write_text('REMARK VINA RESULT: -8 0 0\n')
        log.write_text('Fake partial output.\n')
        return {'returncode': 0}
    monkeypatch.setattr(science, 'run_process', runner)
    result = engines.run_engine('unidock', inp, params, tmp_path / 'output')
    assert result['status'] == 'failed' and result['results'] == []
    assert 'one output file per ligand' in result['error']['message']


def test_unidock_fake_batch_rejects_same_count_with_substituted_filename(tmp_path, monkeypatch):
    inp, params = _request('unidock', tmp_path / 'input')
    _ready_catalog(monkeypatch)
    def runner(argv, cwd, log, timeout, cancelled):
        for name in ('ligand_a_out.pdbqt', 'unrelated_out.pdbqt'):
            (cwd / 'poses' / name).write_text('REMARK VINA RESULT: -8 0 0\n')
        log.write_text('Fake substituted batch output.\n')
        return {'returncode': 0}
    monkeypatch.setattr(science, 'run_process', runner)
    result = engines.run_engine('unidock', inp, params, tmp_path / 'output')
    assert result['status'] == 'failed' and result['results'] == []
    assert 'input-bound filename' in result['error']['message']


@pytest.mark.parametrize('engine', ['gnina', 'gnina_affinity'])
def test_gnina_cpu_and_gpu_requests_have_explicit_hardware_boundary(engine, tmp_path):
    inp, params = _request(engine, tmp_path / 'input')
    cpu = engines.plan_engine_run(engine, inp, params, tmp_path / 'cpu')
    gpu = engines.plan_engine_run(engine, inp, {**params, 'device': 'gpu'}, tmp_path / 'gpu')
    assert '--no_gpu' in cpu['commands'][0] and '--no_gpu' not in gpu['commands'][0]
    assert any('NVIDIA CUDA' in message for message in gpu['blockers'])
    assert ('--score_only' in cpu['commands'][0]) == (engine == 'gnina_affinity')
    assert ('--center_x' in cpu['commands'][0]) == (engine == 'gnina')


@pytest.mark.parametrize('engine', ['mmpbsa', 'mmgbsa'])
@pytest.mark.parametrize('params', [{'start_frame': 5, 'end_frame': 4}, {'end_frame': 2, 'interval': 2}, {'end_frame': 10000, 'interval': 1}])
def test_mm_frame_selection_must_be_bounded_ensemble(engine, params, tmp_path):
    inp, _ = _request(engine, tmp_path / 'input')
    with pytest.raises(ValidationError, match='2–1000'):
        engines.plan_engine_run(engine, inp, params, tmp_path / 'output')


@pytest.mark.parametrize('problem', ['missing_atom_map', 'extra_atom_map', 'unsupported_torsions', 'missing_grid_metadata'])
def test_autodock_native_map_and_torsion_contract(problem, tmp_path):
    inp, params = _request('autodock4', tmp_path / 'input')
    if problem == 'missing_atom_map':
        inp['atom_maps'] = {'N': inp['atom_maps']['C']}
    elif problem == 'extra_atom_map':
        inp['atom_maps']['N'] = inp['atom_maps']['C']
    elif problem == 'unsupported_torsions':
        ligand = Path(inp['ligand'])
        ligand.write_text(ligand.read_text().replace('TORSDOF 0', 'TORSDOF 33'))
    else:
        path = Path(inp['electrostatic_map'])
        path.write_text(path.read_text().replace('SPACING 0.375', 'UNKNOWN 0.375'))
    with pytest.raises(ValueError):
        engines.plan_engine_run('autodock4', inp, params, tmp_path / 'output')


@pytest.mark.parametrize('problem', ['pocket_2d', 'protein_no_alpha_carbon', 'invalid_query_record'])
def test_planet_pocket_and_chemical_query_contract(problem, tmp_path):
    inp, params = _request('planet', tmp_path / 'input')
    if problem == 'pocket_2d':
        _sdf(Path(inp['pocket_reference']), is3d=False)
    elif problem == 'protein_no_alpha_carbon':
        _pdb(Path(inp['receptor']))
    else:
        with Path(inp['ligands']).open('a') as stream:
            stream.write('Broken record\nnot a molblock\n\ninvalid\n$$$$\n')
    with pytest.raises(ValueError):
        engines.plan_engine_run('planet', inp, params, tmp_path / 'output')
