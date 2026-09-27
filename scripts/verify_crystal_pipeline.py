"""Bounded real crystal-pocket redocking and native fixed-crop parity check."""
import csv
import json
import os
from pathlib import Path
from rdkit import Chem
from biocentric import science, discovery_benchmark as bench, discovery_search_runner as search
from biocentric.discovery_affinity import score_groups


def main():
    root = Path(__file__).resolve().parents[1]
    source = root / '.biocentric/evaluation-runs/8d9390d9374a4ef6adb8ddfe11fc4461'
    out = root / 'runs/crystal-pocket-correction'
    out.mkdir(exist_ok=True)
    for filename in ('receptor.pdb', 'reference.sdf'):
        (out / filename).write_bytes((source / filename).read_bytes())
    receptor = out / 'receptor.pdb'
    reference = next(iter(Chem.SDMolSupplier(str(out / 'reference.sdf'), removeHs=False)))
    if (out / 'control.json').exists():
        result = json.loads((out / 'control.json').read_text())
    else:
        pdbqt = Path(science.prepare_receptor(receptor, out / 'receptor-preparation', 'meeko', reviewed=True)['output'])
        pair = {'id': 'crystal-vina-biocentric', 'pocket': {'tool': 'reference_box', 'parameters': {'padding': 4}},
                'docking': {'tool': 'vina', 'parameters': {'exhaustiveness': 2, 'num_modes': 1, 'seed': 42, 'cpu': 4, 'max_evals': 5000}},
                'scoring': {'tool': 'biocentric_affinity', 'parameters': {'device': 'cuda', 'cpu': 2}}}
        record = {'compound_id': 'crystal:300', 'kind': 'reference', 'canonical_smiles': science._canonical(reference)}
        result = search.evaluate(pair, record, 0, out, receptor, pdbqt, reference, {}, {}, {}, {}, lambda: None, lambda: False, lambda event: None)
        science._write_json(out / 'control.json', result)
        if result['status'] != 'completed':
            raise ValueError(str(result))
    docked = next(iter(Chem.SDMolSupplier(str(out / result['pose_path']), removeHs=False)))
    native = out / 'native-shared-parity'
    poses = native / 'poses' / 'complex'
    poses.mkdir(parents=True)
    for name, mol in [('rank1_crystal', reference), ('rank2_docked', docked)]:
        Chem.MolToPDBFile(mol, str(poses / (name + '.pdb')))
    model = Path(os.environ['BIOCENTRIC_MODEL_DIR']) / 'Prediction'
    command = [os.environ['BIOCENTRIC_MODEL_PYTHON'], str(model / 'run_pipeline.py'), '--pocket', str(receptor),
               '--poses', str(poses.parent), '--out', str(native / 'scores.csv'), '--device', 'cuda', '--per-pose']
    science.run_process(command, model, native / 'engine.log', 180)
    with (native / 'scores_per_pose.csv').open() as handle:
        native_rows = list(csv.DictReader(handle))
    adapter = score_groups({'rank1_crystal': [reference], 'rank2_docked': [docked]}, receptor, out / 'adapter-shared-parity',
                           {'device': 'cuda', 'cpu': 2}, lambda: False, lambda event: None, crop_reference=reference)
    differences = []
    for name, group in adapter.items():
        row = next(row for row in native_rows if name in str(row))
        predicted = group['predictions'][0]
        differences.append({'pose': name, 'native_score': float(row['score']), 'app_score': predicted['score'],
                            'absolute_difference': abs(float(row['score']) - predicted['score'])})
        if group['errors'] or differences[-1]['absolute_difference'] > 1e-5:
            raise ValueError('Native shared-crop parity failed: ' + str(differences[-1]))
    report = {'pdb_id': '1O41', 'component': '300', 'scope': 'Crystal control and software parity; not affinity validation',
              'crystal_rmsd_angstrom': result['rmsd_angstrom'], 'physical_pass': result['physical_pass'],
              'docked_biocentric_score': result['score'], 'crystal_biocentric_score': result['reference_score'],
              'native_parity': differences, 'tolerance': 1e-5}
    science._write_json(out / 'verification.json', report)
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
