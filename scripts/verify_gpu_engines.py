"""Bounded native CLI and typed-adapter GPU checks on the bundled 1IEP example."""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import threading
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from biocentric import engines, science


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--engine', choices=['unidock', 'gnina', 'gnina_affinity'], required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    base = Path(__file__).resolve().parents[1]
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    binary = base / '.tools/linux' / ('unidock' if args.engine == 'unidock' else 'gnina')
    os.environ['BIOCENTRIC_' + ('UNIDOCK' if args.engine == 'unidock' else 'GNINA')] = str(binary)
    from rdkit import Chem
    from meeko import MoleculePreparation, PDBQTWriterLegacy
    reference = base / 'examples/1iep/originals/1iep_ligand.sdf'
    prepared = out / 'inputs'; prepared.mkdir()
    molecule = Chem.SDMolSupplier(str(reference), removeHs=False)[0]
    molecule = Chem.AddHs(molecule, addCoords=True)
    setups = MoleculePreparation(charge_model='gasteiger').prepare(molecule)
    if len(setups) != 1:
        raise ValueError('Reference preparation did not produce exactly one reviewed state')
    pdbqt, success, reason = PDBQTWriterLegacy.write_string(setups[0])
    if not success:
        raise ValueError(reason)
    (prepared / 'imatinib.pdbqt').write_text(pdbqt)
    with Chem.SDWriter(str(prepared / 'imatinib.sdf')) as writer:
        writer.write(molecule)
    box = science.derive_box(reference, 4)
    (prepared / 'preparation.json').write_text(json.dumps({'reference_sha256': science.sha256(reference),
        'method': 'Preserve reference coordinates and chemistry; add valence hydrogens; Meeko Gasteiger PDBQT',
        'protonation_enumerated': False, 'box': box}, indent=2) + '\n')
    box = {key: box[key] for key in ['center', 'size']}
    receptor = base / 'examples/1iep/originals/1iep_receptor.pdbqt'
    ligand = prepared / ('imatinib.pdbqt' if args.engine == 'unidock' else 'imatinib.sdf')
    inputs = {'receptor': str(receptor), 'preparation_reviewed': True,
              **({'ligands': [str(ligand)]} if args.engine == 'unidock' else {'ligand': str(ligand)})}
    if args.engine == 'unidock':
        parameters = {'box': box, 'seed': 42, 'num_modes': 1, 'search_mode': 'fast', 'timeout_seconds': 300}
    elif args.engine == 'gnina':
        parameters = {'box': box, 'seed': 42, 'cpu': 1, 'num_modes': 1, 'exhaustiveness': 1, 'device': 'gpu', 'timeout_seconds': 300}
    else:
        parameters = {'cpu': 1, 'device': 'gpu', 'timeout_seconds': 300}
    samples, stop = [], threading.Event()
    def sample():
        while not stop.is_set():
            result = subprocess.run(['nvidia-smi', '--query-gpu=timestamp,name,memory.used,memory.total,utilization.gpu', '--format=csv,noheader,nounits'],
                                    capture_output=True, text=True, timeout=10, env=science._tool_environment())
            processes = subprocess.run(['nvidia-smi', '--query-compute-apps=pid,process_name,used_memory', '--format=csv,noheader,nounits'],
                                       capture_output=True, text=True, timeout=10, env=science._tool_environment())
            samples.append({'elapsed_seconds': time.monotonic() - started, 'returncode': result.returncode, 'gpu': result.stdout.strip(),
                            'compute_processes': processes.stdout.strip()})
            stop.wait(.1)
    started = time.monotonic()
    thread = threading.Thread(target=sample, daemon=True)
    thread.start()
    report = {'engine': args.engine, 'example': 'Bundled official AutoDock Vina 1IEP receptor and cognate imatinib',
              'scope': 'Software/native execution; no independent affinity or enrichment validation', 'inputs': inputs, 'parameters': parameters,
              'source': 'https://github.com/dptech-corp/Uni-Dock/releases/tag/1.1.0' if args.engine == 'unidock' else 'https://github.com/gnina/gnina/releases/tag/v1.3.2'}
    try:
        report['version'] = subprocess.run([str(binary), '--version'], capture_output=True, text=True, timeout=30, env=science._tool_environment()).stdout.strip()
        native = out / 'native'; native.mkdir()
        argv = [str(binary), '--receptor', str(receptor)]
        if args.engine == 'unidock':
            (native / 'poses').mkdir()
            argv += ['--gpu_batch', str(ligand), '--scoring', 'vina', '--search_mode', 'fast', '--seed', '42', '--num_modes', '1', '--dir', str(native / 'poses')]
        else:
            argv += ['--ligand', str(ligand), '--out', str(native / 'poses.sdf'), '--cpu', '1', '--scoring', 'vina', '--cnn_scoring', 'rescore']
            argv += ['--score_only'] if args.engine == 'gnina_affinity' else ['--seed', '42', '--exhaustiveness', '1', '--num_modes', '1']
        if args.engine != 'gnina_affinity':
            for name, values in [('center', box['center']), ('size', box['size'])]:
                for axis, value in zip('xyz', values):
                    argv += [f'--{name}_{axis}', str(value)]
        report['native_argv'] = argv
        report['native_execution'] = science.run_process(argv, native, native / 'engine.log', 300)
        paths = list((native / 'poses').glob('*.pdbqt')) if args.engine == 'unidock' else [native / 'poses.sdf']
        report['native_results'] = [row for path in paths for row in engines.parse_engine_output(args.engine, path)]
        report['native_passed'] = bool(report['native_results'])
        result = engines.run_engine(args.engine, inputs, parameters, out / 'adapter')
        report['adapter_status'] = result['status']
        report['adapter_results'] = result['results']
        report['adapter_error'] = result.get('error')
    except Exception as error:
        report['error'] = {'type': type(error).__name__, 'message': str(error)}
    finally:
        stop.set(); thread.join(timeout=12)
        report['elapsed_seconds'] = time.monotonic() - started
        report['gpu_samples'] = samples
        report['output_files'] = [{'path': str(path.relative_to(out)), 'bytes': path.stat().st_size, 'sha256': science.sha256(path)}
                                  for path in sorted(out.rglob('*')) if path.is_file()]
        (out / 'verification.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({key: value for key, value in report.items() if key not in ['gpu_samples','output_files','native_argv']}, indent=2))
    return 0 if report.get('native_passed') and report.get('adapter_status') == 'completed' else 1


if __name__ == '__main__':
    raise SystemExit(main())
