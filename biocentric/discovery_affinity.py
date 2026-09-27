"""Batch fixed-pose BioCentric inference for discovery benchmarks.

Every pose has an explicit identity. Repeated chemistry is permitted only in
the engine's pose-record mode; a failed pose invalidates its whole pose group.
"""
from pathlib import Path
from typing import Literal

from pydantic import BaseModel, ConfigDict, Field

from . import science


class Parameters(BaseModel):
    model_config = ConfigDict(extra='forbid', allow_inf_nan=False)
    device: Literal['cpu', 'cuda'] = 'cuda'
    cpu: int = Field(default=2, ge=1, le=16)
    batch_size: int = Field(default=16, ge=1, le=32)
    timeout_seconds: int = Field(default=300, ge=30, le=3600)


def score_groups(groups, receptor, folder, parameters, cancelled, emit, crop_reference=None):
    """Return scores and provenance by group key; never discard failed poses."""
    from rdkit import Chem
    from .engines import run_engine
    from .biocentric_model import parse_payload

    par = Parameters.model_validate(parameters)
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    source = folder / 'poses.sdf'
    manifest, expected = [], []
    with Chem.SDWriter(str(source)) as writer:
        for group_index, (key, poses) in enumerate(groups.items()):
            if not poses:
                raise ValueError('A BioCentric pose group is empty')
            for pose_index, pose in enumerate(poses):
                name = f'group-{group_index}-pose-{pose_index}'
                mol = Chem.Mol(pose)
                mol.SetProp('_Name', name)
                writer.write(mol)
                identity = {'record_index': len(expected) + 1, 'compound_id': name,
                            'protein_id': 'evaluation_receptor', 'smiles': science._canonical(mol)}
                expected.append(identity)
                manifest.append({**identity, 'group': key, 'pose_index': pose_index,
                                 'original_name': pose.GetProp('_Name') if pose.HasProp('_Name') else None})
    if not expected:
        return {}
    science._write_json(folder / 'pose-manifest.json', {
        'receptor_sha256': science.sha256(Path(receptor)), 'poses_sha256': science.sha256(source),
        'records': manifest})
    extra_inputs = {}
    if crop_reference is not None:
        reference_path = folder / 'crystal-crop.sdf'
        with Chem.SDWriter(str(reference_path)) as writer:
            writer.write(crop_reference)
        extra_inputs['crop_reference'] = str(reference_path)
    native = run_engine('biocentric_affinity',
                        {'receptor': str(receptor), 'ligands': str(source),
                         'protein_id': 'evaluation_receptor', 'preparation_reviewed': True, **extra_inputs},
                        {'device': par.device, 'cpu': par.cpu, 'timeout_seconds': par.timeout_seconds,
                         'max_pairs': len(expected), 'input_mode': 'pose_records'},
                        folder / 'native', cancelled=cancelled, emit=emit)
    science._cancel(cancelled)
    if native.get('status') not in {'completed', 'partial', 'failed'} or not native.get('results'):
        raise ValueError('; '.join(native.get('blockers', [])) or
                         (native.get('error') or {}).get('message') or 'BioCentric batch failed')
    # Revalidate against our own pose manifest, including recovered native rows.
    rows = parse_payload({'schema_version': 1, 'engine_id': 'biocentric_affinity',
                          'predictions': native['results']}, expected,
                         expected_provenance=native['provenance'], expected_device=par.device)
    result = {key: {'scores': [], 'predictions': [], 'errors': [],
                    'artifact': str(folder / 'native' / 'result.json'),
                    'provenance': native['provenance'], 'device': par.device} for key in groups}
    for mapping, row in zip(manifest, rows):
        group = result[mapping['group']]
        group['predictions'].append(row)
        if row['status'] != 'completed':
            group['errors'].append(f"Pose {mapping['pose_index'] + 1}: {row['failure_reason']}")
        else:
            group['scores'].append(row['score'])
    for group in result.values():
        if group['errors']:
            group['scores'] = []
    return result
