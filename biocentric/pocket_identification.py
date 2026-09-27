"""Fixed, auditable pocket adapters; predicted pockets never use the control ligand."""
from pathlib import Path
import os
import shutil

from . import science

REFERENCE = {'tool': 'reference_box', 'parameters': {'padding': 4}}


def executable():
    local = Path(__file__).resolve().parents[1] / '.tools/linux/fpocket-env/bin/fpocket'
    path = os.environ.get('BIOCENTRIC_FPOCKET') or shutil.which('fpocket') or str(local)
    return str(Path(path).resolve()) if Path(path).is_file() and os.access(path, os.X_OK) else None


def blocker(method):
    tool = method['tool']
    if tool not in {'reference_box', 'fpocket'}:
        return 'Pocket adapter not connected: ' + tool
    if tool == 'fpocket' and not executable():
        return 'fpocket runtime is not installed'
    params = method.get('parameters', {})
    for key, default, low, high in [('padding', 4, 1, 12), ('pocket_rank', 1, 1, 20), ('min_size', 20, 10, 40)]:
        value = params.get(key, default)
        if isinstance(value, bool) or not isinstance(value, int) or not low <= value <= high:
            return f'Invalid pocket parameter: {key}'
    return None


def identify(method, receptor, reference, folder, cancelled):
    import numpy as np
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    reason = blocker(method)
    if reason:
        raise ValueError(reason)
    params = method.get('parameters', {})
    if method['tool'] == 'reference_box':
        box = science.derive_box(reference, padding=params.get('padding', 4))
        if 'min_size' in params:
            box['size'] = [max(side, params['min_size']) for side in box['size']]
            box['minimum_size_angstrom'] = params['min_size']
    else:
        # Use protein atoms only. The selected rank is frozen before validation.
        source = folder / 'protein.pdb'
        shutil.copyfile(receptor, source)
        science.run_process([executable(), '-f', str(source)], folder, folder / 'fpocket.log', 120, cancelled)
        rank = params.get('pocket_rank', 1)
        vertices = folder / 'protein_out' / 'pockets' / f'pocket{rank}_vert.pqr'
        if not vertices.is_file():
            raise ValueError(f'fpocket did not find pocket rank {rank}')
        xyz = np.array([[float(line[30:38]), float(line[38:46]), float(line[46:54])]
                        for line in vertices.read_text().splitlines() if line.startswith(('ATOM  ', 'HETATM'))])
        if not len(xyz) or not np.isfinite(xyz).all():
            raise ValueError('Pocket coordinates are empty or invalid')
        size = np.maximum(xyz.max(axis=0) - xyz.min(axis=0) + 2 * params.get('padding', 4), params.get('min_size', 20))
        if max(size) > 50:
            raise ValueError('Predicted pocket box exceeds 50 Å; review pocket parameters')
        box = {'center': ((xyz.max(axis=0) + xyz.min(axis=0)) / 2).tolist(), 'size': size.tolist(),
               'units': 'angstrom', 'method': 'fpocket_ranked_alpha_spheres', 'pocket_rank': rank,
               'vertices_sha256': science.sha256(vertices), 'receptor_sha256': science.sha256(source),
               'coordinate_frame': 'input receptor, no alignment/transformation'}
    science._write_json(folder / 'pocket.json', {'method': method, 'box': box})
    return box
