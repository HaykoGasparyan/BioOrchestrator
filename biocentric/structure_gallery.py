"""Offline, single-canvas structure gallery with actual geometry thumbnails."""
from __future__ import annotations

import hashlib
from pathlib import Path

from .viewer import VIEWER_SHA256, VIEWER_VERSION, script_json


def gallery_html(items, target):
    scenes = []
    for record, geometry in items:
        if record['id'] != geometry['id'] or hashlib.sha256(geometry['text'].encode()).hexdigest() != geometry['sha256']:
            raise ValueError('Experimental structure and coordinates do not match.')
        ligands = [{'value': r['id'], 'label': r['id'], 'name': r.get('name', '')} for r in record.get('compounds', [])]
        for partner in record.get('partners', []):
            for chain in partner['chains']:
                ligands.append({'value': 'chain:' + chain, 'label': 'Chain ' + chain, 'name': partner['name']})
        scenes.append({'id': record['id'], 'text': geometry['text'], 'format': geometry.get('format', 'cif'),
                       'chains': record['target_chains'], 'models': record.get('models', 1), 'ligands': ligands,
                       'title': record.get('title', ''), 'resolution': record.get('resolution'),
                       'methods': record.get('methods', []), 'mutations': record.get('mutations', [])})
    library = (Path(__file__).parent / 'assets' / f'3Dmol-{VIEWER_VERSION}-min.js').read_bytes()
    if hashlib.sha256(library).hexdigest() != VIEWER_SHA256:
        raise ValueError('3D viewer integrity check failed.')
    template = Path(__file__).with_name('structure_gallery.html').read_text()
    return template.replace('/*__LIBRARY__*/', library.decode().replace('</script', '<\\/script')).replace(
        '/*__SCENES__*/', script_json({'scenes': scenes, 'target': target}))
