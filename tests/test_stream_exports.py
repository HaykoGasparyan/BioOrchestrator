import csv
import hashlib
import json

import pytest

from biocentric.result_paging import page_rows
from biocentric.stream_exports import export_streamed
from test_ligand_results import PID, molecule, pose, register, snapshot


def project(tmp_path, count=5):
    snap = snapshot()
    snap['project']['state'] = 'completed'
    poses = []
    for index in range(count):
        artifact = register(snap, tmp_path, f'{index}.sdf', molecule(str(index)))
        row = pose(compound=f'cmp_{index}', passed=index != 4)
        row.update(pose_sha256=artifact['sha256'], artifact_ids=[artifact['id']])
        poses.append(row)
    snap['data']['screening'] = {'poses': poses}
    snap['candidates'] = [{**poses[0], 'selected': True}]
    return snap


def test_streamed_export_parts_preserve_every_pose_and_hash(tmp_path):
    snap = project(tmp_path)
    out = tmp_path / 'export'
    manifest = export_streamed(snap, tmp_path, PID, out, poses_per_part=2)
    all_parts = [part for part in manifest['sdf_parts'] if part['scope'] == 'all_poses']
    assert [part['records'] for part in all_parts] == [2, 2, 1]
    assert sum((out / part['name']).read_bytes().count(b'$$$$') for part in all_parts) == 5
    for record in manifest['files']:
        assert hashlib.sha256((out / record['name']).read_bytes()).hexdigest() == record['sha256']
    assert len(json.loads((out / 'ranked_ligands.json').read_text())) == 5
    with (out / 'ranked_ligands.csv').open() as handle:
        assert len(list(csv.DictReader(handle))) == 5
    assert sum(part['records'] for part in manifest['sdf_parts'] if part['scope'] == 'shortlist') == 1
    with pytest.raises(ValueError, match='never overwrite'):
        export_streamed(snap, tmp_path, PID, out)


def test_tamper_or_budget_failure_does_not_publish_partial_export(tmp_path):
    snap = project(tmp_path)
    from pathlib import Path
    Path(snap['artifacts'][-1]['path']).write_bytes(molecule('tampered'))
    with pytest.raises(ValueError, match='hash'):
        export_streamed(snap, tmp_path, PID, tmp_path / 'bad')
    assert not (tmp_path / 'bad').exists()
    snap = project(tmp_path / 'second')
    with pytest.raises(ValueError, match='single pose'):
        export_streamed(snap, tmp_path / 'second', PID, tmp_path / 'small', max_part_bytes=1)
    assert not (tmp_path / 'small').exists()
    assert not list(tmp_path.glob('.biocentric-export-*'))


def test_pages_cover_all_rows_without_changing_download_input():
    class UI:
        page = 0
        def selectbox(self, label, values, **kwargs):
            assert list(values) == [0, 1, 2]
            return self.page
        def caption(self, text):
            assert 'Downloads include every recorded row' in text
    ui = UI(); rows = list(range(256)); viewed = []
    for index in range(3):
        ui.page = index
        viewed.extend(page_rows(ui, rows, 'key', 'Results'))
    assert viewed == rows and len(rows) == 256
