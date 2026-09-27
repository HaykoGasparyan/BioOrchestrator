"""Copied immutable history is readable; copied process metadata has no authority."""
import json
import os
import shutil
import sqlite3
import time
from pathlib import Path
from unittest.mock import patch

import pytest

from biocentric.ligand_results import read_registered
from biocentric.relocation import RelocationResolver, configure, file_sha256
from biocentric.service import Service
from biocentric.store import Store


def project(store):
    return store.create({'name': 'Preserved project', 'mode': 'offline_fixture'}, [])


def copied_history(tmp_path):
    old = tmp_path / 'old'
    store = Store(old / '.biocentric')
    pid = project(store)
    content = json.dumps({'path': str(old / 'runs' / 'original.sdf'), 'score': -8.2}).encode()
    artifact = store.artifact(pid, 'original-manifest.json', content)
    store.update(pid, state='completed', data={'historical_path': str(old / 'runs' / 'original.sdf')})
    new = tmp_path / 'new'
    shutil.copytree(old, new)
    return old, new, pid, artifact, content


def test_explicit_mapping_preserves_immutable_bytes_and_database_paths(tmp_path):
    old, new, pid, original, content = copied_history(tmp_path)
    record = configure(new / '.biocentric', [{'old_root': str(old), 'new_root': str(new)}], tmp_path / 'backup')
    service = Service(new / '.biocentric')
    artifact = next(row for row in service.snapshot(pid)['artifacts'] if row['id'] == original['id'])
    assert artifact['path'].startswith(str(new))
    assert artifact['original_path'] == original['path']
    assert record['verified_artifacts'] == 2
    assert service.store.read_artifact(pid, artifact) == content
    assert read_registered(service.snapshot(pid), service.root, pid, artifact['id']) == content
    assert service.store.project(pid)['data']['historical_path'].startswith(str(old))
    with service.store.connect() as db:
        assert db.execute('SELECT path FROM artifacts WHERE id=?', (original['id'],)).fetchone()[0] == original['path']
    backup = json.loads((tmp_path / 'backup' / 'backup.json').read_text())
    assert backup['artifact_records'] == 2
    assert file_sha256(tmp_path / 'backup' / 'evidence.tar.gz') == backup['checksums']['evidence.tar.gz']


def test_mapping_is_required_and_prefix_boundaries_are_respected(tmp_path):
    old, new, pid, original, content = copied_history(tmp_path)
    with pytest.raises(ValueError, match='escapes'):
        Store(new / '.biocentric').read_artifact(pid, original)
    resolver = RelocationResolver([{'old_root': str(old), 'new_root': str(new)}])
    assert resolver.resolve(str(old) + '-different/file') == Path(str(old) + '-different/file')
    assert resolver.resolve(old / 'runs' / 'manifest.json') == new / 'runs' / 'manifest.json'


def test_bad_copy_is_backed_up_but_not_activated(tmp_path):
    old, new, pid, original, content = copied_history(tmp_path)
    copied = new / Path(original['path']).relative_to(old)
    copied.write_bytes(b'corrupted copy')
    with pytest.raises(ValueError, match='SHA256'):
        configure(new / '.biocentric', [{'old_root': str(old), 'new_root': str(new)}], tmp_path / 'backup')
    assert (tmp_path / 'backup' / 'backup.json').exists()
    assert not (new / '.biocentric' / 'relocation.json').exists()


def test_relocation_still_checks_hash_on_each_read(tmp_path):
    old, new, pid, original, content = copied_history(tmp_path)
    configure(new / '.biocentric', [{'old_root': str(old), 'new_root': str(new)}], tmp_path / 'backup')
    store = Store(new / '.biocentric')
    Path(next(row for row in store.artifacts(pid) if row['id'] == original['id'])['path']).write_bytes(b'tampered after configuration')
    with pytest.raises(ValueError, match='integrity'):
        store.read_artifact(pid, original)


def test_mapping_rejects_traversal_and_symlink_escape(tmp_path):
    old, new = tmp_path / 'old', tmp_path / 'new'
    new.mkdir()
    (new / 'escape').symlink_to(tmp_path)
    resolver = RelocationResolver([{'old_root': str(old), 'new_root': str(new)}])
    for path in [old / '..' / 'outside', old / 'escape' / 'outside']:
        with pytest.raises(ValueError):
            resolver.resolve(path)
    with pytest.raises(ValueError):
        RelocationResolver([{'old_root': '/', 'new_root': str(new)}])
    with pytest.raises(ValueError):
        resolver.verified(old / 'file', '')


@pytest.mark.parametrize('host,boot', [(None, None), ('copied-Mac-host', 'copied-boot'), ('same', 'old-boot')])
def test_copied_worker_pid_and_future_lease_are_not_trusted(tmp_path, host, boot):
    store = Store(tmp_path)
    pid = project(store)
    jid = store.enqueue(pid)
    store.claim(jid, owner_token='old-owner', worker_pid=os.getpid())
    with store.connect() as db:
        db.execute('UPDATE jobs SET worker_host=?,worker_boot=?,lease_until=? WHERE id=?',
                   (store.worker_host if host == 'same' else host, boot, time.time() + 86400, jid))
    with patch('biocentric.store.process_identity', side_effect=AssertionError('Do not probe copied PID')), \
         patch('os.kill', side_effect=AssertionError('Do not signal copied PID')):
        assert not store.owns_lease(jid, 'old-owner')
        assert not store.heartbeat(jid, 'old-owner')
        assert not store.finish(jid, 'completed', owner_token='old-owner')
        assert store.enqueue(pid) == jid
    job = store.jobs(pid)[0]
    assert job['status'] == 'queued' and job['worker_pid'] is None and job['lease_until'] is None


def test_reused_pid_with_wrong_process_start_does_not_keep_expired_job_live(tmp_path):
    store = Store(tmp_path)
    pid = project(store)
    jid = store.enqueue(pid)
    store.claim(jid, worker_pid=os.getpid())
    with store.connect() as db:
        db.execute('UPDATE jobs SET worker_start=?,lease_until=? WHERE id=?', ('wrong-incarnation', 0, jid))
    with patch('os.kill', side_effect=AssertionError('PID signal unnecessary')):
        store.enqueue(pid)
    assert store.jobs(pid)[0]['status'] == 'queued'


def test_schema_migration_first_preserves_legacy_database(tmp_path):
    store = Store(tmp_path)
    pid = project(store)
    artifact = store.artifact(pid, 'evidence.json', b'{"observation": 1}')
    with store.connect() as db:
        for column in ['worker_start', 'worker_boot', 'worker_host']:
            db.execute(f'ALTER TABLE jobs DROP COLUMN {column}')
    migrated = Store(tmp_path)
    backups = list((tmp_path / 'metadata-backups').glob('*/metadata.sqlite'))
    assert len(backups) == 1
    with sqlite3.connect(backups[0]) as backup:
        assert 'worker_host' not in {row[1] for row in backup.execute('PRAGMA table_info(jobs)')}
        assert backup.execute('SELECT sha256 FROM artifacts WHERE id=?', (artifact['id'],)).fetchone()[0] == artifact['sha256']
    assert migrated.read_artifact(pid, artifact) == b'{"observation": 1}'


def test_relocated_calculations_use_new_host_scoped_work_directory(tmp_path):
    old, new, pid, original, content = copied_history(tmp_path)
    copied_work = new / '.biocentric' / 'projects' / pid / 'work'
    copied_work.mkdir()
    checkpoint = copied_work / 'completed.json'
    checkpoint.write_text('{"fingerprint": "historical"}')
    configure(new / '.biocentric', [{'old_root': str(old), 'new_root': str(new)}], tmp_path / 'backup')
    store = Store(new / '.biocentric')
    work = store.work_directory(pid)
    assert work != copied_work and work.parent == copied_work
    assert work.name.startswith('host-') and not (work / 'completed.json').exists()
    assert checkpoint.read_text() == '{"fingerprint": "historical"}'
