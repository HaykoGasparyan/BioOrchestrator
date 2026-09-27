"""Content-addressed, verified stage artifacts shared by evaluation workers."""
import fcntl
import json
import os
from pathlib import Path
import shutil
import tempfile
import time

from . import science


class StageCache(dict):
    def __init__(self, archive):
        super().__init__()
        self.archive = archive


class Archive:
    def __init__(self, root, runtime, cancelled, emit):
        self.root = Path(root)
        # Scorer/model changes must not invalidate upstream pocket/docking work.
        self.runtime = runtime
        self.cancelled, self.emit = cancelled, emit
        self.receipts = []

    def materialize(self, stage, identity, destination, produce):
        from .discovery_benchmark import digest
        names = {'pocket': {'pocket_identification.py', 'fpocket'},
                 'ligand': set(), 'docking': {'engines.py', 'engine_workflow.py'},
                 'gnina_scoring': {'engines.py'}}[stage]
        names |= {'science.py', 'discovery_benchmark.py', 'evaluation_archive.py'}
        tool = identity.get('tool')
        executable = None
        if tool == 'fpocket':
            from .pocket_identification import executable as pocket_executable
            executable = pocket_executable()
        elif tool in {'vina', 'gnina', 'unidock'}:
            from .engines import _executable
            executable = _executable(tool)
        files = {path: sha for path, sha in self.runtime.get('files', {}).items()
                 if Path(path).name in names or (executable and Path(path).resolve() == Path(executable).resolve())
                 or (tool == 'gnina' and '/.tools/linux/' in path)}
        key = digest({'schema': 1, 'stage': stage, 'identity': identity,
                      'runtime': files, 'versions': self.runtime.get('versions', {})})
        parent = self.root / stage
        parent.mkdir(parents=True, exist_ok=True)
        entry, destination = parent / key, Path(destination)
        with (parent / (key + '.lock')).open('a') as lock:
            while True:
                science._cancel(self.cancelled)
                try:
                    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
                    break
                except BlockingIOError:
                    time.sleep(.1)
            hit = self._valid(entry, key)
            if hit:
                shutil.copytree(entry / 'artifacts', destination, dirs_exist_ok=True)
            else:
                produce(destination)
                science._cancel(self.cancelled)
                # Publish only successful complete stages. Readers hold the same lock.
                with tempfile.TemporaryDirectory(prefix=key + '-', dir=parent) as temporary:
                    temporary = Path(temporary)
                    shutil.copytree(destination, temporary / 'artifacts')
                    hashes = {str(p.relative_to(temporary / 'artifacts')): science.sha256(p)
                              for p in (temporary / 'artifacts').rglob('*') if p.is_file()}
                    science._write_json(temporary / 'manifest.json', {'key': key, 'identity': identity,
                        'created_at': science._now(), 'files': hashes})
                    if entry.exists():
                        shutil.rmtree(entry)  # Only this invalid generated cache entry.
                    os.replace(temporary, entry)
            receipt = {'stage': stage, 'key': key, 'hit': hit}
            self.receipts.append(receipt)
            self.emit({'tool': stage, 'stage': 'archive',
                       'message': ('Reused archived ' if hit else 'Archived ') + stage + ' · ' + key[:12]})
            return receipt

    @staticmethod
    def _valid(entry, key):
        try:
            manifest = json.loads((entry / 'manifest.json').read_text())
            files = manifest['files']
            if manifest['key'] != key or not files:
                return False
            for name, sha in files.items():
                path = entry / 'artifacts' / name
                if Path(name).is_absolute() or '..' in Path(name).parts or path.is_symlink() or science.sha256(path) != sha:
                    return False
            return True
        except (OSError, ValueError, KeyError, TypeError):
            return False


def identify_pocket(method, receptor, reference, destination, cancelled, archive=None):
    from . import pocket_identification as pockets
    if archive is None:
        return pockets.identify(method, receptor, reference, destination, cancelled)
    identity = {'tool': method['tool'], 'method': method, 'receptor': science.sha256(receptor)}
    if method['tool'] == 'reference_box':
        identity['reference'] = science.sha256(reference)
    archive.materialize('pocket', identity, destination,
                        lambda out: pockets.identify(method, receptor, reference, out, cancelled))
    return json.loads((Path(destination) / 'pocket.json').read_text())['box']
