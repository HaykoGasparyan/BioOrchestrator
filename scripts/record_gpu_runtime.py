"""Record reviewed downloaded GPU executables, embedded weights and runtime files."""
import importlib.metadata
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from biocentric.relocation import file_sha256


def main():
    base = Path(__file__).resolve().parents[1] / '.tools' / 'linux'
    roots = [base / 'cuda12-libs', base / 'cuda12-nvrtc']
    files = [base / 'gnina', base / 'gnina-1.3.2', base / 'unidock-1.1.0',
             *[path for root in roots for path in sorted(root.rglob('*.so*')) if path.is_file()]]
    record = {'packages': {dist.metadata['Name']: dist.version for root in roots
                           for dist in importlib.metadata.distributions(path=[str(root)])},
              'files': [{'path': str(path), 'bytes': path.stat().st_size, 'sha256': file_sha256(path)} for path in files],
              'sources': ['https://github.com/gnina/gnina/releases/tag/v1.3.2',
                          'https://raw.githubusercontent.com/gnina/gnina/v1.3.2/docker/compat/Dockerfile',
                          'https://docs.nvidia.com/deeplearning/cudnn/backend/v9.1.1/reference/support-matrix.html',
                          'https://docs.nvidia.com/deeplearning/cudnn/backend/v9.1.1/installation/linux.html',
                          'https://github.com/dptech-corp/Uni-Dock/releases/tag/1.1.0'],
              'model_weights': 'GNINA default CNN weights are embedded in the hash-recorded native binary. Uni-Dock empirical Vina scoring has no learned checkpoint.'}
    (base / 'gpu-runtime-manifest.json').write_text(json.dumps(record, indent=2) + '\n')
    print(json.dumps({'packages': record['packages'], 'files': len(files)}, indent=2))


if __name__ == '__main__':
    main()
