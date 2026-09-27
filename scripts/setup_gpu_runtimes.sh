#!/usr/bin/env bash
# Project-local runtime only. Does not replace the driver, system CUDA, or ML env.
set -euo pipefail
cd "$(dirname "$0")/.."
if [[ "$(uname -s)" != Linux || "$(uname -m)" != x86_64 ]]; then
  echo 'These pinned binaries require Linux x86_64.' >&2
  exit 1
fi
python=.venv-linux/bin/python
if [[ ! -x "$python" ]]; then
  echo 'Run scripts/setup_ubuntu.sh first.' >&2
  exit 1
fi
mkdir -p .tools/linux
if [[ -f .tools/linux/gpu-runtime-manifest.json ]]; then
  "$python" - <<'PY'
import json
from pathlib import Path
from biocentric.relocation import file_sha256
record=json.loads(Path('.tools/linux/gpu-runtime-manifest.json').read_text())
for item in record['files']:
    if file_sha256(Path(item['path'])) != item['sha256']:
        raise SystemExit('Runtime file changed; review before reinstalling: '+item['path'])
print('Existing GPU runtime hashes verified. No files changed.')
PY
  exit 0
fi
fetch() {
  local url="$1" target="$2" sha="$3"
  if [[ ! -f "$target" ]]; then
    curl --fail --location --retry 3 "$url" --output "$target.download"
    echo "$sha  $target.download" | sha256sum --check
    mv "$target.download" "$target"
  fi
  echo "$sha  $target" | sha256sum --check
  chmod 755 "$target"
}
fetch https://github.com/gnina/gnina/releases/download/v1.3.2/gnina.1.3.2 \
  .tools/linux/gnina-1.3.2 5d33538324b40050a03aa262d51832837e0ea6cc100945abbd2d7b732589690e
fetch https://github.com/dptech-corp/Uni-Dock/releases/download/1.1.0/unidock-1.1.0-cuda120-linux-x86_64 \
  .tools/linux/unidock-1.1.0 8ff3c2394111fc97ef5eed4e0cceac462c1a8f01823c3d7dd937178a98be7838
"$python" -m pip install --target .tools/linux/cuda12-libs --no-cache-dir \
  --report .tools/linux/cuda12-install-report.json -r requirements-gpu-runtime-ubuntu.txt
"$python" -m pip install --target .tools/linux/cuda12-nvrtc --no-cache-dir \
  --report .tools/linux/nvrtc-install-report.json nvidia-cuda-nvrtc-cu12==12.1.105
cat > .tools/linux/gnina <<'SH'
#!/usr/bin/env bash
set -eu
bindir="$(cd -- "$(dirname -- "$0")" && pwd)"
cuda_paths=""
for directory in "$bindir"/cuda12-libs/nvidia/*/lib "$bindir"/cuda12-nvrtc/nvidia/*/lib; do
  cuda_paths="${cuda_paths:+$cuda_paths:}$directory"
done
export LD_LIBRARY_PATH="$cuda_paths${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
exec /lib64/ld-linux-x86-64.so.2 --inhibit-rpath '' --library-path "$LD_LIBRARY_PATH" "$bindir/gnina-1.3.2" "$@"
SH
chmod 755 .tools/linux/gnina
if [[ ! -e .tools/linux/unidock ]]; then
  ln -s unidock-1.1.0 .tools/linux/unidock
fi
"$python" scripts/record_gpu_runtime.py
.tools/linux/gnina --version
.tools/linux/unidock --version
echo 'Runtime installed. Version output is not a GPU execution test; run scripts/verify_gpu_engines.py.'
