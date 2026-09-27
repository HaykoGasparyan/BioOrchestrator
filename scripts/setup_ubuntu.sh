#!/usr/bin/env bash
# Linux application environment. BioCentric and CUDA runtimes are separate.
set -euo pipefail
task_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$task_root"
[[ "$(uname -s)" == Linux && "$(uname -m)" == x86_64 ]] || { echo 'Requires Linux x86_64'; exit 2; }
python3 -c 'import sys; assert (3,11) <= sys.version_info[:2] <= (3,13), "Python 3.11–3.13 required"'
python3 -m venv --without-pip .venv-linux
python3 -m pip --python .venv-linux install -r requirements-lock-ubuntu-x86_64.txt
.venv-linux/bin/python -m pip install --no-deps -e .
mkdir -p .tools/linux
if [[ ! -e .tools/linux/vina-1.2.7 ]]; then
  curl -L --fail --connect-timeout 10 --max-time 180 https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.7/vina_1.2.7_linux_x86_64 -o .tools/linux/vina-1.2.7
fi
echo 'f31f774f723bba7bbe6e9d1c47577020eea9a8da16424284c043d22593570644  .tools/linux/vina-1.2.7' | sha256sum -c -
chmod +x .tools/linux/vina-1.2.7
.tools/linux/vina-1.2.7 --version
.venv-linux/bin/python -m pip check
echo 'Application and Vina ready. See docs/ubuntu-setup.md for the separate model and GPU runtime records.'
