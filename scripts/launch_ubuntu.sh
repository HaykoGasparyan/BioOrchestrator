#!/usr/bin/env bash
set -euo pipefail
task_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$task_root"
export BIOCENTRIC_VINA="$task_root/.tools/linux/vina-1.2.7"
if [[ -x "$task_root/.tools/linux/gnina" ]]; then export BIOCENTRIC_GNINA="$task_root/.tools/linux/gnina"; fi
if [[ -x "$task_root/.tools/linux/unidock" ]]; then export BIOCENTRIC_UNIDOCK="$task_root/.tools/linux/unidock"; fi
export BIOCENTRIC_MODEL_DIR="${BIOCENTRIC_MODEL_DIR:-$HOME/dti}"
export BIOCENTRIC_MODEL_PYTHON="${BIOCENTRIC_MODEL_PYTHON:-$HOME/anaconda3/bin/python}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-2}"
exec "$task_root/.venv-linux/bin/python" -m streamlit run "$task_root/app.py" --server.headless true --server.address 127.0.0.1 --server.port "${BIOCENTRIC_PORT:-8502}" "$@"
