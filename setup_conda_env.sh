#!/usr/bin/env bash
set -euo pipefail

ENV_NAME="${1:-rbh_plus_diamond}"

if command -v mamba >/dev/null 2>&1; then
  CONDA_EXE="mamba"
elif command -v conda >/dev/null 2>&1; then
  CONDA_EXE="conda"
else
  echo "Error: conda or mamba was not found in PATH." >&2
  exit 1
fi

"${CONDA_EXE}" create -y -n "${ENV_NAME}" \
  -c conda-forge \
  -c bioconda \
  python=3.13 \
  diamond=2.1.13 \
  blast \
  numpy=2.3 \
  pandas=2.3

echo
echo "Created conda environment: ${ENV_NAME}"
echo "Activate it with:"
echo "  conda activate ${ENV_NAME}"
