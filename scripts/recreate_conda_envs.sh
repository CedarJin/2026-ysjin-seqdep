#!/usr/bin/env bash
# Recreate dram, qc, and assemble under ~/.conda/envs (symlink may resolve to Quobyte).
# Usage: bash scripts/recreate_conda_envs.sh
# Requires: module load conda/base/latest (or conda in PATH)

set -euo pipefail

if command -v module >/dev/null 2>&1; then
  module load conda/base/latest
fi

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BASE="$(readlink -f "${HOME}/.conda")/envs"
mkdir -p "$BASE"

remove_env() {
  local p="$1"
  if [ -d "$p" ]; then
    conda env remove -p "$p" -y
  fi
}

for name in dram qc assemble; do
  remove_env "$BASE/$name"
done

conda env create -f "$REPO_ROOT/env/env_DRAM.yaml" -p "$BASE/dram"
conda env create -f "$REPO_ROOT/env/env_qc.yaml" -p "$BASE/qc"
conda env create -f "$REPO_ROOT/env/env_assembly_bin.yaml" -p "$BASE/assemble"

echo "Done. Prefixes:"
for name in dram qc assemble; do
  echo "  $BASE/$name"
done
