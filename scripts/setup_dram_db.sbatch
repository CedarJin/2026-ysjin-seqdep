#!/bin/bash

#SBATCH --job-name=dram_db_setup
#SBATCH --output=slurm_dram_db_%j.out
#SBATCH --error=slurm_dram_db_%j.err
#SBATCH --partition=low
#SBATCH --account=publicgrp
#SBATCH --cpus-per-task=16
#SBATCH --mem=512G
#SBATCH --time=12:00:00

### Email notifications
#SBATCH --mail-user=ysjin@health.ucdavis.edu
#SBATCH --mail-type=END,FAIL

# NOTE: DRAM database setup is memory-intensive:
#   - default (with uniref): ~512 GB RAM
#   - --skip_uniref (no KEGG Genes): ~64 GB RAM
#   - --skip_uniref (with KEGG Genes): ~128 GB RAM

# Remove --skip_uniref if you have access to a high-memory node (512G+).

set -euo pipefail

PROJECT_DIR="/home/jys0914/2026-ysjin-seqdep"
DRAM_DB_DIR="${PROJECT_DIR}/reference/DRAM_data"

cd "${PROJECT_DIR}"

export XDG_CACHE_HOME="${PROJECT_DIR}/.cache"
mkdir -p "${XDG_CACHE_HOME}"

# Initialize conda so "conda activate" works in non-interactive shells.
if [[ -z "${CONDA_EXE:-}" ]]; then
  CONDA_BASE="${HOME}/anaconda3"
  if [[ -x "${CONDA_BASE}/bin/conda" ]]; then
    # shellcheck disable=SC1090,SC1091
    source "${CONDA_BASE}/etc/profile.d/conda.sh"
  fi
else
  # shellcheck disable=SC1090,SC1091
  source "$(dirname "${CONDA_EXE}")/../etc/profile.d/conda.sh"
fi

# Activate DRAM env with conda.
conda activate DRAM

echo "[$(date)] Starting DRAM database setup on ${SLURMD_NODENAME}"
echo "CPUs: ${SLURM_CPUS_PER_TASK:-1}, MEM: ${SLURM_MEM_PER_NODE:-unknown}"

mkdir -p "${DRAM_DB_DIR}"

# Disable SSL verification to avoid certificate errors on cluster nodes.
export PYTHONHTTPSVERIFY=0

DRAM-setup.py prepare_databases \
  --output_dir "${DRAM_DB_DIR}" \
  --threads "${SLURM_CPUS_PER_TASK:-1}"
  
echo "[$(date)] DRAM database setup finished. DB saved to: ${DRAM_DB_DIR}"
