#!/bin/bash

#SBATCH --job-name=assembly
#SBATCH --output=slurm_assembly_%A_%a.out
#SBATCH --error=slurm_assembly_%A_%a.err
#SBATCH --partition=low
#SBATCH --account=publicgrp
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --time=5-00:00:00
#SBATCH --array=0-3

### Email notifications
#SBATCH --mail-user=ysjin@health.ucdavis.edu
#SBATCH --mail-type=END,FAIL

set -euo pipefail

PROJECT_DIR="/home/jys0914/2026-ysjin-seqdep"
SMK_FILE="scripts/assembly.smk"
SEEDS=(22 33 44 55)
export SEED_FILTER="${SEEDS[${SLURM_ARRAY_TASK_ID}]}"

cd "${PROJECT_DIR}"

export XDG_CACHE_HOME="${PROJECT_DIR}/.cache"
mkdir -p "${XDG_CACHE_HOME}"

# Initialize conda
if [[ -z "${CONDA_EXE:-}" ]]; then
  CONDA_BASE="${HOME}/anaconda3"
  if [[ -x "${CONDA_BASE}/bin/conda" ]]; then
    source "${CONDA_BASE}/etc/profile.d/conda.sh"
  fi
else
  source "$(dirname "${CONDA_EXE}")/../etc/profile.d/conda.sh"
fi

# Activate assemble env (must be created first: conda env create -f env/env_assemble.yaml -n assemble)
conda activate assemble

# Optional: use project .venv for snakemake if you prefer that Python
if [[ -f ".venv/bin/activate" ]]; then
  source ".venv/bin/activate"
fi

echo "[$(date)] Starting Snakemake assembly job ${SLURM_JOB_ID} array ${SLURM_ARRAY_TASK_ID} (seed=${SEED_FILTER}) on ${SLURMD_NODENAME}"
echo "CPUs: ${SLURM_CPUS_PER_TASK:-1}, MEM: ${SLURM_MEM_PER_NODE:-unknown}"

python -m snakemake -s "${SMK_FILE}" --unlock || true

python -m snakemake \
  -s "${SMK_FILE}" \
  --cores "${SLURM_CPUS_PER_TASK:-1}" \
  --rerun-incomplete \
  --printshellcmds

echo "[$(date)] Snakemake assembly job finished"
