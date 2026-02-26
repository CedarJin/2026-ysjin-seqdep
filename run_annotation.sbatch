#!/bin/bash

#SBATCH --job-name=dram_annotate
#SBATCH --output=slurm_annotation_%j.out
#SBATCH --error=slurm_annotation_%j.err
#SBATCH --partition=low
#SBATCH --account=publicgrp
#SBATCH --cpus-per-task=24
#SBATCH --mem=128G
#SBATCH --time=24:00:00

### Email notifications
#SBATCH --mail-user=ysjin@health.ucdavis.edu
#SBATCH --mail-type=END,FAIL

# NOTE: Memory requirements depend on databases used:
#   - with UniRef90: ~220 GB RAM
#   - KOfam only, no UniRef90: <50 GB RAM
# Adjust --mem above accordingly.

set -euo pipefail

PROJECT_DIR="/home/jys0914/2026-ysjin-seqdep"
SMK_FILE="annotation.smk"

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

# Activate DRAM env (must be created first: conda env create -f DRAM_env.yaml -n DRAM)
conda activate DRAM

# Optional: use project .venv for snakemake
if [[ -f ".venv/bin/activate" ]]; then
  source ".venv/bin/activate"
fi

echo "[$(date)] Starting DRAM annotation job ${SLURM_JOB_ID} on ${SLURMD_NODENAME}"
echo "CPUs: ${SLURM_CPUS_PER_TASK:-1}, MEM: ${SLURM_MEM_PER_NODE:-unknown}"

python -m snakemake -s "${SMK_FILE}" --unlock || true

python -m snakemake \
  -s "${SMK_FILE}" \
  --cores "${SLURM_CPUS_PER_TASK:-1}" \
  --rerun-incomplete \
  --printshellcmds

echo "[$(date)] DRAM annotation job finished"
