#!/bin/bash

#SBATCH --job-name=rmhost_downsample_array
#SBATCH --output=slurm_rmhost_downsample_array_%j.out
#SBATCH --error=slurm_rmhost_downsample_array_%j.err
#SBATCH --partition=low
#SBATCH --account=publicgrp
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --time=24:00:00
#SBATCH --array=0-3

### Email notifications
#SBATCH --mail-user=ysjin@health.ucdavis.edu
#SBATCH --mail-type=END,FAIL

set -euo pipefail

PROJECT_DIR="/home/jys0914/2026-ysjin-seqdep"
SMK_FILE="scripts/remove_host_downsample.smk"

# Array seeds: remaining seeds after pilot seed11.
SEEDS=(22 33 44 55)
export SEED_FILTER="${SEEDS[${SLURM_ARRAY_TASK_ID}]}"

cd "${PROJECT_DIR}"

# Keep caches local to project (avoid home permission issues).
export XDG_CACHE_HOME="${PROJECT_DIR}/.cache"
mkdir -p "${XDG_CACHE_HOME}"

# Initialize conda in non-interactive shell.
if [[ -z "${CONDA_EXE:-}" ]]; then
  CONDA_BASE="${HOME}/anaconda3"
  if [[ -x "${CONDA_BASE}/bin/conda" ]]; then
    # shellcheck disable=SC1090,SC1091
    source "${CONDA_BASE}/etc/profile.d/conda.sh"
  else
    module load conda/base/latest
    # shellcheck disable=SC1091
    source "$(conda info --base)/etc/profile.d/conda.sh"
  fi
else
  # shellcheck disable=SC1090,SC1091
  source "$(dirname "${CONDA_EXE}")/../etc/profile.d/conda.sh"
fi

# assemble env contains bowtie2/samtools.
conda activate assemble

# Use project .venv Python/Snakemake if present.
if [[ -f ".venv/bin/activate" ]]; then
  # shellcheck disable=SC1091
  source ".venv/bin/activate"
fi

echo "[$(date)] Starting host removal: job ${SLURM_JOB_ID} array ${SLURM_ARRAY_TASK_ID} seed ${SEED_FILTER} on ${SLURMD_NODENAME}"
echo "CPUs: ${SLURM_CPUS_PER_TASK:-1}, MEM: ${SLURM_MEM_PER_NODE:-unknown}"

python -m snakemake -s "${SMK_FILE}" --unlock || true

python -m snakemake \
  -s "${SMK_FILE}" \
  --cores "${SLURM_CPUS_PER_TASK:-1}" \
  --rerun-incomplete \
  --printshellcmds

echo "[$(date)] Host removal workflow finished"
