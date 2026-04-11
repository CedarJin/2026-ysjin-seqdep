#!/bin/bash

#SBATCH --job-name=downsample_50M
#SBATCH --output=slurm_downsample_%j.out
#SBATCH --error=slurm_downsample_%j.err
#SBATCH --partition=low
#SBATCH --account=publicgrp
#SBATCH --cpus-per-task=4
#SBATCH --mem=256G
#SBATCH --time=72:00:00

### Email notifications
#SBATCH --mail-user=ysjin@health.ucdavis.edu
#SBATCH --mail-type=END,FAIL

set -euo pipefail


# Load required modules
module load seqtk
module load conda/base/latest

PROJECT_DIR="/home/jys0914/2026-ysjin-seqdep"
SMK_FILE="scripts/downsample_50M.smk"

cd "${PROJECT_DIR}"

# Keep caches local to the project to avoid home permission issues.
export XDG_CACHE_HOME="${PROJECT_DIR}/.cache"
mkdir -p "${XDG_CACHE_HOME}"

# Initialize and activate conda env with python + snakemake.
# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate qc

# Optional: activate project venv if present.
if [[ -f ".venv/bin/activate" ]]; then
  # shellcheck disable=SC1091
  source ".venv/bin/activate"
fi

echo "[$(date)] Starting Snakemake job ${SLURM_JOB_ID} on ${SLURMD_NODENAME}"
echo "CPUs: ${SLURM_CPUS_PER_TASK:-1}, MEM: ${SLURM_MEM_PER_NODE:-unknown}"

python -m snakemake -s "${SMK_FILE}" --unlock || true

python -m snakemake \
  -s "${SMK_FILE}" \
  --cores "${SLURM_CPUS_PER_TASK:-1}" \
  --rerun-incomplete \
  --printshellcmds

echo "[$(date)] Snakemake job finished"
