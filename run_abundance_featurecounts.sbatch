#!/bin/bash

#SBATCH --job-name=abundance_fc
#SBATCH --output=slurm_abundance_featurecounts_%j.out
#SBATCH --error=slurm_abundance_featurecounts_%j.err
#SBATCH --partition=low
#SBATCH --account=publicgrp
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00

### Email notifications
#SBATCH --mail-user=ysjin@health.ucdavis.edu
#SBATCH --mail-type=END,FAIL

set -euo pipefail

PROJECT_DIR="/home/jys0914/2026-ysjin-seqdep"
SMK_FILE="abundance_featurecounts.smk"

cd "${PROJECT_DIR}"

export XDG_CACHE_HOME="${PROJECT_DIR}/.cache"
mkdir -p "${XDG_CACHE_HOME}"

# Initialize conda
if [[ -z "${CONDA_EXE:-}" ]]; then
  CONDA_BASE="$(conda info --base)"
  if [[ -x "${CONDA_BASE}/bin/conda" ]]; then
    source "${CONDA_BASE}/etc/profile.d/conda.sh"
  fi
else
  source "$(dirname "${CONDA_EXE}")/../etc/profile.d/conda.sh"
fi

conda activate abundance

echo "[$(date)] Starting abundance featureCounts job ${SLURM_JOB_ID} on ${SLURMD_NODENAME}"
echo "CPUs: ${SLURM_CPUS_PER_TASK:-1}, MEM: ${SLURM_MEM_PER_NODE:-unknown}"

snakemake -s "${SMK_FILE}" --unlock || true

snakemake \
  -s "${SMK_FILE}" \
  --cores "${SLURM_CPUS_PER_TASK:-1}" \
  --rerun-incomplete \
  --printshellcmds

echo "[$(date)] Abundance featureCounts job finished"
