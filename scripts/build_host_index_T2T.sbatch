#!/bin/bash

#SBATCH --job-name=build_T2T_bt2
#SBATCH --output=slurm_build_T2T_bt2_%j.out
#SBATCH --error=slurm_build_T2T_bt2_%j.err
#SBATCH --partition=high
#SBATCH --account=publicgrp
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=8:00:00

### Email notifications
#SBATCH --mail-user=ysjin@health.ucdavis.edu
#SBATCH --mail-type=ALL

set -euo pipefail

PROJECT_DIR="/home/jys0914/2026-ysjin-seqdep"
SMK_FILE="scripts/remove_host_downsample_T2T.smk"

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

echo "[$(date)] Starting T2T bowtie2 index build: job ${SLURM_JOB_ID} on ${SLURMD_NODENAME}"
echo "CPUs: ${SLURM_CPUS_PER_TASK:-1}, MEM: ${SLURM_MEM_PER_NODE:-unknown}"

python -m snakemake \
  -s "${SMK_FILE}" \
  build_host_index \
  --cores "${SLURM_CPUS_PER_TASK:-1}" \
  --rerun-incomplete \
  --printshellcmds

echo "[$(date)] T2T bowtie2 index build finished"
