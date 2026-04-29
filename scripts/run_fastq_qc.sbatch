#!/bin/bash

#SBATCH --job-name=fastq_qc
#SBATCH --output=slurm_fastq_qc_%j.out
#SBATCH --error=slurm_fastq_qc_%j.err
#SBATCH --partition=high
#SBATCH --account=publicgrp
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G      
   
#SBATCH --time=48:00:00

### Email notifications
#SBATCH --mail-user=ysjin@health.ucdavis.edu
#SBATCH --mail-type=END,FAIL


#### 16CPU and 256G mem for all runs



set -euo pipefail

# Use both: conda qc (fastqc, fastp) and project .venv (snakemake, multiqc).
# Order: activate qc first, then .venv so .venv's Python/snakemake take precedence;
# qc stays on PATH so fastqc and fastp are found when rules run.

PROJECT_DIR="/home/jys0914/2026-ysjin-seqdep"
SMK_FILE="scripts/fastq_qc.smk"

cd "${PROJECT_DIR}"

# Keep caches local to the project to avoid home permission issues.
export XDG_CACHE_HOME="${PROJECT_DIR}/.cache"
mkdir -p "${XDG_CACHE_HOME}"

# Initialize conda so "conda activate" works in this script (required in non-interactive shells).
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

# Activate conda qc env (fastqc, fastp).
conda activate qc

# Then activate project .venv (snakemake, multiqc). This prepends .venv to PATH
# so we use .venv's Python and snakemake; qc's bin remains on PATH for fastqc/fastp.
if [[ -f ".venv/bin/activate" ]]; then
  # shellcheck disable=SC1091
  source ".venv/bin/activate"
fi

echo "[$(date)] Starting Snakemake job ${SLURM_JOB_ID} on ${SLURMD_NODENAME}"
echo "CPUs: ${SLURM_CPUS_PER_TASK:-1}, MEM: ${SLURM_MEM_PER_NODE:-unknown}"

# Remove stale lock from a previous killed/crashed run (no-op if no lock).
python -m snakemake -s "${SMK_FILE}" --unlock || true

python -m snakemake \
  -s "${SMK_FILE}" \
  --cores "${SLURM_CPUS_PER_TASK:-1}" \
  --rerun-incomplete \
  --printshellcmds

echo "[$(date)] Snakemake job finished"
