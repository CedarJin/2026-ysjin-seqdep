#!/bin/bash
#SBATCH --job-name=rawdata_fqc_mqc
#SBATCH --output=slurm_rawdata_fastqc_multiqc_%j.out
#SBATCH --error=slurm_rawdata_fastqc_multiqc_%j.err
#SBATCH --partition=low
#SBATCH --account=publicgrp
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --mail-user=ysjin@health.ucdavis.edu
#SBATCH --mail-type=END,FAIL

set -euo pipefail

PROJECT_DIR="/home/jys0914/2026-ysjin-seqdep"
SMK_FILE="scripts/rawdata_fastqc_multiqc.smk"

cd "${PROJECT_DIR}"
export XDG_CACHE_HOME="${PROJECT_DIR}/.cache"
mkdir -p "${XDG_CACHE_HOME}"

# Load conda in batch shell.
if [[ -f "${HOME}/anaconda3/etc/profile.d/conda.sh" ]]; then
  source "${HOME}/anaconda3/etc/profile.d/conda.sh"
else
  module load conda/base/latest
  # module load alone is not enough for `conda activate` in non-interactive shells
  if [[ -n "${CONDA_EXE:-}" ]]; then
    # shellcheck disable=SC1090
    source "$(dirname "${CONDA_EXE}")/../etc/profile.d/conda.sh"
  else
    # shellcheck disable=SC1091
    source "$(conda info --base)/etc/profile.d/conda.sh"
  fi
fi

conda activate qc
[[ -f .venv/bin/activate ]] && source .venv/bin/activate

command -v snakemake >/dev/null || { echo "snakemake not found in PATH"; exit 1; }
command -v fastqc >/dev/null || { echo "fastqc not found in PATH"; exit 1; }
command -v multiqc >/dev/null || { echo "multiqc not found in PATH"; exit 1; }

python -m snakemake -s "${SMK_FILE}" --unlock || true
python -m snakemake -s "${SMK_FILE}" --cores "${SLURM_CPUS_PER_TASK:-1}" --rerun-incomplete --printshellcmds
