#!/bin/bash
#SBATCH --job-name=curation_ec
#SBATCH --output=slurm_logs/curation_modules_id_ec_%j.out
#SBATCH --error=slurm_logs/curation_modules_id_ec_%j.err
# Run id_ec curated-module abundance for all subsamples + id_ec-only long table.
# Does NOT touch id_only outputs or all_modules_per_module_long_id_only.tsv.
#
# Outputs:
#   annotation/.../{prefix}/custom_modules/id_ec/{module}_per_{feature,module}.tsv
#   annotation/dram_contigs_T2T/curation_module_summary/all_modules_per_module_long_id_ec.tsv
#SBATCH --partition=high
#SBATCH --account=publicgrp
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=2-00:00:00
#SBATCH --mail-user=ysjin@health.ucdavis.edu
#SBATCH --mail-type=ALL

set -euo pipefail

PROJECT_DIR="/home/jys0914/2026-ysjin-seqdep"
SMK_FILE="scripts/curation_modules.smk"
ID_EC_LONG="annotation/dram_contigs_T2T/curation_module_summary/all_modules_per_module_long_id_ec.tsv"

cd "${PROJECT_DIR}"
export XDG_CACHE_HOME="${PROJECT_DIR}/.cache"
mkdir -p "${XDG_CACHE_HOME}" slurm_logs logs/curation_modules

echo "[$(date)] curation_modules.smk — id_ec (rule all -> ${ID_EC_LONG})"

if command -v module >/dev/null 2>&1; then
  module load conda/base/latest
fi
if [[ -n "${CONDA_EXE:-}" ]]; then
  source "$(dirname "${CONDA_EXE}")/../etc/profile.d/conda.sh"
elif command -v conda >/dev/null 2>&1; then
  CONDA_EXE="$(command -v conda)"
  source "$(dirname "${CONDA_EXE}")/../etc/profile.d/conda.sh"
else
  echo "ERROR: conda not found" >&2
  exit 1
fi
conda activate dram

THREADS="${SLURM_CPUS_PER_TASK:-8}"

# Target the concrete long-table output (rule all), not the wildcard rule name.
python -m snakemake \
  -s "${SMK_FILE}" \
  --cores "${THREADS}" \
  --config curation_id_variants='["id_ec"]' curation_summary="${ID_EC_LONG}" \
  --rerun-incomplete \
  --nolock \
  --printshellcmds \
  "${ID_EC_LONG}"

echo "[$(date)] DONE"
