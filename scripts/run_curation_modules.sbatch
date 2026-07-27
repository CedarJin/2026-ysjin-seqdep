#!/bin/bash
#SBATCH --job-name=curation_mod
#SBATCH --output=slurm_logs/curation_modules_%j.out
#SBATCH --error=slurm_logs/curation_modules_%j.err
# Run SAP curated-module abundance (id_only by default; see env/config.yaml).
#
# Prerequisites per subsample (from annotation_contigs_abundance.smk):
#   annotation/dram_contigs_T2T/{omic}/{sample}/{prefix}/distill_input_fixed/annotations.tsv
#   abundance/dram_contigs_T2T/{omic}/{sample}/{prefix}/gene_counts.tsv
#
# Outputs:
#   annotation/.../{prefix}/custom_modules/id_only/{module}_per_{feature,module}.tsv
#   annotation/dram_contigs_T2T/curation_module_summary/all_modules_per_module_long_id_only.tsv
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

cd "${PROJECT_DIR}"
export XDG_CACHE_HOME="${PROJECT_DIR}/.cache"
mkdir -p "${XDG_CACHE_HOME}" slurm_logs logs/curation_modules

echo "[$(date)] curation_modules.smk (id_only)"

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

THREADS="${SLURM_CPUS_PER_TASK:-16}"

python -m snakemake \
  -s "${SMK_FILE}" \
  --cores "${THREADS}" \
  --rerun-incomplete \
  --nolock \
  --printshellcmds

echo "[$(date)] DONE"
