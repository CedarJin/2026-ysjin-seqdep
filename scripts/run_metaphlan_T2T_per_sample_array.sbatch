#!/bin/bash
#SBATCH --job-name=metaphlan_sample
#SBATCH --output=slurm_logs/metaphlan_sample_%A_%a.out
#SBATCH --error=slurm_logs/metaphlan_sample_%A_%a.err
# One array task per (omic, sample): MetaPhlAn profiling for every paired FASTQ
# prefix under cleandata_T2T/{omic}/{sample}/ (all depths × all seeds on disk).
#
# Array indices 0–15: MG0001–MG0008 (metaG) + MT0001–MT0008 (metaT).
# Submit: sbatch scripts/run_metaphlan_T2T_per_sample_array.sbatch
#SBATCH --array=9
#SBATCH --partition=high
#SBATCH --account=publicgrp
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=3-00:00:00
#SBATCH --mail-user=ysjin@health.ucdavis.edu
#SBATCH --mail-type=ALL

set -euo pipefail

PROJECT_DIR="/home/jys0914/2026-ysjin-seqdep"
SMK_FILE="scripts/metaphlan_T2T.smk"

# (omic, sample) list — fixed order for reproducible array indices.
SAMPLES=(
  "metaG MG0001" "metaG MG0002" "metaG MG0003" "metaG MG0004"
  "metaG MG0005" "metaG MG0006" "metaG MG0007" "metaG MG0008"
  "metaT MT0001" "metaT MT0002" "metaT MT0003" "metaT MT0004"
  "metaT MT0005" "metaT MT0006" "metaT MT0007" "metaT MT0008"
)

ENTRY="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"
OMIC="${ENTRY%% *}"
SAMPLE="${ENTRY##* }"

cd "${PROJECT_DIR}"
export XDG_CACHE_HOME="${PROJECT_DIR}/.cache"
mkdir -p "${XDG_CACHE_HOME}" slurm_logs logs/metaphlan_T2T

echo "[$(date)] Task ${SLURM_ARRAY_TASK_ID}: OMIC=${OMIC} SAMPLE=${SAMPLE}"

CLEANDATA_SAMPLE_DIR="cleandata_T2T/${OMIC}/${SAMPLE}"
if [[ ! -d "${CLEANDATA_SAMPLE_DIR}" ]]; then
  echo "ERROR: cleandata dir missing: ${CLEANDATA_SAMPLE_DIR}" >&2
  exit 1
fi

TARGETS=()
N_PREFIX=0
while IFS= read -r R1; do
  [[ -s "${R1}" ]] || continue
  R2="${R1/_R1.fastq/_R2.fastq}"
  if [[ ! -s "${R2}" ]]; then
    echo "  WARN: R2 missing for ${R1}" >&2
    continue
  fi
  PREFIX="$(basename "${R1}" _R1.fastq)"
  TARGETS+=("taxonomy/metaphlan_T2T/${OMIC}/${SAMPLE}/${PREFIX}/profile.tsv")
  N_PREFIX=$((N_PREFIX + 1))
done < <(find "${CLEANDATA_SAMPLE_DIR}" -maxdepth 1 -name '*_R1.fastq' 2>/dev/null | sort)

if [[ ${N_PREFIX} -eq 0 ]]; then
  echo "ERROR: no paired FASTQ under ${CLEANDATA_SAMPLE_DIR}" >&2
  exit 1
fi
echo "Found ${N_PREFIX} FASTQ prefixes (depth×seed combos), ${#TARGETS[@]} profile targets"

if command -v module >/dev/null 2>&1; then
  module load conda/base/latest
fi
if [[ -n "${CONDA_EXE:-}" ]]; then
  # shellcheck disable=SC1090,SC1091
  source "$(dirname "${CONDA_EXE}")/../etc/profile.d/conda.sh"
elif command -v conda >/dev/null 2>&1; then
  CONDA_EXE="$(command -v conda)"
  # shellcheck disable=SC1090,SC1091
  source "$(dirname "${CONDA_EXE}")/../etc/profile.d/conda.sh"
else
  echo "ERROR: conda not found" >&2
  exit 1
fi
conda activate metaphlan

THREADS="${SLURM_CPUS_PER_TASK:-8}"
export METAPHLAN_THREADS="${THREADS}"

python -m snakemake \
  -s "${SMK_FILE}" \
  --cores "${THREADS}" \
  --rerun-incomplete \
  --nolock \
  --printshellcmds \
  "${TARGETS[@]}"

echo "[$(date)] Task ${SLURM_ARRAY_TASK_ID} DONE: ${OMIC}/${SAMPLE} (${N_PREFIX} prefixes)"
