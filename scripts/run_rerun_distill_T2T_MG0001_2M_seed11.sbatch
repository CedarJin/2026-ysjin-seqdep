#!/bin/bash

#SBATCH --job-name=dram_distill_rerun_s11
#SBATCH --output=slurm_dram_distill_rerun_MG0001_2M_seed11_%j.out
#SBATCH --error=slurm_dram_distill_rerun_MG0001_2M_seed11_%j.err
#SBATCH --partition=high
#SBATCH --account=publicgrp
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=0-08:00:00

### Email notifications
#SBATCH --mail-user=ysjin@health.ucdavis.edu
#SBATCH --mail-type=ALL

set -euo pipefail

PROJECT_DIR="/home/jys0914/2026-ysjin-seqdep"
SMK_FILE="scripts/annotation_contigs.smk"

TEST_OMIC="metaG"
TEST_SAMPLE="MG0001"
TEST_PREFIX="MG0001_2M_seed11"
SAMPLE_ROOT="annotation/dram_contigs_T2T/${TEST_OMIC}/${TEST_SAMPLE}/${TEST_PREFIX}"
TARGET="${SAMPLE_ROOT}/distillate/genome_stats.tsv"

cd "${PROJECT_DIR}"

export XDG_CACHE_HOME="${PROJECT_DIR}/.cache"
mkdir -p "${XDG_CACHE_HOME}" logs/dram_contigs_T2T

if [[ ! -s "${SAMPLE_ROOT}/annotations.tsv" ]]; then
  echo "ERROR: missing annotation input: ${SAMPLE_ROOT}/annotations.tsv" >&2
  exit 1
fi

if [[ ! -s "${SAMPLE_ROOT}/trnas.tsv" || ! -s "${SAMPLE_ROOT}/rrnas.tsv" ]]; then
  echo "ERROR: missing tRNA/rRNA input under ${SAMPLE_ROOT}" >&2
  exit 1
fi

# Initialize conda in the non-interactive Slurm shell.
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
  echo "ERROR: conda was not found; module load conda/base/latest did not provide conda" >&2
  exit 1
fi

conda activate dram

echo "[$(date)] Starting DRAM distill rerun ${SLURM_JOB_ID} on ${SLURMD_NODENAME}"
echo "Sample: ${TEST_OMIC}/${TEST_SAMPLE}/${TEST_PREFIX}"
echo "Target: ${TARGET}"
echo "CPUs: ${SLURM_CPUS_PER_TASK:-1}, MEM: ${SLURM_MEM_PER_NODE:-unknown}"

python -m snakemake \
  -s "${SMK_FILE}" \
  --cores "${SLURM_CPUS_PER_TASK:-1}" \
  --config dram_threads="${SLURM_CPUS_PER_TASK:-8}" \
  --forcerun fix_dram_contig_inputs_for_distill prepare_contig_distill_config_with_camper dram_contig_distill \
  --rerun-incomplete \
  --nolock \
  --printshellcmds \
  "${TARGET}"

echo "[$(date)] DRAM distill rerun finished"
