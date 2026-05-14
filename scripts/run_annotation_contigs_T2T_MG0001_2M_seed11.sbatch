#!/bin/bash

#SBATCH --job-name=dram_T2T_MG0001_2M_s11
#SBATCH --output=slurm_dram_T2T_MG0001_2M_seed11_%j.out
#SBATCH --error=slurm_dram_T2T_MG0001_2M_seed11_%j.err
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

CONTIGS="assembly_T2T/megahit/${TEST_OMIC}/${TEST_SAMPLE}/${TEST_PREFIX}/final.contigs.fa"
TARGET="annotation/dram_contigs_T2T/${TEST_OMIC}/${TEST_SAMPLE}/${TEST_PREFIX}/distillate/genome_stats.tsv"

cd "${PROJECT_DIR}"

export XDG_CACHE_HOME="${PROJECT_DIR}/.cache"
mkdir -p "${XDG_CACHE_HOME}" logs/dram_contigs_T2T

if [[ ! -s "${CONTIGS}" ]]; then
  echo "ERROR: expected contigs not found or empty: ${CONTIGS}" >&2
  echo "Check TEST_OMIC/TEST_SAMPLE/TEST_PREFIX in this sbatch script." >&2
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

# Activate DRAM env (must be created first: conda env create -f env/env_DRAM.yaml -n dram)
conda activate dram

echo "[$(date)] Starting T2T contig DRAM test job ${SLURM_JOB_ID} on ${SLURMD_NODENAME}"
echo "Test contigs: ${CONTIGS}"
echo "Target output: ${TARGET}"
echo "CPUs: ${SLURM_CPUS_PER_TASK:-1}, MEM: ${SLURM_MEM_PER_NODE:-unknown}"

python -m snakemake \
  -s "${SMK_FILE}" \
  --cores "${SLURM_CPUS_PER_TASK:-1}" \
  --config dram_threads="${SLURM_CPUS_PER_TASK:-8}" \
  --rerun-incomplete \
  --nolock \
  --printshellcmds \
  "${TARGET}"

echo "[$(date)] T2T contig DRAM test job finished"
