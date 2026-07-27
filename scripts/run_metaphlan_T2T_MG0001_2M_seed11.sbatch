#!/bin/bash

#SBATCH --job-name=metaphlan_T2T_test
#SBATCH --output=slurm_metaphlan_T2T_MG0001_2M_seed11_%j.out
#SBATCH --error=slurm_metaphlan_T2T_MG0001_2M_seed11_%j.err
#SBATCH --partition=high
#SBATCH --account=publicgrp
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=0-08:00:00

### Email notifications
#SBATCH --mail-user=ysjin@health.ucdavis.edu
#SBATCH --mail-type=ALL

# Test MetaPhlAn profiling on the first metaG sample/prefix:
#   cleandata_T2T/metaG/MG0001/MG0001_2M_seed11_{R1,R2}.fastq
# Submit: sbatch scripts/run_metaphlan_T2T_MG0001_2M_seed11.sbatch

set -euo pipefail

PROJECT_DIR="/home/jys0914/2026-ysjin-seqdep"
SMK_FILE="scripts/metaphlan_T2T.smk"

TEST_OMIC="metaG"
TEST_SAMPLE="MG0001"
TEST_PREFIX="MG0001_2M_seed11"

R1="cleandata_T2T/${TEST_OMIC}/${TEST_SAMPLE}/${TEST_PREFIX}_R1.fastq"
R2="cleandata_T2T/${TEST_OMIC}/${TEST_SAMPLE}/${TEST_PREFIX}_R2.fastq"
TARGET="taxonomy/metaphlan_T2T/${TEST_OMIC}/${TEST_SAMPLE}/${TEST_PREFIX}/profile.tsv"

cd "${PROJECT_DIR}"

export XDG_CACHE_HOME="${PROJECT_DIR}/.cache"
mkdir -p "${XDG_CACHE_HOME}" logs/metaphlan_T2T

if [[ ! -s "${R1}" || ! -s "${R2}" ]]; then
  echo "ERROR: expected paired FASTQ not found or empty:" >&2
  echo "  ${R1}" >&2
  echo "  ${R2}" >&2
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

# metaphlan env must include snakemake-minimal (env/env_metaphlan.yaml).
conda activate metaphlan

echo "[$(date)] Starting MetaPhlAn T2T test job ${SLURM_JOB_ID} on ${SLURMD_NODENAME}"
echo "Input R1: ${R1}"
echo "Input R2: ${R2}"
echo "Target:   ${TARGET}"
echo "CPUs: ${SLURM_CPUS_PER_TASK:-1}, MEM: ${SLURM_MEM_PER_NODE:-unknown}"

python -m snakemake \
  -s "${SMK_FILE}" \
  --cores "${SLURM_CPUS_PER_TASK:-1}" \
  --rerun-incomplete \
  --nolock \
  --printshellcmds \
  "${TARGET}"

echo "[$(date)] MetaPhlAn T2T test job finished"
