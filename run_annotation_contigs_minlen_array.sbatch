#!/bin/bash

#SBATCH --job-name=dram_contigs_minlen
#SBATCH --output=slurm_annotation_contigs_minlen_%A_%a.out
#SBATCH --error=slurm_annotation_contigs_minlen_%A_%a.err
#SBATCH --partition=low
#SBATCH --account=publicgrp
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=72:00:00
#SBATCH --array=0

set -euo pipefail

PROJECT_DIR="/home/jys0914/2026-ysjin-seqdep"
SMK_FILE="annotation_contigs_minlen_scan.smk"
MINLENS=(200 500 1000 10000)
MINLEN="${MINLENS[$SLURM_ARRAY_TASK_ID]}"
RUN="SRR10692699"
DEPTH="50M"
SEED="11"
CONTIGS="assembly/megahit/${RUN}_${DEPTH}_seed${SEED}/final.contigs.fa"

cd "${PROJECT_DIR}"
export XDG_CACHE_HOME="${PROJECT_DIR}/.cache"
mkdir -p "${XDG_CACHE_HOME}"

# Initialize conda
if [[ -z "${CONDA_EXE:-}" ]]; then
  CONDA_BASE="${HOME}/anaconda3"
  if [[ -x "${CONDA_BASE}/bin/conda" ]]; then
    source "${CONDA_BASE}/etc/profile.d/conda.sh"
  fi
else
  source "$(dirname "${CONDA_EXE}")/../etc/profile.d/conda.sh"
fi

conda activate DRAM

echo "[$(date)] Starting DRAM minlen scan job ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID} (minlen=${MINLEN}) on ${SLURMD_NODENAME}"
echo "Input contigs: ${CONTIGS}"
echo "CPUs: ${SLURM_CPUS_PER_TASK:-1}, MEM: ${SLURM_MEM_PER_NODE:-unknown}"

if [[ ! -s "${CONTIGS}" ]]; then
  echo "ERROR: Missing contigs file: ${CONTIGS}" >&2
  exit 1
fi

python -m snakemake -s "${SMK_FILE}" --unlock || true

python -m snakemake \
  -s "${SMK_FILE}" \
  --cores "${SLURM_CPUS_PER_TASK:-1}" \
  --rerun-incomplete \
  --printshellcmds \
  --config \
    min_contig_len="${MINLEN}" \
    runtime_scan_run="${RUN}" \
    runtime_scan_depth="${DEPTH}" \
    runtime_scan_seed="${SEED}" \
    runtime_scan_threads="${SLURM_CPUS_PER_TASK:-1}" \
    runtime_scan_contigs="${CONTIGS}"

echo "[$(date)] DRAM minlen scan job finished for minlen=${MINLEN}"
