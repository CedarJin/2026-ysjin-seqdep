#!/bin/bash

#SBATCH --job-name=dram_contigs_annotate
#SBATCH --output=slurm_annotation_contigs_%j.out
#SBATCH --error=slurm_annotation_contigs_%j.err
#SBATCH --partition=low
#SBATCH --account=publicgrp
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --time=24:00:00

### Email notifications
#SBATCH --mail-user=ysjin@health.ucdavis.edu
#SBATCH --mail-type=END,FAIL

# NOTE:
#   - This workflow runs DRAM on assembled contigs from MEGAHIT:
#       assembly/megahit/{run}_{depth}_seed{seed}/final.contigs.fa
#   - To disable CAMPER annotation, set in config.yaml:
#       dram_use_camper: false

set -euo pipefail

PROJECT_DIR="/home/jys0914/2026-ysjin-seqdep"
SMK_FILE="annotation_contigs.smk"

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

# Activate DRAM env
conda activate DRAM

echo "[$(date)] Starting DRAM contig annotation job ${SLURM_JOB_ID} on ${SLURMD_NODENAME}"
echo "CPUs: ${SLURM_CPUS_PER_TASK:-1}, MEM: ${SLURM_MEM_PER_NODE:-unknown}"

python -m snakemake -s "${SMK_FILE}" --unlock || true

python -m snakemake \
  -s "${SMK_FILE}" \
  --cores "${SLURM_CPUS_PER_TASK:-1}" \
  --rerun-incomplete \
  --printshellcmds

echo "[$(date)] DRAM contig annotation job finished"
