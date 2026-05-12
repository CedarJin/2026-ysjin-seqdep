#!/bin/bash

#SBATCH --job-name=assembly_T2T_seed11
#SBATCH --output=slurm_assembly_T2T_seed11_%j.out
#SBATCH --error=slurm_assembly_T2T_seed11_%j.err
#SBATCH --partition=low
#SBATCH --account=publicgrp
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --time=7-00:00:00

### Email notifications
#SBATCH --mail-user=ysjin@health.ucdavis.edu
#SBATCH --mail-type=ALL

set -euo pipefail

PROJECT_DIR="/home/jys0914/2026-ysjin-seqdep"
SMK_FILE="scripts/assembly_T2T.smk"
export SEED_FILTER="11"

cd "${PROJECT_DIR}"

export XDG_CACHE_HOME="${PROJECT_DIR}/.cache"
mkdir -p "${XDG_CACHE_HOME}"
mkdir -p assembly_T2T/megahit assembly_T2T/alignments logs/megahit_T2T logs/bowtie2_T2T

# Initialize conda
module load conda/base/latest
source "$(conda info --base)/etc/profile.d/conda.sh"

# Activate assemble env (must be created first: conda env create -f env/env_assemble.yaml -n assemble)
conda activate assemble

# Pre-create seed11 output directories to reduce metadata races on the shared filesystem.
python - <<'PY'
import glob
import os

for r1 in glob.glob("cleandata_T2T/*/*/*_seed11_R1.fastq"):
    omic = r1.split(os.sep)[1]
    sample = r1.split(os.sep)[2]
    prefix = os.path.basename(r1).replace("_R1.fastq", "")
    os.makedirs(f"assembly_T2T/megahit/{omic}/{sample}/{prefix}", exist_ok=True)
    os.makedirs(f"assembly_T2T/alignments/{omic}/{sample}", exist_ok=True)
PY

# Optional: use project .venv for snakemake if you prefer that Python
if [[ -f ".venv/bin/activate" ]]; then
  source ".venv/bin/activate"
fi

echo "[$(date)] Starting Snakemake assembly seed11 job ${SLURM_JOB_ID} on ${SLURMD_NODENAME}"
echo "CPUs: ${SLURM_CPUS_PER_TASK:-1}, MEM: ${SLURM_MEM_PER_NODE:-unknown}"

python -m snakemake \
  -s "${SMK_FILE}" \
  --cores "${SLURM_CPUS_PER_TASK:-1}" \
  --rerun-incomplete \
  --nolock \
  --printshellcmds

echo "[$(date)] Snakemake assembly seed11 job finished"
