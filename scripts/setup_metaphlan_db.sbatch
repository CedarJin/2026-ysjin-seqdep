#!/bin/bash
#SBATCH --job-name=metaphlan_db
#SBATCH --output=slurm_metaphlan_db_%j.out
#SBATCH --error=slurm_metaphlan_db_%j.err
#SBATCH --partition=low
#SBATCH --account=publicgrp
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=72:00:00

#SBATCH --mail-user=ysjin@health.ucdavis.edu
#SBATCH --mail-type=END,FAIL

# Remove any partial MetaPhlAn DB, then run metaphlan --install.
# Submit: sbatch scripts/setup_metaphlan_db.sbatch

set -euo pipefail

PROJECT_DIR="/home/jys0914/2026-ysjin-seqdep"
DB_DIR="$(readlink -f "${HOME}/.conda")/envs/metaphlan/lib/python3.13/site-packages/metaphlan/metaphlan_databases"

cd "${PROJECT_DIR}"

module load conda/base/latest
# shellcheck disable=SC1091
source "$(dirname "$(command -v conda)")/../etc/profile.d/conda.sh"
conda activate metaphlan

echo "[$(date)] Removing old database at ${DB_DIR}"
rm -rf "${DB_DIR}"

echo "[$(date)] Running metaphlan --install"
metaphlan --install --verbose

echo "[$(date)] Done"
metaphlan --version
