#!/bin/bash

#SBATCH --job-name=rclone_trimmed
#SBATCH --output=slurm_logs/rclone_trimmed_%j.out
#SBATCH --error=slurm_logs/rclone_trimmed_%j.err
#SBATCH --partition=low
#SBATCH --account=publicgrp
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=7-00:00:00
#SBATCH --mail-user=ysjin@health.ucdavis.edu
#SBATCH --mail-type=ALL

set -euo pipefail

# rclone is usually provided via environment modules on compute nodes.
if command -v module >/dev/null 2>&1; then
  # Change if `module spider rclone` shows a different name on your cluster.
  module load rclone
fi

# ── edit these if needed ─────────────────────────────────────────────
PROJECT_DIR="${PROJECT_DIR:-/home/jys0914/2026-ysjin-seqdep}"
SRC="${PROJECT_DIR}/trimmed"
# Same remote path layout as your manual runs:
REMOTE_ROOT="${RCLONE_REMOTE_ROOT:-mybox:/2026-ysjin-seqdep-backup}"
DEST="${REMOTE_ROOT%/}/trimmed"

mkdir -p "${PROJECT_DIR}/slurm_logs"

if [[ ! -d "${SRC}" ]]; then
  echo "ERROR: source dir missing: ${SRC}" >&2
  exit 1
fi

echo "[$(date)] SRC=${SRC}"
echo "[$(date)] DEST=${DEST}"

if ! command -v rclone >/dev/null 2>&1; then
  echo "ERROR: rclone not in PATH (module load rclone / conda install rclone?)" >&2
  exit 1
fi

# --stats-log-level NOTICE: progress summary in the Slurm .out file
# --copy-links (-L): follow symlinks (same issue as reference/DRAM tmp/latest)
# Drop --copy-links if trimmed has no symlinks and you want to preserve links as-is.
rclone copy "${SRC}" "${DEST}" \
  --stats 1m \
  --stats-log-level NOTICE \
  --transfers 8 \
  --checkers 16 \
  --buffer-size 64M

echo "[$(date)] DONE"
