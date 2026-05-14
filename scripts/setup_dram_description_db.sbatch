#!/bin/bash

#SBATCH --job-name=dram_desc_db
#SBATCH --output=slurm_dram_desc_db_%j.out
#SBATCH --error=slurm_dram_desc_db_%j.err
#SBATCH --partition=low
#SBATCH --account=publicgrp
#SBATCH --cpus-per-task=4
#SBATCH --mem=512G
#SBATCH --time=48:00:00

### Email notifications
#SBATCH --mail-user=ysjin@health.ucdavis.edu
#SBATCH --mail-type=END,FAIL

set -euo pipefail

PROJECT_DIR="/home/jys0914/2026-ysjin-seqdep"
FINAL_DB_DIR="${PROJECT_DIR}/reference/DRAM_data"
FINAL_DESC_DB="${FINAL_DB_DIR}/description_db.sqlite"
DRAM_ENV="/quobyte/angelazgrp/ysjin/.conda/envs/dram"
DRAM_CONFIG="/quobyte/angelazgrp/ysjin/.conda/envs/dram/lib/python3.10/site-packages/mag_annotator/CONFIG"
DRAM_SETUP="${DRAM_ENV}/bin/DRAM-setup.py"

cd "${PROJECT_DIR}"

export XDG_CACHE_HOME="${PROJECT_DIR}/.cache"
mkdir -p "${XDG_CACHE_HOME}"

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

conda activate "${DRAM_ENV}"

if [[ ! -x "${DRAM_SETUP}" ]]; then
  echo "ERROR: DRAM-setup.py is not executable at ${DRAM_SETUP}" >&2
  exit 1
fi

check_description_db() {
  local db_path="$1"
  echo "[$(date)] Checking SQLite integrity for ${db_path}"
  "${DRAM_ENV}/bin/python" - "${db_path}" <<'PY'
import sqlite3
import sys
from pathlib import Path

db_path = Path(sys.argv[1])
if not db_path.is_file() or db_path.stat().st_size == 0:
    raise SystemExit(f"ERROR: SQLite DB is missing or empty: {db_path}")

conn = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
try:
    integrity = [row[0] for row in conn.execute("PRAGMA integrity_check;")]
    if integrity != ["ok"]:
        print("\n".join(integrity), file=sys.stderr)
        raise SystemExit(f"ERROR: SQLite integrity_check failed for {db_path}")

    for table in ("peptidase_description", "pfam_description", "dbcan_description"):
        count = conn.execute(f"SELECT COUNT(*) FROM {table};").fetchone()[0]
        if count <= 0:
            raise SystemExit(f"ERROR: table {table} is empty in {db_path}")
        print(f"{table}: {count}")
finally:
    conn.close()
PY
}

echo "[$(date)] Starting DRAM description DB setup on ${SLURMD_NODENAME:-unknown}"
echo "CPUs: ${SLURM_CPUS_PER_TASK:-1}, MEM: ${SLURM_MEM_PER_NODE:-unknown}"

# Prefer node-local scratch because SQLite writes are very slow on Quobyte.
# If the cluster does not provide local scratch, fall back to a Quobyte temp dir.
if [[ -n "${SLURM_TMPDIR:-}" && -d "${SLURM_TMPDIR}" && -w "${SLURM_TMPDIR}" ]]; then
  WORK_DIR="${SLURM_TMPDIR}/dram_description_db"
elif [[ -n "${TMPDIR:-}" && -d "${TMPDIR}" && -w "${TMPDIR}" ]]; then
  WORK_DIR="${TMPDIR}/dram_description_db_${SLURM_JOB_ID}"
elif [[ -d /tmp && -w /tmp ]]; then
  WORK_DIR="/tmp/${USER}/dram_description_db_${SLURM_JOB_ID}"
else
  WORK_DIR="${PROJECT_DIR}/reference/DRAM_data/tmp_description_db_${SLURM_JOB_ID}"
fi

mkdir -p "${WORK_DIR}"
TMP_DESC_DB="${WORK_DIR}/description_db.sqlite"
rm -f "${TMP_DESC_DB}"

echo "[$(date)] Working directory: ${WORK_DIR}"
df -h "${WORK_DIR}" "${FINAL_DB_DIR}" || true

echo "[$(date)] Populating DRAM description DB at ${TMP_DESC_DB}"
"${DRAM_SETUP}" update_description_db \
  --output_loc "${TMP_DESC_DB}" \
  --config_loc "${DRAM_CONFIG}"

check_description_db "${TMP_DESC_DB}"

if [[ -s "${FINAL_DESC_DB}" ]]; then
  BACKUP_DESC_DB="${FINAL_DESC_DB}.bak_${SLURM_JOB_ID:-manual}"
  echo "[$(date)] Backing up current description DB to ${BACKUP_DESC_DB}"
  cp -f "${FINAL_DESC_DB}" "${BACKUP_DESC_DB}"
fi

echo "[$(date)] Copying checked description DB to ${FINAL_DESC_DB}"
cp -f "${TMP_DESC_DB}" "${FINAL_DESC_DB}.new"
check_description_db "${FINAL_DESC_DB}.new"
mv -f "${FINAL_DESC_DB}.new" "${FINAL_DESC_DB}"
check_description_db "${FINAL_DESC_DB}"

echo "[$(date)] Updating DRAM CONFIG to final description DB path"
"${DRAM_ENV}/bin/python" - <<PY
import json
from pathlib import Path

config_path = Path("${DRAM_CONFIG}")
config = json.loads(config_path.read_text())
config["description_db"] = "${FINAL_DESC_DB}"
config_path.write_text(json.dumps(config, indent=2))
PY

echo "[$(date)] DRAM description DB setup finished: ${FINAL_DESC_DB}"
