#!/bin/bash

#SBATCH --job-name=check_scratch
#SBATCH --output=slurm_check_scratch_%j.out
#SBATCH --error=slurm_check_scratch_%j.err
#SBATCH --partition=low
#SBATCH --account=publicgrp
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:05:00

set -euo pipefail

echo "[$(date)] node=${SLURMD_NODENAME:-unknown}"
echo "SLURM_JOB_ID=${SLURM_JOB_ID:-unset}"
echo "SLURM_TMPDIR=${SLURM_TMPDIR:-unset}"
echo "TMPDIR=${TMPDIR:-unset}"
echo "TMP=${TMP:-unset}"
echo "TEMP=${TEMP:-unset}"

echo
echo "Candidate scratch directories:"
for d in "${SLURM_TMPDIR:-}" "${TMPDIR:-}" /scratch /lscratch /localscratch /tmp; do
  [[ -n "${d}" ]] || continue
  if [[ -e "${d}" ]]; then
    echo "--- ${d}"
    ls -ld "${d}" || true
    df -h "${d}" || true
    if [[ -w "${d}" ]]; then
      test_file="${d%/}/scratch_write_test_${SLURM_JOB_ID:-$$}"
      echo "write_test=attempting ${test_file}"
      echo "ok" > "${test_file}"
      rm -f "${test_file}"
      echo "write_test=ok"
    else
      echo "write_test=not_writable"
    fi
  fi
done

echo
echo "Mounts mentioning scratch/tmp/local:"
mount | grep -Ei 'scratch|tmp|local' || true
