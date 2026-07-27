#!/bin/bash
#SBATCH --job-name=dram_sample
#SBATCH --output=slurm_logs/dram_sample_%A_%a.out
#SBATCH --error=slurm_logs/dram_sample_%A_%a.err
# One array task per (omic, sample): runs annotate→distill→SCFA→featureCounts→
# aggregate for every assembled prefix under that sample (all depths × all
# seeds present on disk — no hard-coded seed).
#
# Array indices 0–15: MG0001–MG0008 (metaG) + MT0001–MT0008 (metaT).
# Adjust --cpus-per-task / partition / time to match your QOS.
#SBATCH --array=1
#SBATCH --partition=low
#SBATCH --account=publicgrp
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --time=7-00:00:00
#SBATCH --mail-user=ysjin@health.ucdavis.edu
#SBATCH --mail-type=ALL

set -euo pipefail

PROJECT_DIR="/home/jys0914/2026-ysjin-seqdep"
SMK_FILE="scripts/annotation_contigs_abundance.smk"

# (omic, sample) list — fixed order for reproducible array indices.
SAMPLES=(
  "metaG MG0001" "metaG MG0002" "metaG MG0003" "metaG MG0004"
  "metaG MG0005" "metaG MG0006" "metaG MG0007" "metaG MG0008"
  "metaT MT0001" "metaT MT0002" "metaT MT0003" "metaT MT0004"
  "metaT MT0005" "metaT MT0006" "metaT MT0007" "metaT MT0008"
)

ENTRY="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"
OMIC="${ENTRY%% *}"
SAMPLE="${ENTRY##* }"

cd "${PROJECT_DIR}"
export XDG_CACHE_HOME="${PROJECT_DIR}/.cache"
mkdir -p "${XDG_CACHE_HOME}" slurm_logs logs/dram_contigs_T2T

echo "[$(date)] Task ${SLURM_ARRAY_TASK_ID}: OMIC=${OMIC} SAMPLE=${SAMPLE}"

ASSEMBLY_SAMPLE_DIR="assembly_T2T/megahit/${OMIC}/${SAMPLE}"
if [[ ! -d "${ASSEMBLY_SAMPLE_DIR}" ]]; then
  echo "ERROR: assembly dir missing: ${ASSEMBLY_SAMPLE_DIR}" >&2
  exit 1
fi

TARGETS=()
N_PREFIX=0
while IFS= read -r CONTIGS; do
  [[ -f "${CONTIGS}" ]] || continue
  PREFIX="$(basename "$(dirname "${CONTIGS}")")"   # e.g. MG0001_2M_seed11
  BAM="assembly_T2T/alignments/${OMIC}/${SAMPLE}/${PREFIX}.sorted.bam"
  if [[ ! -f "${BAM}" ]]; then
    echo "  WARN: BAM missing for ${PREFIX}, abundance step may fail without BAM" >&2
  fi
  SAMPLE_DIR="annotation/dram_contigs_T2T/${OMIC}/${SAMPLE}/${PREFIX}"
  ABUND_DIR="abundance/dram_contigs_T2T/${OMIC}/${SAMPLE}/${PREFIX}"
  TARGETS+=(
    "${SAMPLE_DIR}/distillate/genome_stats.tsv"
    "${SAMPLE_DIR}/scfa_extra/scfa_per_scfa.tsv"
    "${ABUND_DIR}/kegg_module_abundance.tsv"
  )
  N_PREFIX=$((N_PREFIX + 1))
done < <(find "${ASSEMBLY_SAMPLE_DIR}" -mindepth 2 -maxdepth 2 -name final.contigs.fa 2>/dev/null | sort)

if [[ ${N_PREFIX} -eq 0 ]]; then
  echo "ERROR: no final.contigs.fa under ${ASSEMBLY_SAMPLE_DIR}" >&2
  exit 1
fi
echo "Found ${N_PREFIX} assemblies (depth×seed combos), ${#TARGETS[@]} leaf targets"

if command -v module >/dev/null 2>&1; then
  module load conda/base/latest
fi
if [[ -n "${CONDA_EXE:-}" ]]; then
  source "$(dirname "${CONDA_EXE}")/../etc/profile.d/conda.sh"
elif command -v conda >/dev/null 2>&1; then
  CONDA_EXE="$(command -v conda)"
  source "$(dirname "${CONDA_EXE}")/../etc/profile.d/conda.sh"
else
  echo "ERROR: conda not found" >&2
  exit 1
fi
conda activate dram

THREADS="${SLURM_CPUS_PER_TASK:-8}"

python -m snakemake \
  -s "${SMK_FILE}" \
  --cores "${THREADS}" \
  --config dram_threads="${THREADS}" \
            featurecounts_threads="${THREADS}" \
  --rerun-incomplete \
  --nolock \
  --printshellcmds \
  "${TARGETS[@]}"

echo "[$(date)] Task ${SLURM_ARRAY_TASK_ID} DONE: ${OMIC}/${SAMPLE} (${N_PREFIX} prefixes)"
