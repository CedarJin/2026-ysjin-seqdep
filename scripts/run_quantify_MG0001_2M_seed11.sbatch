#!/usr/bin/env bash
#SBATCH --job-name=quant_MG0001_2M_seed11
#SBATCH --output=slurm_logs/slurm_quant_MG0001_2M_seed11_%j.out
#SBATCH --error=slurm_logs/slurm_quant_MG0001_2M_seed11_%j.err
#SBATCH --partition=low
#SBATCH --account=publicgrp
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --mail-user=ysjin@health.ucdavis.edu
#SBATCH --mail-type=ALL

# Single-sample quantification test: MG0001_2M_seed11
# Steps:
#   1) featureCounts on sorted BAM against the prodigal GFF produced by DRAM
#   2) aggregate gene counts -> TPM -> ID TPM -> SCFA + KEGG module abundance
#
# Outputs go to:
#   abundance/dram_contigs_T2T/metaG/MG0001/MG0001_2M_seed11/
#
# NOTE: This standalone script is superseded by the snakemake workflow in
#   scripts/annotation_contigs_abundance.smk
# which runs annotate -> distill -> verify -> SCFA -> featureCounts ->
# aggregate as a single DAG. Keep this around only for ad-hoc quantification
# tests against an existing DRAM annotation.

set -euo pipefail

PROJECT_DIR="/home/jys0914/2026-ysjin-seqdep"
cd "${PROJECT_DIR}"
mkdir -p slurm_logs

OMIC="metaG"
SAMPLE="MG0001"
PREFIX="MG0001_2M_seed11"
THREADS="${SLURM_CPUS_PER_TASK:-8}"

BAM="${PROJECT_DIR}/assembly_T2T/alignments/${OMIC}/${SAMPLE}/${PREFIX}.sorted.bam"
GFF="${PROJECT_DIR}/annotation/dram_contigs_T2T/${OMIC}/${SAMPLE}/${PREFIX}/genes.gff"
ANNOTATIONS="${PROJECT_DIR}/annotation/dram_contigs_T2T/${OMIC}/${SAMPLE}/${PREFIX}/distill_input_fixed/annotations.tsv"
DISTILL_CONFIG="${PROJECT_DIR}/annotation/dram_contigs_T2T/${OMIC}/${SAMPLE}/${PREFIX}/distill_input_fixed/distill_with_camper.config.json"
SCFA_REF="${PROJECT_DIR}/scripts/scfa_reference.tsv"
OUT_DIR="${PROJECT_DIR}/abundance/dram_contigs_T2T/${OMIC}/${SAMPLE}/${PREFIX}"
mkdir -p "${OUT_DIR}"

for f in "${BAM}" "${GFF}" "${ANNOTATIONS}" "${DISTILL_CONFIG}" "${SCFA_REF}"; do
  if [[ ! -f "${f}" ]]; then
    echo "ERROR: required input missing: ${f}" >&2
    exit 1
  fi
done

# ---------- conda init ----------
if command -v module >/dev/null 2>&1; then
  module load conda/base/latest || true
fi

if [[ -n "${CONDA_EXE:-}" ]]; then
  # shellcheck disable=SC1090,SC1091
  source "$(dirname "${CONDA_EXE}")/../etc/profile.d/conda.sh"
elif command -v conda >/dev/null 2>&1; then
  CONDA_EXE="$(command -v conda)"
  # shellcheck disable=SC1090,SC1091
  source "$(dirname "${CONDA_EXE}")/../etc/profile.d/conda.sh"
else
  echo "ERROR: conda not found" >&2
  exit 1
fi

# ---------- ensure quant env exists ----------
QUANT_ENV="quant"
QUANT_ENV_YAML="${PROJECT_DIR}/env/env_quant.yaml"
if ! conda env list | awk '{print $1}' | grep -qx "${QUANT_ENV}"; then
  echo "[$(date)] Creating conda env '${QUANT_ENV}' from ${QUANT_ENV_YAML}"
  conda env create -f "${QUANT_ENV_YAML}" -n "${QUANT_ENV}"
fi

# ---------- step 1: featureCounts ----------
# DRAM prodigal output prefixes every contig name with "final.contigs_" but
# the bowtie2 index used to build the BAM keeps the bare MEGAHIT names
# (k141_<n>). Generate a normalised GFF on the fly so featureCounts can
# match chromosomes; we keep the gene ID column unchanged so downstream
# aggregation still joins on the DRAM annotations.tsv gene_id.
NORM_GFF="${OUT_DIR}/genes.normalized.gff"
echo "[$(date)] Normalising GFF chromosome names -> ${NORM_GFF}"
awk -F'\t' -v OFS='\t' '
    /^#/ {print; next}
    { sub(/^final\.contigs_/, "", $1); print }
' "${GFF}" > "${NORM_GFF}"

# Sanity check the rename worked and now matches the BAM SQ names.
# Use awk-internal limit instead of `| head -N` so we don't trip SIGPIPE
# (141) under `set -euo pipefail`.
echo "[$(date)] GFF first seqids:"
awk -F'\t' '$1 !~ /^#/ { print $1; n++; if (n==3) exit }' "${NORM_GFF}"

GENE_COUNTS="${OUT_DIR}/gene_counts.tsv"
echo "[$(date)] Running featureCounts (threads=${THREADS})"

conda activate "${QUANT_ENV}"
featureCounts \
  -a "${NORM_GFF}" \
  -F GFF \
  -t CDS \
  -g ID \
  -T "${THREADS}" \
  -p --countReadPairs \
  -O --fraction \
  --primary \
  -o "${GENE_COUNTS}" \
  "${BAM}"

head -2 "${GENE_COUNTS}"
echo "[$(date)] featureCounts done -> ${GENE_COUNTS}"
echo "[$(date)] featureCounts summary:"
cat "${GENE_COUNTS}.summary"

# ---------- step 2: aggregate (TPM -> ID -> module) ----------
# Pull the three DRAM sheet paths from the distill config so we
# always use the same DBs that produced product.tsv.
read MODULE_STEP_FORM ETC_MODULE_DB FUNCTION_HEATMAP <<EOF
$(/quobyte/angelazgrp/ysjin/.conda/envs/dram/bin/python -c '
import json, sys
with open(sys.argv[1]) as f:
    cfg = json.load(f)
s = cfg["dram_sheets"]
print(s.get("module_step_form", ""), s.get("etc_module_database", ""), s.get("function_heatmap_form", ""))
' "${DISTILL_CONFIG}")
EOF

for kv in "module_step_form:${MODULE_STEP_FORM}" \
          "etc_module_database:${ETC_MODULE_DB}" \
          "function_heatmap_form:${FUNCTION_HEATMAP}"; do
  k="${kv%%:*}"; v="${kv#*:}"
  if [[ -z "${v}" || ! -f "${v}" ]]; then
    echo "WARN: ${k} not found at '${v}'; that section will be skipped." >&2
  fi
done

echo "[$(date)] Aggregating gene counts -> TPM -> ID -> module abundance"

# mag_annotator (used by aggregate_ko_module_abundance.py) lives in the
# DRAM env, so we call its python explicitly.
conda deactivate
DRAM_PY="/quobyte/angelazgrp/ysjin/.conda/envs/dram/bin/python"

AGG_ARGS=(
  "${PROJECT_DIR}/scripts/aggregate_ko_module_abundance.py"
  --counts "${GENE_COUNTS}"
  --annotations "${ANNOTATIONS}"
  --scfa-reference "${SCFA_REF}"
  --out-dir "${OUT_DIR}"
  --sample-label "${PREFIX}"
)
[[ -f "${MODULE_STEP_FORM}" ]] && AGG_ARGS+=(--module-step-form "${MODULE_STEP_FORM}")
[[ -f "${ETC_MODULE_DB}" ]]    && AGG_ARGS+=(--etc-module-database "${ETC_MODULE_DB}")
[[ -f "${FUNCTION_HEATMAP}" ]] && AGG_ARGS+=(--function-heatmap-form "${FUNCTION_HEATMAP}")

"${DRAM_PY}" "${AGG_ARGS[@]}"

echo "[$(date)] DONE. Outputs in ${OUT_DIR}"
ls -1 "${OUT_DIR}"
