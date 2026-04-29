#!/usr/bin/env bash
# Download metatranscriptomic raw FASTQ files from a metaT rawdata links TSV,
# renaming them to MT{NNNN}_R1.fastq.gz / MT{NNNN}_R2.fastq.gz.
#
# Usage:
#   bash scripts/download_metaT_rawdata.sh                         # rawdata/metaT
#   bash scripts/download_metaT_rawdata.sh rawdata/metaT_batch2    # rawdata/metaT_batch2
#   bash scripts/download_metaT_rawdata.sh /your/path              # legacy: default links, custom output
#
# If the target directory contains metaT_rawdata_links.tsv or one *_rawdata_links*.tsv
# file, that file is used. Otherwise, the default rawdata/metaT links TSV is used.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
OUTDIR="${1:-${PROJECT_ROOT}/rawdata/metaT}"

if [[ "$OUTDIR" != /* ]]; then
    OUTDIR="${PROJECT_ROOT}/${OUTDIR}"
fi

DEFAULT_LINKS_TSV="${PROJECT_ROOT}/rawdata/metaT/metaT_rawdata_links.tsv"
LINKS_TSV=""
if [[ -f "${OUTDIR}/metaT_rawdata_links.tsv" ]]; then
    LINKS_TSV="${OUTDIR}/metaT_rawdata_links.tsv"
else
    shopt -s nullglob
    link_candidates=("${OUTDIR}"/*_rawdata_links*.tsv)
    shopt -u nullglob
    if [[ "${#link_candidates[@]}" -eq 1 ]]; then
        LINKS_TSV="${link_candidates[0]}"
    elif [[ "${#link_candidates[@]}" -eq 0 ]]; then
        LINKS_TSV="${DEFAULT_LINKS_TSV}"
    else
        echo "ERROR: Multiple rawdata links TSV files found in: $OUTDIR" >&2
        printf '  %s\n' "${link_candidates[@]}" >&2
        exit 1
    fi
fi

if [[ ! -f "$LINKS_TSV" ]]; then
    echo "ERROR: Rawdata links TSV not found: $LINKS_TSV" >&2
    exit 1
fi

mkdir -p "$OUTDIR"

echo "Downloading metatranscriptomic FASTQ files to: $OUTDIR"
echo "Using rawdata links TSV: $LINKS_TSV"
echo "----------------------------------------------------"

download_fastq() {
    local label="$1"
    local url="$2"
    local output="$3"
    local partial="${output}.part"

    if [[ -s "$output" ]]; then
        echo "  ${label} already exists, skipping: $output"
        return
    fi

    if [[ -f "$output" ]]; then
        echo "  ${label} exists but is empty/incomplete, re-downloading: $output"
        rm -f "$output"
    fi

    rm -f "$partial"
    echo "  Downloading ${label} -> $output"
    wget -q --show-progress -O "$partial" "$url"
    mv "$partial" "$output"
}

# Skip header line; iterate rows and assign sequential MT IDs
idx=0
while IFS=$'\t' read -r internal_id customer_label url_r1 url_r2; do
    # skip header
    [[ "$internal_id" == "Internal_id" ]] && continue
    # skip empty lines
    [[ -z "$internal_id" ]] && continue

    idx=$(( idx + 1 ))
    mt_id=$(printf "MT%04d" "$idx")

    out_r1="${OUTDIR}/${mt_id}_R1.fastq.gz"
    out_r2="${OUTDIR}/${mt_id}_R2.fastq.gz"

    echo "[${mt_id}] ${customer_label} (${internal_id})"

    download_fastq "R1" "$url_r1" "$out_r1"
    download_fastq "R2" "$url_r2" "$out_r2"

done < "$LINKS_TSV"

echo "----------------------------------------------------"
echo "Done. Files saved in: $OUTDIR"

echo ""
echo "Decompressing .fastq.gz files in: $OUTDIR"
echo "----------------------------------------------------"

for gz_file in "${OUTDIR}"/*.fastq.gz; do
    [[ -f "$gz_file" ]] || { echo "  No .fastq.gz files found."; break; }
    fastq_file="${gz_file%.gz}"
    if [[ -f "$fastq_file" ]]; then
        echo "  Already decompressed, skipping: $fastq_file"
    else
        echo "  Decompressing: $gz_file"
        gunzip -k "$gz_file"
    fi
done

echo "----------------------------------------------------"
echo "Decompression complete."
