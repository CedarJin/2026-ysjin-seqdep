#!/usr/bin/env bash
# Download metagenomic raw FASTQ files from metaG_rawdata_links.tsv,
# renaming them to MG{NNNN}_R1.fastq.gz / MG{NNNN}_R2.fastq.gz.
#
# Usage:
#   bash download_metaG_rawdata.sh [output_dir]
#
# Default output directory: rawdata/metaG

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LINKS_TSV="${SCRIPT_DIR}/metaG_rawdata_links.tsv"
OUTDIR="${1:-rawdata/metaG}"

mkdir -p "$OUTDIR"

echo "Downloading metagenomic FASTQ files to: $OUTDIR"
echo "----------------------------------------------------"

# Skip header line; iterate rows and assign sequential MG IDs
idx=0
while IFS=$'\t' read -r internal_id customer_label url_r1 url_r2; do
    # skip header
    [[ "$internal_id" == "Internal_id" ]] && continue
    # skip empty lines
    [[ -z "$internal_id" ]] && continue

    idx=$(( idx + 1 ))
    mg_id=$(printf "MG%04d" "$idx")

    out_r1="${OUTDIR}/${mg_id}_R1.fastq.gz"
    out_r2="${OUTDIR}/${mg_id}_R2.fastq.gz"

    echo "[${mg_id}] ${customer_label} (${internal_id})"

    if [[ -f "$out_r1" ]]; then
        echo "  R1 already exists, skipping: $out_r1"
    else
        echo "  Downloading R1 -> $out_r1"
        wget -q --show-progress -O "$out_r1" "$url_r1"
    fi

    if [[ -f "$out_r2" ]]; then
        echo "  R2 already exists, skipping: $out_r2"
    else
        echo "  Downloading R2 -> $out_r2"
        wget -q --show-progress -O "$out_r2" "$url_r2"
    fi

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
