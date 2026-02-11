#!/usr/bin/env bash

set -euo pipefail

# Download two test metagenome runs from NCBI SRA using fasterq-dump:
#   SRR10692699
#   SRR10692860
#
# Usage:
#   ./download_test_metagenome.sh
#   ./download_test_metagenome.sh /absolute/or/relative/output_dir

RUNS=("SRR10692699" "SRR10692860")
OUTDIR="${1:-test-metagenome}"
THREADS="${THREADS:-8}"


mkdir -p "${OUTDIR}"
cd "${OUTDIR}"

echo "Output directory: $(pwd)"
echo "Downloading and converting SRA runs to FASTQ: ${RUNS[*]}"
echo "Using ${THREADS} threads..."
echo "Started at: $(date)"
echo ""

for run in "${RUNS[@]}"; do
  echo "=========================================="
  echo "Processing ${run} at $(date)..."
  echo "=========================================="
  prefetch "${run}" --max-size 30G -O $(pwd)
  fasterq-dump "${run}" -e "${THREADS}" -O $(pwd)/${run} --progress
  echo "Completed ${run} at $(date)"
  echo ""
done

echo "Compressing FASTQ files..."
gzip -f SRR10692699*.fastq SRR10692860*.fastq 2>/dev/null || true

echo ""
echo "=========================================="
echo "Download completed at: $(date)"
echo "=========================================="
echo ""
echo "Final files:"
ls -lh SRR10692699* SRR10692860* 2>/dev/null || echo "No files found matching pattern"
echo "Total disk usage:"
du -sh .
