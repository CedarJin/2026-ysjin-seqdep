#!/bin/bash

#SBATCH --job-name=contig_len_scan
#SBATCH --output=slurm_contig_len_scan_%j.out
#SBATCH --error=slurm_contig_len_scan_%j.err
#SBATCH --partition=low
#SBATCH --account=publicgrp
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=24:00:00

set -euo pipefail

PROJECT_DIR="/home/jys0914/2026-ysjin-seqdep"
cd "${PROJECT_DIR}"

OUTDIR="analysis/contig_len_alignment_scan/SRR10692699_50M_seed11"
mkdir -p "${OUTDIR}"

python scripts/scan_contig_minlen_alignment.py \
  --contigs assembly/megahit/SRR10692699_50M_seed11/final.contigs.fa \
  --r1 host_removed/SRR10692699_50M_seed11_R1.fastq \
  --r2 host_removed/SRR10692699_50M_seed11_R2.fastq \
  --outdir "${OUTDIR}" \
  --thresholds 200,500,1000,2500 \
  --threads "${SLURM_CPUS_PER_TASK:-32}" \
  --reuse-log-200 logs/bowtie2/map_contigs_SRR10692699_50M_seed11.log

echo "Done. Results in ${OUTDIR}"
