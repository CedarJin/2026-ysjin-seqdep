#!/bin/bash
#SBATCH --job-name=download_sra
#SBATCH --output=download_sra_%j.out
#SBATCH --error=download_sra_%j.err
#SBATCH --time=6:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --partition=low

### Email notifications
#SBATCH --mail-user=ysjin@health.ucdavis.edu
#SBATCH --mail-type=END,FAIL

# Load required modules
module load sratoolkit

# Set working directory
cd /home/jys0914/2026-ysjin-seqdep

# Set number of threads for fasterq-dump
export THREADS=8

# Run the download script
echo "Starting SRA download job at $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Running on node: $SLURM_NODELIST"
echo ""

./download_test_metagenome.sh

echo ""
echo "Download job completed at $(date)"
