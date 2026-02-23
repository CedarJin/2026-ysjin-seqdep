# 2026-ysjin-seqdep
Microbiome metagenome and metatranscriptome sequencing depth

## Test workflow using SRA data
### Set up sratoolkit

Please refer to https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit to download and install sratoolkit.

This project uses SRA test data to benchmark metagenome depth effects. The pipeline is:
1. download SRA reads,
2. downsample to fixed depths,
3. perform read QC and trimming,
4. remove host reads,
5. assemble metagenomes.

Before running workflows, confirm tools/environments are available:
- `run_testdata_downsample.sbatch` and `run_fastq_qc.sbatch`: QC/downsample environment.
- `run_assembly.sbatch`: `assemble` environment (megahit, bowtie2, samtools, etc.).


### Download and Downsample the SRA data to 10M, 20M, 30M, 40M, 50M

Use the helper script for a quick test run:
```bash
uv run ./download_test_metagenome.sh # alternatively, use the following snakemake workflow. 
```

Expected output:
- raw/downsampled FASTQs in `test-metagenome/downsample/<run>/`
- paired files named like `{run}_{depth}_seed{seed}_R1.fastq` and `_R2.fastq`

Snakemake workflow:
```bash
module load seqtk
uv run snakemake -s testdata_download_downsample.smk -n # -n : Dry run, just test if DAG can be build.
uv run snakemake -s testdata_download_downsample.smk -j 1 # -j : specify the # of cores used in the workflow. This workflow is memory intensive and time-consuming. You can submit a sbatch job and see if 4 cores and 64G mem are faster. (Only run this when you want to run snakemake interactively)
```

Recommended use:
- run `-n` first to validate DAG and file paths,
- then run with `sbatch` for long jobs on HPC.

Submit a batch job:
```bash
sbatch run_testdata_downsample.sbatch # use conda qc environment for this task
squeue -u $USER # see the progress
scancel <jobid> # cancel if needed # or use scancel -u $USER
```
Once downloading and downsampling are finished, spot-check several FASTQ files to confirm:
- R1/R2 files exist for each run-depth-seed combination,
- read counts match expected downsampling depth,
- read pairing is consistent.


### Fastq data quality control
`fastq_qc.smk` performs:
1. FastQC on raw downsampled reads,
2. MultiQC summary of raw read quality,
3. fastp trimming + deduplication,
4. FastQC on trimmed reads,
5. MultiQC summary of trimmed reads.

Outputs are primarily in:
- `trimmed/fastp/` (trimmed paired FASTQs and fastp reports),
- `qc/` (FastQC and MultiQC reports),
- `logs/` (per-rule logs).

Run:
```bash
sbatch run_fastq_qc.sbatch
```

### Remove host genome

Download the latest human reference genome from NCBI RefSeq
(`GCF_000001405.40_GRCh38.p14`):
```bash
mkdir -p reference/human
cd reference/human
wget -O GCF_000001405.40_GRCh38.p14_genomic.fna.gz \
  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz
gunzip -f GCF_000001405.40_GRCh38.p14_genomic.fna.gz
cd ../..
```

The assembly workflow in `assembly.smk` now:
1. builds a Bowtie2 index for the host reference,
2. removes host-mapped reads from `trimmed/fastp`,
3. assembles host-filtered read pairs with MEGAHIT.

Default host reference paths are:
- `reference/human/GCF_000001405.40_GRCh38.p14_genomic.fna`
- Bowtie2 index prefix: `reference/human/GCF_000001405.40_GRCh38.p14_genomic`

You can override these in `config.yaml`:
```yaml
host_reference_fasta: reference/human/GCF_000001405.40_GRCh38.p14_genomic.fna
host_index_prefix: reference/human/GCF_000001405.40_GRCh38.p14_genomic
```

Run assembly:
```bash
sbatch run_assembly.sbatch
```
