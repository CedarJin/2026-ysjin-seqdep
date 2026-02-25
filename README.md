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
5. assemble metagenomes,
6. bin contigs and annotate bins.

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
3. assembles host-filtered read pairs with MEGAHIT,
4. maps host-filtered reads to contigs and estimates contig depth,
5. generates microbial bins with MetaBAT2,
6. runs CheckM2 quality prediction on bins,
7. runs GTDB-Tk taxonomy classification for bins.

Default host reference paths are:
- `reference/human/GCF_000001405.40_GRCh38.p14_genomic.fna`
- Bowtie2 index prefix: `reference/human/GCF_000001405.40_GRCh38.p14_genomic`

You can override these in `config.yaml`:
```yaml
host_reference_fasta: reference/human/GCF_000001405.40_GRCh38.p14_genomic.fna
host_index_prefix: reference/human/GCF_000001405.40_GRCh38.p14_genomic
```

Database setup (required for CheckM2 and GTDB-Tk):
```bash
# 1) Make sure environment contains required tools
conda activate assemble
conda env update -f env_assembly_bin.yaml

# 2) Configure CheckM2 database
# Download/extract the CheckM2 database according to your cluster instructions.
# Then either set in config.yaml (recommended) or pass via environment/module setup.
# Example:
# checkm2_database_path: /path/to/CheckM2_database.dmnd
checkm2 database --download --path reference/



# 3) Download and configure GTDB-Tk database (save in reference/)
mkdir -p reference/gtdbtk
cd reference/gtdbtk

# Primary source
wget https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/auxillary_files/gtdbtk_package/full_package/gtdbtk_data.tar.gz

# Mirror for Australia (use this if primary is unavailable)
# wget https://data.gtdb.ecogenomic.org/releases/latest/auxillary_files/gtdbtk_package/full_package/gtdbtk_data.tar.gz

# Unarchive
tar xvzf gtdbtk_data.tar.gz

# Set GTDB-Tk database path (replace release* with your extracted release folder)
export GTDBTK_DATA_PATH=$PWD/release*
cd ../..
```

Optional `config.yaml` setting for CheckM2 DB path:
```yaml
checkm2_database_path: reference/CheckM2_database/uniref100.KO.1.dmnd
```

Key outputs after `assembly.smk`:
- host-filtered reads: `host_removed/`
- assemblies: `assembly/megahit/.../final.contigs.fa`
- bins: `bins/metabat2/`
- CheckM2 reports: `qc/checkm2/.../quality_report.tsv`
- GTDB-Tk results: `taxonomy/gtdbtk/`

Run assembly:
```bash
sbatch run_assembly.sbatch
```

### DRAM Workflow for annotation
Read the publication here: https://academic.oup.com/nar/article/48/16/8883/5884738?login=true
