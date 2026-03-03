# 2026-ysjin-seqdep
Microbiome metagenome and metatranscriptome sequencing depth for real data

## Real metagenomic data

### Sample metadata

Raw FASTQ files are downloaded and renamed using sequential IDs (`MG0001`–`MG0008`).
The mapping between these IDs and the customer labels is in `metaG_metadata.tsv`:

| sample_id | customer_label           | internal_id | R1                    | R2                    |
|-----------|--------------------------|-------------|-----------------------|-----------------------|
| MG0001    | GB2006.Day0.metagenome   | zr27836_1   | MG0001_R1.fastq.gz    | MG0001_R2.fastq.gz    |
| MG0002    | GB2006.Day180.metagenome | zr27836_2   | MG0002_R1.fastq.gz    | MG0002_R2.fastq.gz    |
| MG0003    | GB2033.Day0.metagenome   | zr27836_3   | MG0003_R1.fastq.gz    | MG0003_R2.fastq.gz    |
| MG0004    | GB2033.Day180.metagenome | zr27836_4   | MG0004_R1.fastq.gz    | MG0004_R2.fastq.gz    |
| MG0005    | GB2003.Day0.metagenome   | zr27836_5   | MG0005_R1.fastq.gz    | MG0005_R2.fastq.gz    |
| MG0006    | GB2003.Day180.metagenome | zr27836_6   | MG0006_R1.fastq.gz    | MG0006_R2.fastq.gz    |
| MG0007    | GB2032.Day0.metagenome   | zr27836_7   | MG0007_R1.fastq.gz    | MG0007_R2.fastq.gz    |
| MG0008    | GB2032.Day180.metagenome | zr27836_8   | MG0008_R1.fastq.gz    | MG0008_R2.fastq.gz    |

Download URLs are stored in `metaG_rawdata_links.tsv`.

### Download raw FASTQ files

```bash
bash download_metaG_rawdata.sh            # saves to rawdata/metaG/ (default)
```

The script reads `metaG_rawdata_links.tsv`, downloads each pair of reads, and saves them as
`MG{NNNN}_R1.fastq.gz` / `MG{NNNN}_R2.fastq.gz` in the output directory.
Files that already exist are skipped automatically (safe to re-run).

---

## SRA data download & downsample
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

##  Conda environment preperation
### qc environment
```bash
# Export to a yaml file
conda env export -n qc > env_qc.yaml

# Recreate from it (on same or another machine)
conda env create -f env_qc.yaml -n qc
```
### assembly environment
```bash
conda env create -f env_assembly_bin.yaml -n assemble
# if you already have the environment and want to update:
conda env update -f env_assembly_bin.yaml
```


### DRAM environment
```bash
wget https://raw.githubusercontent.com/WrightonLabCSU/DRAM/master/environment.yaml 
# replace the last row with: git+https://github.com/WrightonLabCSU/DRAM.git to fix syntax issues in DRAM-setup.py
# add pip setuptools<74
mv environment.yaml env_DRAM.yaml  
conda env create -f env_DRAM.yaml -n DRAM
```

## Fastq data quality control
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

## Remove host genome, assembly, binning

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
# conda env update -f env_assembly_bin.yaml

# 2) Configure CheckM2 database
# Download/extract the CheckM2 database according to your cluster instructions.
# Then either set in config.yaml (recommended) or pass via environment/module setup.
# Example:
# checkm2_database_path: /path/to/CheckM2_database.dmnd
checkm2 database --download --path reference/ # ONLY parent path!!



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
export GTDBTK_DATA_PATH=/home/jys0914/2026-ysjin-seqdep/reference/gtdbtk/release226 # MUST DO this!

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

## DRAM Workflow for annotation
Read the publication here: https://academic.oup.com/nar/article/48/16/8883/5884738?login=true

### activate the environment
```bash
conda activate DRAM
```

### Set up the DRAM database

NOTE: Setting up DRAM can take a long time (up to 5 hours) and uses a large amount of memory (512 gb) by default. To use less memory you can use the --skip_uniref flag which will reduce memory usage to ~64 gb if you do not provide KEGG Genes and 128 gb if you do. Depending on the number of processors which you tell it to use (using the --threads argument) and the speed of your internet connection. On a less than 5 year old server with 10 processors it takes about 2 hours to process the data when databases do not need to be downloaded.

```bash
sbatch setup_dram_db.sbatch # BUGGY!!! 
# Mofidied the script: 
# /home/jys0914/.conda/envs/DRAM/lib/python3.10/site-packages/mag_annotator/database_processing.py
# 1) Modified the 3 download dbCAN database URLs in the functions because the old ones don't work. 
# 2) Modified dbCAN version
# 3) Modified process_vogdb, merge_files(glob(path.join(hmm_dir, 'hmm/VOG*.hmm')), vog_hmms)
```
### Set up CAMPER database
```bash
git clone https://github.com/WrightonLabCSU/CAMPER.git
cd CAMPER
conda env update --name DRAM -f CAMPER_DRAMKit/environment.yaml
pip install CAMPER_DRAMKit/dist/camper_dramkit-1.0.13.tar.gz
```

## get CAMPER database from release:
https://github.com/WrightonLabCSU/CAMPER/releases


## DEBUG
site-packages/mag_annotator/database_handler.py
	Fix typo: camper_fa_db_cotoffs → camper_fa_db_cutoffs.
	•	Add handling for camper_tar_gz_loc: if detected in the wrapper, call process_camper_tar_gz to extract and process the .tar.gz, then pass the resulting paths to set_database_paths.
	•	Update the database_handler wrapper to support the CAMPER-1.0.0 structure (instead of CAMPER-1.0.0-beta.1) and accommodate varying subdirectory layouts.
	•	Implement more robust extraction logic in database_handler to ensure compatibility with CAMPER-1.0.0 and alternative directory structures.

```bash

DRAM-setup.py set_database_locations \
  --camper_tar_gz_loc /home/jys0914/2026-ysjin-seqdep/reference/DRAM_data/CAMPER_v1.0.0.tar.gz

# Verify
DRAM-setup.py print_config 2>&1 | grep -i camper
```

