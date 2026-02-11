# Snakefile for downloading and processing SRA metagenome data
# Based on download_test_metagenome.sh

# Configuration
configfile: "config.yaml"

# SRA run IDs
RUNS = config.get("runs", ["SRR10692699", "SRR10692860"])
OUTDIR = config.get("outdir", "test-metagenome")
THREADS = config.get("threads", 8)

# Define final output files (compressed FASTQ)
FASTQ_FILES = expand(
    f"{OUTDIR}/{{run}}/{{run}}_{{read}}.fastq",
    run=RUNS,
    read=["1", "2"]
)

rule all:
    """
    Default rule: download and process all SRA runs.
    """
    input:
        FASTQ_FILES

rule prefetch_sra:
    """
    Download SRA files using prefetch.
    Creates the SRA file in the output directory.
    """
    output:
        sra_file = f"{OUTDIR}/{{run}}/{{run}}.sra"
    params:
        run = "{run}",
        outdir = OUTDIR,
        max_size = "30G"
    log:
        f"{OUTDIR}/{{run}}/prefetch.log"
    shell:
        """
        prefetch {params.run} --max-size {params.max_size} -O {params.outdir} 2>&1 | tee {log}
        """

rule fasterq_dump:
    """
    Convert SRA files to FASTQ format using fasterq-dump.
    Produces paired-end FASTQ files (_1.fastq and _2.fastq).
    """
    input:
        sra_file = f"{OUTDIR}/{{run}}/{{run}}.sra"
    output:
        fastq1 = f"{OUTDIR}/{{run}}/{{run}}_1.fastq",
        fastq2 = f"{OUTDIR}/{{run}}/{{run}}_2.fastq"
    params:
        run = "{run}",
        outdir = OUTDIR,
        threads = THREADS
    log:
        f"{OUTDIR}/{{run}}/fasterq_dump.log"
    shell:
        """
        fasterq-dump {params.run} -e {params.threads} -O {params.outdir}/{params.run} --progress 
        """

