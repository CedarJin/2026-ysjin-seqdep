# Snakefile for downloading and processing SRA metagenome data
# Based on download_test_metagenome.sh

# Configuration
configfile: "config.yaml"

# SRA run IDs
RUNS = config.get("runs", ["SRR10692699"])
OUTDIR = config.get("outdir", "test-metagenome")
THREADS = config.get("threads", 8)
SEEDS = [str(s) for s in config.get("downsample_seeds", [11, 22, 33, 44, 55])]
DEPTH_LABELS = config.get("downsample_depths", ["10M", "20M", "30M", "40M", "50M"])
DEPTH_TO_READS = {
    "10M": 10000000,
    "20M": 20000000,
    "30M": 30000000,
    "40M": 40000000,
    "50M": 50000000,
}

# Define final output files (compressed FASTQ)
FASTQ_FILES = expand(
    f"{OUTDIR}/{{run}}/{{run}}_{{read}}.fastq",
    run=RUNS,
    read=["1", "2"]
)

DOWNSAMPLED_FASTQ_FILES = expand(
    f"{OUTDIR}/downsample/{{run}}/{{run}}_{{depth}}_seed{{seed}}_R{{read}}.fastq",
    run=RUNS,
    depth=DEPTH_LABELS,
    seed=SEEDS,
    read=["1", "2"]
)

rule all:
    """
    Default rule: download and process all SRA runs.
    """
    input:
        FASTQ_FILES,
        DOWNSAMPLED_FASTQ_FILES

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

rule downsample_fastq:
    """
    Downsample paired FASTQ files with seqtk to fixed read-pair depths.
    """
    input:
        r1=f"{OUTDIR}/{{run}}/{{run}}_1.fastq",
        r2=f"{OUTDIR}/{{run}}/{{run}}_2.fastq"
    output:
        out_r1=f"{OUTDIR}/downsample/{{run}}/{{run}}_{{depth}}_seed{{seed}}_R1.fastq",
        out_r2=f"{OUTDIR}/downsample/{{run}}/{{run}}_{{depth}}_seed{{seed}}_R2.fastq"
    params:
        n_reads=lambda wildcards: DEPTH_TO_READS[wildcards.depth],
        seed=lambda wildcards: wildcards.seed
    shell:
        """
        mkdir -p $(dirname {output.out_r1})
        seqtk sample -s{params.seed} {input.r1} {params.n_reads} > {output.out_r1}
        seqtk sample -s{params.seed} {input.r2} {params.n_reads} > {output.out_r2}
        """
