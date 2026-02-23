# Snakemake workflow: host read removal (Bowtie2/Samtools) -> MEGAHIT assembly
# Input: trimmed/fastp/{run}_{depth}_seed{seed}_R1.fastq, R2.fastq (from fastq_qc.smk)
# Output: assembly/megahit/{run}_{depth}_seed{seed}/final.contigs.fa

configfile: "config.yaml"

RUNS = config.get("runs", ["SRR10692699"])
SEEDS = [str(s) for s in config.get("downsample_seeds", [11, 22, 33])]
DEPTH_LABELS = config.get("downsample_depths", ["10M", "20M", "30M", "40M", "50M"])
THREADS = config.get("threads", 8)
HOST_REFERENCE_FASTA = config.get(
    "host_reference_fasta",
    "reference/human/GCF_000001405.40_GRCh38.p14_genomic.fna",
)
HOST_INDEX_PREFIX = config.get(
    "host_index_prefix",
    "reference/human/GCF_000001405.40_GRCh38.p14_genomic",
)

HOST_INDEX_FILES = expand(
    "{prefix}.{ext}",
    prefix=HOST_INDEX_PREFIX,
    ext=["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"],
)

# All MEGAHIT contig outputs (one per sample)
MEGAHIT_CONTIGS = expand(
    "assembly/megahit/{run}_{depth}_seed{seed}/final.contigs.fa",
    run=RUNS,
    depth=DEPTH_LABELS,
    seed=SEEDS,
)

rule all:
    """Request all host-filtered MEGAHIT assemblies."""
    input:
        MEGAHIT_CONTIGS,

# ---------------------------------------------------------------------------
# Build Bowtie2 index for the host reference genome
# ---------------------------------------------------------------------------
rule build_host_index:
    """Build Bowtie2 index for host filtering."""
    input:
        ref=HOST_REFERENCE_FASTA,
    output:
        HOST_INDEX_FILES,
    params:
        prefix=HOST_INDEX_PREFIX,
    log:
        "logs/bowtie2/build_host_index.log",
    threads: THREADS
    shell:
        "mkdir -p $(dirname {params.prefix}) logs/bowtie2 && "
        "bowtie2-build --threads {threads} {input.ref} {params.prefix} >> {log} 2>&1"

# ---------------------------------------------------------------------------
# Remove host reads: keep only read pairs that do not map to the host genome
# ---------------------------------------------------------------------------
rule remove_host_reads:
    """Filter host-matching reads with Bowtie2 and keep non-host pairs."""
    input:
        r1="trimmed/fastp/{run}_{depth}_seed{seed}_R1.fastq",
        r2="trimmed/fastp/{run}_{depth}_seed{seed}_R2.fastq",
        idx=HOST_INDEX_FILES,
    output:
        r1="host_removed/{run}_{depth}_seed{seed}_R1.fastq",
        r2="host_removed/{run}_{depth}_seed{seed}_R2.fastq",
    params:
        idx_prefix=HOST_INDEX_PREFIX,
    log:
        "logs/bowtie2/{run}_{depth}_seed{seed}.log",
    threads: THREADS
    shell:
        "mkdir -p host_removed logs/bowtie2 && "
        "bowtie2 -x {params.idx_prefix} -1 {input.r1} -2 {input.r2} "
        "--threads {threads} --very-sensitive "
        "--un-conc host_removed/{wildcards.run}_{wildcards.depth}_seed{wildcards.seed}_R%.fastq "
        "-S /dev/null >> {log} 2>&1"

# ---------------------------------------------------------------------------
# MEGAHIT: assemble host-filtered paired-end FASTQs (default parameters)
# ---------------------------------------------------------------------------
rule megahit:
    """Run MEGAHIT on host-filtered R1/R2 with default parameters."""
    input:
        r1="host_removed/{run}_{depth}_seed{seed}_R1.fastq",
        r2="host_removed/{run}_{depth}_seed{seed}_R2.fastq",
    output:
        contigs="assembly/megahit/{run}_{depth}_seed{seed}/final.contigs.fa",
    log:
        "logs/megahit/{run}_{depth}_seed{seed}.log",
    threads: THREADS
    shell:
        "megahit -1 {input.r1} -2 {input.r2} "
        "-f "
        "-o $(dirname {output.contigs}) "
        "-t {threads} "
        ">> {log} 2>&1"
