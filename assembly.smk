# Snakemake workflow: MEGAHIT assembly of trimmed metagenome FASTQs
# Input: trimmed/fastp/{run}_{depth}_seed{seed}_R1.fastq, R2.fastq (from fastq_qc.smk)
# Output: assembly/megahit/{run}_{depth}_seed{seed}/final.contigs.fa

configfile: "config.yaml"

RUNS = config.get("runs", ["SRR10692699"])
OUTDIR = config.get("outdir", "test-metagenome")
SEEDS = [str(s) for s in config.get("downsample_seeds", [11, 22, 33])]
DEPTH_LABELS = config.get("downsample_depths", ["10M", "20M", "30M", "40M", "50M"])
THREADS = config.get("threads", 8)

# All MEGAHIT contig outputs (one per sample)
MEGAHIT_CONTIGS = expand(
    "assembly/megahit/{run}_{depth}_seed{seed}/final.contigs.fa",
    run=RUNS,
    depth=DEPTH_LABELS,
    seed=SEEDS,
)

rule all:
    """Request all MEGAHIT assemblies."""
    input:
        MEGAHIT_CONTIGS,

# ---------------------------------------------------------------------------
# MEGAHIT: assemble paired-end trimmed FASTQs (default parameters)
# ---------------------------------------------------------------------------
rule megahit:
    """Run MEGAHIT on trimmed fastp R1/R2 with default parameters."""
    input:
        r1="trimmed/fastp/{run}_{depth}_seed{seed}_R1.fastq",
        r2="trimmed/fastp/{run}_{depth}_seed{seed}_R2.fastq",
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
