# Snakemake workflow: FastQC → MultiQC (raw) → fastp (trim+dedup) → FastQC (trimmed) → MultiQC (trimmed)
# Input: FASTQs in test-metagenome/downsample/SRR10692699 (or config runs)

configfile: "config.yaml"

RUNS = config.get("runs", ["SRR10692699"])
OUTDIR = config.get("outdir", "test-metagenome")
SEEDS = [str(s) for s in config.get("downsample_seeds", [11, 22, 33, 44, 55])]
DEPTH_LABELS = config.get("downsample_depths", ["10M", "20M", "30M", "40M", "50M"])

# All raw FASTQs (one per read)
RAW_FASTQS = expand(
    f"{OUTDIR}/downsample/{{run}}/{{run}}_{{depth}}_seed{{seed}}_R{{read}}.fastq",
    run=RUNS,
    depth=DEPTH_LABELS,
    seed=SEEDS,
    read=["1", "2"],
)

# FastQC reports (one per FASTQ)
FASTQC_RAW = expand(
    "qc/fastqc_raw/{run}_{depth}_seed{seed}_R{read}_fastqc.html",
    run=RUNS,
    depth=DEPTH_LABELS,
    seed=SEEDS,
    read=["1", "2"],
)

# FastQC reports for trimmed+deduped FASTQs
FASTQC_TRIMMED = expand(
    "qc/fastqc_trimmed/{run}_{depth}_seed{seed}_R{read}_fastqc.html",
    run=RUNS,
    depth=DEPTH_LABELS,
    seed=SEEDS,
    read=["1", "2"],
)

rule all:
    """Request raw MultiQC and trimmed+dedup MultiQC reports."""
    input:
        "qc/multiqc_raw/multiqc_report.html",
        "qc/multiqc_trimmed/multiqc_report.html",

# ---------------------------------------------------------------------------
# 1) FastQC on all raw FASTQs
# ---------------------------------------------------------------------------
rule fastqc_raw:
    """Run FastQC on each raw FASTQ."""
    input:
        fastq=f"{OUTDIR}/downsample/{{run}}/{{run}}_{{depth}}_seed{{seed}}_R{{read}}.fastq",
    output:
        html="qc/fastqc_raw/{run}_{depth}_seed{seed}_R{read}_fastqc.html",
        zip="qc/fastqc_raw/{run}_{depth}_seed{seed}_R{read}_fastqc.zip",
    log:
        "logs/fastqc/{run}_{depth}_seed{seed}_R{read}.log",
    shell:
        "fastqc -o qc/fastqc_raw {input.fastq} "
        ">> {log} 2>&1"

# ---------------------------------------------------------------------------
# 2) MultiQC on raw FastQC reports
# ---------------------------------------------------------------------------
rule multiqc_raw:
    """Aggregate raw FastQC reports with MultiQC."""
    input:
        FASTQC_RAW,
    output:
        report="qc/multiqc_raw/multiqc_report.html",
    log:
        "logs/multiqc_raw.log",
    shell:
        "multiqc qc/fastqc_raw -o qc/multiqc_raw --force >> {log} 2>&1"

# ---------------------------------------------------------------------------
# 3) fastp: adapter trimming, deduplication, and QC in one step
# ---------------------------------------------------------------------------
rule fastp:
    """Trim adapters, deduplicate, and run QC with fastp (one step)."""
    input:
        r1=f"{OUTDIR}/downsample/{{run}}/{{run}}_{{depth}}_seed{{seed}}_R1.fastq",
        r2=f"{OUTDIR}/downsample/{{run}}/{{run}}_{{depth}}_seed{{seed}}_R2.fastq",
    output:
        r1="trimmed/fastp/{run}_{depth}_seed{seed}_R1.fastq",
        r2="trimmed/fastp/{run}_{depth}_seed{seed}_R2.fastq",
        html="trimmed/fastp/{run}_{depth}_seed{seed}.fastp.html",
        json="trimmed/fastp/{run}_{depth}_seed{seed}.fastp.json",
    log:
        "logs/fastp/{run}_{depth}_seed{seed}.log",
    shell:
        "fastp -i {input.r1} -I {input.r2} "
        "--detect_adapter_for_pe --dedup "
        "-o {output.r1} -O {output.r2} "
        "-h {output.html} -j {output.json} "
        ">> {log} 2>&1"

# ---------------------------------------------------------------------------
# 4) FastQC on trimmed+deduped FASTQs
# ---------------------------------------------------------------------------
rule fastqc_trimmed:
    """Run FastQC on each trimmed+deduped FASTQ."""
    input:
        fastq="trimmed/fastp/{run}_{depth}_seed{seed}_R{read}.fastq",
    output:
        html="qc/fastqc_trimmed/{run}_{depth}_seed{seed}_R{read}_fastqc.html",
        zip="qc/fastqc_trimmed/{run}_{depth}_seed{seed}_R{read}_fastqc.zip",
    log:
        "logs/fastqc_trimmed/{run}_{depth}_seed{seed}_R{read}.log",
    shell:
        "fastqc -o qc/fastqc_trimmed {input.fastq} >> {log} 2>&1"

# ---------------------------------------------------------------------------
# 5) MultiQC on FastQC reports for trimmed+deduped FASTQs
# ---------------------------------------------------------------------------
rule multiqc_trimmed:
    """Aggregate FastQC (trimmed+deduped) reports with MultiQC."""
    input:
        FASTQC_TRIMMED,
    output:
        report="qc/multiqc_trimmed/multiqc_report.html",
    log:
        "logs/multiqc_trimmed.log",
    shell:
        "multiqc qc/fastqc_trimmed -o qc/multiqc_trimmed --force >> {log} 2>&1"
