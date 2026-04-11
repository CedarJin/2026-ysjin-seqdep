# Workflow:
# FastQC (raw) -> MultiQC (raw) -> fastp -> FastQC (post-fastp) -> MultiQC (post-fastp)

import glob
import os

configfile: "env/config.yaml"

OMICS = config.get("omics", ["metaG", "metaT"])
DEPTHS = config.get("downsample_depths", ["10M", "20M", "30M", "40M"])
SEEDS = [str(s) for s in config.get("downsample_seeds", [11, 22, 33, 44, 55])]
FASTQC_THREADS = int(config.get("fastqc_threads", 2))
FASTP_THREADS = int(config.get("fastp_threads", 4))
QC_ROOT = "qc/downsample"
TRIM_ROOT = "trimmed/downsample/fastp"
LOG_ROOT = "logs/downsample"


def sample_ids(omic):
    dirs = [d for d in glob.glob(f"rawdata/downsample/{omic}/*") if os.path.isdir(d)]
    return sorted(s for s in (os.path.basename(d) for d in dirs) if s != "MT0002")
## delete if s != MT0002 when you want to add it back

SAMPLES = {omic: sample_ids(omic) for omic in OMICS}
ACTIVE_OMICS = [omic for omic in OMICS if SAMPLES[omic]]


def raw_fastqc_htmls(omic):
    return [
        f"{QC_ROOT}/fastqc_raw/{omic}/{s}_{d}_seed{k}_R{r}_fastqc.html"
        for s in SAMPLES[omic] for d in DEPTHS for k in SEEDS for r in ("1", "2")
    ]


def post_fastqc_htmls(omic):
    return [
        f"{QC_ROOT}/fastqc_post/{omic}/{s}_{d}_seed{k}_R{r}_fastqc.html"
        for s in SAMPLES[omic] for d in DEPTHS for k in SEEDS for r in ("1", "2")
    ]


rule all:
    input:
        expand(f"{QC_ROOT}/multiqc_raw/{{omic}}/multiqc_raw_report.html", omic=ACTIVE_OMICS),
        expand(f"{QC_ROOT}/multiqc_post/{{omic}}/multiqc_post_report.html", omic=ACTIVE_OMICS),


rule fastqc_raw:
    input:
        "rawdata/downsample/{omic}/{sample}/{sample}_{depth}_seed{seed}_R{read}.fastq"
    output:
        html=f"{QC_ROOT}/fastqc_raw/{{omic}}/{{sample}}_{{depth}}_seed{{seed}}_R{{read}}_fastqc.html",
        zip=f"{QC_ROOT}/fastqc_raw/{{omic}}/{{sample}}_{{depth}}_seed{{seed}}_R{{read}}_fastqc.zip"
    log:
        f"{LOG_ROOT}/fastqc_raw/{{omic}}/{{sample}}_{{depth}}_seed{{seed}}_R{{read}}.log"
    threads:
        FASTQC_THREADS
    shell:
        f"mkdir -p {QC_ROOT}/fastqc_raw/{{wildcards.omic}} && "
        f"fastqc -t {{threads}} -o {QC_ROOT}/fastqc_raw/{{wildcards.omic}} {{input}} >> {{log}} 2>&1"


rule multiqc_raw:
    input:
        lambda wc: raw_fastqc_htmls(wc.omic)
    output:
        report=f"{QC_ROOT}/multiqc_raw/{{omic}}/multiqc_raw_report.html"
    log:
        f"{LOG_ROOT}/multiqc_raw/{{omic}}.log"
    shell:
        f"mkdir -p {QC_ROOT}/multiqc_raw/{{wildcards.omic}} && "
        f"multiqc {QC_ROOT}/fastqc_raw/{{wildcards.omic}} -o {QC_ROOT}/multiqc_raw/{{wildcards.omic}} "
        "-n multiqc_raw_report --force >> {log} 2>&1"


rule fastp:
    input:
        r1="rawdata/downsample/{omic}/{sample}/{sample}_{depth}_seed{seed}_R1.fastq",
        r2="rawdata/downsample/{omic}/{sample}/{sample}_{depth}_seed{seed}_R2.fastq"
    output:
        r1=f"{TRIM_ROOT}/{{omic}}/{{sample}}_{{depth}}_seed{{seed}}_R1.fastq",
        r2=f"{TRIM_ROOT}/{{omic}}/{{sample}}_{{depth}}_seed{{seed}}_R2.fastq",
        html=f"{TRIM_ROOT}/{{omic}}/{{sample}}_{{depth}}_seed{{seed}}.fastp.html",
        json=f"{TRIM_ROOT}/{{omic}}/{{sample}}_{{depth}}_seed{{seed}}.fastp.json"
    log:
        f"{LOG_ROOT}/fastp/{{omic}}/{{sample}}_{{depth}}_seed{{seed}}.log"
    params:
        dedup_opt=lambda wc: "--dedup" if wc.omic == "metaG" else ""
    threads:
        FASTP_THREADS
    shell:
        f"mkdir -p {TRIM_ROOT}/{{wildcards.omic}} && "
        "fastp -w {threads} -i {input.r1} -I {input.r2} --detect_adapter_for_pe {params.dedup_opt} "
        "-o {output.r1} -O {output.r2} -h {output.html} -j {output.json} >> {log} 2>&1"


rule fastqc_post:
    input:
        f"{TRIM_ROOT}/{{omic}}/{{sample}}_{{depth}}_seed{{seed}}_R{{read}}.fastq"
    output:
        html=f"{QC_ROOT}/fastqc_post/{{omic}}/{{sample}}_{{depth}}_seed{{seed}}_R{{read}}_fastqc.html",
        zip=f"{QC_ROOT}/fastqc_post/{{omic}}/{{sample}}_{{depth}}_seed{{seed}}_R{{read}}_fastqc.zip"
    log:
        f"{LOG_ROOT}/fastqc_post/{{omic}}/{{sample}}_{{depth}}_seed{{seed}}_R{{read}}.log"
    threads:
        FASTQC_THREADS
    shell:
        f"mkdir -p {QC_ROOT}/fastqc_post/{{wildcards.omic}} && "
        f"fastqc -t {{threads}} -o {QC_ROOT}/fastqc_post/{{wildcards.omic}} {{input}} >> {{log}} 2>&1"


rule multiqc_post:
    input:
        lambda wc: post_fastqc_htmls(wc.omic)
    output:
        report=f"{QC_ROOT}/multiqc_post/{{omic}}/multiqc_post_report.html"
    log:
        f"{LOG_ROOT}/multiqc_post/{{omic}}.log"
    shell:
        f"mkdir -p {QC_ROOT}/multiqc_post/{{wildcards.omic}} && "
        f"multiqc {QC_ROOT}/fastqc_post/{{wildcards.omic}} -o {QC_ROOT}/multiqc_post/{{wildcards.omic}} "
        "-n multiqc_post_report --force >> {log} 2>&1"
