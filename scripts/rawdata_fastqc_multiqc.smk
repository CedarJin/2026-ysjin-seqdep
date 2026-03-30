# gunzip .fastq.gz in rawdata/{metaG,metaT}/ (removes .gz), then FastQC and per-omic MultiQC.
# Run from repo root: snakemake -s scripts/rawdata_fastqc_multiqc.smk --cores N

import glob
import os

THREADS_FASTQC = 2
OMICS = ["metaG", "metaT"]


def _samples(omic):
    """Return sample IDs detected from *_R1.fastq.gz and/or *_R1.fastq."""
    gz = [
        os.path.basename(p).replace("_R1.fastq.gz", "")
        for p in glob.glob(f"rawdata/{omic}/*_R1.fastq.gz")
    ]
    fq = [
        os.path.basename(p).replace("_R1.fastq", "")
        for p in glob.glob(f"rawdata/{omic}/*_R1.fastq")
    ]
    return sorted(set(gz) | set(fq))


SAMPLES = {omic: _samples(omic) for omic in OMICS}
ACTIVE_OMICS = [omic for omic in OMICS if SAMPLES[omic]]


rule all:
    input:
        expand("qc/multiqc/{omic}/multiqc_raw_report.html", omic=ACTIVE_OMICS),


rule decompress:
    input:
        "rawdata/{omic}/{sample}_R{read}.fastq.gz"
    output:
        "rawdata/{omic}/{sample}_R{read}.fastq"
    wildcard_constraints:
        omic="metaG|metaT",
        read="1|2"
    shell:
        "gunzip -f {input}"


rule fastqc:
    input:
        fq="rawdata/{omic}/{sample}_R{read}.fastq"
    output:
        html="qc/fastqc/raw/{omic}/{sample}_R{read}_fastqc.html",
        zip="qc/fastqc/raw/{omic}/{sample}_R{read}_fastqc.zip"
    log:
        "logs/fastqc/{omic}/{sample}_R{read}.log"
    wildcard_constraints:
        omic="metaG|metaT",
        read="1|2"
    threads:
        THREADS_FASTQC
    shell:
        "mkdir -p qc/fastqc/raw/{wildcards.omic} && "
        "fastqc {input.fq} -o qc/fastqc/raw/{wildcards.omic} -t {threads} >> {log} 2>&1"


rule multiqc:
    input:
        lambda wc: expand(
            "qc/fastqc/raw/{omic}/{sample}_R{read}_fastqc.html",
            omic=wc.omic,
            sample=SAMPLES[wc.omic],
            read=["1", "2"],
        )
    output:
        report="qc/multiqc/{omic}/multiqc_raw_report.html"
    log:
        "logs/multiqc/{omic}.log"
    wildcard_constraints:
        omic="metaG|metaT"
    shell:
        "mkdir -p qc/multiqc/{wildcards.omic} && "
        "multiqc qc/fastqc/raw/{wildcards.omic} -o qc/multiqc/{wildcards.omic} "
        "-n multiqc_raw_report --force >> {log} 2>&1"
