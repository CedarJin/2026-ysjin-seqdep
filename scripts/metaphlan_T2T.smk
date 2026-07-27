import glob
import os
import re

# Taxonomic profiling of T2T host-removed reads with MetaPhlAn 4.
# Input:  cleandata_T2T/{omic}/{sample}/{prefix}_R1.fastq and *_R2.fastq
# Output: taxonomy/metaphlan_T2T/{omic}/{sample}/{prefix}/profile.tsv
#         taxonomy/metaphlan_T2T/{omic}/{sample}/{prefix}/mapout.txt
#
# Optional filters (comma-separated env vars):
#   OMIC_FILTER    e.g. metaG
#   SAMPLE_FILTER  e.g. MG0001
#   PREFIX_FILTER  e.g. MG0001_2M_seed11
#   SEED_FILTER    e.g. 11,22
#
# Usage:
#   snakemake -s scripts/metaphlan_T2T.smk --cores 8
#   snakemake -s scripts/metaphlan_T2T.smk --cores 8 \
#     taxonomy/metaphlan_T2T/metaG/MG0001/MG0001_2M_seed11/profile.tsv

CLEANDATA_ROOT = "cleandata_T2T"
OUTPUT_ROOT = "taxonomy/metaphlan_T2T"
OMICS = ["metaG", "metaT"]
THREADS = int(os.environ.get("METAPHLAN_THREADS", "8"))
OMIC_FILTER = {s.strip() for s in os.environ.get("OMIC_FILTER", "").split(",") if s.strip()}
SAMPLE_FILTER = {s.strip() for s in os.environ.get("SAMPLE_FILTER", "").split(",") if s.strip()}
PREFIX_FILTER = {s.strip() for s in os.environ.get("PREFIX_FILTER", "").split(",") if s.strip()}
SEED_FILTER = {s.strip() for s in os.environ.get("SEED_FILTER", "").split(",") if s.strip()}


def list_pairs(omic):
    if OMIC_FILTER and omic not in OMIC_FILTER:
        return []

    pattern = os.path.join(CLEANDATA_ROOT, omic, "*", "*_R1.fastq")
    pairs = []
    for r1 in glob.glob(pattern):
        r2 = r1.replace("_R1.fastq", "_R2.fastq")
        if not os.path.exists(r2):
            continue

        sample = os.path.basename(os.path.dirname(r1))
        prefix = os.path.basename(r1).replace("_R1.fastq", "")

        if SAMPLE_FILTER and sample not in SAMPLE_FILTER:
            continue
        if PREFIX_FILTER and prefix not in PREFIX_FILTER:
            continue
        if SEED_FILTER:
            match = re.search(r"_seed(\d+)$", prefix)
            if not match or match.group(1) not in SEED_FILTER:
                continue

        pairs.append((sample, prefix))
    return sorted(set(pairs))


PAIRS = {omic: list_pairs(omic) for omic in OMICS}
ALL_PROFILES = [
    f"{OUTPUT_ROOT}/{omic}/{sample}/{prefix}/profile.tsv"
    for omic in OMICS
    for sample, prefix in PAIRS[omic]
]


rule all:
    input:
        ALL_PROFILES


rule metaphlan_profile:
    input:
        r1=f"{CLEANDATA_ROOT}/{{omic}}/{{sample}}/{{prefix}}_R1.fastq",
        r2=f"{CLEANDATA_ROOT}/{{omic}}/{{sample}}/{{prefix}}_R2.fastq",
    output:
        profile=f"{OUTPUT_ROOT}/{{omic}}/{{sample}}/{{prefix}}/profile.tsv",
        mapout=f"{OUTPUT_ROOT}/{{omic}}/{{sample}}/{{prefix}}/mapout.txt",
    log:
        "logs/metaphlan_T2T/{omic}_{sample}_{prefix}.log",
    threads:
        THREADS
    shell:
        """
        mkdir -p $(dirname {output.profile}) logs/metaphlan_T2T
        metaphlan \
            {input.r1},{input.r2} \
            --input_type fastq \
            --mapout {output.mapout} \
            --nproc {threads} \
            --sample_id {wildcards.prefix} \
            -o {output.profile} \
            >> {log} 2>&1
        """
