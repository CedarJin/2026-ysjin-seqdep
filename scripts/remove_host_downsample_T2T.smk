import glob
import os
import re

# Remove host reads from all paired FASTQ files after fastp trimming.
# Input:  trimmed/downsample/fastp/{omic}/{sample}_{depth}_seed{seed}_R1.fastq and *_R2.fastq
# Output: cleandata_T2T/{omic}/{sample}/{prefix}_R1.fastq and {prefix}_R2.fastq

TRIM_ROOT = "trimmed/downsample/fastp"
OMICS = ["metaG", "metaT"]
HOST_REFERENCE_FASTA = "reference/human_T2T/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna"
HOST_INDEX_PREFIX = "reference/human_T2T/GCF_009914755.1_T2T-CHM13v2.0_genomic"
BOWTIE2_THREADS = 8
SEED_FILTER = {s.strip() for s in os.environ.get("SEED_FILTER", "").split(",") if s.strip()}

INDEX_EXTENSIONS = ["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"]
HOST_INDEX_FILES = expand("{prefix}.{ext}", prefix=HOST_INDEX_PREFIX, ext=INDEX_EXTENSIONS)


def list_prefixes(omic):
    pattern = os.path.join(TRIM_ROOT, omic, "*_R1.fastq")
    prefixes = []
    for r1 in glob.glob(pattern):
        r2 = r1.replace("_R1.fastq", "_R2.fastq")
        if not os.path.exists(r2):
            continue
        prefix = os.path.basename(r1).replace("_R1.fastq", "")
        sample_match = re.match(r"(.+?)_[^_]+_seed\d+$", prefix)
        if not sample_match:
            continue
        sample = sample_match.group(1)
        if SEED_FILTER:
            match = re.search(r"_seed(\d+)$", prefix)
            if not match or match.group(1) not in SEED_FILTER:
                continue
        prefixes.append((sample, prefix))
    return sorted(set(prefixes))


PAIRS = {omic: list_prefixes(omic) for omic in OMICS}
ALL_OUTPUTS = [
    f"cleandata_T2T/{omic}/{sample}/{prefix}_R{read}.fastq"
    for omic in OMICS
    for sample, prefix in PAIRS[omic]
    for read in ("1", "2")
]


rule all:
    input:
        ALL_OUTPUTS


rule build_host_index:
    input:
        ref=HOST_REFERENCE_FASTA
    output:
        HOST_INDEX_FILES
    params:
        prefix=HOST_INDEX_PREFIX
    log:
        "logs/bowtie2/build_host_index_downsample_T2T.log"
    threads:
        BOWTIE2_THREADS
    shell:
        """
        mkdir -p $(dirname {params.prefix}) logs/bowtie2
        bowtie2-build --threads {threads} {input.ref} {params.prefix} >> {log} 2>&1
        """


rule remove_host_reads:
    input:
        r1=lambda wc: f"{TRIM_ROOT}/{wc.omic}/{wc.prefix}_R1.fastq",
        r2=lambda wc: f"{TRIM_ROOT}/{wc.omic}/{wc.prefix}_R2.fastq",
        idx=HOST_INDEX_FILES
    output:
        r1="cleandata_T2T/{omic}/{sample}/{prefix}_R1.fastq",
        r2="cleandata_T2T/{omic}/{sample}/{prefix}_R2.fastq"
    params:
        idx_prefix=HOST_INDEX_PREFIX,
        out_prefix="cleandata_T2T/{omic}/{sample}/{prefix}_R"
    log:
        "logs/bowtie2/remove_host_T2T_{omic}_{sample}_{prefix}.log"
    threads:
        BOWTIE2_THREADS
    shell:
        """
        mkdir -p cleandata_T2T/{wildcards.omic}/{wildcards.sample} logs/bowtie2
        bowtie2 -x {params.idx_prefix} -1 {input.r1} -2 {input.r2} \
            --threads {threads} --very-sensitive-local \
            --un-conc {params.out_prefix}%.fastq \
            -S /dev/null >> {log} 2>&1
        """
