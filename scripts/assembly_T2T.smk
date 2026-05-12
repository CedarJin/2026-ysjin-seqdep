import glob
import os
import re

# Minimal assembly workflow for reads cleaned against the T2T host reference:
# cleandata_T2T paired reads -> megahit -> contig index -> read remapping

CLEANDATA_ROOT = "cleandata_T2T"
ASSEMBLY_ROOT = "assembly_T2T"
OMICS = ["metaG", "metaT"]
THREADS = int(os.environ.get("ASSEMBLY_THREADS", "8"))
MIN_CONTIG_LEN = int(os.environ.get("MEGAHIT_MIN_CONTIG_LEN", "1000"))
SEED_FILTER = {s.strip() for s in os.environ.get("SEED_FILTER", "").split(",") if s.strip()}
CONTIG_INDEX_EXT = ["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"]


def list_pairs(omic):
    pattern = os.path.join(CLEANDATA_ROOT, omic, "*", "*_R1.fastq")
    pairs = []
    for r1 in glob.glob(pattern):
        r2 = r1.replace("_R1.fastq", "_R2.fastq")
        if not os.path.exists(r2):
            continue

        sample = os.path.basename(os.path.dirname(r1))
        prefix = os.path.basename(r1).replace("_R1.fastq", "")

        if SEED_FILTER:
            match = re.search(r"_seed(\d+)$", prefix)
            if not match or match.group(1) not in SEED_FILTER:
                continue

        pairs.append((sample, prefix))
    return sorted(set(pairs))


PAIRS = {omic: list_pairs(omic) for omic in OMICS}
ALL_BAMS = [
    f"{ASSEMBLY_ROOT}/alignments/{omic}/{sample}/{prefix}.sorted.bam"
    for omic in OMICS
    for sample, prefix in PAIRS[omic]
]
ALL_BAIS = [f"{bam}.bai" for bam in ALL_BAMS]


rule all:
    input:
        ALL_BAMS,
        ALL_BAIS


rule megahit:
    input:
        r1=f"{CLEANDATA_ROOT}/{{omic}}/{{sample}}/{{prefix}}_R1.fastq",
        r2=f"{CLEANDATA_ROOT}/{{omic}}/{{sample}}/{{prefix}}_R2.fastq"
    output:
        contigs=f"{ASSEMBLY_ROOT}/megahit/{{omic}}/{{sample}}/{{prefix}}/final.contigs.fa"
    log:
        "logs/megahit_T2T/{omic}_{sample}_{prefix}.log"
    threads:
        THREADS
    params:
        min_contig_len=MIN_CONTIG_LEN
    shell:
        """
        mkdir -p $(dirname {output.contigs}) logs/megahit_T2T
        megahit -1 {input.r1} -2 {input.r2} -f -o $(dirname {output.contigs}) -t {threads} --min-contig-len {params.min_contig_len} >> {log} 2>&1
        """


rule build_contig_index:
    input:
        contigs=f"{ASSEMBLY_ROOT}/megahit/{{omic}}/{{sample}}/{{prefix}}/final.contigs.fa"
    output:
        expand(
            f"{ASSEMBLY_ROOT}/megahit/{{{{omic}}}}/{{{{sample}}}}/{{{{prefix}}}}/contigs.{{ext}}",
            ext=CONTIG_INDEX_EXT
        )
    params:
        index_prefix=f"{ASSEMBLY_ROOT}/megahit/{{omic}}/{{sample}}/{{prefix}}/contigs"
    log:
        "logs/bowtie2_T2T/contig_index_{omic}_{sample}_{prefix}.log"
    threads:
        THREADS
    shell:
        """
        mkdir -p logs/bowtie2_T2T
        bowtie2-build --threads {threads} {input.contigs} {params.index_prefix} >> {log} 2>&1
        """


rule map_reads_to_contigs:
    input:
        r1=f"{CLEANDATA_ROOT}/{{omic}}/{{sample}}/{{prefix}}_R1.fastq",
        r2=f"{CLEANDATA_ROOT}/{{omic}}/{{sample}}/{{prefix}}_R2.fastq",
        idx=expand(
            f"{ASSEMBLY_ROOT}/megahit/{{{{omic}}}}/{{{{sample}}}}/{{{{prefix}}}}/contigs.{{ext}}",
            ext=CONTIG_INDEX_EXT
        )
    output:
        bam=f"{ASSEMBLY_ROOT}/alignments/{{omic}}/{{sample}}/{{prefix}}.sorted.bam",
        bai=f"{ASSEMBLY_ROOT}/alignments/{{omic}}/{{sample}}/{{prefix}}.sorted.bam.bai"
    params:
        index_prefix=f"{ASSEMBLY_ROOT}/megahit/{{omic}}/{{sample}}/{{prefix}}/contigs"
    log:
        "logs/bowtie2_T2T/map_contigs_{omic}_{sample}_{prefix}.log"
    threads:
        THREADS
    shell:
        """
        mkdir -p "$(dirname '{output.bam}')" logs/bowtie2_T2T
        SAM_TMP="$(dirname '{output.bam}')/.st_sort_tmp_{wildcards.prefix}"
        rm -rf "$SAM_TMP" && mkdir -p "$SAM_TMP"
        bowtie2 -x {params.index_prefix} -1 {input.r1} -2 {input.r2} --threads {threads} --very-sensitive 2>> {log} | \
            samtools view -@ {threads} -b - | \
            samtools sort -@ {threads} -T "$SAM_TMP/samtools_sort" -o '{output.bam}' - >> {log} 2>&1
        samtools index -@ {threads} '{output.bam}' '{output.bai}' >> {log} 2>&1
        rm -rf "$SAM_TMP"
        """
