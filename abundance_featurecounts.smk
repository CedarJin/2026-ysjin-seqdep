# Snakemake workflow: gene abundance counting across DRAM contig min-length annotations
#
# This workflow reuses one coordinate-sorted BAM and counts reads/fragments against
# DRAM-predicted CDS features for multiple annotation sets:
#   - annotation/dram_contigs_minlen/*_min200
#   - annotation/dram_contigs_minlen/*_min500
#   - annotation/dram_contigs_minlen/*_min1000
#   - annotation/dram_contigs/*            (treated as min contig length 2500)

configfile: "config.yaml"

RUN = config.get("abundance_run", config.get("runs", ["SRR10692699"])[0])
DEPTH = config.get("abundance_depth", config.get("downsample_depths", ["50M"])[0])
SEED = str(config.get("abundance_seed", config.get("downsample_seeds", [11])[0]))
THREADS = int(config.get("abundance_threads", 8))
MINLENS = [int(x) for x in config.get("abundance_min_contig_lengths", [200, 500, 1000, 2500])]

SAMPLE = f"{RUN}_{DEPTH}_seed{SEED}"
BAM = config.get("abundance_bam", f"assembly/alignments/{SAMPLE}.sorted.bam")
BAM_BAI = config.get("abundance_bai", f"{BAM}.bai")
FEATURECOUNTS = config.get(
    "featurecounts_bin",
    "/group/datalabgrp/ctbrown/ysjin/.conda/envs/abundance/bin/featureCounts",
)


def annotation_dir(minlen):
    if int(minlen) == 2500:
        return f"annotation/dram_contigs/{SAMPLE}"
    return f"annotation/dram_contigs_minlen/{SAMPLE}_min{minlen}"


def gff_path(minlen):
    return f"{annotation_dir(minlen)}/genes.gff"


def annotations_tsv_path(minlen):
    return f"{annotation_dir(minlen)}/annotations.tsv"


RAW_COUNTS = expand(
    f"abundance/featurecounts/{SAMPLE}/min{{minlen}}/featurecounts.txt",
    minlen=MINLENS,
)
ABUNDANCE_TABLES = expand(
    f"abundance/featurecounts/{SAMPLE}/min{{minlen}}/gene_abundance.tsv",
    minlen=MINLENS,
)


rule all:
    input:
        RAW_COUNTS,
        ABUNDANCE_TABLES,
        f"abundance/featurecounts/{SAMPLE}/combined_gene_abundance.tsv",
        f"abundance/featurecounts/{SAMPLE}/combined_gene_abundance_matrix.tsv",
        f"abundance/featurecounts/{SAMPLE}/min_contig_len_count_summary.tsv",
        f"abundance/featurecounts/{SAMPLE}/cazy_cpm.tsv",
        f"abundance/featurecounts/{SAMPLE}/camper_cpm.tsv",
        f"abundance/featurecounts/{SAMPLE}/cazy_cpm_distribution.tsv",
        f"abundance/featurecounts/{SAMPLE}/camper_cpm_distribution.tsv",
        f"abundance/featurecounts/{SAMPLE}/cazy_cpm_distribution.svg",
        f"abundance/featurecounts/{SAMPLE}/camper_cpm_distribution.svg",


rule normalize_dram_gff_for_featurecounts:
    input:
        gff=lambda wc: gff_path(wc.minlen),
    output:
        gff=f"abundance/featurecounts/{SAMPLE}/min{{minlen}}/genes.featurecounts.gff",
    log:
        f"logs/abundance_featurecounts/{SAMPLE}/min{{minlen}}/normalize_gff.log",
    shell:
        "mkdir -p $(dirname {output.gff}) $(dirname {log}) && "
        "python scripts/normalize_dram_gff_for_featurecounts.py "
        "{input.gff} {output.gff} "
        ">> {log} 2>&1"


rule featurecounts_gene_abundance:
    input:
        bam=BAM,
        bai=BAM_BAI,
        gff=f"abundance/featurecounts/{SAMPLE}/min{{minlen}}/genes.featurecounts.gff",
    output:
        counts=f"abundance/featurecounts/{SAMPLE}/min{{minlen}}/featurecounts.txt",
        summary=f"abundance/featurecounts/{SAMPLE}/min{{minlen}}/featurecounts.txt.summary",
    params:
        featurecounts=FEATURECOUNTS,
    log:
        f"logs/abundance_featurecounts/{SAMPLE}/min{{minlen}}/featurecounts.log",
    threads: THREADS
    shell:
        "mkdir -p $(dirname {output.counts}) $(dirname {log}) && "
        "{params.featurecounts} "
        "-F GFF "
        "-t CDS "
        "-g ID "
        "-p "
        "--countReadPairs "
        "-B "
        "-C "
        "-T {threads} "
        "-a {input.gff} "
        "-o {output.counts} "
        "{input.bam} "
        ">> {log} 2>&1"


rule build_gene_abundance_table:
    input:
        counts=f"abundance/featurecounts/{SAMPLE}/min{{minlen}}/featurecounts.txt",
        annotations=lambda wc: annotations_tsv_path(wc.minlen),
    output:
        table=f"abundance/featurecounts/{SAMPLE}/min{{minlen}}/gene_abundance.tsv",
    log:
        f"logs/abundance_featurecounts/{SAMPLE}/min{{minlen}}/merge_counts.log",
    shell:
        "python scripts/merge_featurecounts_with_dram.py "
        "--annotations {input.annotations} "
        "--counts {input.counts} "
        "--min-contig-len {wildcards.minlen} "
        "--output {output.table} "
        ">> {log} 2>&1"


rule summarize_featurecounts_across_minlen:
    input:
        tables=ABUNDANCE_TABLES,
        summaries=expand(
            f"abundance/featurecounts/{SAMPLE}/min{{minlen}}/featurecounts.txt.summary",
            minlen=MINLENS,
        ),
    output:
        long=f"abundance/featurecounts/{SAMPLE}/combined_gene_abundance.tsv",
        matrix=f"abundance/featurecounts/{SAMPLE}/combined_gene_abundance_matrix.tsv",
        summary=f"abundance/featurecounts/{SAMPLE}/min_contig_len_count_summary.tsv",
    log:
        f"logs/abundance_featurecounts/{SAMPLE}/summarize.log",
    shell:
        "python scripts/summarize_featurecounts_thresholds.py "
        "--tables {input.tables} "
        "--summaries {input.summaries} "
        "--output-long {output.long} "
        "--output-matrix {output.matrix} "
        "--output-summary {output.summary} "
        ">> {log} 2>&1"


rule summarize_cazy_camper_cpm:
    input:
        tables=ABUNDANCE_TABLES,
    output:
        cazy_table=f"abundance/featurecounts/{SAMPLE}/cazy_cpm.tsv",
        camper_table=f"abundance/featurecounts/{SAMPLE}/camper_cpm.tsv",
        cazy_dist=f"abundance/featurecounts/{SAMPLE}/cazy_cpm_distribution.tsv",
        camper_dist=f"abundance/featurecounts/{SAMPLE}/camper_cpm_distribution.tsv",
        cazy_svg=f"abundance/featurecounts/{SAMPLE}/cazy_cpm_distribution.svg",
        camper_svg=f"abundance/featurecounts/{SAMPLE}/camper_cpm_distribution.svg",
    log:
        f"logs/abundance_featurecounts/{SAMPLE}/functional_cpm.log",
    shell:
        "python scripts/summarize_functional_cpm.py "
        "--tables {input.tables} "
        "--cazy-output {output.cazy_table} "
        "--camper-output {output.camper_table} "
        "--cazy-dist-output {output.cazy_dist} "
        "--camper-dist-output {output.camper_dist} "
        "--cazy-plot-output {output.cazy_svg} "
        "--camper-plot-output {output.camper_svg} "
        ">> {log} 2>&1"
