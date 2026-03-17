# Snakemake workflow: DRAM annotation directly from assembled contigs
# Steps:
#   1) DRAM.py annotate on MEGAHIT final.contigs.fa
#   2) Normalize IDs for distill compatibility
#   3) Build per-sample distill config (with optional CAMPER heatmap injection)
#   4) DRAM.py distill
#
# Usage:
#   snakemake -s annotation_contigs.smk --cores 32

configfile: "config.yaml"

RUNS = config.get("runs", ["SRR10692699"])
SEEDS = [str(s) for s in config.get("downsample_seeds", [11, 22, 33])]
DEPTH_LABELS = config.get("downsample_depths", ["10M", "20M", "30M", "40M", "50M"])
THREADS = 32
USE_CAMPER = config.get("dram_use_camper", True)
DRAM_CONFIG = config.get(
    "dram_config",
    "/home/jys0914/.conda/envs/DRAM/lib/python3.10/site-packages/mag_annotator/CONFIG",
)
DRAM_PY = config.get("dram_py", "/home/jys0914/.conda/envs/DRAM/bin/DRAM.py")

DRAM_CONTIG_ANNOTATIONS = expand(
    "annotation/dram_contigs/{run}_{depth}_seed{seed}/annotations.tsv",
    run=RUNS,
    depth=DEPTH_LABELS,
    seed=SEEDS,
)
DRAM_CONTIG_DISTILL = expand(
    "annotation/dram_contigs/{run}_{depth}_seed{seed}/distillate/genome_stats.tsv",
    run=RUNS,
    depth=DEPTH_LABELS,
    seed=SEEDS,
)


rule all:
    """Request contig-based DRAM annotations and distillate for all samples."""
    input:
        DRAM_CONTIG_ANNOTATIONS,
        DRAM_CONTIG_DISTILL,


rule dram_annotate_contigs:
    """Run DRAM.py annotate on assembled contigs for one sample."""
    input:
        contigs="assembly/megahit/{run}_{depth}_seed{seed}/final.contigs.fa",
    output:
        annotations="annotation/dram_contigs/{run}_{depth}_seed{seed}/annotations.tsv",
        trnas="annotation/dram_contigs/{run}_{depth}_seed{seed}/trnas.tsv",
        rrnas="annotation/dram_contigs/{run}_{depth}_seed{seed}/rrnas.tsv",
    params:
        outdir="annotation/dram_contigs/{run}_{depth}_seed{seed}",
        camper_flag=lambda wc: "--use_camper" if USE_CAMPER else "",
    log:
        "logs/dram_contigs/annotate_{run}_{depth}_seed{seed}.log",
    threads: THREADS
    shell:
        "rm -rf {params.outdir} && "
        "{DRAM_PY} annotate "
        "-i {input.contigs} "
        "-o {params.outdir} "
        "--threads {threads} "
        "{params.camper_flag} "
        "--verbose "
        ">> {log} 2>&1"


rule fix_dram_contig_inputs_for_distill:
    """Normalize numeric-only fasta IDs across DRAM outputs before distill."""
    input:
        annotations="annotation/dram_contigs/{run}_{depth}_seed{seed}/annotations.tsv",
        trnas="annotation/dram_contigs/{run}_{depth}_seed{seed}/trnas.tsv",
        rrnas="annotation/dram_contigs/{run}_{depth}_seed{seed}/rrnas.tsv",
    output:
        annotations_fixed="annotation/dram_contigs/{run}_{depth}_seed{seed}/distill_input_fixed/annotations.tsv",
        trnas_fixed="annotation/dram_contigs/{run}_{depth}_seed{seed}/distill_input_fixed/trnas.tsv",
        rrnas_fixed="annotation/dram_contigs/{run}_{depth}_seed{seed}/distill_input_fixed/rrnas.tsv",
    log:
        "logs/dram_contigs/fix_distill_input_{run}_{depth}_seed{seed}.log",
    shell:
        "python scripts/normalize_dram_fasta_ids.py {input.annotations} {output.annotations_fixed} >> {log} 2>&1 && "
        "python scripts/normalize_dram_fasta_ids.py {input.trnas} {output.trnas_fixed} >> {log} 2>&1 && "
        "python scripts/normalize_dram_fasta_ids.py {input.rrnas} {output.rrnas_fixed} >> {log} 2>&1"


rule prepare_contig_distill_config_with_camper:
    """Build a sample-specific distill config that injects CAMPER into heatmaps."""
    input:
        annotations_fixed="annotation/dram_contigs/{run}_{depth}_seed{seed}/distill_input_fixed/annotations.tsv",
    output:
        config_loc="annotation/dram_contigs/{run}_{depth}_seed{seed}/distill_input_fixed/distill_with_camper.config.json",
        heatmap_form="annotation/dram_contigs/{run}_{depth}_seed{seed}/distill_input_fixed/function_heatmap_with_camper.tsv",
    log:
        "logs/dram_contigs/prepare_distill_config_{run}_{depth}_seed{seed}.log",
    shell:
        "python scripts/prepare_distill_config_with_camper.py "
        "--base-config {DRAM_CONFIG} "
        "--annotations {input.annotations_fixed} "
        "--out-config {output.config_loc} "
        "--out-heatmap {output.heatmap_form} "
        ">> {log} 2>&1"


rule dram_contig_distill:
    """Run DRAM.py distill to generate functional category summaries."""
    input:
        annotations="annotation/dram_contigs/{run}_{depth}_seed{seed}/distill_input_fixed/annotations.tsv",
        trnas="annotation/dram_contigs/{run}_{depth}_seed{seed}/distill_input_fixed/trnas.tsv",
        rrnas="annotation/dram_contigs/{run}_{depth}_seed{seed}/distill_input_fixed/rrnas.tsv",
        config_loc="annotation/dram_contigs/{run}_{depth}_seed{seed}/distill_input_fixed/distill_with_camper.config.json",
    output:
        genome_stats="annotation/dram_contigs/{run}_{depth}_seed{seed}/distillate/genome_stats.tsv",
    params:
        outdir="annotation/dram_contigs/{run}_{depth}_seed{seed}/distillate",
    log:
        "logs/dram_contigs/distill_{run}_{depth}_seed{seed}.log",
    shell:
        "rm -rf {params.outdir} && "
        "{DRAM_PY} distill "
        "-i {input.annotations} "
        "-o {params.outdir} "
        "--trna_path {input.trnas} "
        "--rrna_path {input.rrnas} "
        "--config_loc {input.config_loc} "
        ">> {log} 2>&1"
