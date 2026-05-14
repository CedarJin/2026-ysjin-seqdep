# Snakemake workflow: DRAM annotation directly from T2T-cleaned assembled contigs
# Steps:
#   1) DRAM.py annotate on assembly_T2T/megahit final.contigs.fa
#   2) Normalize IDs for distill compatibility
#   3) Build per-sample distill config (with optional CAMPER heatmap injection)
#   4) DRAM.py distill
#
# Usage:
#   snakemake -s scripts/annotation_contigs.smk --cores 8

import glob
import os

configfile: "env/config.yaml"

ASSEMBLY_ROOT = config.get("assembly_root", "assembly_T2T")
ANNOTATION_ROOT = config.get("dram_contig_annotation_root", "annotation/dram_contigs_T2T")
THREADS = int(config.get("dram_threads", 8))
MIN_CONTIG_SIZE = int(config.get("dram_min_contig_size", 1000))
USE_CAMPER = config.get("dram_use_camper", True)
DRAM_CONFIG = config.get(
    "dram_config",
    "/home/jys0914/.conda/envs/dram/lib/python3.10/site-packages/mag_annotator/CONFIG",
)
DRAM_PY = config.get("dram_py", "/home/jys0914/.conda/envs/dram/bin/DRAM.py")


def list_t2t_contigs():
    pattern = os.path.join(ASSEMBLY_ROOT, "megahit", "*", "*", "*", "final.contigs.fa")
    contigs = []
    for path in glob.glob(pattern):
        parts = path.split(os.sep)
        try:
            megahit_idx = parts.index("megahit")
        except ValueError:
            continue
        if len(parts) < megahit_idx + 5:
            continue
        omic, sample, prefix = parts[megahit_idx + 1 : megahit_idx + 4]
        contigs.append((omic, sample, prefix))
    return sorted(set(contigs))


CONTIGS = list_t2t_contigs()
DRAM_CONTIG_ANNOTATIONS = expand(
    f"{ANNOTATION_ROOT}/{{omic}}/{{sample}}/{{prefix}}/annotations.tsv",
    zip,
    omic=[c[0] for c in CONTIGS],
    sample=[c[1] for c in CONTIGS],
    prefix=[c[2] for c in CONTIGS],
)
DRAM_CONTIG_DISTILL = expand(
    f"{ANNOTATION_ROOT}/{{omic}}/{{sample}}/{{prefix}}/distillate/genome_stats.tsv",
    zip,
    omic=[c[0] for c in CONTIGS],
    sample=[c[1] for c in CONTIGS],
    prefix=[c[2] for c in CONTIGS],
)


rule all:
    """Request contig-based DRAM annotations and distillate for all samples."""
    input:
        DRAM_CONTIG_ANNOTATIONS,
        DRAM_CONTIG_DISTILL,


rule dram_annotate_contigs:
    """Run DRAM.py annotate on assembled contigs for one sample."""
    input:
        contigs=f"{ASSEMBLY_ROOT}/megahit/{{omic}}/{{sample}}/{{prefix}}/final.contigs.fa",
    output:
        annotations=f"{ANNOTATION_ROOT}/{{omic}}/{{sample}}/{{prefix}}/annotations.tsv",
        trnas=f"{ANNOTATION_ROOT}/{{omic}}/{{sample}}/{{prefix}}/trnas.tsv",
        rrnas=f"{ANNOTATION_ROOT}/{{omic}}/{{sample}}/{{prefix}}/rrnas.tsv",
    params:
        outdir=f"{ANNOTATION_ROOT}/{{omic}}/{{sample}}/{{prefix}}",
        camper_flag=lambda wc: "--use_camper" if USE_CAMPER else "",
        min_contig_size=MIN_CONTIG_SIZE,
    log:
        "logs/dram_contigs_T2T/annotate_{omic}_{sample}_{prefix}.log",
    threads: THREADS
    shell:
        "mkdir -p $(dirname {log}) && "
        "rm -rf {params.outdir} && "
        "{DRAM_PY} annotate "
        "-i {input.contigs} "
        "-o {params.outdir} "
        "--min_contig_size {params.min_contig_size} "
        "--threads {threads} "
        "{params.camper_flag} "
        "--verbose "
        ">> {log} 2>&1"


rule fix_dram_contig_inputs_for_distill:
    """Normalize numeric-only fasta IDs across DRAM outputs before distill."""
    input:
        annotations=f"{ANNOTATION_ROOT}/{{omic}}/{{sample}}/{{prefix}}/annotations.tsv",
        trnas=f"{ANNOTATION_ROOT}/{{omic}}/{{sample}}/{{prefix}}/trnas.tsv",
        rrnas=f"{ANNOTATION_ROOT}/{{omic}}/{{sample}}/{{prefix}}/rrnas.tsv",
    output:
        annotations_fixed=f"{ANNOTATION_ROOT}/{{omic}}/{{sample}}/{{prefix}}/distill_input_fixed/annotations.tsv",
        trnas_fixed=f"{ANNOTATION_ROOT}/{{omic}}/{{sample}}/{{prefix}}/distill_input_fixed/trnas.tsv",
        rrnas_fixed=f"{ANNOTATION_ROOT}/{{omic}}/{{sample}}/{{prefix}}/distill_input_fixed/rrnas.tsv",
    log:
        "logs/dram_contigs_T2T/fix_distill_input_{omic}_{sample}_{prefix}.log",
    shell:
        "mkdir -p $(dirname {log}) && "
        "python scripts/normalize_dram_fasta_ids.py {input.annotations} {output.annotations_fixed} >> {log} 2>&1 && "
        "python scripts/normalize_dram_fasta_ids.py {input.trnas} {output.trnas_fixed} >> {log} 2>&1 && "
        "python scripts/normalize_dram_fasta_ids.py {input.rrnas} {output.rrnas_fixed} >> {log} 2>&1"


rule prepare_contig_distill_config_with_camper:
    """Build a sample-specific distill config that injects CAMPER into heatmaps."""
    input:
        annotations_fixed=f"{ANNOTATION_ROOT}/{{omic}}/{{sample}}/{{prefix}}/distill_input_fixed/annotations.tsv",
    output:
        config_loc=f"{ANNOTATION_ROOT}/{{omic}}/{{sample}}/{{prefix}}/distill_input_fixed/distill_with_camper.config.json",
        heatmap_form=f"{ANNOTATION_ROOT}/{{omic}}/{{sample}}/{{prefix}}/distill_input_fixed/function_heatmap_with_camper.tsv",
    log:
        "logs/dram_contigs_T2T/prepare_distill_config_{omic}_{sample}_{prefix}.log",
    shell:
        "mkdir -p $(dirname {log}) && "
        "python scripts/prepare_distill_config_with_camper.py "
        "--base-config {DRAM_CONFIG} "
        "--annotations {input.annotations_fixed} "
        "--out-config {output.config_loc} "
        "--out-heatmap {output.heatmap_form} "
        ">> {log} 2>&1"


rule dram_contig_distill:
    """Run DRAM.py distill to generate functional category summaries."""
    input:
        annotations=f"{ANNOTATION_ROOT}/{{omic}}/{{sample}}/{{prefix}}/distill_input_fixed/annotations.tsv",
        trnas=f"{ANNOTATION_ROOT}/{{omic}}/{{sample}}/{{prefix}}/distill_input_fixed/trnas.tsv",
        rrnas=f"{ANNOTATION_ROOT}/{{omic}}/{{sample}}/{{prefix}}/distill_input_fixed/rrnas.tsv",
        config_loc=f"{ANNOTATION_ROOT}/{{omic}}/{{sample}}/{{prefix}}/distill_input_fixed/distill_with_camper.config.json",
    output:
        genome_stats=f"{ANNOTATION_ROOT}/{{omic}}/{{sample}}/{{prefix}}/distillate/genome_stats.tsv",
    params:
        outdir=f"{ANNOTATION_ROOT}/{{omic}}/{{sample}}/{{prefix}}/distillate",
    log:
        "logs/dram_contigs_T2T/distill_{omic}_{sample}_{prefix}.log",
    shell:
        "mkdir -p $(dirname {log}) && "
        "rm -rf {params.outdir} && "
        "{DRAM_PY} distill "
        "-i {input.annotations} "
        "-o {params.outdir} "
        "--trna_path {input.trnas} "
        "--rrna_path {input.rrnas} "
        "--config_loc {input.config_loc} "
        ">> {log} 2>&1"
