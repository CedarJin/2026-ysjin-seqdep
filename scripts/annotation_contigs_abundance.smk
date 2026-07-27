# Snakemake workflow: DRAM annotation + abundance quantification
# Two output trees with mirrored {omic}/{sample}/{prefix} layout:
#   annotation/dram_contigs_T2T/...  -- annotate / distill / scfa_extra
#   abundance/dram_contigs_T2T/...   -- featureCounts + module TPM tables
#
# Steps:
#   1) DRAM.py annotate on assembly_T2T/megahit final.contigs.fa
#   2) Normalize IDs for distill compatibility
#   3) Build per-sample distill config (with optional CAMPER heatmap injection)
#   4) DRAM.py distill
#   5) extra_scfa_check.py extends SCFA pathway coverage against scfa_reference
#   6) featureCounts -> gene-level read counts from the sorted BAM     (-> abundance tree)
#   7) aggregate_ko_module_abundance.py -> gene/ID/module/SCFA/ETC TPM (-> abundance tree)
#
# Usage:
#   snakemake -s scripts/annotation_contigs_abundance.smk --cores 8

import glob
import json
import os

configfile: "env/config.yaml"

ASSEMBLY_ROOT = config.get("assembly_root", "assembly_T2T")
ANNOTATION_ROOT = config.get("dram_contig_annotation_root", "annotation/dram_contigs_T2T")
ABUNDANCE_ROOT = config.get("dram_contig_abundance_root", "abundance/dram_contigs_T2T")
ALIGN_ROOT = config.get("align_root", "assembly_T2T/alignments")
THREADS = int(config.get("dram_threads", 8))
MIN_CONTIG_SIZE = int(config.get("dram_min_contig_size", 1000))
USE_CAMPER = config.get("dram_use_camper", True)
DRAM_CONFIG = config.get(
    "dram_config",
    "/home/jys0914/.conda/envs/dram/lib/python3.10/site-packages/mag_annotator/CONFIG",
)
DRAM_PY = config.get("dram_py", "/home/jys0914/.conda/envs/dram/bin/DRAM.py")
DRAM_PYTHON = config.get("dram_python", "/home/jys0914/.conda/envs/dram/bin/python")
FEATURECOUNTS_BIN = config.get(
    "featurecounts_bin",
    "/quobyte/angelazgrp/ysjin/.conda/envs/quant/bin/featureCounts",
)
FEATURECOUNTS_THREADS = int(config.get("featurecounts_threads", 8))
SCFA_REFERENCE = config.get("scfa_reference", "scripts/scfa_reference.tsv")


def _load_dram_sheets(dram_config_path: str) -> dict:
    """Pull module_step_form / etc_module_database paths out of DRAM CONFIG.
    function_heatmap_form is generated per-sample, so we read it from the
    distill_with_camper.config.json instead."""
    try:
        with open(dram_config_path) as fh:
            cfg = json.load(fh)
    except Exception:
        return {}
    return cfg.get("dram_sheets", {}) or {}


_DRAM_SHEETS = _load_dram_sheets(DRAM_CONFIG)
MODULE_STEP_FORM = _DRAM_SHEETS.get("module_step_form", "")
ETC_MODULE_DATABASE = _DRAM_SHEETS.get("etc_module_database", "")


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
_OMICS = [c[0] for c in CONTIGS]
_SAMPLES = [c[1] for c in CONTIGS]
_PREFIXES = [c[2] for c in CONTIGS]

DRAM_CONTIG_ANNOTATIONS = expand(
    f"{ANNOTATION_ROOT}/{{omic}}/{{sample}}/{{prefix}}/annotations.tsv",
    zip, omic=_OMICS, sample=_SAMPLES, prefix=_PREFIXES,
)
DRAM_CONTIG_DISTILL = expand(
    f"{ANNOTATION_ROOT}/{{omic}}/{{sample}}/{{prefix}}/distillate/genome_stats.tsv",
    zip, omic=_OMICS, sample=_SAMPLES, prefix=_PREFIXES,
)
DRAM_SCFA_EXTRA = expand(
    f"{ANNOTATION_ROOT}/{{omic}}/{{sample}}/{{prefix}}/scfa_extra/scfa_per_scfa.tsv",
    zip, omic=_OMICS, sample=_SAMPLES, prefix=_PREFIXES,
)
DRAM_ABUNDANCE = expand(
    f"{ABUNDANCE_ROOT}/{{omic}}/{{sample}}/{{prefix}}/kegg_module_abundance.tsv",
    zip, omic=_OMICS, sample=_SAMPLES, prefix=_PREFIXES,
)


rule all:
    """End-to-end contig pipeline: annotate -> distill -> SCFA -> abundance."""
    input:
        DRAM_CONTIG_ANNOTATIONS,
        DRAM_CONTIG_DISTILL,
        DRAM_SCFA_EXTRA,
        DRAM_ABUNDANCE,


rule dram_annotate_contigs:
    """Run DRAM.py annotate on assembled contigs for one sample."""
    input:
        contigs=f"{ASSEMBLY_ROOT}/megahit/{{omic}}/{{sample}}/{{prefix}}/final.contigs.fa",
    output:
        annotations=f"{ANNOTATION_ROOT}/{{omic}}/{{sample}}/{{prefix}}/annotations.tsv",
        trnas=f"{ANNOTATION_ROOT}/{{omic}}/{{sample}}/{{prefix}}/trnas.tsv",
        rrnas=f"{ANNOTATION_ROOT}/{{omic}}/{{sample}}/{{prefix}}/rrnas.tsv",
        genes_gff=f"{ANNOTATION_ROOT}/{{omic}}/{{sample}}/{{prefix}}/genes.gff",
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


rule extra_scfa_check:
    """Re-annotate SCFA pathways against the curated scfa_reference.tsv.

    Outputs per-step / per-pathway / per-SCFA tables and a side-by-side
    comparison vs DRAM product.tsv columns."""
    input:
        product=f"{ANNOTATION_ROOT}/{{omic}}/{{sample}}/{{prefix}}/distillate/genome_stats.tsv",
        annotations_fixed=f"{ANNOTATION_ROOT}/{{omic}}/{{sample}}/{{prefix}}/distill_input_fixed/annotations.tsv",
        reference=SCFA_REFERENCE,
    output:
        per_scfa=f"{ANNOTATION_ROOT}/{{omic}}/{{sample}}/{{prefix}}/scfa_extra/scfa_per_scfa.tsv",
        per_pathway=f"{ANNOTATION_ROOT}/{{omic}}/{{sample}}/{{prefix}}/scfa_extra/scfa_per_pathway.tsv",
        per_step=f"{ANNOTATION_ROOT}/{{omic}}/{{sample}}/{{prefix}}/scfa_extra/scfa_per_step.tsv",
        diff=f"{ANNOTATION_ROOT}/{{omic}}/{{sample}}/{{prefix}}/scfa_extra/scfa_vs_dram_product.tsv",
    params:
        sample_dir=f"{ANNOTATION_ROOT}/{{omic}}/{{sample}}/{{prefix}}",
        out_dir=f"{ANNOTATION_ROOT}/{{omic}}/{{sample}}/{{prefix}}/scfa_extra",
    log:
        "logs/dram_contigs_T2T/extra_scfa_{omic}_{sample}_{prefix}.log",
    shell:
        "mkdir -p $(dirname {log}) {params.out_dir} && "
        "{DRAM_PYTHON} scripts/extra_scfa_check.py "
        "--sample-dir {params.sample_dir} "
        "--reference {input.reference} "
        "--out-dir {params.out_dir} "
        ">> {log} 2>&1"


rule featurecounts:
    """Gene-level read counts from the sorted BAM.

    DRAM prodigal prefixes every contig name with `final.contigs_`, while
    the bowtie2 index used for mapping kept the bare MEGAHIT names. We
    therefore generate a normalised GFF on the fly that strips the prefix
    from the first column so featureCounts can match chromosomes; the
    gene ID (GFF attribute `ID=`) is left untouched so downstream joins
    on DRAM annotations.tsv still work."""
    input:
        gff=f"{ANNOTATION_ROOT}/{{omic}}/{{sample}}/{{prefix}}/genes.gff",
        bam=f"{ALIGN_ROOT}/{{omic}}/{{sample}}/{{prefix}}.sorted.bam",
    output:
        counts=f"{ABUNDANCE_ROOT}/{{omic}}/{{sample}}/{{prefix}}/gene_counts.tsv",
        summary=f"{ABUNDANCE_ROOT}/{{omic}}/{{sample}}/{{prefix}}/gene_counts.tsv.summary",
        normalised_gff=f"{ABUNDANCE_ROOT}/{{omic}}/{{sample}}/{{prefix}}/genes.normalized.gff",
    params:
        outdir=f"{ABUNDANCE_ROOT}/{{omic}}/{{sample}}/{{prefix}}",
    log:
        "logs/dram_contigs_T2T/featurecounts_{omic}_{sample}_{prefix}.log",
    threads: FEATURECOUNTS_THREADS
    shell:
        "mkdir -p $(dirname {log}) {params.outdir} && "
        "awk -F'\\t' -v OFS='\\t' '"
        "/^#/ {{print; next}} "
        "{{ sub(/^final\\.contigs_/, \"\", $1); print }}' "
        "{input.gff} > {output.normalised_gff} && "
        "{FEATURECOUNTS_BIN} "
        "-a {output.normalised_gff} "
        "-F GFF "
        "-t CDS "
        "-g ID "
        "-T {threads} "
        "-p --countReadPairs "
        "-O --fraction "
        "--primary "
        "-o {output.counts} "
        "{input.bam} "
        ">> {log} 2>&1 && "
        "cat {output.summary} >> {log}"


rule aggregate_abundance:
    """Gene counts -> TPM -> ID TPM -> KEGG / ETC / SCFA / functional abundance.

    Pulls the KEGG `module_step_form` and `etc_module_database` paths from
    the DRAM CONFIG and the sample-specific `function_heatmap_with_camper`
    from the distill_input_fixed directory."""
    input:
        counts=f"{ABUNDANCE_ROOT}/{{omic}}/{{sample}}/{{prefix}}/gene_counts.tsv",
        annotations_fixed=f"{ANNOTATION_ROOT}/{{omic}}/{{sample}}/{{prefix}}/distill_input_fixed/annotations.tsv",
        function_heatmap=f"{ANNOTATION_ROOT}/{{omic}}/{{sample}}/{{prefix}}/distill_input_fixed/function_heatmap_with_camper.tsv",
        scfa_reference=SCFA_REFERENCE,
    output:
        gene_tpm=f"{ABUNDANCE_ROOT}/{{omic}}/{{sample}}/{{prefix}}/gene_tpm.tsv",
        id_tpm=f"{ABUNDANCE_ROOT}/{{omic}}/{{sample}}/{{prefix}}/id_tpm.tsv",
        kegg=f"{ABUNDANCE_ROOT}/{{omic}}/{{sample}}/{{prefix}}/kegg_module_abundance.tsv",
        etc=f"{ABUNDANCE_ROOT}/{{omic}}/{{sample}}/{{prefix}}/etc_complex_abundance.tsv",
        functional=f"{ABUNDANCE_ROOT}/{{omic}}/{{sample}}/{{prefix}}/functional_abundance.tsv",
        scfa_pathway=f"{ABUNDANCE_ROOT}/{{omic}}/{{sample}}/{{prefix}}/scfa_per_pathway_abundance.tsv",
        scfa_per_scfa=f"{ABUNDANCE_ROOT}/{{omic}}/{{sample}}/{{prefix}}/scfa_per_scfa_abundance.tsv",
        scfa_step=f"{ABUNDANCE_ROOT}/{{omic}}/{{sample}}/{{prefix}}/scfa_per_step_abundance.tsv",
    params:
        outdir=f"{ABUNDANCE_ROOT}/{{omic}}/{{sample}}/{{prefix}}",
        sample_label="{prefix}",
        module_step_form=MODULE_STEP_FORM,
        etc_module_database=ETC_MODULE_DATABASE,
    log:
        "logs/dram_contigs_T2T/aggregate_abundance_{omic}_{sample}_{prefix}.log",
    shell:
        "mkdir -p $(dirname {log}) && "
        "{DRAM_PYTHON} scripts/aggregate_ko_module_abundance.py "
        "--counts {input.counts} "
        "--annotations {input.annotations_fixed} "
        "--scfa-reference {input.scfa_reference} "
        "--module-step-form {params.module_step_form} "
        "--etc-module-database {params.etc_module_database} "
        "--function-heatmap-form {input.function_heatmap} "
        "--out-dir {params.outdir} "
        "--sample-label {params.sample_label} "
        ">> {log} 2>&1"
