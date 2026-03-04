# Snakemake workflow: DRAM assembly-free annotation from host-removed reads
# Steps:
#   1) Convert paired FASTQ reads to FASTA
#   2) Predict genes with FragGeneScanRs
#   3) DRAM.py annotate_genes on predicted proteins
#   4) DRAM.py distill summaries from annotations
#
# Usage:
#   snakemake -s annotation_assembly_free.smk --cores 24

configfile: "config.yaml"

RUNS = config.get("runs", ["SRR10692699"])
SEEDS = [str(s) for s in config.get("downsample_seeds", [11, 22, 33])]
DEPTH_LABELS = config.get("downsample_depths", ["10M", "20M", "30M", "40M", "50M"])
THREADS = 32
DRAM_CONFIG = config.get(
    "dram_config",
    "/home/jys0914/.conda/envs/DRAM/lib/python3.10/site-packages/mag_annotator/CONFIG",
)
DRAM_PY = config.get("dram_py", "/home/jys0914/.conda/envs/DRAM/bin/DRAM.py")
FGS_MODEL = config.get("fraggenescanrs_model", "illumina_10")

DRAM_AF_ANNOTATIONS = expand(
    "annotation/dram_assembly_free/{run}_{depth}_seed{seed}/annotations.tsv",
    run=RUNS,
    depth=DEPTH_LABELS,
    seed=SEEDS,
)
DRAM_AF_DISTILL = expand(
    "annotation/dram_assembly_free/{run}_{depth}_seed{seed}/distillate/genome_stats.tsv",
    run=RUNS,
    depth=DEPTH_LABELS,
    seed=SEEDS,
)


rule all:
    """Request assembly-free DRAM annotations and distillate for all samples."""
    input:
        DRAM_AF_ANNOTATIONS,
        DRAM_AF_DISTILL,


rule fastq_to_reads_fasta:
    """Convert paired host-removed FASTQ to a single FASTA for FragGeneScanRs."""
    input:
        r1="host_removed/{run}_{depth}_seed{seed}_R1.fastq",
        r2="host_removed/{run}_{depth}_seed{seed}_R2.fastq",
    output:
        reads_fa="annotation/dram_assembly_free/{run}_{depth}_seed{seed}/reads.fa",
    log:
        "logs/dram_assembly_free/reads_to_fasta_{run}_{depth}_seed{seed}.log",
    shell:
        "mkdir -p $(dirname {output.reads_fa}) && "
        "awk 'NR%4==1{{print \">\" substr($0,2)}} NR%4==2{{print}}' {input.r1} > {output.reads_fa} && "
        "awk 'NR%4==1{{print \">\" substr($0,2)}} NR%4==2{{print}}' {input.r2} >> {output.reads_fa} "
        ">> {log} 2>&1"


rule predict_genes_fraggenescanrs:
    """Predict proteins directly from read-derived FASTA with FragGeneScanRs."""
    input:
        reads_fa="annotation/dram_assembly_free/{run}_{depth}_seed{seed}/reads.fa",
    output:
        genes_faa="annotation/dram_assembly_free/{run}_{depth}_seed{seed}/genes_predicted.faa",
    log:
        "logs/dram_assembly_free/fraggenescanrs_{run}_{depth}_seed{seed}.log",
    shell:
        "python scripts/run_fraggenescanrs.py {input.reads_fa} {output.genes_faa} {FGS_MODEL} "
        ">> {log} 2>&1"


rule dram_annotate_genes:
    """Run DRAM.py annotate_genes for assembly-free predicted proteins."""
    input:
        genes_faa="annotation/dram_assembly_free/{run}_{depth}_seed{seed}/genes_predicted.faa",
    output:
        annotations="annotation/dram_assembly_free/{run}_{depth}_seed{seed}/annotations.tsv",
        genes_annotated="annotation/dram_assembly_free/{run}_{depth}_seed{seed}/genes.faa",
    params:
        outdir="annotation/dram_assembly_free/{run}_{depth}_seed{seed}",
    log:
        "logs/dram_assembly_free/annotate_genes_{run}_{depth}_seed{seed}.log",
    threads: THREADS
    shell:
        "rm -rf {params.outdir}/working_dir && "
        "{DRAM_PY} annotate_genes "
        "-i {input.genes_faa} "
        "-o {params.outdir} "
        "--threads {threads} "
        "--use_camper "
        "--verbose "
        ">> {log} 2>&1"


rule fix_dram_af_annotations_for_distill:
    """Normalize numeric-only fasta IDs in assembly-free annotations before distill."""
    input:
        annotations="annotation/dram_assembly_free/{run}_{depth}_seed{seed}/annotations.tsv",
    output:
        annotations_fixed="annotation/dram_assembly_free/{run}_{depth}_seed{seed}/distill_input_fixed/annotations.tsv",
    log:
        "logs/dram_assembly_free/fix_distill_input_{run}_{depth}_seed{seed}.log",
    shell:
        "python scripts/normalize_dram_fasta_ids.py {input.annotations} {output.annotations_fixed} "
        ">> {log} 2>&1"


rule prepare_af_distill_config_with_camper:
    """Build assembly-free distill config with CAMPER functions injected to heatmaps."""
    input:
        annotations_fixed="annotation/dram_assembly_free/{run}_{depth}_seed{seed}/distill_input_fixed/annotations.tsv",
    output:
        config_loc="annotation/dram_assembly_free/{run}_{depth}_seed{seed}/distill_input_fixed/distill_with_camper.config.json",
        heatmap_form="annotation/dram_assembly_free/{run}_{depth}_seed{seed}/distill_input_fixed/function_heatmap_with_camper.tsv",
    log:
        "logs/dram_assembly_free/prepare_distill_config_{run}_{depth}_seed{seed}.log",
    shell:
        "python scripts/prepare_distill_config_with_camper.py "
        "--base-config {DRAM_CONFIG} "
        "--annotations {input.annotations_fixed} "
        "--out-config {output.config_loc} "
        "--out-heatmap {output.heatmap_form} "
        ">> {log} 2>&1"


rule dram_af_distill:
    """Run DRAM.py distill on assembly-free annotation output."""
    input:
        annotations="annotation/dram_assembly_free/{run}_{depth}_seed{seed}/distill_input_fixed/annotations.tsv",
        config_loc="annotation/dram_assembly_free/{run}_{depth}_seed{seed}/distill_input_fixed/distill_with_camper.config.json",
    output:
        genome_stats="annotation/dram_assembly_free/{run}_{depth}_seed{seed}/distillate/genome_stats.tsv",
    params:
        outdir="annotation/dram_assembly_free/{run}_{depth}_seed{seed}/distillate",
    log:
        "logs/dram_assembly_free/distill_{run}_{depth}_seed{seed}.log",
    shell:
        "rm -rf {params.outdir} && "
        "{DRAM_PY} distill "
        "-i {input.annotations} "
        "-o {params.outdir} "
        "--config_loc {input.config_loc} "
        ">> {log} 2>&1"
