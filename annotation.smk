# Snakemake workflow: DRAM annotation of MetaBAT2 bins
# Steps:
#   1) DRAM.py annotate  -- annotate all bins per sample
#   2) DRAM.py distill   -- summarize annotations into functional categories
#
# Requires:
#   - bins from assembly.smk under bins/metabat2/{run}_{depth}_seed{seed}/
#   - DRAM env activated (conda activate DRAM)
#   - DRAM databases set up (DRAM-setup.py prepare_databases)
#
# Usage:
#   snakemake -s annotation.smk --cores 24

configfile: "config.yaml"

RUNS = config.get("runs", ["SRR10692699"])
SEEDS = [str(s) for s in config.get("downsample_seeds", [11, 22, 33])]
DEPTH_LABELS = config.get("downsample_depths", ["10M", "20M", "30M", "40M", "50M"])
THREADS = 32
# Extension used for DAS_Tool bins (must match dastool rule in assembly.smk)
BIN_EXT = config.get("bin_extension", "fa")
DRAM_CONFIG = config.get(
    "dram_config",
    "/home/jys0914/.conda/envs/DRAM/lib/python3.10/site-packages/mag_annotator/CONFIG",
)
DRAM_PY = config.get("dram_py", "/home/jys0914/.conda/envs/DRAM/bin/DRAM.py")

DRAM_ANNOTATIONS = expand(
    "annotation/dram/{run}_{depth}_seed{seed}/annotations.tsv",
    run=RUNS,
    depth=DEPTH_LABELS,
    seed=SEEDS,
)
DRAM_DISTILL = expand(
    "annotation/dram/{run}_{depth}_seed{seed}/distillate/genome_stats.tsv",
    run=RUNS,
    depth=DEPTH_LABELS,
    seed=SEEDS,
)

rule all:
    """Request DRAM annotation and distillate for all samples."""
    input:
        DRAM_ANNOTATIONS,
        DRAM_DISTILL,


# ---------------------------------------------------------------------------
# Pre-annotation: validate bins vs quality report and normalize quality report # for debug
# ---------------------------------------------------------------------------
rule validate_annotation_inputs:
    """Fail immediately if bin names and quality report Name column do not match (or type mismatch)."""
    input:
        bins_dir="bins/dastool/{run}_{depth}_seed{seed}/dastool_DASTool_bins",
        quality_report="qc/checkm2/{run}_{depth}_seed{seed}/quality_report.tsv",
    output:
        touch("annotation/dram/{run}_{depth}_seed{seed}/.validated"),
    log:
        "logs/dram/validate_{run}_{depth}_seed{seed}.log",
    shell:
        "mkdir -p $(dirname {output}) && "
        "python scripts/validate_annotation_inputs.py "
        "{input.bins_dir} {input.quality_report} --ext {BIN_EXT} "
        ">> {log} 2>&1 && touch {output}"

rule normalize_quality_report:
    """Force Name column to string so DRAM/distill do not see mixed types (e.g. 4 vs '4')."""
    input:
        "qc/checkm2/{run}_{depth}_seed{seed}/quality_report.tsv",
    output:
        "qc/checkm2/{run}_{depth}_seed{seed}/quality_report_for_dram.tsv",
    log:
        "logs/dram/normalize_qc_{run}_{depth}_seed{seed}.log",
    shell:
        "python scripts/normalize_quality_report.py "
        "{input} {output} >> {log} 2>&1"

# ---------------------------------------------------------------------------
# Annotate bins with DRAM.py annotate
# ---------------------------------------------------------------------------
rule dram_annotate:
    """Run DRAM.py annotate on all bins for one sample. Uses wrapper to fail early if glycan subfamily description is missing."""
    input:
        bins_dir="bins/dastool/{run}_{depth}_seed{seed}/dastool_DASTool_bins",
        gtdb_taxonomy="taxonomy/gtdbtk/{run}_{depth}_seed{seed}/gtdbtk.bac120.summary.tsv",
        checkm_quality="qc/checkm2/{run}_{depth}_seed{seed}/quality_report_for_dram.tsv",
        validated="annotation/dram/{run}_{depth}_seed{seed}/.validated",
    output:
        annotations="annotation/dram/{run}_{depth}_seed{seed}/annotations.tsv",
        trnas="annotation/dram/{run}_{depth}_seed{seed}/trnas.tsv",
        rrnas="annotation/dram/{run}_{depth}_seed{seed}/rrnas.tsv",
    params:
        outdir="annotation/dram/{run}_{depth}_seed{seed}",
    log:
        "logs/dram/annotate_{run}_{depth}_seed{seed}.log",
    threads: THREADS
    shell:
        "rm -rf {params.outdir} && "
        "python scripts/run_dram_annotate.py "
        "-i '{input.bins_dir}/*.{BIN_EXT}' "
        "-o {params.outdir} "
        "--threads {threads} "
        "--gtdb_taxonomy {input.gtdb_taxonomy} "
        "--checkm_quality {input.checkm_quality} "
        "--use_camper "
        "--verbose "
        ">> {log} 2>&1"


# ---------------------------------------------------------------------------
# Distill DRAM annotations into functional summaries
# ---------------------------------------------------------------------------
rule fix_dram_inputs_for_distill:
    """Normalize numeric-only fasta IDs across DRAM outputs before distill."""
    input:
        annotations="annotation/dram/{run}_{depth}_seed{seed}/annotations.tsv",
        trnas="annotation/dram/{run}_{depth}_seed{seed}/trnas.tsv",
        rrnas="annotation/dram/{run}_{depth}_seed{seed}/rrnas.tsv",
    output:
        annotations_fixed="annotation/dram/{run}_{depth}_seed{seed}/distill_input_fixed/annotations.tsv",
        trnas_fixed="annotation/dram/{run}_{depth}_seed{seed}/distill_input_fixed/trnas.tsv",
        rrnas_fixed="annotation/dram/{run}_{depth}_seed{seed}/distill_input_fixed/rrnas.tsv",
    log:
        "logs/dram/fix_distill_input_{run}_{depth}_seed{seed}.log",
    shell:
        "python scripts/normalize_dram_fasta_ids.py {input.annotations} {output.annotations_fixed} >> {log} 2>&1 && "
        "python scripts/normalize_dram_fasta_ids.py {input.trnas} {output.trnas_fixed} >> {log} 2>&1 && "
        "python scripts/normalize_dram_fasta_ids.py {input.rrnas} {output.rrnas_fixed} >> {log} 2>&1"

rule prepare_distill_config_with_camper:
    """Build a sample-specific distill config that injects CAMPER into product heatmaps."""
    input:
        annotations_fixed="annotation/dram/{run}_{depth}_seed{seed}/distill_input_fixed/annotations.tsv",
    output:
        config_loc="annotation/dram/{run}_{depth}_seed{seed}/distill_input_fixed/distill_with_camper.config.json",
        heatmap_form="annotation/dram/{run}_{depth}_seed{seed}/distill_input_fixed/function_heatmap_with_camper.tsv",
    log:
        "logs/dram/prepare_distill_config_{run}_{depth}_seed{seed}.log",
    shell:
        "python scripts/prepare_distill_config_with_camper.py "
        "--base-config {DRAM_CONFIG} "
        "--annotations {input.annotations_fixed} "
        "--out-config {output.config_loc} "
        "--out-heatmap {output.heatmap_form} "
        ">> {log} 2>&1"

rule dram_distill:
    """Run DRAM.py distill to generate functional category summaries."""
    input:
        annotations="annotation/dram/{run}_{depth}_seed{seed}/distill_input_fixed/annotations.tsv",
        trnas="annotation/dram/{run}_{depth}_seed{seed}/distill_input_fixed/trnas.tsv",
        rrnas="annotation/dram/{run}_{depth}_seed{seed}/distill_input_fixed/rrnas.tsv",
        config_loc="annotation/dram/{run}_{depth}_seed{seed}/distill_input_fixed/distill_with_camper.config.json",
    output:
        genome_stats="annotation/dram/{run}_{depth}_seed{seed}/distillate/genome_stats.tsv",
    params:
        outdir="annotation/dram/{run}_{depth}_seed{seed}/distillate",
    log:
        "logs/dram/distill_{run}_{depth}_seed{seed}.log",
    shell:
        "rm -rf {params.outdir} && "
        "{DRAM_PY} distill "
        "-i {input.annotations} "
        "-o {params.outdir} "
        "--trna_path {input.trnas} "
        "--rrna_path {input.rrnas} "
        "--config_loc {input.config_loc} "
        ">> {log} 2>&1"
