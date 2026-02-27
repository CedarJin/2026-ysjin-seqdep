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
THREADS = 24
# Extension used for DAS_Tool bins (must match dastool rule in assembly.smk)
BIN_EXT = config.get("bin_extension", "fa")

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
# Annotate bins with DRAM.py annotate
# ---------------------------------------------------------------------------
rule dram_annotate:
    """Run DRAM.py annotate on all bins for one sample."""
    input:
        bins_dir="bins/dastool/{run}_{depth}_seed{seed}/dastool_DASTool_bins",
        gtdb_taxonomy="taxonomy/gtdbtk/{run}_{depth}_seed{seed}/gtdbtk.bac120.summary.tsv",
        checkm_quality="qc/checkm2/{run}_{depth}_seed{seed}/quality_report.tsv",
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
        "DRAM.py annotate "
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
rule dram_distill:
    """Run DRAM.py distill to generate functional category summaries."""
    input:
        annotations="annotation/dram/{run}_{depth}_seed{seed}/annotations.tsv",
        trnas="annotation/dram/{run}_{depth}_seed{seed}/trnas.tsv",
        rrnas="annotation/dram/{run}_{depth}_seed{seed}/rrnas.tsv",
    output:
        genome_stats="annotation/dram/{run}_{depth}_seed{seed}/distillate/genome_stats.tsv",
    params:
        outdir="annotation/dram/{run}_{depth}_seed{seed}/distillate",
    log:
        "logs/dram/distill_{run}_{depth}_seed{seed}.log",
    shell:
        "rm -rf {params.outdir} && "
        "DRAM.py distill "
        "-i {input.annotations} "
        "-o {params.outdir} "
        "--trna_path {input.trnas} "
        "--rrna_path {input.rrnas} "
        ">> {log} 2>&1"
