# Snakemake workflow: FragGeneScanRs score scan on assembled contigs
# Steps:
#   1) Run FragGeneScanRs on MEGAHIT final.contigs.fa (output raw FAA + metadata)
#   2) Scan multiple score thresholds from metadata
#
# Usage:
#   snakemake -s fgs_contigs_score_scan.smk --cores 32

configfile: "config.yaml"

RUNS = config.get("runs", ["SRR10692699"])
SEEDS = [str(s) for s in config.get("downsample_seeds", [11, 22, 33])]
DEPTH_LABELS = config.get("downsample_depths", ["10M", "20M", "30M", "40M", "50M"])
THREADS = 16
FGS_CONTIGS_MODEL = config.get("fraggenescanrs_contigs_model", "complete")
FGS_SCAN_SCORES = [float(x) for x in config.get("fraggenescanrs_scan_scores", [1.20, 1.25, 1.30, 1.35])]

FGS_CONTIG_SCORE_SCAN = expand(
    "annotation/fgs_contigs/{run}_{depth}_seed{seed}/genes_predicted.score_scan.tsv",
    run=RUNS,
    depth=DEPTH_LABELS,
    seed=SEEDS,
)


rule all:
    """Request FragGeneScanRs score-scan summaries for all contig datasets."""
    input:
        FGS_CONTIG_SCORE_SCAN,


rule predict_genes_fraggenescanrs_on_contigs:
    """Predict proteins and metadata from MEGAHIT final contigs with FragGeneScanRs."""
    input:
        contigs="assembly/megahit/{run}_{depth}_seed{seed}/final.contigs.fa",
    output:
        genes_faa_raw="annotation/fgs_contigs/{run}_{depth}_seed{seed}/genes_predicted.raw.faa",
        genes_meta="annotation/fgs_contigs/{run}_{depth}_seed{seed}/genes_predicted.meta.tsv",
    log:
        "logs/fgs_contigs/fraggenescanrs_{run}_{depth}_seed{seed}.log",
    threads: THREADS
    shell:
        "mkdir -p $(dirname {output.genes_faa_raw}) logs/fgs_contigs && "
        "python scripts/run_fraggenescanrs.py {input.contigs} {output.genes_faa_raw} {FGS_CONTIGS_MODEL} "
        "--meta-out {output.genes_meta} --threads {threads} "
        ">> {log} 2>&1"


rule scan_fraggenescanrs_score_thresholds_on_contigs:
    """Scan FragGeneScanRs score thresholds for contig-based predictions."""
    input:
        genes_meta="annotation/fgs_contigs/{run}_{depth}_seed{seed}/genes_predicted.meta.tsv",
    output:
        summary="annotation/fgs_contigs/{run}_{depth}_seed{seed}/genes_predicted.score_scan.tsv",
    params:
        thresholds=",".join(str(x) for x in FGS_SCAN_SCORES),
    log:
        "logs/fgs_contigs/score_scan_{run}_{depth}_seed{seed}.log",
    shell:
        "python scripts/scan_fraggenescanrs_scores.py "
        "--meta-in {input.genes_meta} "
        "--summary-out {output.summary} "
        "--thresholds {params.thresholds} "
        ">> {log} 2>&1"
