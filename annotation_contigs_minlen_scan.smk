# Snakemake workflow: DRAM annotation runtime scan across contig min-length values
# Target sample is fixed to preserve comparability with the previous default run.

configfile: "config.yaml"

RUN = config.get("runtime_scan_run", "SRR10692699")
DEPTH = config.get("runtime_scan_depth", "50M")
SEED = str(config.get("runtime_scan_seed", 11))
MIN_CONTIG_LEN = int(config.get("min_contig_len", 200))
THREADS = int(config.get("runtime_scan_threads", 8))

DRAM_CONFIG = config.get(
    "dram_config",
    "/home/jys0914/.conda/envs/DRAM/lib/python3.10/site-packages/mag_annotator/CONFIG",
)
DRAM_PY = config.get("dram_py", "/home/jys0914/.conda/envs/DRAM/bin/DRAM.py")
USE_CAMPER = config.get("dram_use_camper", True)

CONTIGS = config.get(
    "runtime_scan_contigs",
    f"assembly/megahit/{RUN}_{DEPTH}_seed{SEED}/final.contigs.fa",
)
OUTDIR = f"annotation/dram_contigs_minlen/{RUN}_{DEPTH}_seed{SEED}_min{MIN_CONTIG_LEN}"
LOGDIR = f"logs/dram_contigs_minlen/min{MIN_CONTIG_LEN}"


rule all:
    input:
        f"{OUTDIR}/annotations.tsv",
        f"{OUTDIR}/distillate/genome_stats.tsv",


rule dram_annotate_contigs_minlen:
    input:
        contigs=CONTIGS,
    output:
        annotations=f"{OUTDIR}/annotations.tsv",
        trnas=f"{OUTDIR}/trnas.tsv",
        rrnas=f"{OUTDIR}/rrnas.tsv",
    params:
        outdir=OUTDIR,
        camper_flag="--use_camper" if USE_CAMPER else "",
    log:
        f"{LOGDIR}/annotate.log",
    threads: THREADS
    shell:
        "mkdir -p {LOGDIR} && "
        "rm -rf {params.outdir} && "
        "{DRAM_PY} annotate "
        "-i {input.contigs} "
        "-o {params.outdir} "
        "--min_contig_size {MIN_CONTIG_LEN} "
        "--threads {threads} "
        "{params.camper_flag} "
        "--verbose "
        ">> {log} 2>&1"


rule fix_dram_contig_inputs_for_distill_minlen:
    input:
        annotations=f"{OUTDIR}/annotations.tsv",
        trnas=f"{OUTDIR}/trnas.tsv",
        rrnas=f"{OUTDIR}/rrnas.tsv",
    output:
        annotations_fixed=f"{OUTDIR}/distill_input_fixed/annotations.tsv",
        trnas_fixed=f"{OUTDIR}/distill_input_fixed/trnas.tsv",
        rrnas_fixed=f"{OUTDIR}/distill_input_fixed/rrnas.tsv",
    log:
        f"{LOGDIR}/fix_distill_input.log",
    shell:
        "python scripts/normalize_dram_fasta_ids.py {input.annotations} {output.annotations_fixed} >> {log} 2>&1 && "
        "python scripts/normalize_dram_fasta_ids.py {input.trnas} {output.trnas_fixed} >> {log} 2>&1 && "
        "python scripts/normalize_dram_fasta_ids.py {input.rrnas} {output.rrnas_fixed} >> {log} 2>&1"


rule prepare_contig_distill_config_with_camper_minlen:
    input:
        annotations_fixed=f"{OUTDIR}/distill_input_fixed/annotations.tsv",
    output:
        config_loc=f"{OUTDIR}/distill_input_fixed/distill_with_camper.config.json",
        heatmap_form=f"{OUTDIR}/distill_input_fixed/function_heatmap_with_camper.tsv",
    log:
        f"{LOGDIR}/prepare_distill_config.log",
    shell:
        "python scripts/prepare_distill_config_with_camper.py "
        "--base-config {DRAM_CONFIG} "
        "--annotations {input.annotations_fixed} "
        "--out-config {output.config_loc} "
        "--out-heatmap {output.heatmap_form} "
        ">> {log} 2>&1"


rule dram_contig_distill_minlen:
    input:
        annotations=f"{OUTDIR}/distill_input_fixed/annotations.tsv",
        trnas=f"{OUTDIR}/distill_input_fixed/trnas.tsv",
        rrnas=f"{OUTDIR}/distill_input_fixed/rrnas.tsv",
        config_loc=f"{OUTDIR}/distill_input_fixed/distill_with_camper.config.json",
    output:
        genome_stats=f"{OUTDIR}/distillate/genome_stats.tsv",
    params:
        outdir=f"{OUTDIR}/distillate",
    log:
        f"{LOGDIR}/distill.log",
    shell:
        "rm -rf {params.outdir} && "
        "{DRAM_PY} distill "
        "-i {input.annotations} "
        "-o {params.outdir} "
        "--trna_path {input.trnas} "
        "--rrna_path {input.rrnas} "
        "--config_loc {input.config_loc} "
        ">> {log} 2>&1"
