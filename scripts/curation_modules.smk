# Snakemake workflow: SAP curated gene-module abundance (STANDALONE).
#
# This workflow is intentionally SEPARATE from
# scripts/annotation_contigs_abundance.smk and never touches its outputs.
# It only *reads* two files the main pipeline already produced per subsample:
#   annotation/dram_contigs_T2T/{omic}/{sample}/{prefix}/distill_input_fixed/annotations.tsv
#   abundance/dram_contigs_T2T/{omic}/{sample}/{prefix}/gene_counts.tsv
#
# and writes the curated-module tables into a dedicated per-subsample subfolder
#   annotation/dram_contigs_T2T/{omic}/{sample}/{prefix}/custom_modules/
# so the DRAM-native SCFA / polyphenol(CAMPER) / CAZy tables are kept as-is
# (you keep BOTH the DRAM version and your curated version of each).
#
# Steps:
#   1) build curated reference TSVs from SAP/gene_curation_*.xlsx  (once)
#   2) per subsample: curation_module_abundance.py -> <name>_per_{feature,module}.tsv
#   3) concatenate *_per_module.tsv into one long table PER id_variant (separate files)
#
# Usage:
#   snakemake -s scripts/curation_modules.smk --cores 8
#   snakemake -s scripts/curation_modules.smk --cores 1 curation_long   # just the summary

import glob
import os

configfile: "env/config.yaml"

ANNOTATION_ROOT = config.get("dram_contig_annotation_root", "annotation/dram_contigs_T2T")
ABUNDANCE_ROOT = config.get("dram_contig_abundance_root", "abundance/dram_contigs_T2T")
DRAM_PYTHON = config.get("dram_python", "/home/jys0914/.conda/envs/dram/bin/python")

SAP_DIR = config.get("sap_dir", "SAP")
CURATION_REFS_DIR = config.get("curation_refs_dir", "scripts/curation_refs")
CURATION_ID_VARIANTS = config.get("curation_id_variants", ["id_only", "id_ec"])
CURATION_NAMES = config.get(
    "curation_names",
    ["carotenoids", "lbp", "lps", "polyphenol", "bile_acid", "scfa48"],
)
# Per-subsample output layout:
#   custom_modules/{id_only|id_ec}/{name}_per_{feature,module}.tsv
CUSTOM_SUBDIR = config.get("curation_subdir", "custom_modules")
# Cross-subsample long table for depth-benchmarking (override via --config curation_summary=...).
#   id_only: .../all_modules_per_module_long_id_only.tsv
#   id_ec:   .../all_modules_per_module_long_id_ec.tsv
CURATION_SUMMARY = config.get(
    "curation_summary",
    f"{ANNOTATION_ROOT}/curation_module_summary/all_modules_per_module_long_id_only.tsv",
)
# Python used only for the xlsx->TSV reference build (needs pandas + openpyxl).
CURATION_BUILD_PY = config.get("curation_build_python", "python3")

# Source workbooks for the reference-build step.
CURATION_XLSX = [
    f"{SAP_DIR}/gene_curation_Carotenoids.xlsx",
    f"{SAP_DIR}/gene_curation_LBP.xlsx",
    f"{SAP_DIR}/gene_curation_LPS_biosynthesis_BRITE_ko01005.xlsx",
    f"{SAP_DIR}/gene_curation_Polyphenol.xlsx",
    f"{SAP_DIR}/gene_curation_SCFA_48_enzyme.xlsx",
    f"{SAP_DIR}/gene_curation_microbial_bile_acid.xlsx",
]

CURATION_REFS = expand(
    f"{CURATION_REFS_DIR}/{{id_variant}}/{{name}}.tsv",
    id_variant=CURATION_ID_VARIANTS,
    name=CURATION_NAMES,
)


wildcard_constraints:
    omic=r"[^/]+",
    sample=r"[^/]+",
    prefix=r"[^/]+",


def discover_subsamples():
    """Subsamples with BOTH a fixed annotations.tsv and gene_counts.tsv.

    Discovering by presence means we never emit a MissingInput error for
    partially processed subsamples, and we skip *.bak.* backup dirs."""
    pattern = os.path.join(
        ANNOTATION_ROOT, "*", "*", "*", "distill_input_fixed", "annotations.tsv"
    )
    out = []
    for ann in glob.glob(pattern):
        parts = ann.split(os.sep)
        # .../{omic}/{sample}/{prefix}/distill_input_fixed/annotations.tsv
        prefix = parts[-3]
        sample = parts[-4]
        omic = parts[-5]
        if ".bak" in prefix:
            continue
        counts = os.path.join(ABUNDANCE_ROOT, omic, sample, prefix, "gene_counts.tsv")
        if not os.path.isfile(counts):
            continue
        out.append((omic, sample, prefix))
    return sorted(set(out))


SUBSAMPLES = discover_subsamples()
_OMICS = [s[0] for s in SUBSAMPLES]
_SAMPLES = [s[1] for s in SUBSAMPLES]
_PREFIXES = [s[2] for s in SUBSAMPLES]

PER_MODULE_ALL = [
    f"{ANNOTATION_ROOT}/{o}/{s}/{p}/{CUSTOM_SUBDIR}/{variant}/{name}_per_module.tsv"
    for (o, s, p) in SUBSAMPLES
    for variant in CURATION_ID_VARIANTS
    for name in CURATION_NAMES
]


rule all:
    input:
        CURATION_REFS,
        CURATION_SUMMARY,


rule build_curation_refs:
    """Build normalized curated reference TSVs from the SAP workbooks."""
    input:
        xlsx=CURATION_XLSX,
    output:
        refs=CURATION_REFS,
    params:
        sap_dir=SAP_DIR,
        out_dir=CURATION_REFS_DIR,
    log:
        "logs/curation_modules/build_curation_refs.log",
    shell:
        "mkdir -p $(dirname {log}) && "
        "{CURATION_BUILD_PY} scripts/build_curation_references.py "
        "--sap-dir {params.sap_dir} "
        "--out-dir {params.out_dir} "
        ">> {log} 2>&1"


rule curation_module_abundance:
    """Per-subsample curated-module abundance (id_only by default)."""
    input:
        counts=f"{ABUNDANCE_ROOT}/{{omic}}/{{sample}}/{{prefix}}/gene_counts.tsv",
        annotations=f"{ANNOTATION_ROOT}/{{omic}}/{{sample}}/{{prefix}}/distill_input_fixed/annotations.tsv",
        refs=CURATION_REFS,
    output:
        per_feature=expand(
            f"{ANNOTATION_ROOT}/{{{{omic}}}}/{{{{sample}}}}/{{{{prefix}}}}/{CUSTOM_SUBDIR}/{{id_variant}}/{{name}}_per_feature.tsv",
            id_variant=CURATION_ID_VARIANTS,
            name=CURATION_NAMES,
        ),
        per_module=expand(
            f"{ANNOTATION_ROOT}/{{{{omic}}}}/{{{{sample}}}}/{{{{prefix}}}}/{CUSTOM_SUBDIR}/{{id_variant}}/{{name}}_per_module.tsv",
            id_variant=CURATION_ID_VARIANTS,
            name=CURATION_NAMES,
        ),
    params:
        out_dir=f"{ANNOTATION_ROOT}/{{omic}}/{{sample}}/{{prefix}}/{CUSTOM_SUBDIR}",
        refs_dir=CURATION_REFS_DIR,
        sample_label="{prefix}",
        id_variants=" ".join(CURATION_ID_VARIANTS),
    log:
        "logs/curation_modules/abundance_{omic}_{sample}_{prefix}.log",
    shell:
        "mkdir -p $(dirname {log}) {params.out_dir} && "
        "for variant in {params.id_variants}; do "
        "  mkdir -p {params.out_dir}/$variant && "
        "  {DRAM_PYTHON} scripts/curation_module_abundance.py "
        "  --counts {input.counts} "
        "  --annotations {input.annotations} "
        "  --curation-refs-dir {params.refs_dir}/$variant "
        "  --out-dir {params.out_dir}/$variant "
        "  --sample-label {params.sample_label}; "
        "done "
        ">> {log} 2>&1"


rule curation_long:
    """Concatenate every subsample's *_per_module.tsv into one long table."""
    input:
        per_module=PER_MODULE_ALL,
    output:
        long=CURATION_SUMMARY,
    log:
        "logs/curation_modules/curation_long.log",
    shell:
        "mkdir -p $(dirname {log}) $(dirname {output.long}) && "
        "{DRAM_PYTHON} scripts/aggregate_curation_long.py "
        "--out {output.long} "
        "{input.per_module} "
        ">> {log} 2>&1"
