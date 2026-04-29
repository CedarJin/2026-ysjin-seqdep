# Downsample paired reads in rawdata/metaG and rawdata/metaT.
# Config: env/config.yaml

import glob
import os

configfile: "env/config.yaml"

OMICS = config.get("omics", ["metaG", "metaT"])
SEEDS = [str(s) for s in config.get("downsample_seeds", [11, 22, 33, 44, 55])]
DEPTHS = config.get("downsample_depths", ["10M", "20M", "30M", "40M", "50M"])
THREADS = int(config.get("downsample_threads", 1))


def depth_to_reads(depth):
    s = str(depth).strip().upper()
    if s.endswith("M"):
        return int(float(s[:-1]) * 1_000_000)
    return int(s)


def list_samples(omic):
    r1 = glob.glob(f"rawdata/{omic}/*_R1.fastq")
    return sorted(os.path.basename(p).replace("_R1.fastq", "") for p in r1)


SAMPLES = {omic: list_samples(omic) for omic in OMICS}
ALL_OUTPUTS = [
    f"rawdata/downsample/{omic}/{sample}/{sample}_{depth}_seed{seed}_R{read}.fastq"
    for omic in OMICS
    for sample in SAMPLES[omic]
    for depth in DEPTHS
    for seed in SEEDS
    for read in ["1", "2"]
]


rule all:
    input:
        ALL_OUTPUTS,


rule downsample:
    input:
        r1="rawdata/{omic}/{sample}_R1.fastq",
        r2="rawdata/{omic}/{sample}_R2.fastq",
    output:
        r1="rawdata/downsample/{omic}/{sample}/{sample}_{depth}_seed{seed}_R1.fastq",
        r2="rawdata/downsample/{omic}/{sample}/{sample}_{depth}_seed{seed}_R2.fastq",
    params:
        n=lambda wc: depth_to_reads(wc.depth),
    threads:
        THREADS
    shell:
        """
        mkdir -p $(dirname {output.r1})
        seqtk sample -s{wildcards.seed} {input.r1} {params.n} > {output.r1}
        seqtk sample -s{wildcards.seed} {input.r2} {params.n} > {output.r2}
        """
