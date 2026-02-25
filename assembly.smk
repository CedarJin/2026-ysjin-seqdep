# Snakemake workflow:
# host removal -> assembly -> read remapping -> binning -> bin QC -> taxonomy

configfile: "config.yaml"

RUNS = config.get("runs", ["SRR10692699"])
SEEDS = [str(s) for s in config.get("downsample_seeds", [11, 22, 33])]
DEPTH_LABELS = config.get("downsample_depths", ["10M", "20M", "30M", "40M", "50M"])
THREADS = config.get("threads", 8)
HOST_REFERENCE_FASTA = config.get(
    "host_reference_fasta",
    "reference/human/GCF_000001405.40_GRCh38.p14_genomic.fna",
)
HOST_INDEX_PREFIX = config.get(
    "host_index_prefix",
    "reference/human/GCF_000001405.40_GRCh38.p14_genomic",
)
CHECKM2_DB_PATH = config.get("checkm2_database_path", "reference/CheckM2_database/uniref100.KO.1.dmnd")
METABAT_MIN_CONTIG = config.get("metabat_min_contig", 1500)
CONTIG_INDEX_EXT = ["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"]

HOST_INDEX_FILES = expand(
    "{prefix}.{ext}",
    prefix=HOST_INDEX_PREFIX,
    ext=CONTIG_INDEX_EXT,
)

SAMPLE_CONTIGS = expand(
    "assembly/megahit/{run}_{depth}_seed{seed}/final.contigs.fa",
    run=RUNS,
    depth=DEPTH_LABELS,
    seed=SEEDS,
)
SAMPLE_CHECKM2 = expand(
    "qc/checkm2/{run}_{depth}_seed{seed}/quality_report.tsv",
    run=RUNS,
    depth=DEPTH_LABELS,
    seed=SEEDS,
)
SAMPLE_GTDB_DONE = expand(
    "taxonomy/gtdbtk/{run}_{depth}_seed{seed}/.done",
    run=RUNS,
    depth=DEPTH_LABELS,
    seed=SEEDS,
)

rule all:
    """Final targets for assembly, binning QC, and taxonomy."""
    input:
        SAMPLE_CONTIGS,
        SAMPLE_CHECKM2,
        SAMPLE_GTDB_DONE,


rule build_host_index:
    """Build Bowtie2 index for host read depletion."""
    input:
        ref=HOST_REFERENCE_FASTA,
    output:
        HOST_INDEX_FILES,
    params:
        prefix=HOST_INDEX_PREFIX,
    log:
        "logs/bowtie2/build_host_index.log",
    threads: THREADS
    shell:
        "mkdir -p $(dirname {params.prefix}) logs/bowtie2 && "
        "bowtie2-build --threads {threads} {input.ref} {params.prefix} >> {log} 2>&1"


rule remove_host_reads:
    """Keep paired reads that do not concordantly map to host reference."""
    input:
        r1="trimmed/fastp/{run}_{depth}_seed{seed}_R1.fastq",
        r2="trimmed/fastp/{run}_{depth}_seed{seed}_R2.fastq",
        idx=HOST_INDEX_FILES,
    output:
        r1="host_removed/{run}_{depth}_seed{seed}_R1.fastq",
        r2="host_removed/{run}_{depth}_seed{seed}_R2.fastq",
    params:
        idx_prefix=HOST_INDEX_PREFIX,
    log:
        "logs/bowtie2/{run}_{depth}_seed{seed}.log",
    threads: THREADS
    shell:
        "mkdir -p host_removed logs/bowtie2 && "
        "bowtie2 -x {params.idx_prefix} -1 {input.r1} -2 {input.r2} "
        "--threads {threads} --very-sensitive "
        "--un-conc host_removed/{wildcards.run}_{wildcards.depth}_seed{wildcards.seed}_R%.fastq "
        "-S /dev/null >> {log} 2>&1"


rule megahit:
    """Assemble host-depleted read pairs with default MEGAHIT settings."""
    input:
        r1="host_removed/{run}_{depth}_seed{seed}_R1.fastq",
        r2="host_removed/{run}_{depth}_seed{seed}_R2.fastq",
    output:
        contigs="assembly/megahit/{run}_{depth}_seed{seed}/final.contigs.fa",
    log:
        "logs/megahit/{run}_{depth}_seed{seed}.log",
    threads: THREADS
    shell:
        "megahit -1 {input.r1} -2 {input.r2} "
        "-f "
        "-o $(dirname {output.contigs}) "
        "-t {threads} "
        ">> {log} 2>&1"

rule build_contig_index:
    """Build Bowtie2 index for assembled contigs (used for coverage mapping)."""
    input:
        contigs="assembly/megahit/{run}_{depth}_seed{seed}/final.contigs.fa",
    output:
        expand(
            "assembly/megahit/{{run}}_{{depth}}_seed{{seed}}/contigs.{ext}",
            ext=CONTIG_INDEX_EXT,
        ),
    params:
        prefix="assembly/megahit/{run}_{depth}_seed{seed}/contigs",
    log:
        "logs/bowtie2/contig_index_{run}_{depth}_seed{seed}.log",
    threads: THREADS
    shell:
        "bowtie2-build --threads {threads} {input.contigs} {params.prefix} >> {log} 2>&1"


rule map_reads_to_contigs:
    """Map host-depleted reads back to contigs and generate sorted/indexed BAM."""
    input:
        r1="host_removed/{run}_{depth}_seed{seed}_R1.fastq",
        r2="host_removed/{run}_{depth}_seed{seed}_R2.fastq",
        idx=expand(
            "assembly/megahit/{{run}}_{{depth}}_seed{{seed}}/contigs.{ext}",
            ext=CONTIG_INDEX_EXT,
        ),
    output:
        bam="assembly/alignments/{run}_{depth}_seed{seed}.sorted.bam",
        bai="assembly/alignments/{run}_{depth}_seed{seed}.sorted.bam.bai",
    params:
        idx_prefix="assembly/megahit/{run}_{depth}_seed{seed}/contigs",
    log:
        "logs/bowtie2/map_contigs_{run}_{depth}_seed{seed}.log",
    threads: THREADS
    shell:
        "mkdir -p assembly/alignments logs/bowtie2 && "
        "bowtie2 -x {params.idx_prefix} -1 {input.r1} -2 {input.r2} "
        "--threads {threads} --very-sensitive 2>> {log} | "
        "samtools view -@ {threads} -b - | "
        "samtools sort -@ {threads} -o {output.bam} - >> {log} 2>&1 && "
        "samtools index -@ {threads} {output.bam} {output.bai} >> {log} 2>&1"


rule contig_depth:
    """Compute per-contig depth table needed by MetaBAT2."""
    input:
        bam="assembly/alignments/{run}_{depth}_seed{seed}.sorted.bam",
        bai="assembly/alignments/{run}_{depth}_seed{seed}.sorted.bam.bai",
    output:
        depth="assembly/depth/{run}_{depth}_seed{seed}.depth.txt",
    log:
        "logs/metabat2/depth_{run}_{depth}_seed{seed}.log",
    shell:
        "mkdir -p assembly/depth logs/metabat2 && "
        "jgi_summarize_bam_contig_depths --outputDepth {output.depth} {input.bam} >> {log} 2>&1"


rule metabat2:
    """Generate metagenome bins from contigs and depth profiles."""
    input:
        contigs="assembly/megahit/{run}_{depth}_seed{seed}/final.contigs.fa",
        depth="assembly/depth/{run}_{depth}_seed{seed}.depth.txt",
    output:
        bins_dir=directory("bins/metabat2/{run}_{depth}_seed{seed}"),
    log:
        "logs/metabat2/binning_{run}_{depth}_seed{seed}.log",
    threads: THREADS
    shell:
        "mkdir -p {output.bins_dir} logs/metabat2 && "
        "metabat2 -i {input.contigs} -a {input.depth} -o {output.bins_dir}/bin "
        "-m {METABAT_MIN_CONTIG} -t {threads} >> {log} 2>&1"


rule checkm2:
    """Run CheckM2 quality assessment for generated bins."""
    input:
        bins_dir="bins/metabat2/{run}_{depth}_seed{seed}",
    output:
        report="qc/checkm2/{run}_{depth}_seed{seed}/quality_report.tsv",
    params:
        outdir="qc/checkm2/{run}_{depth}_seed{seed}",
        db_path=CHECKM2_DB_PATH,
    log:
        "logs/checkm2/{run}_{depth}_seed{seed}.log",
    threads: THREADS
    shell:
        "mkdir -p {params.outdir} logs/checkm2 && "
        "if [[ -n '{params.db_path}' ]]; then db_arg=\"--database_path {params.db_path}\"; else db_arg=\"\"; fi && "
        "checkm2 predict --threads {threads} --input {input.bins_dir} -x fa "
        "--output-directory {params.outdir} --force $db_arg >> {log} 2>&1 && "
        "test -s {output.report}"


rule gtdbtk:
    """Classify metagenome bins with GTDB-Tk (requires GTDBTK_DATA_PATH)."""
    input:
        bins_dir="bins/metabat2/{run}_{depth}_seed{seed}",
    output:
        done="taxonomy/gtdbtk/{run}_{depth}_seed{seed}/.done",
    params:
        outdir="taxonomy/gtdbtk/{run}_{depth}_seed{seed}",
    log:
        "logs/gtdbtk/{run}_{depth}_seed{seed}.log",
    threads: THREADS
    shell:
        "mkdir -p {params.outdir} logs/gtdbtk && "
        "if [[ -z \"${{GTDBTK_DATA_PATH:-}}\" ]]; then "
        "echo 'GTDBTK_DATA_PATH is not set. Please configure GTDB-Tk database path.' >> {log}; "
        "exit 1; "
        "fi && "
        "gtdbtk classify_wf --genome_dir {input.bins_dir} --out_dir {params.outdir} "
        "--cpus {threads} --extension fa >> {log} 2>&1 && "
        "touch {output.done}"
