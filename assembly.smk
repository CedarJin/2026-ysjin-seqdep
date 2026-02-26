# Snakemake workflow:
# host removal -> assembly -> read remapping
# -> binning (MetaBAT2, MaxBin2, CONCOCT)
# -> bin refinement (DAS_Tool)
# -> bin QC (CheckM2)
# -> taxonomy (GTDB-Tk)

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
# CheckM2 now runs on DAS_Tool-refined bins
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
    """Final targets: assembly, DAS_Tool-refined bin QC, and taxonomy."""
    input:
        SAMPLE_CONTIGS,
        SAMPLE_CHECKM2,
        SAMPLE_GTDB_DONE,


# ---------------------------------------------------------------------------
# Host read removal
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# Assembly
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# Read-to-contig mapping and coverage
# ---------------------------------------------------------------------------

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
    """Compute per-contig depth table needed by MetaBAT2 and MaxBin2."""
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


# ---------------------------------------------------------------------------
# Binning: MetaBAT2, MaxBin2, CONCOCT
# ---------------------------------------------------------------------------

rule metabat2:
    """Bin contigs using MetaBAT2 with depth-based coverage profiles."""
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


rule maxbin2:
    """Bin contigs using MaxBin2 with paired read coverage."""
    input:
        contigs="assembly/megahit/{run}_{depth}_seed{seed}/final.contigs.fa",
        r1="host_removed/{run}_{depth}_seed{seed}_R1.fastq",
        r2="host_removed/{run}_{depth}_seed{seed}_R2.fastq",
    output:
        bins_dir=directory("bins/maxbin2/{run}_{depth}_seed{seed}"),
    log:
        "logs/maxbin2/{run}_{depth}_seed{seed}.log",
    threads: THREADS
    shell:
        "mkdir -p {output.bins_dir} logs/maxbin2 && "
        "run_MaxBin.pl "
        "-contig {input.contigs} "
        "-reads {input.r1} "
        "-reads2 {input.r2} "
        "-out {output.bins_dir}/bin "
        "-thread {threads} "
        ">> {log} 2>&1"


rule concoct:
    """Bin contigs using CONCOCT with coverage-based clustering (5 steps)."""
    input:
        contigs="assembly/megahit/{run}_{depth}_seed{seed}/final.contigs.fa",
        bam="assembly/alignments/{run}_{depth}_seed{seed}.sorted.bam",
        bai="assembly/alignments/{run}_{depth}_seed{seed}.sorted.bam.bai",
    output:
        bins_dir=directory("bins/concoct/{run}_{depth}_seed{seed}"),
    params:
        workdir="bins/concoct/{run}_{depth}_seed{seed}",
    log:
        "logs/concoct/{run}_{depth}_seed{seed}.log",
    threads: THREADS
    shell:
        "mkdir -p {params.workdir} logs/concoct && "
        # Step 1: cut contigs into 10 kb chunks for uniform coverage estimation
        "cut_up_fasta.py {input.contigs} -c 10000 -o 0 --merge_last "
        "-b {params.workdir}/contigs_10k.bed > {params.workdir}/contigs_10k.fa 2>> {log} && "
        # Step 2: compute per-chunk coverage from BAM
        "concoct_coverage_table.py {params.workdir}/contigs_10k.bed {input.bam} "
        "> {params.workdir}/coverage_table.tsv 2>> {log} && "
        # Step 3: cluster contigs
        "concoct --composition_file {params.workdir}/contigs_10k.fa "
        "--coverage_file {params.workdir}/coverage_table.tsv "
        "-b {params.workdir}/concoct_out/ "
        "--threads {threads} >> {log} 2>&1 && "
        # Step 4: merge chunk-level clusters back to original contigs
        "merge_cutup_clustering.py {params.workdir}/concoct_out/clustering_gt1000.csv "
        "> {params.workdir}/clustering_merged.csv 2>> {log} && "
        # Step 5: extract per-bin FASTA files
        "extract_fasta_bins.py {input.contigs} {params.workdir}/clustering_merged.csv "
        "--output_path {output.bins_dir}/ >> {log} 2>&1"


# ---------------------------------------------------------------------------
# Bin refinement: DAS_Tool aggregates MetaBAT2 + MaxBin2 + CONCOCT
# ---------------------------------------------------------------------------

rule dastool:
    """Aggregate and dereplicate bins from all three binners using DAS_Tool."""
    input:
        contigs="assembly/megahit/{run}_{depth}_seed{seed}/final.contigs.fa",
        metabat2_dir="bins/metabat2/{run}_{depth}_seed{seed}",
        maxbin2_dir="bins/maxbin2/{run}_{depth}_seed{seed}",
        concoct_dir="bins/concoct/{run}_{depth}_seed{seed}",
    output:
        # DAS_Tool writes refined bins to <prefix>_DASTool_bins/
        bins_dir=directory("bins/dastool/{run}_{depth}_seed{seed}/dastool_DASTool_bins"),
        done="bins/dastool/{run}_{depth}_seed{seed}/.done",
    params:
        outprefix="bins/dastool/{run}_{depth}_seed{seed}/dastool",
    log:
        "logs/dastool/{run}_{depth}_seed{seed}.log",
    threads: THREADS
    shell:
        "mkdir -p bins/dastool/{wildcards.run}_{wildcards.depth}_seed{wildcards.seed} logs/dastool && "
        # Convert each binner's output directory to a scaffold-to-bin TSV
        # MetaBAT2 outputs .fa, MaxBin2 outputs .fasta, CONCOCT outputs .fa
        "Fasta_to_Contig2Bin.sh -i {input.metabat2_dir} -e fa "
        "> {params.outprefix}_metabat2.tsv 2>> {log} && "
        "Fasta_to_Contig2Bin.sh -i {input.maxbin2_dir} -e fasta "  
        "> {params.outprefix}_maxbin2.tsv 2>> {log} && "
        "Fasta_to_Contig2Bin.sh -i {input.concoct_dir} -e fa "
        "> {params.outprefix}_concoct.tsv 2>> {log} && "
        # Run DAS_Tool: select best bins and write refined FASTA files
        "DAS_Tool "
        "-i {params.outprefix}_metabat2.tsv,"
        "{params.outprefix}_maxbin2.tsv,"
        "{params.outprefix}_concoct.tsv "
        "-l metabat2,maxbin2,concoct "
        "-c {input.contigs} "
        "-o {params.outprefix} "
        "--threads {threads} "
        "--write_bins >> {log} 2>&1 && "
        "touch {output.done}"


# ---------------------------------------------------------------------------
# Bin quality and taxonomy on DAS_Tool-refined bins
# ---------------------------------------------------------------------------

rule checkm2:
    """Run CheckM2 quality assessment on DAS_Tool-refined bins."""
    input:
        bins_dir="bins/dastool/{run}_{depth}_seed{seed}/dastool_DASTool_bins",
        done="bins/dastool/{run}_{depth}_seed{seed}/.done",
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
    """Classify DAS_Tool-refined bins with GTDB-Tk (requires GTDBTK_DATA_PATH)."""
    input:
        bins_dir="bins/dastool/{run}_{depth}_seed{seed}/dastool_DASTool_bins",
        done="bins/dastool/{run}_{depth}_seed{seed}/.done",
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
