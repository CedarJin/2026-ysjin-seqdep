# 2026-ysjin-seqdep
Microbiome metagenome and metatranscriptome sequencing depth for real data

## Real metagenomic data

### Sample metadata

Raw FASTQ files are downloaded and renamed using sequential IDs (`MG0001`–`MG0008`).
The mapping between these IDs and the customer labels is in `rawdata/metaG/metaG_metadata.tsv`:

| sample_id | customer_label           | internal_id | R1                    | R2                    |
|-----------|--------------------------|-------------|-----------------------|-----------------------|
| MG0001    | GB2006.Day0.metagenome   | zr27836_1   | MG0001_R1.fastq.gz    | MG0001_R2.fastq.gz    |
| MG0002    | GB2006.Day180.metagenome | zr27836_2   | MG0002_R1.fastq.gz    | MG0002_R2.fastq.gz    |
| MG0003    | GB2033.Day0.metagenome   | zr27836_3   | MG0003_R1.fastq.gz    | MG0003_R2.fastq.gz    |
| MG0004    | GB2033.Day180.metagenome | zr27836_4   | MG0004_R1.fastq.gz    | MG0004_R2.fastq.gz    |
| MG0005    | GB2003.Day0.metagenome   | zr27836_5   | MG0005_R1.fastq.gz    | MG0005_R2.fastq.gz    |
| MG0006    | GB2003.Day180.metagenome | zr27836_6   | MG0006_R1.fastq.gz    | MG0006_R2.fastq.gz    |
| MG0007    | GB2032.Day0.metagenome   | zr27836_7   | MG0007_R1.fastq.gz    | MG0007_R2.fastq.gz    |
| MG0008    | GB2032.Day180.metagenome | zr27836_8   | MG0008_R1.fastq.gz    | MG0008_R2.fastq.gz    |

Download URLs are stored in `rawdata/metaG/metaG_rawdata_links.tsv`.

### Download raw FASTQ files

```bash
bash scripts/download_metaG_rawdata.sh            # saves to rawdata/metaG/ (default)
```

The script reads `metaG_rawdata_links.tsv`, downloads each pair of reads, and saves them as
`MG{NNNN}_R1.fastq.gz` / `MG{NNNN}_R2.fastq.gz` in the output directory.
Files that already exist are skipped automatically (safe to re-run).

## Real metatranscriptomic data

### Sample metadata

Raw FASTQ files are downloaded and renamed using sequential IDs (`MT0001`–`MT0008`).
The mapping between these IDs and the customer labels is in `rawdata/metaT/metaT_metadata.tsv`:

| sample_id | customer_label                 | internal_id  | R1                    | R2                    |
|-----------|-------------------------------|--------------|-----------------------|-----------------------|
| MT0001    | GB2006.Day0.metatranscriptome  | zr27836_9T   | MT0001_R1.fastq.gz    | MT0001_R2.fastq.gz    |
| MT0002    | GB2006.Day180.metatranscriptome| zr27836_10T  | MT0002_R1.fastq.gz    | MT0002_R2.fastq.gz    |
| MT0003    | GB2033.Day0.metatranscriptome  | zr27836_11T  | MT0003_R1.fastq.gz    | MT0003_R2.fastq.gz    |
| MT0004    | GB2033.Day180.metatranscriptome| zr27836_12T  | MT0004_R1.fastq.gz    | MT0004_R2.fastq.gz    |
| MT0005    | GB2003.Day0.metatranscriptome  | zr27836_13T  | MT0005_R1.fastq.gz    | MT0005_R2.fastq.gz    |
| MT0006    | GB2003.Day180.metatranscriptome| zr27836_14T  | MT0006_R1.fastq.gz    | MT0006_R2.fastq.gz    |
| MT0007    | GB2032.Day0.metatranscriptome  | zr27836_15T  | MT0007_R1.fastq.gz    | MT0007_R2.fastq.gz    |
| MT0008    | GB2032.Day180.metatranscriptome| zr27836_16T  | MT0008_R1.fastq.gz    | MT0008_R2.fastq.gz    |

Download URLs are stored in `rawdata/metaT/metaT_rawdata_links.tsv`.

### Download raw FASTQ files

```bash
bash scripts/download_metaT_rawdata.sh              # saves to rawdata/metaT/ (default)
bash scripts/download_metaT_rawdata.sh /your/path   # custom output directory
```

The script reads `rawdata/metaT/metaT_rawdata_links.tsv`, downloads each pair of reads, and saves them as
`MT{NNNN}_R1.fastq.gz` / `MT{NNNN}_R2.fastq.gz` in the output directory.
Files that already exist are skipped automatically (safe to re-run).

---
