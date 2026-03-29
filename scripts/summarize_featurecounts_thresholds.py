#!/usr/bin/env python3
"""Combine per-threshold featureCounts abundance tables into summary outputs."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


MATRIX_METADATA_FIELDS = [
    "scaffold",
    "gene_position",
    "start_position",
    "end_position",
    "strandedness",
    "rank",
    "ko_id",
    "kegg_hit",
    "camper_id",
    "camper_definition",
    "peptidase_id",
    "cazy_ids",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tables", nargs="+", required=True, help="Per-threshold gene_abundance.tsv files")
    parser.add_argument("--summaries", nargs="+", required=True, help="Per-threshold featureCounts .summary files")
    parser.add_argument("--output-long", required=True, help="Concatenated long-format abundance table")
    parser.add_argument("--output-matrix", required=True, help="Wide abundance matrix across thresholds")
    parser.add_argument("--output-summary", required=True, help="Compact threshold summary table")
    return parser.parse_args()


def read_table(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader)


def parse_minlen_from_table(rows: list[dict[str, str]], path: Path) -> int:
    if not rows:
        raise ValueError(f"No rows found in {path}")
    return int(rows[0]["min_contig_len"])


def parse_summary(path: Path) -> tuple[str, dict[str, int]]:
    with path.open() as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        sample_column = header[-1]
        counts: dict[str, int] = {}
        for row in reader:
            if not row:
                continue
            counts[row[0]] = int(row[-1])
    return sample_column, counts


def write_long_table(rows: list[dict[str, str]], path: Path) -> None:
    fieldnames = list(rows[0].keys())
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_matrix(rows_by_gene: dict[str, dict[str, str]], thresholds: list[int], path: Path) -> None:
    fieldnames = ["gene_id"] + MATRIX_METADATA_FIELDS + [f"assigned_fragments_min{thr}" for thr in thresholds]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for gene_id in sorted(rows_by_gene):
            writer.writerow(rows_by_gene[gene_id])


def write_threshold_summary(rows: list[dict[str, str]], path: Path) -> None:
    fieldnames = [
        "min_contig_len",
        "total_fragments",
        "assigned_fragments",
        "assigned_pct",
        "genes_total",
        "genes_with_counts",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    long_rows: list[dict[str, str]] = []
    matrix_rows: dict[str, dict[str, str]] = {}
    threshold_summary_rows: list[dict[str, str]] = []

    table_info = []
    for table_path in args.tables:
        rows = read_table(Path(table_path))
        minlen = parse_minlen_from_table(rows, Path(table_path))
        table_info.append((minlen, Path(table_path), rows))

    summary_info = {}
    for summary_path in args.summaries:
        sample_column, counts = parse_summary(Path(summary_path))
        summary_info[Path(summary_path).parent.name] = (sample_column, counts)

    thresholds = sorted(minlen for minlen, _, _ in table_info)

    for minlen, table_path, rows in sorted(table_info, key=lambda x: x[0]):
        long_rows.extend(rows)

        genes_total = len(rows)
        genes_with_counts = sum(float(row.get("assigned_fragments", "0") or 0) > 0 for row in rows)

        summary_key = f"min{minlen}"
        sample_column, counts = summary_info[summary_key]
        assigned = counts.get("Assigned", 0)
        total = sum(counts.values())
        assigned_pct = 0.0 if total == 0 else (assigned / total) * 100.0

        threshold_summary_rows.append(
            {
                "min_contig_len": str(minlen),
                "total_fragments": str(total),
                "assigned_fragments": str(assigned),
                "assigned_pct": f"{assigned_pct:.2f}",
                "genes_total": str(genes_total),
                "genes_with_counts": str(genes_with_counts),
            }
        )

        for row in rows:
            gene_id = row["gene_id"]
            if gene_id not in matrix_rows:
                matrix_rows[gene_id] = {"gene_id": gene_id}
                for field in MATRIX_METADATA_FIELDS:
                    matrix_rows[gene_id][field] = row.get(field, "")

            matrix_rows[gene_id][f"assigned_fragments_min{minlen}"] = row.get("assigned_fragments", "0")

    for gene_row in matrix_rows.values():
        for minlen in thresholds:
            gene_row.setdefault(f"assigned_fragments_min{minlen}", "0")

    out_long = Path(args.output_long)
    out_matrix = Path(args.output_matrix)
    out_summary = Path(args.output_summary)
    out_long.parent.mkdir(parents=True, exist_ok=True)

    write_long_table(long_rows, out_long)
    write_matrix(matrix_rows, thresholds, out_matrix)
    write_threshold_summary(threshold_summary_rows, out_summary)

    print(f"Wrote combined long-format abundance table to {out_long}")
    print(f"Wrote abundance matrix to {out_matrix}")
    print(f"Wrote min-contig threshold summary to {out_summary}")


if __name__ == "__main__":
    main()
