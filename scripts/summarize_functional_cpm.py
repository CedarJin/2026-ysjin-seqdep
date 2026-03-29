#!/usr/bin/env python3
"""Summarize CAZy and CAMPER CPM across contig-length thresholds."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


TOTAL_LABEL = "__TOTAL__"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tables", nargs="+", required=True, help="Per-threshold gene_abundance.tsv files")
    parser.add_argument("--cazy-output", required=True, help="Output CAZy CPM table")
    parser.add_argument("--camper-output", required=True, help="Output CAMPER CPM table")
    parser.add_argument("--cazy-dist-output", required=True, help="Output CAZy CPM distribution summary")
    parser.add_argument("--camper-dist-output", required=True, help="Output CAMPER CPM distribution summary")
    parser.add_argument("--cazy-plot-output", required=True, help="Output CAZy CPM distribution SVG")
    parser.add_argument("--camper-plot-output", required=True, help="Output CAMPER CPM distribution SVG")
    return parser.parse_args()


def read_table(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader)


def split_hits(raw_value: str) -> list[str]:
    if not raw_value:
        return []
    hits = []
    seen = set()
    for part in raw_value.split(";"):
        token = part.strip()
        if token and token not in seen:
            hits.append(token)
            seen.add(token)
    return hits


def percentile(sorted_values: list[float], p: float) -> float:
    if not sorted_values:
        return 0.0
    if len(sorted_values) == 1:
        return sorted_values[0]
    pos = (len(sorted_values) - 1) * p
    lower = math.floor(pos)
    upper = math.ceil(pos)
    if lower == upper:
        return sorted_values[lower]
    frac = pos - lower
    return sorted_values[lower] * (1.0 - frac) + sorted_values[upper] * frac


def summarize_feature(rows: list[dict[str, str]], minlen: int, field: str) -> tuple[dict[str, float], float]:
    denominator = 0.0
    counts_by_feature: dict[str, float] = {}

    for row in rows:
        hits = split_hits(row.get(field, ""))
        if not hits:
            continue

        assigned = float(row.get("assigned_fragments", "0") or 0.0)
        denominator += assigned
        for hit in hits:
            counts_by_feature[hit] = counts_by_feature.get(hit, 0.0) + assigned

    cpm_by_feature = {}
    total_count = 0.0
    for hit, count in counts_by_feature.items():
        total_count += count
        cpm_by_feature[hit] = 0.0 if denominator == 0.0 else (count / denominator) * 1_000_000.0

    cpm_by_feature[TOTAL_LABEL] = 0.0 if denominator == 0.0 else (total_count / denominator) * 1_000_000.0
    return cpm_by_feature, denominator


def write_cpm_table(path: Path, label_name: str, thresholds: list[int], values: dict[str, dict[int, float]]) -> None:
    fieldnames = [label_name] + [f"min{threshold}_cpm" for threshold in thresholds]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        labels = sorted(x for x in values if x != TOTAL_LABEL) + ([TOTAL_LABEL] if TOTAL_LABEL in values else [])
        for label in labels:
            row = {label_name: label}
            for threshold in thresholds:
                row[f"min{threshold}_cpm"] = f"{values[label].get(threshold, 0.0):.6f}"
            writer.writerow(row)


def build_distribution_rows(
    thresholds: list[int],
    values: dict[str, dict[int, float]],
    denominators: dict[int, float],
) -> list[dict[str, str]]:
    rows = []
    labels = [label for label in values if label != TOTAL_LABEL]
    for threshold in thresholds:
        cpm_values = sorted(values[label].get(threshold, 0.0) for label in labels)
        rows.append(
            {
                "min_contig_len": str(threshold),
                "denominator_assigned_fragments": f"{denominators.get(threshold, 0.0):.6f}",
                "n_features": str(len(cpm_values)),
                "min_cpm": f"{(cpm_values[0] if cpm_values else 0.0):.6f}",
                "q25_cpm": f"{percentile(cpm_values, 0.25):.6f}",
                "median_cpm": f"{percentile(cpm_values, 0.50):.6f}",
                "q75_cpm": f"{percentile(cpm_values, 0.75):.6f}",
                "max_cpm": f"{(cpm_values[-1] if cpm_values else 0.0):.6f}",
            }
        )
    return rows


def write_distribution_table(path: Path, rows: list[dict[str, str]]) -> None:
    fieldnames = [
        "min_contig_len",
        "denominator_assigned_fragments",
        "n_features",
        "min_cpm",
        "q25_cpm",
        "median_cpm",
        "q75_cpm",
        "max_cpm",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_quantile_trend_svg(path: Path, rows: list[dict[str, str]], title: str) -> None:
    width, height = 900, 520
    left, right, top, bottom = 90, 40, 40, 90
    plot_w = width - left - right
    plot_h = height - top - bottom
    thresholds = [int(row["min_contig_len"]) for row in rows]
    q25_values = [float(row["q25_cpm"]) for row in rows]
    median_values = [float(row["median_cpm"]) for row in rows]
    q75_values = [float(row["q75_cpm"]) for row in rows]
    max_y = max(q75_values + median_values + q25_values) if rows else 1.0
    if max_y <= 0:
        max_y = 1.0

    def x_to_px(value: int) -> float:
        if len(thresholds) == 1 or min(thresholds) == max(thresholds):
            return left + plot_w / 2
        return left + (value - min(thresholds)) / (max(thresholds) - min(thresholds)) * plot_w

    def y_to_px(value: float) -> float:
        return top + (max_y - value) / max_y * plot_h

    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="{width/2:.1f}" y="24" text-anchor="middle" font-size="18" font-family="Arial">{title}</text>',
    ]

    for i in range(6):
        value = max_y * i / 5.0
        py = y_to_px(value)
        lines.append(f'<line x1="{left}" y1="{py:.1f}" x2="{width-right}" y2="{py:.1f}" stroke="#e5e7eb" stroke-width="1"/>')
        lines.append(f'<text x="{left-10}" y="{py+4:.1f}" text-anchor="end" font-size="12" font-family="Arial">{value:.0f}</text>')

    lines.append(f'<line x1="{left}" y1="{top}" x2="{left}" y2="{height-bottom}" stroke="black" stroke-width="1.5"/>')
    lines.append(f'<line x1="{left}" y1="{height-bottom}" x2="{width-right}" y2="{height-bottom}" stroke="black" stroke-width="1.5"/>')

    series = [
        ("25% quantile", q25_values, "#0f766e"),
        ("Median", median_values, "#1d4ed8"),
        ("75% quantile", q75_values, "#b45309"),
    ]

    for _, values, color in series:
        points = [(x_to_px(threshold), y_to_px(value)) for threshold, value in zip(thresholds, values)]
        polyline = " ".join(f"{x:.1f},{y:.1f}" for x, y in points)
        lines.append(f'<polyline points="{polyline}" fill="none" stroke="{color}" stroke-width="2.5"/>')
        for (x, y), value in zip(points, values):
            lines.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="4.5" fill="{color}"/>')
            lines.append(
                f'<text x="{x:.1f}" y="{max(y - 10, top + 12):.1f}" text-anchor="middle" font-size="11" font-family="Arial" fill="{color}">{value:.1f}</text>'
            )

    for row in rows:
        x = x_to_px(int(row["min_contig_len"]))
        lines.append(f'<text x="{x:.1f}" y="{height-bottom+24}" text-anchor="middle" font-size="12" font-family="Arial">{row["min_contig_len"]}</text>')
        lines.append(f'<line x1="{x:.1f}" y1="{height-bottom}" x2="{x:.1f}" y2="{height-bottom+6}" stroke="black"/>')

    lines.append(f'<text x="{left + plot_w/2:.1f}" y="{height-28}" text-anchor="middle" font-size="13" font-family="Arial">Min Contig Length (bp)</text>')
    lines.append(
        f'<text x="22" y="{top + plot_h/2:.1f}" transform="rotate(-90 22,{top + plot_h/2:.1f})" text-anchor="middle" font-size="13" font-family="Arial">Counts Per Million (CPM)</text>'
    )

    legend_x = width - right - 180
    legend_y = top + 10
    for idx, (label, _, color) in enumerate(series):
        y = legend_y + idx * 24
        lines.append(f'<line x1="{legend_x}" y1="{y}" x2="{legend_x+28}" y2="{y}" stroke="{color}" stroke-width="2.5"/>')
        lines.append(f'<circle cx="{legend_x+14}" cy="{y}" r="4.5" fill="{color}"/>')
        lines.append(f'<text x="{legend_x+38}" y="{y+4}" font-size="12" font-family="Arial">{label}</text>')

    lines.append("</svg>")
    path.write_text("\n".join(lines))


def collect_feature_tables(
    table_paths: list[str],
    field: str,
) -> tuple[list[int], dict[str, dict[int, float]], dict[int, float]]:
    thresholds = []
    values: dict[str, dict[int, float]] = {}
    denominators: dict[int, float] = {}

    for table_path in table_paths:
        rows = read_table(Path(table_path))
        if not rows:
            continue
        minlen = int(rows[0]["min_contig_len"])
        thresholds.append(minlen)
        cpm_by_feature, denominator = summarize_feature(rows, minlen, field)
        denominators[minlen] = denominator

        for feature, cpm in cpm_by_feature.items():
            values.setdefault(feature, {})
            values[feature][minlen] = cpm

    thresholds = sorted(set(thresholds))
    for feature in values.values():
        for threshold in thresholds:
            feature.setdefault(threshold, 0.0)
    return thresholds, values, denominators


def main() -> None:
    args = parse_args()

    cazy_thresholds, cazy_values, cazy_denominators = collect_feature_tables(args.tables, "cazy_ids")
    camper_thresholds, camper_values, camper_denominators = collect_feature_tables(args.tables, "camper_id")

    cazy_dist_rows = build_distribution_rows(cazy_thresholds, cazy_values, cazy_denominators)
    camper_dist_rows = build_distribution_rows(camper_thresholds, camper_values, camper_denominators)

    cazy_output = Path(args.cazy_output)
    cazy_output.parent.mkdir(parents=True, exist_ok=True)

    write_cpm_table(cazy_output, "cazy_id", cazy_thresholds, cazy_values)
    write_cpm_table(Path(args.camper_output), "camper_id", camper_thresholds, camper_values)
    write_distribution_table(Path(args.cazy_dist_output), cazy_dist_rows)
    write_distribution_table(Path(args.camper_dist_output), camper_dist_rows)
    write_quantile_trend_svg(Path(args.cazy_plot_output), cazy_dist_rows, "CAZy CPM Quantiles Across Min Contig Lengths")
    write_quantile_trend_svg(Path(args.camper_plot_output), camper_dist_rows, "CAMPER CPM Quantiles Across Min Contig Lengths")

    print(f"Wrote CAZy CPM table to {args.cazy_output}")
    print(f"Wrote CAMPER CPM table to {args.camper_output}")
    print(f"Wrote CAZy CPM distribution summary to {args.cazy_dist_output}")
    print(f"Wrote CAMPER CPM distribution summary to {args.camper_dist_output}")


if __name__ == "__main__":
    main()
