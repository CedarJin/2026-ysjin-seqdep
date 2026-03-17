#!/usr/bin/env python3
"""
Scan read alignment rate across contig minimum-length thresholds.

For each threshold:
1) Filter contigs by length
2) Build Bowtie2 index
3) Align paired reads to filtered contigs
4) Parse Bowtie2 overall alignment rate

Outputs:
- alignment_scan.tsv
- alignment_scan.svg

No third-party Python libraries required.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import subprocess
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Contig min-length alignment scan with Bowtie2.")
    parser.add_argument("--contigs", required=True, help="Input contigs FASTA")
    parser.add_argument("--r1", required=True, help="R1 FASTQ")
    parser.add_argument("--r2", required=True, help="R2 FASTQ")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument(
        "--thresholds",
        default="200,500,1000,2500",
        help="Comma-separated min contig lengths (default: 200,500,1000,2500)",
    )
    parser.add_argument("--threads", type=int, default=16, help="Bowtie2 threads (default: 16)")
    parser.add_argument("--bowtie2", default="/home/jys0914/.conda/envs/assemble/bin/bowtie2")
    parser.add_argument("--bowtie2-build", default="/home/jys0914/.conda/envs/assemble/bin/bowtie2-build")
    parser.add_argument(
        "--reuse-log-200",
        default="",
        help="Optional existing bowtie2 log for threshold 200 to reuse overall alignment rate",
    )
    parser.add_argument(
        "--skip-existing",
        action="store_true",
        help="Reuse existing filtered fasta/index/log if present",
    )
    return parser.parse_args()


def run(cmd: list[str], log_path: Path | None = None) -> None:
    if log_path is None:
        subprocess.run(cmd, check=True)
        return
    with log_path.open("w") as log:
        subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=True)


def n50_from_lengths(lengths: list[int]) -> int:
    if not lengths:
        return 0
    total = sum(lengths)
    half = total / 2.0
    running = 0
    for length in sorted(lengths, reverse=True):
        running += length
        if running >= half:
            return length
    return 0


def fasta_stats(fa_path: Path) -> tuple[int, int, int, int]:
    lengths: list[int] = []
    with fa_path.open() as fh:
        seq_len = 0
        saw_header = False
        for line in fh:
            if line.startswith(">"):
                if saw_header:
                    lengths.append(seq_len)
                saw_header = True
                seq_len = 0
            else:
                seq_len += len(line.strip())
        if saw_header:
            lengths.append(seq_len)
    n = len(lengths)
    total_bp = sum(lengths)
    max_bp = max(lengths) if lengths else 0
    n50_bp = n50_from_lengths(lengths)
    return n, total_bp, max_bp, n50_bp


def filter_fasta_by_len(in_fa: Path, out_fa: Path, min_len: int) -> tuple[int, int, int, int]:
    out_fa.parent.mkdir(parents=True, exist_ok=True)
    lengths: list[int] = []
    header = None
    seq_lines: list[str] = []
    with in_fa.open() as inp, out_fa.open("w") as out:
        for raw in inp:
            line = raw.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    seq = "".join(seq_lines)
                    if len(seq) >= min_len:
                        out.write(f"{header}\n")
                        if seq:
                            for i in range(0, len(seq), 80):
                                out.write(seq[i : i + 80] + "\n")
                        lengths.append(len(seq))
                header = line
                seq_lines = []
            else:
                seq_lines.append(line.strip())
        if header is not None:
            seq = "".join(seq_lines)
            if len(seq) >= min_len:
                out.write(f"{header}\n")
                if seq:
                    for i in range(0, len(seq), 80):
                        out.write(seq[i : i + 80] + "\n")
                lengths.append(len(seq))
    n = len(lengths)
    total_bp = sum(lengths)
    max_bp = max(lengths) if lengths else 0
    n50_bp = n50_from_lengths(lengths)
    return n, total_bp, max_bp, n50_bp


def parse_overall_alignment_rate(log_path: Path) -> float:
    pat = re.compile(r"([0-9]+(?:\.[0-9]+)?)%\s+overall alignment rate")
    for line in log_path.read_text().splitlines():
        m = pat.search(line)
        if m:
            return float(m.group(1))
    raise RuntimeError(f"Could not find overall alignment rate in log: {log_path}")


def write_svg(results: list[dict], out_svg: Path) -> None:
    width, height = 900, 520
    left, right, top, bottom = 90, 40, 40, 90
    plot_w = width - left - right
    plot_h = height - top - bottom
    xs = [r["min_contig_len"] for r in results]
    ys = [r["overall_alignment_rate_pct"] for r in results]
    x_min, x_max = min(xs), max(xs)
    y_min = 0.0
    y_max = max(100.0, max(ys) + 1.0)

    def x_to_px(x: float) -> float:
        if x_max == x_min:
            return left + plot_w / 2
        return left + (x - x_min) / (x_max - x_min) * plot_w

    def y_to_px(y: float) -> float:
        return top + (y_max - y) / (y_max - y_min) * plot_h

    grid_vals = [0, 20, 40, 60, 80, 100]
    points = [(x_to_px(x), y_to_px(y), x, y) for x, y in zip(xs, ys)]
    polyline = " ".join(f"{px:.1f},{py:.1f}" for px, py, _, _ in points)

    lines = []
    lines.append(f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">')
    lines.append('<rect width="100%" height="100%" fill="white"/>')
    lines.append('<text x="450" y="24" text-anchor="middle" font-size="18" font-family="Arial">Read Alignment Rate vs Min Contig Length</text>')

    # Grid and Y ticks
    for gv in grid_vals:
        py = y_to_px(gv)
        lines.append(f'<line x1="{left}" y1="{py:.1f}" x2="{width-right}" y2="{py:.1f}" stroke="#e5e7eb" stroke-width="1"/>')
        lines.append(f'<text x="{left-10}" y="{py+4:.1f}" text-anchor="end" font-size="12" font-family="Arial">{gv}</text>')

    # Axes
    lines.append(f'<line x1="{left}" y1="{top}" x2="{left}" y2="{height-bottom}" stroke="black" stroke-width="1.5"/>')
    lines.append(f'<line x1="{left}" y1="{height-bottom}" x2="{width-right}" y2="{height-bottom}" stroke="black" stroke-width="1.5"/>')

    # X ticks
    for x in xs:
        px = x_to_px(x)
        lines.append(f'<line x1="{px:.1f}" y1="{height-bottom}" x2="{px:.1f}" y2="{height-bottom+6}" stroke="black"/>')
        lines.append(f'<text x="{px:.1f}" y="{height-bottom+24}" text-anchor="middle" font-size="12" font-family="Arial">{x}</text>')

    # Axis labels
    lines.append(f'<text x="{left + plot_w/2:.1f}" y="{height-28}" text-anchor="middle" font-size="13" font-family="Arial">Min Contig Length (bp)</text>')
    lines.append(
        f'<text x="22" y="{top + plot_h/2:.1f}" transform="rotate(-90 22,{top + plot_h/2:.1f})" text-anchor="middle" font-size="13" font-family="Arial">Overall Alignment Rate (%)</text>'
    )

    # Data line
    lines.append(f'<polyline points="{polyline}" fill="none" stroke="#1d4ed8" stroke-width="2.5"/>')
    for px, py, x, y in points:
        lines.append(f'<circle cx="{px:.1f}" cy="{py:.1f}" r="4.5" fill="#1d4ed8"/>')
        lines.append(f'<text x="{px:.1f}" y="{py-10:.1f}" text-anchor="middle" font-size="11" font-family="Arial">{y:.2f}%</text>')

    lines.append("</svg>")
    out_svg.write_text("\n".join(lines))


def write_metric_svg(
    results: list[dict],
    out_svg: Path,
    y_key: str,
    title: str,
    y_label: str,
    line_color: str = "#0f766e",
) -> None:
    width, height = 900, 520
    left, right, top, bottom = 90, 40, 40, 90
    plot_w = width - left - right
    plot_h = height - top - bottom
    xs = [r["min_contig_len"] for r in results]
    ys = [float(r[y_key]) for r in results]
    x_min, x_max = min(xs), max(xs)
    y_min = 0.0
    y_max = max(ys) if ys else 1.0
    if y_max <= y_min:
        y_max = y_min + 1.0

    def x_to_px(x: float) -> float:
        if x_max == x_min:
            return left + plot_w / 2
        return left + (x - x_min) / (x_max - x_min) * plot_w

    def y_to_px(y: float) -> float:
        return top + (y_max - y) / (y_max - y_min) * plot_h

    grid_vals = [y_max * i / 5 for i in range(6)]
    points = [(x_to_px(x), y_to_px(y), x, y) for x, y in zip(xs, ys)]
    polyline = " ".join(f"{px:.1f},{py:.1f}" for px, py, _, _ in points)

    lines = []
    lines.append(f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">')
    lines.append('<rect width="100%" height="100%" fill="white"/>')
    lines.append(f'<text x="450" y="24" text-anchor="middle" font-size="18" font-family="Arial">{title}</text>')

    for gv in grid_vals:
        py = y_to_px(gv)
        lines.append(f'<line x1="{left}" y1="{py:.1f}" x2="{width-right}" y2="{py:.1f}" stroke="#e5e7eb" stroke-width="1"/>')
        lines.append(f'<text x="{left-10}" y="{py+4:.1f}" text-anchor="end" font-size="12" font-family="Arial">{gv:.0f}</text>')

    lines.append(f'<line x1="{left}" y1="{top}" x2="{left}" y2="{height-bottom}" stroke="black" stroke-width="1.5"/>')
    lines.append(f'<line x1="{left}" y1="{height-bottom}" x2="{width-right}" y2="{height-bottom}" stroke="black" stroke-width="1.5"/>')

    for x in xs:
        px = x_to_px(x)
        lines.append(f'<line x1="{px:.1f}" y1="{height-bottom}" x2="{px:.1f}" y2="{height-bottom+6}" stroke="black"/>')
        lines.append(f'<text x="{px:.1f}" y="{height-bottom+24}" text-anchor="middle" font-size="12" font-family="Arial">{x}</text>')

    lines.append(f'<text x="{left + plot_w/2:.1f}" y="{height-28}" text-anchor="middle" font-size="13" font-family="Arial">Min Contig Length (bp)</text>')
    lines.append(
        f'<text x="22" y="{top + plot_h/2:.1f}" transform="rotate(-90 22,{top + plot_h/2:.1f})" text-anchor="middle" font-size="13" font-family="Arial">{y_label}</text>'
    )

    lines.append(f'<polyline points="{polyline}" fill="none" stroke="{line_color}" stroke-width="2.5"/>')
    for px, py, _, y in points:
        lines.append(f'<circle cx="{px:.1f}" cy="{py:.1f}" r="4.5" fill="{line_color}"/>')
        lines.append(f'<text x="{px:.1f}" y="{py-10:.1f}" text-anchor="middle" font-size="11" font-family="Arial">{y:.0f}</text>')

    lines.append("</svg>")
    out_svg.write_text("\n".join(lines))


def main() -> None:
    args = parse_args()
    contigs = Path(args.contigs)
    r1 = Path(args.r1)
    r2 = Path(args.r2)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    thresholds = [int(x.strip()) for x in args.thresholds.split(",") if x.strip()]
    thresholds = sorted(set(thresholds))

    results = []

    rate_200_from_log = None
    if args.reuse_log_200:
        reuse_log = Path(args.reuse_log_200)
        if reuse_log.exists():
            rate_200_from_log = parse_overall_alignment_rate(reuse_log)

    for t in thresholds:
        filtered_fa = outdir / f"contigs.min{t}.fa"
        idx_prefix = outdir / f"idx.min{t}"
        align_log = outdir / f"bowtie2.min{t}.log"
        build_log = outdir / f"bowtie2_build.min{t}.log"

        if args.skip_existing and filtered_fa.exists():
            n, total_bp, max_bp, n50_bp = fasta_stats(filtered_fa)
        else:
            n, total_bp, max_bp, n50_bp = filter_fasta_by_len(contigs, filtered_fa, t)

        if n == 0:
            raise RuntimeError(f"No contigs remain after min length {t}")

        if t == 200 and rate_200_from_log is not None:
            rate = rate_200_from_log
        elif args.skip_existing and align_log.exists():
            rate = parse_overall_alignment_rate(align_log)
        else:
            bt2_idx_files = [f"{idx_prefix}.{ext}" for ext in ["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"]]
            if not (args.skip_existing and all(Path(p).exists() for p in bt2_idx_files)):
                run([args.bowtie2_build, str(filtered_fa), str(idx_prefix), "--threads", str(args.threads)], build_log)
            run(
                [
                    args.bowtie2,
                    "-x",
                    str(idx_prefix),
                    "-1",
                    str(r1),
                    "-2",
                    str(r2),
                    "--threads",
                    str(args.threads),
                    "--very-sensitive",
                    "-S",
                    "/dev/null",
                ],
                align_log,
            )
            rate = parse_overall_alignment_rate(align_log)

        results.append(
            {
                "min_contig_len": t,
                "n_contigs": n,
                "assembly_bp": total_bp,
                "avg_bp": (total_bp / n) if n else 0.0,
                "max_bp": max_bp,
                "n50_bp": n50_bp,
                "overall_alignment_rate_pct": rate,
            }
        )

    tsv_path = outdir / "alignment_scan.tsv"
    with tsv_path.open("w", newline="") as out:
        writer = csv.DictWriter(
            out,
            fieldnames=["min_contig_len", "n_contigs", "assembly_bp", "avg_bp", "max_bp", "n50_bp", "overall_alignment_rate_pct"],
            delimiter="\t",
        )
        writer.writeheader()
        for row in results:
            writer.writerow(row)

    svg_path = outdir / "alignment_scan.svg"
    write_svg(results, svg_path)
    write_metric_svg(results, outdir / "contig_count_scan.svg", "n_contigs", "Contig Count vs Min Contig Length", "Contig Count")
    write_metric_svg(results, outdir / "total_bp_scan.svg", "assembly_bp", "Total Assembly bp vs Min Contig Length", "Total bp")
    write_metric_svg(results, outdir / "avg_bp_scan.svg", "avg_bp", "Average Contig Length vs Min Contig Length", "Average bp")
    write_metric_svg(results, outdir / "max_bp_scan.svg", "max_bp", "Max Contig Length vs Min Contig Length", "Max bp")
    write_metric_svg(results, outdir / "n50_bp_scan.svg", "n50_bp", "N50 vs Min Contig Length", "N50 bp")

    print(f"Wrote {tsv_path}")
    print(f"Wrote {svg_path}")


if __name__ == "__main__":
    main()
