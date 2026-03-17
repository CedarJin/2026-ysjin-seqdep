#!/usr/bin/env python3
"""
Summarize DRAM contig-minlen scans and draw simple SVG plots.

Metrics:
- raw CAZy/CAMPER hit counts from annotations.tsv
- unique CAZy/CAMPER ID counts from annotations.tsv
- runtime from selected SLURM stdout files
- functional module count from metabolism_summary.xlsx
  defined as the number of non-zero rows across functional sheets
  (all sheets except MISC, rRNA, tRNA)
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import zipfile
from pathlib import Path
import xml.etree.ElementTree as ET


NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
SHEET_COLORS = {
    "cazy": "#0f766e",
    "camper": "#b45309",
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--outdir", required=True)
    p.add_argument("--default-annotation", required=True, help="Default DRAM contig annotation dir (used as min2500)")
    p.add_argument(
        "--scan-annotation-root",
        required=True,
        help="Root dir containing dram_contigs_minlen outputs",
    )
    p.add_argument(
        "--runtime",
        action="append",
        default=[],
        help="Runtime mapping in the form minlen=slurm_stdout_path",
    )
    p.add_argument(
        "--minlens",
        default="200,500,1000,2500",
        help="Comma-separated min contig lengths to summarize",
    )
    return p.parse_args()


def parse_runtime_map(items: list[str]) -> dict[int, list[Path]]:
    out: dict[int, list[Path]] = {}
    for item in items:
        key, value = item.split("=", 1)
        parts = [Path(x) for x in value.split(",") if x]
        out[int(key)] = parts
    return out


def parse_walltime_to_hours(text: str) -> float:
    text = text.strip()
    days = 0
    if "-" in text:
        day_part, text = text.split("-", 1)
        days = int(day_part)
    h, m, s = [int(x) for x in text.split(":")]
    return days * 24 + h + m / 60.0 + s / 3600.0


def runtime_hours_from_slurm(paths: list[Path]) -> float:
    from datetime import datetime, timedelta

    walltime_pat = re.compile(r"Used walltime\s*:\s*([0-9-:]+)")
    stamp_pat = re.compile(r"\[([A-Z][a-z]{2} [A-Z][a-z]{2}\s+\d+ \d{2}:\d{2}:\d{2}) (PST|PDT) (\d{4})\]")

    texts = []
    for path in paths:
        if path.exists():
            texts.append((path, path.read_text()))

    for path, text in texts:
        matches = walltime_pat.findall(text)
        if matches:
            return parse_walltime_to_hours(matches[-1])

    def parse_stamp(parts: tuple[str, str, str]) -> datetime:
        body, tz_name, year = parts
        dt = datetime.strptime(f"{body} {year}", "%a %b %d %H:%M:%S %Y")
        offset_hours = -8 if tz_name == "PST" else -7
        return dt - timedelta(hours=offset_hours)

    for path, text in texts:
        stamped = stamp_pat.findall(text)
        if len(stamped) >= 2:
            start = parse_stamp(stamped[0])
            end = parse_stamp(stamped[-1])
            return (end - start).total_seconds() / 3600.0

    joined = ', '.join(str(p) for p in paths)
    raise RuntimeError(f"Could not parse runtime from any of: {joined}")


def annotation_metrics(path: Path) -> dict[str, int]:
    raw_camper = 0
    raw_cazy = 0
    unique_camper: set[str] = set()
    unique_cazy: set[str] = set()
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="	")
        for row in reader:
            camper_hits = (row.get("camper_hits") or "").strip()
            cazy_hits = (row.get("cazy_hits") or "").strip()
            camper_id = (row.get("camper_id") or "").strip()
            cazy_ids = (row.get("cazy_ids") or "").strip()
            if camper_hits:
                raw_camper += 1
            if cazy_hits:
                raw_cazy += 1
            if camper_id:
                unique_camper.update(x.strip() for x in camper_id.split(";") if x.strip())
            if cazy_ids:
                unique_cazy.update(x.strip() for x in cazy_ids.split(";") if x.strip())
    return {
        "raw_camper_hits": raw_camper,
        "raw_cazy_hits": raw_cazy,
        "unique_camper_ids": len(unique_camper),
        "unique_cazy_ids": len(unique_cazy),
    }


def cell_text(cell: ET.Element) -> str:
    cell_type = cell.attrib.get("t")
    if cell_type == "inlineStr":
        return "".join(node.text or "" for node in cell.iter(NS + "t"))
    value = cell.find(NS + "v")
    return "" if value is None else (value.text or "")


def functional_module_count(xlsx_path: Path) -> int:
    excluded = {"MISC", "rRNA", "tRNA"}
    rel_ns = "{http://schemas.openxmlformats.org/package/2006/relationships}"
    col_re = re.compile(r"([A-Z]+)")
    total_nonzero = 0

    with zipfile.ZipFile(xlsx_path) as zf:
        workbook = ET.fromstring(zf.read("xl/workbook.xml"))
        rels = ET.fromstring(zf.read("xl/_rels/workbook.xml.rels"))
        rel_map = {rel.attrib["Id"]: rel.attrib["Target"].lstrip("/") for rel in rels.findall(rel_ns + "Relationship")}

        for sheet in workbook.find(NS + "sheets"):
            name = sheet.attrib["name"]
            if name in excluded:
                continue
            rid = sheet.attrib["{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id"]
            target = rel_map[rid]
            if not target.startswith("xl/"):
                target = "xl/" + target
            root = ET.fromstring(zf.read(target))
            rows = root.find(NS + "sheetData").findall(NS + "row")
            if not rows:
                continue
            header = {
                col_re.match(cell.attrib["r"]).group(1): cell_text(cell)
                for cell in rows[0].findall(NS + "c")
            }
            sample_cols = [col for col, val in header.items() if val == "final.contigs"]
            if not sample_cols:
                continue
            sample_col = sample_cols[0]
            for row in rows[1:]:
                row_map = {
                    col_re.match(cell.attrib["r"]).group(1): cell_text(cell)
                    for cell in row.findall(NS + "c")
                }
                value = row_map.get(sample_col, "")
                if not value:
                    continue
                try:
                    number = float(value)
                except ValueError:
                    continue
                if not math.isnan(number) and number > 0:
                    total_nonzero += 1
    return total_nonzero


def write_tsv(rows: list[dict[str, object]], out_path: Path) -> None:
    fieldnames = [
        "min_contig_len",
        "raw_cazy_hits",
        "raw_camper_hits",
        "unique_cazy_ids",
        "unique_camper_ids",
        "runtime_hours",
        "functional_module_count",
    ]
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def svg_header(width: int, height: int) -> list[str]:
    return [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">',
        '<rect width="100%" height="100%" fill="white"/>',
    ]


def write_single_series_svg(
    rows: list[dict[str, object]],
    x_key: str,
    y_key: str,
    title: str,
    y_label: str,
    color: str,
    out_path: Path,
    integer_labels: bool = False,
) -> None:
    width, height = 920, 540
    left, right, top, bottom = 90, 40, 45, 90
    plot_w = width - left - right
    plot_h = height - top - bottom
    xs = [float(r[x_key]) for r in rows]
    ys = [float(r[y_key]) for r in rows]
    x_min, x_max = min(xs), max(xs)
    y_min = 0.0
    y_max = max(ys) * 1.08 if max(ys) > 0 else 1.0

    def x_to_px(x: float) -> float:
        return left + (x - x_min) / (x_max - x_min) * plot_w if x_max != x_min else left + plot_w / 2

    def y_to_px(y: float) -> float:
        return top + (y_max - y) / (y_max - y_min) * plot_h

    points = [(x_to_px(x), y_to_px(y), x, y) for x, y in zip(xs, ys)]
    polyline = " ".join(f"{px:.1f},{py:.1f}" for px, py, _, _ in points)
    grid_vals = [y_max * i / 5 for i in range(6)]

    lines = svg_header(width, height)
    lines.append(f'<text x="{width/2:.1f}" y="24" text-anchor="middle" font-size="18" font-family="Arial">{title}</text>')
    for gv in grid_vals:
        py = y_to_px(gv)
        lines.append(f'<line x1="{left}" y1="{py:.1f}" x2="{width-right}" y2="{py:.1f}" stroke="#e5e7eb" stroke-width="1"/>')
        tick_label = f"{int(round(gv))}" if integer_labels else f"{gv:.1f}"
        lines.append(f'<text x="{left-10}" y="{py+4:.1f}" text-anchor="end" font-size="12" font-family="Arial">{tick_label}</text>')
    lines.append(f'<line x1="{left}" y1="{top}" x2="{left}" y2="{height-bottom}" stroke="black" stroke-width="1.5"/>')
    lines.append(f'<line x1="{left}" y1="{height-bottom}" x2="{width-right}" y2="{height-bottom}" stroke="black" stroke-width="1.5"/>')
    for x in xs:
        px = x_to_px(x)
        lines.append(f'<line x1="{px:.1f}" y1="{height-bottom}" x2="{px:.1f}" y2="{height-bottom+6}" stroke="black"/>')
        lines.append(f'<text x="{px:.1f}" y="{height-bottom+24}" text-anchor="middle" font-size="12" font-family="Arial">{int(x)}</text>')
    lines.append(f'<text x="{left + plot_w/2:.1f}" y="{height-28}" text-anchor="middle" font-size="13" font-family="Arial">Min Contig Length (bp)</text>')
    lines.append(f'<text x="22" y="{top + plot_h/2:.1f}" transform="rotate(-90 22,{top + plot_h/2:.1f})" text-anchor="middle" font-size="13" font-family="Arial">{y_label}</text>')
    lines.append(f'<polyline points="{polyline}" fill="none" stroke="{color}" stroke-width="2.5"/>')
    for px, py, _, y in points:
        lines.append(f'<circle cx="{px:.1f}" cy="{py:.1f}" r="4.5" fill="{color}"/>')
        point_label = f"{int(round(y))}" if integer_labels else f"{y:.1f}"
        lines.append(f'<text x="{px:.1f}" y="{py-10:.1f}" text-anchor="middle" font-size="11" font-family="Arial">{point_label}</text>')
    lines.append("</svg>")
    out_path.write_text("\n".join(lines))


def write_dual_series_svg(
    rows: list[dict[str, object]],
    x_key: str,
    y1_key: str,
    y2_key: str,
    title: str,
    y_label: str,
    legend1: str,
    legend2: str,
    out_path: Path,
    integer_labels: bool = False,
) -> None:
    width, height = 920, 560
    left, right, top, bottom = 90, 40, 45, 100
    plot_w = width - left - right
    plot_h = height - top - bottom
    xs = [float(r[x_key]) for r in rows]
    y1s = [float(r[y1_key]) for r in rows]
    y2s = [float(r[y2_key]) for r in rows]
    x_min, x_max = min(xs), max(xs)
    y_min = 0.0
    y_max = max(max(y1s), max(y2s)) * 1.08 if max(max(y1s), max(y2s)) > 0 else 1.0

    def x_to_px(x: float) -> float:
        return left + (x - x_min) / (x_max - x_min) * plot_w if x_max != x_min else left + plot_w / 2

    def y_to_px(y: float) -> float:
        return top + (y_max - y) / (y_max - y_min) * plot_h

    p1 = [(x_to_px(x), y_to_px(y), x, y) for x, y in zip(xs, y1s)]
    p2 = [(x_to_px(x), y_to_px(y), x, y) for x, y in zip(xs, y2s)]
    grid_vals = [y_max * i / 5 for i in range(6)]

    lines = svg_header(width, height)
    lines.append(f'<text x="{width/2:.1f}" y="24" text-anchor="middle" font-size="18" font-family="Arial">{title}</text>')
    for gv in grid_vals:
        py = y_to_px(gv)
        lines.append(f'<line x1="{left}" y1="{py:.1f}" x2="{width-right}" y2="{py:.1f}" stroke="#e5e7eb" stroke-width="1"/>')
        tick_label = f"{int(round(gv))}" if integer_labels else f"{gv:.1f}"
        lines.append(f'<text x="{left-10}" y="{py+4:.1f}" text-anchor="end" font-size="12" font-family="Arial">{tick_label}</text>')
    lines.append(f'<line x1="{left}" y1="{top}" x2="{left}" y2="{height-bottom}" stroke="black" stroke-width="1.5"/>')
    lines.append(f'<line x1="{left}" y1="{height-bottom}" x2="{width-right}" y2="{height-bottom}" stroke="black" stroke-width="1.5"/>')
    for x in xs:
        px = x_to_px(x)
        lines.append(f'<line x1="{px:.1f}" y1="{height-bottom}" x2="{px:.1f}" y2="{height-bottom+6}" stroke="black"/>')
        lines.append(f'<text x="{px:.1f}" y="{height-bottom+24}" text-anchor="middle" font-size="12" font-family="Arial">{int(x)}</text>')
    lines.append(f'<text x="{left + plot_w/2:.1f}" y="{height-28}" text-anchor="middle" font-size="13" font-family="Arial">Min Contig Length (bp)</text>')
    lines.append(f'<text x="22" y="{top + plot_h/2:.1f}" transform="rotate(-90 22,{top + plot_h/2:.1f})" text-anchor="middle" font-size="13" font-family="Arial">{y_label}</text>')

    poly1 = " ".join(f"{px:.1f},{py:.1f}" for px, py, _, _ in p1)
    poly2 = " ".join(f"{px:.1f},{py:.1f}" for px, py, _, _ in p2)
    lines.append(f'<polyline points="{poly1}" fill="none" stroke="{SHEET_COLORS["cazy"]}" stroke-width="2.5"/>')
    lines.append(f'<polyline points="{poly2}" fill="none" stroke="{SHEET_COLORS["camper"]}" stroke-width="2.5"/>')
    for px, py, _, y in p1:
        lines.append(f'<circle cx="{px:.1f}" cy="{py:.1f}" r="4.5" fill="{SHEET_COLORS["cazy"]}"/>')
        point_label = f"{int(round(y))}" if integer_labels else f"{y:.1f}"
        lines.append(f'<text x="{px:.1f}" y="{py-10:.1f}" text-anchor="middle" font-size="11" font-family="Arial">{point_label}</text>')
    for px, py, _, y in p2:
        lines.append(f'<rect x="{px-4:.1f}" y="{py-4:.1f}" width="8" height="8" fill="{SHEET_COLORS["camper"]}"/>')
        point_label = f"{int(round(y))}" if integer_labels else f"{y:.1f}"
        lines.append(f'<text x="{px:.1f}" y="{py+18:.1f}" text-anchor="middle" font-size="11" font-family="Arial">{point_label}</text>')

    legend_y = height - 56
    lines.append(f'<line x1="{left}" y1="{legend_y}" x2="{left+28}" y2="{legend_y}" stroke="{SHEET_COLORS["cazy"]}" stroke-width="2.5"/>')
    lines.append(f'<circle cx="{left+14}" cy="{legend_y}" r="4.5" fill="{SHEET_COLORS["cazy"]}"/>')
    lines.append(f'<text x="{left+36}" y="{legend_y+4}" font-size="12" font-family="Arial">{legend1}</text>')
    lines.append(f'<line x1="{left+190}" y1="{legend_y}" x2="{left+218}" y2="{legend_y}" stroke="{SHEET_COLORS["camper"]}" stroke-width="2.5"/>')
    lines.append(f'<rect x="{left+200:.1f}" y="{legend_y-4:.1f}" width="8" height="8" fill="{SHEET_COLORS["camper"]}"/>')
    lines.append(f'<text x="{left+226}" y="{legend_y+4}" font-size="12" font-family="Arial">{legend2}</text>')

    lines.append("</svg>")
    out_path.write_text("\n".join(lines))


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    minlens = [int(x) for x in args.minlens.split(",") if x.strip()]
    runtime_map = parse_runtime_map(args.runtime)

    rows: list[dict[str, object]] = []
    for minlen in minlens:
        if minlen == 2500:
            ann_dir = Path(args.default_annotation)
        else:
            ann_dir = Path(args.scan_annotation_root) / f"SRR10692699_50M_seed11_min{minlen}"
        metrics = annotation_metrics(ann_dir / "annotations.tsv")
        module_count = functional_module_count(ann_dir / "distillate" / "metabolism_summary.xlsx")
        runtime_hours = runtime_hours_from_slurm(runtime_map[minlen])
        row = {
            "min_contig_len": minlen,
            "runtime_hours": runtime_hours,
            "functional_module_count": module_count,
            **metrics,
        }
        rows.append(row)

    rows.sort(key=lambda r: int(r["min_contig_len"]))
    write_tsv(rows, outdir / "dram_minlen_summary.tsv")
    write_dual_series_svg(
        rows,
        "min_contig_len",
        "raw_cazy_hits",
        "raw_camper_hits",
        "Raw CAZy/CAMPER Hits vs Min Contig Length",
        "Raw hit count",
        "CAZy",
        "CAMPER",
        outdir / "dram_minlen_raw_hits.svg",
        integer_labels=True,
    )
    write_dual_series_svg(
        rows,
        "min_contig_len",
        "unique_cazy_ids",
        "unique_camper_ids",
        "Unique CAZy/CAMPER IDs vs Min Contig Length",
        "Unique ID count",
        "CAZy",
        "CAMPER",
        outdir / "dram_minlen_unique_ids.svg",
        integer_labels=True,
    )
    write_single_series_svg(
        rows,
        "min_contig_len",
        "runtime_hours",
        "Annotation Runtime vs Min Contig Length",
        "Runtime (hours)",
        "#7c3aed",
        outdir / "dram_minlen_runtime.svg",
    )
    write_single_series_svg(
        rows,
        "min_contig_len",
        "functional_module_count",
        "Functional Module Count vs Min Contig Length",
        "Non-zero functional modules",
        "#be123c",
        outdir / "dram_minlen_functional_modules.svg",
        integer_labels=True,
    )
    print(f"Wrote summary table to {outdir / 'dram_minlen_summary.tsv'}")
    print(f"Wrote plots to {outdir}")


if __name__ == "__main__":
    main()
