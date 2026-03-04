#!/usr/bin/env python3
"""
Prepare a DRAM distill config that injects CAMPER functions into product heatmap.

This script writes:
1) a merged function_heatmap_form TSV
2) a copied CONFIG JSON with dram_sheets.function_heatmap_form redirected to the merged TSV

Usage:
  python scripts/prepare_distill_config_with_camper.py \
    --base-config <DRAM CONFIG> \
    --annotations <annotations.tsv> \
    --out-config <config.json> \
    --out-heatmap <function_heatmap.tsv>
"""

import argparse
import csv
import json
from pathlib import Path


HEATMAP_COLUMNS = [
    "category",
    "subcategory",
    "function_name",
    "function_ids",
    "long_function_name",
    "gene_symbol",
]


def load_tsv(path: Path):
    with path.open("r", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Empty TSV: {path}")
        return reader.fieldnames, list(reader)


def write_tsv(path: Path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=HEATMAP_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in HEATMAP_COLUMNS})


def has_camper_column(annotations_tsv: Path) -> bool:
    with annotations_tsv.open("r", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None:
            return False
        return "camper_id" in reader.fieldnames


def camper_to_heatmap_rows(camper_rows):
    grouped = {}

    def polyphenol_family(header: str, module: str) -> str:
        parts = [p.strip() for p in header.split(";") if p.strip()]
        if len(parts) >= 2:
            return f"{parts[0]};{parts[1]}"
        if len(parts) == 1:
            return parts[0]
        return module

    for row in camper_rows:
        gene_ids = (row.get("gene_id") or "").strip()
        module = (row.get("module") or "").strip()
        if not gene_ids:
            continue
        header = (row.get("header") or "").strip()
        group_name = polyphenol_family(header, module)
        group = grouped.setdefault(
            group_name,
            {
                "ids": set(),
                "modules": set(),
                "headers": set(),
            },
        )
        for gid in [g.strip() for g in gene_ids.split(",") if g.strip()]:
            group["ids"].add(gid)
        if module:
            group["modules"].add(module)
        if header:
            group["headers"].add(header)

    out = []
    for group_name, data in sorted(grouped.items()):
        full_headers = sorted(data["headers"])
        modules = sorted(data["modules"])
        long_desc_parts = []
        if full_headers:
            long_desc_parts.append("Headers: " + " | ".join(full_headers))
        if modules:
            long_desc_parts.append("Modules: " + " | ".join(modules))
        out.append(
            {
                "category": "CAMPER",
                "subcategory": "CAMPER",
                "function_name": group_name,
                "function_ids": ", ".join(sorted(data["ids"])),
                "long_function_name": " || ".join(long_desc_parts),
                "gene_symbol": "",
            }
        )
    return out


def dedup_rows(rows):
    seen = set()
    out = []
    for row in rows:
        key = tuple((row.get(k, "") or "").strip() for k in HEATMAP_COLUMNS)
        if key in seen:
            continue
        seen.add(key)
        out.append({k: key[i] for i, k in enumerate(HEATMAP_COLUMNS)})
    return out


def main():
    ap = argparse.ArgumentParser(description="Build CAMPER-aware distill config for DRAM product outputs.")
    ap.add_argument("--base-config", required=True, help="Path to DRAM CONFIG JSON")
    ap.add_argument("--annotations", required=True, help="Sample annotations.tsv path")
    ap.add_argument("--out-config", required=True, help="Output config path")
    ap.add_argument("--out-heatmap", required=True, help="Output merged function_heatmap_form TSV")
    args = ap.parse_args()

    base_config = Path(args.base_config)
    annotations = Path(args.annotations)
    out_config = Path(args.out_config)
    out_heatmap = Path(args.out_heatmap)

    config = json.loads(base_config.read_text())
    dram_sheets = config.get("dram_sheets", {})
    heatmap_path = Path(dram_sheets["function_heatmap_form"])

    heatmap_columns, base_rows = load_tsv(heatmap_path)
    missing = [c for c in HEATMAP_COLUMNS if c not in heatmap_columns]
    if missing:
        raise ValueError(f"function_heatmap_form missing required columns: {missing}")

    merged_rows = [{k: (r.get(k, "") or "").strip() for k in HEATMAP_COLUMNS} for r in base_rows]
    camper_added = 0

    if has_camper_column(annotations) and "camper_distillate" in dram_sheets:
        camper_path = Path(dram_sheets["camper_distillate"])
        _, camper_rows = load_tsv(camper_path)
        camper_heatmap_rows = camper_to_heatmap_rows(camper_rows)
        camper_added = len(camper_heatmap_rows)
        merged_rows.extend(camper_heatmap_rows)

    merged_rows = dedup_rows(merged_rows)
    write_tsv(out_heatmap, merged_rows)

    config["dram_sheets"]["function_heatmap_form"] = str(out_heatmap.resolve())
    out_config.parent.mkdir(parents=True, exist_ok=True)
    out_config.write_text(json.dumps(config, indent=2))

    print(
        f"Wrote merged heatmap form to {out_heatmap} with {len(merged_rows)} rows "
        f"(added CAMPER rows before dedup: {camper_added})"
    )
    print(f"Wrote DRAM config for distill to {out_config}")


if __name__ == "__main__":
    main()
