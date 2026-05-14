#!/usr/bin/env python3
"""
Verify DRAM distillate outputs (product.tsv) for one sample.

For a given DRAM contig-annotation sample directory, this script:
  1. Recomputes the three product.tsv blocks (KEGG module coverage,
     ETC complex coverage, functional True/False) using DRAM's own
     internal functions on the same annotations.tsv that was fed to
     distill.
  2. Compares the recomputed values against the values currently in
     product.tsv and prints any mismatches.
  3. For functional True/False mismatches, reports which CAZy/KO/PFAM
     IDs from the annotation matched (or failed to match) each row in
     the function_heatmap form.

Usage:
  /path/to/dram/python scripts/verify_dram_distillate.py \\
    --sample-dir annotation/dram_contigs_T2T/metaG/MG0001/MG0001_2M_seed11

Run this with the same DRAM Python that produced the annotations so we
can import mag_annotator and re-use DRAM's algorithms.
"""

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def import_dram():
    try:
        from mag_annotator.summarize_genomes import (
            HEATMAP_MODULES,
            build_module_net,
            make_etc_coverage_df,
            make_functional_df,
            make_module_coverage_frame,
        )
    except Exception as exc:  # pragma: no cover - import-time failure path
        print(
            f"ERROR: failed to import mag_annotator: {exc}\n"
            "Run this script with the DRAM conda env's python, e.g.\n"
            "/quobyte/angelazgrp/ysjin/.conda/envs/dram/bin/python "
            "scripts/verify_dram_distillate.py ...",
            file=sys.stderr,
        )
        raise SystemExit(2)
    return {
        "HEATMAP_MODULES": HEATMAP_MODULES,
        "build_module_net": build_module_net,
        "make_etc_coverage_df": make_etc_coverage_df,
        "make_functional_df": make_functional_df,
        "make_module_coverage_frame": make_module_coverage_frame,
    }


def load_config_sheets(config_loc: Path):
    config = json.loads(config_loc.read_text())
    sheets = config.get("dram_sheets", {})
    required = ["module_step_form", "etc_module_database", "function_heatmap_form"]
    missing = [k for k in required if k not in sheets]
    if missing:
        raise SystemExit(f"ERROR: missing dram_sheets in config: {missing}")
    return {k: Path(sheets[k]) for k in required}


def recompute_product(sample_dir: Path):
    distill_input = sample_dir / "distill_input_fixed"
    annotations_path = distill_input / "annotations.tsv"
    config_path = distill_input / "distill_with_camper.config.json"
    product_path = sample_dir / "distillate" / "product.tsv"
    for p in (annotations_path, config_path, product_path):
        if not p.is_file():
            raise SystemExit(f"ERROR: missing required file {p}")

    dram = import_dram()
    sheets = load_config_sheets(config_path)

    annotations = pd.read_csv(annotations_path, sep="\t", index_col=0)
    module_steps_form = pd.read_csv(sheets["module_step_form"], sep="\t")
    etc_module_df = pd.read_csv(sheets["etc_module_database"], sep="\t")
    function_heatmap_form = pd.read_csv(sheets["function_heatmap_form"], sep="\t")

    module_nets = {
        module: dram["build_module_net"](module_df)
        for module, module_df in module_steps_form.groupby("module")
        if module in dram["HEATMAP_MODULES"]
    }

    module_coverage = dram["make_module_coverage_frame"](
        annotations, module_nets, groupby_column="fasta"
    )
    etc_coverage = dram["make_etc_coverage_df"](
        etc_module_df, annotations, groupby_column="fasta"
    )

    class _SilentLogger:
        def info(self, *args, **kwargs):
            pass

        warning = info
        error = info

    function_df = dram["make_functional_df"](
        annotations, function_heatmap_form, _SilentLogger(), groupby_column="fasta"
    )

    return {
        "annotations": annotations,
        "function_heatmap_form": function_heatmap_form,
        "module_coverage": module_coverage,
        "etc_coverage": etc_coverage,
        "function_df": function_df,
        "product_path": product_path,
    }


def compare_module_coverage(module_coverage: pd.DataFrame, product: pd.DataFrame, atol: float = 1e-9):
    expected = module_coverage.pivot(
        index="genome", columns="module_name", values="step_coverage"
    )
    mismatches = []
    for module_name in expected.columns:
        if module_name not in product.columns:
            mismatches.append((module_name, None, expected[module_name].to_dict(), "missing in product.tsv"))
            continue
        for genome in expected.index:
            exp_val = expected.at[genome, module_name]
            got_val = product.at[genome, module_name]
            if pd.isna(exp_val) and pd.isna(got_val):
                continue
            try:
                if math.isclose(float(exp_val), float(got_val), abs_tol=atol):
                    continue
            except Exception:
                pass
            mismatches.append((module_name, genome, exp_val, got_val))
    return mismatches


def compare_etc_coverage(etc_coverage: pd.DataFrame, product: pd.DataFrame, atol: float = 1e-9):
    expected = etc_coverage.pivot(
        index="genome", columns="complex_module_name", values="percent_coverage"
    )
    mismatches = []
    for col in expected.columns:
        if col not in product.columns:
            mismatches.append((col, None, expected[col].to_dict(), "missing in product.tsv"))
            continue
        for genome in expected.index:
            exp_val = expected.at[genome, col]
            got_val = product.at[genome, col]
            if pd.isna(exp_val) and pd.isna(got_val):
                continue
            try:
                if math.isclose(float(exp_val), float(got_val), abs_tol=atol):
                    continue
            except Exception:
                pass
            mismatches.append((col, genome, exp_val, got_val))
    return mismatches


def compare_function_presence(function_df: pd.DataFrame, product: pd.DataFrame):
    expected = function_df.pivot(
        index="genome", columns="category_function_name", values="present"
    )
    mismatches = []
    for col in expected.columns:
        if col not in product.columns:
            mismatches.append((col, None, expected[col].to_dict(), "missing in product.tsv"))
            continue
        for genome in expected.index:
            exp_val = bool(expected.at[genome, col])
            raw = product.at[genome, col]
            if isinstance(raw, str):
                got_val = raw.strip().lower() == "true"
            else:
                got_val = bool(raw)
            if exp_val == got_val:
                continue
            mismatches.append((col, genome, exp_val, got_val))
    return mismatches


def annotate_function_mismatch(
    function_heatmap_form: pd.DataFrame,
    annotations: pd.DataFrame,
    function_full_name: str,
):
    from mag_annotator.summarize_genomes import get_ids_from_annotations_all

    if ": " not in function_full_name:
        return ""
    _, function_name = function_full_name.split(": ", 1)
    sub = function_heatmap_form[function_heatmap_form["function_name"] == function_name]
    if sub.empty:
        return f"  (no rule rows found for function_name={function_name!r})"
    detail_lines = [f"  function_name={function_name!r} ({len(sub)} rule row(s))"]
    for _, frame in sub.groupby("fasta", sort=False) if "fasta" in sub.columns else [(None, sub)]:
        for _, row in sub.iterrows():
            ids = {i.strip() for i in str(row["function_ids"]).split(",") if i.strip()}
            # observed IDs per row are computed across the whole sample
            obs = set(get_ids_from_annotations_all(annotations).keys())
            hit = sorted(ids & obs)
            miss = sorted(ids - obs)
            detail_lines.append(
                f"    rule[{row.get('long_function_name', '')!r}] hits={len(hit)} (e.g. {hit[:5]}); "
                f"missing={len(miss)} (e.g. {miss[:5]})"
            )
        break
    return "\n".join(detail_lines)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument(
        "--sample-dir",
        required=True,
        help="Path to sample annotation directory, e.g. "
        "annotation/dram_contigs_T2T/metaG/MG0001/MG0001_2M_seed11",
    )
    ap.add_argument(
        "--atol",
        type=float,
        default=1e-9,
        help="Absolute tolerance for numeric comparisons (default 1e-9).",
    )
    ap.add_argument(
        "--report-dir",
        default=None,
        help="If set, write per-column TSV reports (modules, ETC complexes, "
        "functional True/False) into this directory for human review.",
    )
    args = ap.parse_args()

    sample_dir = Path(args.sample_dir).resolve()
    result = recompute_product(sample_dir)
    product = pd.read_csv(result["product_path"], sep="\t", index_col=0)

    print(f"Sample: {sample_dir}")
    print(f"product.tsv genomes={list(product.index)}")
    print()

    mod_mismatch = compare_module_coverage(result["module_coverage"], product, atol=args.atol)
    etc_mismatch = compare_etc_coverage(result["etc_coverage"], product, atol=args.atol)
    fun_mismatch = compare_function_presence(result["function_df"], product)

    print("== KEGG module step_coverage ==")
    if not mod_mismatch:
        print("  OK: all module step coverage values match.")
    else:
        for col, genome, exp_val, got_val in mod_mismatch:
            print(f"  MISMATCH {col!r} [{genome}] expected={exp_val} got={got_val}")

    print()
    print("== ETC complex percent_coverage ==")
    if not etc_mismatch:
        print("  OK: all ETC complex coverage values match.")
    else:
        for col, genome, exp_val, got_val in etc_mismatch:
            print(f"  MISMATCH {col!r} [{genome}] expected={exp_val} got={got_val}")

    print()
    print("== Functional True/False ==")
    if not fun_mismatch:
        print("  OK: all functional presence flags match.")
    else:
        for col, genome, exp_val, got_val in fun_mismatch:
            print(f"  MISMATCH {col!r} [{genome}] expected={exp_val} got={got_val}")
            detail = annotate_function_mismatch(
                result["function_heatmap_form"], result["annotations"], col
            )
            if detail:
                print(detail)

    total = len(mod_mismatch) + len(etc_mismatch) + len(fun_mismatch)
    print()
    print(f"Total mismatches: {total}")

    if args.report_dir:
        report_dir = Path(args.report_dir).resolve()
        report_dir.mkdir(parents=True, exist_ok=True)
        write_evidence_reports(
            report_dir=report_dir,
            sample_dir=sample_dir,
            annotations=result["annotations"],
            module_coverage=result["module_coverage"],
            etc_coverage=result["etc_coverage"],
            function_df=result["function_df"],
            function_heatmap_form=result["function_heatmap_form"],
            product=product,
        )
        print(f"Evidence reports written under {report_dir}")

    raise SystemExit(0 if total == 0 else 1)


def write_evidence_reports(
    *,
    report_dir: Path,
    sample_dir: Path,
    annotations: pd.DataFrame,
    module_coverage: pd.DataFrame,
    etc_coverage: pd.DataFrame,
    function_df: pd.DataFrame,
    function_heatmap_form: pd.DataFrame,
    product: pd.DataFrame,
):
    from mag_annotator.summarize_genomes import get_ids_from_annotations_all

    observed_ids_by_genome = {
        genome: set(get_ids_from_annotations_all(frame).keys())
        for genome, frame in annotations.groupby("fasta")
    }

    module_rows = []
    for _, row in module_coverage.iterrows():
        genome = row["genome"]
        module_name = row["module_name"]
        steps_present = row["steps_present"]
        steps = row["steps"]
        coverage = row["step_coverage"]
        product_val = product.at[genome, module_name] if module_name in product.columns else ""
        module_rows.append({
            "genome": genome,
            "category": "KEGG module",
            "name": module_name,
            "product_value": product_val,
            "recomputed_value": coverage,
            "steps_present": steps_present,
            "steps": steps,
            "ko_count": row.get("ko_count", ""),
            "kos_present": row.get("kos_present", ""),
            "genes_present": row.get("genes_present", ""),
        })
    pd.DataFrame(module_rows).to_csv(
        report_dir / "module_coverage_evidence.tsv", sep="\t", index=False
    )

    etc_rows = []
    for _, row in etc_coverage.iterrows():
        genome = row["genome"]
        col = row["complex_module_name"]
        product_val = product.at[genome, col] if col in product.columns else ""
        etc_rows.append({
            "genome": genome,
            "category": "ETC complex",
            "name": col,
            "product_value": product_val,
            "recomputed_value": row["percent_coverage"],
            "path_length_coverage": row.get("path_length_coverage", ""),
            "path_length": row.get("path_length", ""),
            "genes_present": row.get("genes", ""),
            "genes_missing": row.get("missing_genes", ""),
        })
    pd.DataFrame(etc_rows).to_csv(
        report_dir / "etc_coverage_evidence.tsv", sep="\t", index=False
    )

    func_rows = []
    grouped_rules = function_heatmap_form.groupby(
        ["category", "subcategory", "function_name"], sort=False
    )
    for genome, observed_ids in observed_ids_by_genome.items():
        for (cat, subcat, fname), group in grouped_rules:
            col_name = f"{cat}: {fname}"
            rule_hits_summary = []
            all_hits = []
            per_rule_pass = []
            for _, rule in group.iterrows():
                rule_ids = {i.strip() for i in str(rule["function_ids"]).split(",") if i.strip()}
                hits = sorted(rule_ids & observed_ids)
                missing = sorted(rule_ids - observed_ids)
                per_rule_pass.append(len(hits) > 0)
                rule_hits_summary.append(
                    f"[{rule.get('long_function_name', fname)}] "
                    f"hits={len(hits)}/{len(rule_ids)}; "
                    f"hit_ids={','.join(hits) if hits else '-'}; "
                    f"missing={','.join(missing) if missing else '-'}"
                )
                all_hits.extend(hits)
            expected_present = bool(per_rule_pass) and all(per_rule_pass)
            product_val = product.at[genome, col_name] if col_name in product.columns else ""
            func_rows.append({
                "genome": genome,
                "category": cat,
                "subcategory": subcat,
                "function_name": fname,
                "product_column": col_name,
                "product_value": product_val,
                "recomputed_value": expected_present,
                "num_rules": len(group),
                "rules_passed": sum(per_rule_pass),
                "all_hits_union": ",".join(sorted(set(all_hits))) if all_hits else "",
                "per_rule_detail": " || ".join(rule_hits_summary),
            })
    pd.DataFrame(func_rows).to_csv(
        report_dir / "functional_presence_evidence.tsv", sep="\t", index=False
    )


if __name__ == "__main__":
    main()
