#!/usr/bin/env python3
"""
Convert the manually curated gene-curation workbooks in SAP/ into normalized
reference TSVs that plug into the abundance pipeline exactly like
scripts/scfa_reference.tsv does.

Each output row lists, for one curated feature (a gene / enzyme / reaction),
the set of *alternative* annotation IDs that identify it, formatted to match
the IDs produced in abundance/.../id_tpm.tsv:

    KO      -> bare        e.g. K00929
    EC      -> 'EC:' prefix e.g. EC:2.7.2.7
    CAZy    -> bare family  e.g. GH5, GT19, PL1, CBM32
    CAMPER  -> bare D-id    e.g. D00001

Two reference variants are written under --out-dir (for all curated modules):

    id_only/   KO + CAZy + CAMPER  (no EC — stricter, less broad-match FP)
    id_ec/     KO + EC + CAZy + CAMPER

Component columns (ko_ids, ec_ids, cazy_ids, camper_ids) are identical in both
variants; only the composed `ids` / `n_ids` columns differ.

Unified output schema (one row per curated feature):
    feature_class  module  feature  gene_symbol  ids
    n_ids  ko_ids  ec_ids  cazy_ids  camper_ids  confidence  description  reference

Usage:
    python scripts/build_curation_references.py \
        --sap-dir SAP --out-dir scripts/curation_refs
"""
import argparse
import re
from pathlib import Path

import pandas as pd

# --- ID extraction --------------------------------------------------------
KO_RE = re.compile(r"\bK\d{5}\b")
# EC: allow full (1.2.3.4) and partial (1.2.3.-) codes. No trailing \b because
# '-' is not a word char and would break the boundary. Leading guard avoids
# matching inside longer dotted numbers.
EC_RE = re.compile(r"(?<![\d.])\d+\.\d+\.\d+\.(?:\d+|-)")
CAZY_RE = re.compile(r"\b(?:GH|GT|PL|CE|CBM|AA)\d+(?:_\d+)?\b")
CAMPER_RE = re.compile(r"\bD\d{5}\b")

ID_VARIANTS = ("id_only", "id_ec")

OUT_COLUMNS = [
    "feature_class", "module", "feature", "gene_symbol", "ids",
    "n_ids", "ko_ids", "ec_ids", "cazy_ids", "camper_ids",
    "confidence", "description", "reference",
]


def _txt(row, col):
    """Safe cell -> stripped string ('' for NaN / missing)."""
    if col is None or col not in row:
        return ""
    v = row[col]
    if v is None or (isinstance(v, float) and pd.isna(v)):
        return ""
    return str(v).strip()


def extract_ids(
    text,
    use_ko=True,
    use_ec=True,
    use_cazy=True,
    use_camper=True,
):
    """Return (ko, ec, cazy, camper) lists found in `text`, deduped+sorted."""
    ko = sorted(set(KO_RE.findall(text))) if use_ko else []
    ec = sorted({"EC:" + e for e in EC_RE.findall(text)}) if use_ec else []
    cazy = sorted(set(CAZY_RE.findall(text))) if use_cazy else []
    camper = sorted(set(CAMPER_RE.findall(text))) if use_camper else []
    return ko, ec, cazy, camper


def _split_col(s):
    if not s or (isinstance(s, float) and pd.isna(s)):
        return []
    return [p for p in str(s).split(",") if p]


def compose_ids(ko, ec, cazy, camper, variant):
    """Build the `ids` string for id_only vs id_ec reference variants."""
    if variant == "id_only":
        return ko + cazy + camper
    if variant == "id_ec":
        return ko + ec + cazy + camper
    raise ValueError(f"unknown id variant: {variant}")


def apply_id_variant(df: pd.DataFrame, variant: str) -> pd.DataFrame:
    """Set composed `ids` / `n_ids` from component ID columns."""
    rows = []
    for _, r in df.iterrows():
        ko = _split_col(r.get("ko_ids", ""))
        ec = _split_col(r.get("ec_ids", ""))
        cazy = _split_col(r.get("cazy_ids", ""))
        camper = _split_col(r.get("camper_ids", ""))
        ids = compose_ids(ko, ec, cazy, camper, variant)
        row = r.to_dict()
        row["ids"] = ",".join(ids)
        row["n_ids"] = len(ids)
        rows.append(row)
    return pd.DataFrame(rows, columns=OUT_COLUMNS)


def build_generic(df, cfg):
    """Build unified rows from a simple 'one feature per row' workbook."""
    rows = []
    for _, r in df.iterrows():
        id_text = " ; ".join(_txt(r, c) for c in cfg["id_cols"])
        ko, ec, cazy, camper = extract_ids(
            id_text,
            use_ko=cfg.get("use_ko", True),
            use_ec=cfg.get("use_ec", True),
            use_cazy=cfg.get("use_cazy", True),
            use_camper=cfg.get("use_camper", True),
        )
        module = _txt(r, cfg["module_col"]) or cfg.get("module_default", "unassigned")
        gene_symbol = _txt(r, cfg.get("gene_symbol_col"))
        desc = " | ".join(x for x in (_txt(r, c) for c in cfg.get("desc_cols", [])) if x)
        rows.append({
            "feature_class": cfg["feature_class"],
            "module": module,
            "feature": gene_symbol or (
                ko[0] if ko else (camper[0] if camper else (cazy[0] if cazy else desc[:40]))
            ),
            "gene_symbol": gene_symbol,
            "ids": "",
            "n_ids": 0,
            "ko_ids": ",".join(ko),
            "ec_ids": ",".join(ec),
            "cazy_ids": ",".join(cazy),
            "camper_ids": ",".join(camper),
            "confidence": _txt(r, cfg.get("confidence_col")),
            "description": desc,
            "reference": _txt(r, cfg.get("reference_col")),
        })
    return pd.DataFrame(rows, columns=OUT_COLUMNS)


def build_scfa48(sap_dir):
    """SCFA 48-enzyme workbook: join enzyme systems with expanded KO mapping."""
    path = sap_dir / "gene_curation_SCFA_48_enzyme.xlsx"
    systems = pd.read_excel(path, sheet_name="48_Enzyme_Systems")
    komap = pd.read_excel(path, sheet_name="KO_Mapping_Expanded")

    ko_by_enzyme = {}
    for _, r in komap.iterrows():
        eid = _txt(r, "Enzyme ID")
        ko_text = _txt(r, "KO identifier")
        for k in KO_RE.findall(ko_text):
            ko_by_enzyme.setdefault(eid, set()).add(k)

    rows = []
    for _, r in systems.iterrows():
        eid = _txt(r, "Enzyme ID")
        product = _txt(r, "Product").lower()
        variant = _txt(r, "Pathway variant")
        ec_text = _txt(r, "EC number")
        _, ec, _, _ = extract_ids(
            ec_text, use_ko=False, use_ec=True, use_cazy=False, use_camper=False
        )
        ko = sorted(ko_by_enzyme.get(eid, set()))
        rows.append({
            "feature_class": "SCFA48",
            "module": f"{product}:{variant}",
            "feature": eid,
            "gene_symbol": _txt(r, "Gene / enzyme label"),
            "ids": "",
            "n_ids": 0,
            "ko_ids": ",".join(ko),
            "ec_ids": ",".join(ec),
            "cazy_ids": "",
            "camper_ids": "",
            "confidence": _txt(r, "KO status"),
            "description": _txt(r, "Enzyme name"),
            "reference": _txt(r, "Supporting references"),
        })
    return pd.DataFrame(rows, columns=OUT_COLUMNS)


# Per-workbook configuration (simple one-row-per-feature files).
CONFIGS = {
    "carotenoids": {
        "file": "gene_curation_Carotenoids.xlsx",
        "feature_class": "Carotenoids",
        "module_col": "Carotenoid_derivative",
        "gene_symbol_col": "Symbol",
        "id_cols": ["KO identifier", "EC number"],
        "desc_cols": ["Name", "Metabolic transformation"],
        "confidence_col": "Substrate specificity for listed compound",
        "reference_col": "Compound reference",
        "use_cazy": False,
    },
    "lbp": {
        "file": "gene_curation_LBP.xlsx",
        "feature_class": "LBP_glycan",
        "module_col": "Name",
        "gene_symbol_col": "Putative microbial enzyme gene / CAZyme activity",
        "id_cols": ["EC Number", "Likely CAZy family / class"],
        "desc_cols": ["Structure location", "Residue / linkage from Table 1"],
        "confidence_col": "Confidence",
        "reference_col": "Source URL",
        "use_cazy": True,
    },
    "lps": {
        "file": "gene_curation_LPS_biosynthesis_BRITE_ko01005.xlsx",
        "feature_class": "LPS_biosynthesis",
        "module_col": "Analysis set",
        "gene_symbol_col": "Gene symbol",
        "id_cols": ["KO identifier", "EC number"],
        "desc_cols": ["KEGG KO symbol/name", "Functional role"],
        "confidence_col": None,
        "reference_col": "KEGG source URL",
        "use_cazy": False,
    },
    "polyphenol": {
        "file": "gene_curation_Polyphenol.xlsx",
        "feature_class": "Polyphenol",
        "module_col": "Raw compound",
        "gene_symbol_col": "CAMPER gene_description",
        "id_cols": [
            "CAMPER raw gene_id",
            "CAMPER identifier",
            "KEGG Orthology identifier",
            "EC number from raw gene_id",
        ],
        "desc_cols": ["Polyphenol class", "Polyphenol subclass", "Specific reaction"],
        "confidence_col": "Confidence",
        "reference_col": "Reference(s)",
        "use_cazy": False,
    },
    "bile_acid": {
        "file": "gene_curation_microbial_bile_acid.xlsx",
        "feature_class": "Bile_acid",
        "module_col": "Functional module",
        "gene_symbol_col": "Gene symbol",
        "id_cols": ["KEGG Orthology", "EC number"],
        "desc_cols": ["Enzyme / protein", "Microbial reaction or role"],
        "confidence_col": None,
        "reference_col": "Paper reference(s)",
        "use_cazy": False,
    },
}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sap-dir", default="SAP")
    ap.add_argument("--out-dir", default="scripts/curation_refs")
    args = ap.parse_args()

    sap_dir = Path(args.sap_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    bases = {}
    for name, cfg in CONFIGS.items():
        src = sap_dir / cfg["file"]
        df = pd.read_excel(src)  # first sheet
        bases[name] = build_generic(df, cfg)
    bases["scfa48"] = build_scfa48(sap_dir)

    summary = []
    for variant in ID_VARIANTS:
        variant_dir = out_dir / variant
        variant_dir.mkdir(parents=True, exist_ok=True)
        print(f"\n--- variant: {variant} -> {variant_dir} ---")
        for name, base in bases.items():
            out = apply_id_variant(base, variant)
            dest = variant_dir / f"{name}.tsv"
            out.to_csv(dest, sep="\t", index=False)
            n_with_ids = int((out["n_ids"].astype(int) > 0).sum())
            summary.append((variant, name, len(out), n_with_ids, out["module"].nunique()))
            print(f"[{name}] {len(out)} rows ({n_with_ids} with ids) -> {dest}")

    print("\n=== summary (variant / module: rows / rows_with_ids / modules) ===")
    for variant, name, n, n_ids, n_mod in summary:
        print(f"  {variant:8s} {name:16s} {n:4d} rows | {n_ids:4d} with IDs | {n_mod:3d} modules")


if __name__ == "__main__":
    main()
