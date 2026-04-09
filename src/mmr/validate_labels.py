"""Validate MMR/MSI labels in the TCGA-COAD GDC phenotype file.

Original author: Ayan Choudhury (validate_coad_mmr_labels.py).
Adapted by Keith Coughlin to run from the project config and to fit the
src/ module conventions used by the rest of the pipeline.

Reads the GDC phenotype TSV, locates the MMR/MSI column, normalizes the
raw values into a binary MSI/MSS label, intersects against patients
present in the expression matrix, and writes a small validation report.
This is run before training so we can confirm class counts before any
modeling work begins.
"""

from __future__ import annotations

import json
import os
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd

from src.utils import load_config

TAG = "mmr_validate"


def _read_table_auto(path: str) -> pd.DataFrame:
    return pd.read_csv(path, sep=None, engine="python")


def _normalize_tcga_barcode(value: str, level: str = "patient") -> Optional[str]:
    if pd.isna(value):
        return None

    s = str(value).strip().upper().replace(".", "-")
    s = re.sub(r"[^A-Z0-9\-]", "", s)

    if not s.startswith("TCGA-"):
        return None

    parts = s.split("-")
    if len(parts) < 3:
        return None

    if level == "patient":
        return "-".join(parts[:3])

    if level == "sample":
        if len(parts) >= 4:
            return "-".join(parts[:4])
        return None

    raise ValueError("level must be 'patient' or 'sample'")


def _patient_from_sample(sample_barcode: str) -> Optional[str]:
    return _normalize_tcga_barcode(sample_barcode, level="patient")


def _find_barcode_column(df: pd.DataFrame) -> Optional[str]:
    candidates = [
        c for c in df.columns
        if any(k in c.lower() for k in ["barcode", "submitter", "sample", "patient"])
    ]

    for col in candidates:
        vals = df[col].dropna().astype(str).head(100).tolist()
        if any("TCGA-" in v.upper().replace(".", "-") for v in vals):
            return col

    return candidates[0] if candidates else None


def _find_msi_candidate_columns(df: pd.DataFrame) -> List[str]:
    keywords = ["msi", "microsatellite", "mmr", "mismatch"]
    return [c for c in df.columns if any(k in c.lower() for k in keywords)]


def _normalize_msi_value(value: str) -> Optional[str]:
    """Map raw phenotype values to binary MSI/MSS labels.

    The IHC-based MMR loss field uses YES/NO; other fields use MSI-H/MSS/etc.
    Both encodings are mapped to the same MSI/MSS binary used downstream.
    """
    if pd.isna(value):
        return None

    s = str(value).strip().upper()
    s = s.replace("-", "_").replace(" ", "_")
    s = re.sub(r"_+", "_", s)

    if s in {"", "NA", "NAN", "NONE", "NULL", "NOT_REPORTED", "UNKNOWN"}:
        return None

    if s == "YES":
        return "MSI"
    if s == "NO":
        return "MSS"

    positive_patterns = {"MSI_H", "MSIHIGH", "MSI_HIGH", "HIGH", "UNSTABLE", "INSTABLE"}
    negative_patterns = {"MSS", "STABLE", "MSI_L", "MSILOW", "MSI_LOW", "LOW"}

    if s in positive_patterns:
        return "MSI"
    if s in negative_patterns:
        return "MSS"

    if "MSS" in s or "STABLE" in s:
        return "MSS"
    if "MSI" in s and ("HIGH" in s or s.endswith("_H")):
        return "MSI"
    if "MSI" in s and ("LOW" in s or s.endswith("_L")):
        return "MSS"

    return None


def _summarize_candidate_column(df: pd.DataFrame, col: str) -> Dict:
    raw_counts = (
        df[col]
        .fillna("<NA>")
        .astype(str)
        .value_counts(dropna=False)
        .head(20)
        .to_dict()
    )

    mapped = df[col].apply(_normalize_msi_value)
    mapped_counts = mapped.value_counts(dropna=False).to_dict()

    return {
        "column": col,
        "raw_top_values": raw_counts,
        "mapped_counts": {str(k): int(v) for k, v in mapped_counts.items()},
        "usable_rows": int(mapped.notna().sum()),
        "n_mapped_classes": int(mapped.dropna().nunique()),
    }


def _select_best_msi_column(df: pd.DataFrame, candidates: List[str]) -> Optional[str]:
    best_col = None
    best_score = -1

    for col in candidates:
        mapped = df[col].apply(_normalize_msi_value)
        usable = mapped.notna().sum()
        n_classes = mapped.dropna().nunique()

        # Prefer columns with 2 classes, then by number of usable rows.
        score = usable + (100000 if n_classes == 2 else 0)

        if n_classes >= 2 and score > best_score:
            best_col = col
            best_score = score

    return best_col


def _build_expression_patient_set(expr_path: str) -> Tuple[int, int, List[str]]:
    expr_df = _read_table_auto(expr_path)
    sample_cols = [c for c in expr_df.columns if c != expr_df.columns[0]]

    patient_ids = []
    for s in sample_cols:
        norm_sample = _normalize_tcga_barcode(s, level="sample")
        if norm_sample is not None:
            patient = _patient_from_sample(norm_sample)
            if patient is not None:
                patient_ids.append(patient)

    return len(sample_cols), len(set(patient_ids)), sorted(set(patient_ids))


def run_mmr_label_validation(config_path: str) -> str:
    """Validate MMR phenotype labels against the expression matrix.

    Args:
        config_path: Path to the project YAML config.

    Returns:
        Path to the JSON validation report.
    """
    config = load_config(config_path)
    mmr_config = config.get("mmr", {})

    phenotype_path = mmr_config.get("phenotype_path")
    if not phenotype_path:
        raise ValueError(
            "config['mmr']['phenotype_path'] is required to validate labels."
        )

    expression_path = mmr_config.get(
        "expression_path", config["omics"]["counts_path"]
    )
    output_dir = Path(
        mmr_config.get("output_dir", "outputs/coad_mmr_classifier")
    )
    min_class_size = int(mmr_config.get("min_class_size", 15))

    output_dir.mkdir(parents=True, exist_ok=True)

    if not os.path.exists(phenotype_path):
        raise FileNotFoundError(
            f"Phenotype file not found: {phenotype_path}. "
            "Download TCGA-COAD.GDC_phenotype.tsv from the GDC/Xena hub "
            "and place it in the data/ directory."
        )

    if not os.path.exists(expression_path):
        raise FileNotFoundError(f"Expression file not found: {expression_path}")

    print(f"[{TAG}] Loading phenotype: {phenotype_path}")
    pheno_df = _read_table_auto(phenotype_path)
    print(f"[{TAG}] Phenotype shape: {pheno_df.shape}")

    barcode_col = _find_barcode_column(pheno_df)
    if barcode_col is None:
        raise ValueError("Could not identify barcode column in phenotype file.")

    print(f"[{TAG}] Barcode column: {barcode_col}")

    pheno_df["raw_barcode"] = pheno_df[barcode_col].astype(str)
    pheno_df["patient_barcode"] = pheno_df["raw_barcode"].apply(
        lambda x: _normalize_tcga_barcode(x, level="patient")
    )

    candidate_cols = _find_msi_candidate_columns(pheno_df)
    print(f"[{TAG}] Found {len(candidate_cols)} MSI/MMR candidate columns.")

    candidate_summaries = [
        _summarize_candidate_column(pheno_df, col) for col in candidate_cols
    ]

    with open(output_dir / "msi_candidate_columns_summary.json", "w", encoding="utf-8") as f:
        json.dump(candidate_summaries, f, indent=2)

    best_col = _select_best_msi_column(pheno_df, candidate_cols)
    if best_col is None:
        report = {
            "status": "FAILED",
            "reason": "No MSI/MMR phenotype column produced usable binary labels.",
            "phenotype_barcode_column": barcode_col,
            "candidate_columns": candidate_cols,
        }
        report_path = output_dir / "msi_validation_report.json"
        with open(report_path, "w", encoding="utf-8") as f:
            json.dump(report, f, indent=2)
        print(f"[{TAG}] FAILED: no usable MSI column. Report at {report_path}")
        return str(report_path)

    print(f"[{TAG}] Selected MMR/MSI column: {best_col}")

    pheno_df["msi_label"] = pheno_df[best_col].apply(_normalize_msi_value)

    label_df = (
        pheno_df[["patient_barcode", best_col, "msi_label"]]
        .dropna(subset=["patient_barcode"])
        .drop_duplicates()
        .copy()
    )

    raw_label_counts = label_df["msi_label"].value_counts(dropna=False).to_dict()

    print(f"[{TAG}] Building expression patient set from {expression_path}")
    n_expr_samples, n_expr_patients, expr_patients = _build_expression_patient_set(
        expression_path
    )
    expr_patient_set = set(expr_patients)

    merged_df = label_df[label_df["patient_barcode"].isin(expr_patient_set)].copy()
    merged_df = merged_df.dropna(subset=["msi_label"]).copy()
    merged_counts = merged_df["msi_label"].value_counts().to_dict()

    status = "PASS"
    warnings: List[str] = []
    if merged_df["msi_label"].nunique() != 2:
        status = "FAILED"
        warnings.append(
            "Merged expression/phenotype set does not contain both MSI and MSS classes."
        )
    if len(merged_counts) == 2 and min(merged_counts.values()) < min_class_size:
        status = "WARNING"
        warnings.append(
            f"At least one class has fewer than {min_class_size} matched samples: {merged_counts}"
        )

    report = {
        "status": status,
        "phenotype_shape": list(pheno_df.shape),
        "phenotype_barcode_column": barcode_col,
        "selected_msi_column": best_col,
        "raw_candidate_columns": candidate_cols,
        "raw_label_counts_before_expression_merge": {
            str(k): int(v) for k, v in raw_label_counts.items()
        },
        "expression_n_sample_columns": int(n_expr_samples),
        "expression_n_unique_patients": int(n_expr_patients),
        "n_matched_patients_with_valid_labels": int(merged_df["patient_barcode"].nunique()),
        "matched_label_counts": {str(k): int(v) for k, v in merged_counts.items()},
        "expression_patient_overlap_rate": float(
            len(set(merged_df["patient_barcode"])) / max(1, n_expr_patients)
        ),
        "warnings": warnings,
        "notes": [
            "YES in IHC-based MMR loss field maps to MSI.",
            "NO in IHC-based MMR loss field maps to MSS.",
            "MSI-H maps to MSI; MSS and MSI-L map to MSS.",
        ],
    }

    report_path = output_dir / "msi_validation_report.json"
    with open(report_path, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2)

    merged_df.to_csv(
        output_dir / "validated_msi_labels_matched_to_expression.csv", index=False
    )

    print(f"[{TAG}] Status: {status}")
    print(f"[{TAG}] Matched label counts: {merged_counts}")
    print(f"[{TAG}] Wrote validation report to {report_path}")
    if warnings:
        for w in warnings:
            print(f"[{TAG}] WARNING: {w}")

    return str(report_path)


if __name__ == "__main__":
    import sys
    cfg = sys.argv[1] if len(sys.argv) > 1 else "config/config.yaml"
    run_mmr_label_validation(cfg)
