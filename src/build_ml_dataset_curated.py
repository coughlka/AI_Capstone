"""
build_ml_dataset_curated.py

Module-8D
Create a curated supervised ML dataset using an external CRC-positive gene list.

Approach
--------
1. Load classifier_dataset.csv
2. Normalize malformed Ensembl IDs in classifier dataset
3. Map Ensembl IDs to gene symbols using ensembl_to_symbol_cache.tsv
4. Label genes as positive if gene_symbol is in curated CRC list
5. Sample negatives from genes not in curated list
6. Save ml_training_dataset_curated.csv

Inputs
------
data/processed/classifier_dataset.csv
data/ensembl_to_symbol_cache.tsv

Outputs
-------
data/processed/ml_training_dataset_curated.csv
data/processed/ml_dataset_curated_summary.txt
"""

from __future__ import annotations

from pathlib import Path
import re
import pandas as pd


DEFAULT_PROJECT_ROOT = Path(r"D:\Ayan\PSU Capstone Project\AI_Capstone")

# Starter curated CRC list
CURATED_CRC_GENES = {
    "APC", "TP53", "KRAS", "PIK3CA", "SMAD4", "FBXW7", "BRAF", "NRAS",
    "PTEN", "TCF7L2", "MLH1", "MSH2", "MSH6", "PMS2", "MUTYH", "BMPR1A",
    "AXIN2", "POLE", "POLD1", "TGFBR2", "CTNNB1", "RNF43", "ARID1A",
    "SOX9", "B2M", "ACVR2A", "FAT4", "ATM", "ERBB2", "MAP2K1", "NTRK1",
    "NTRK2", "NTRK3", "EGFR", "VEGFA", "MYC", "CCND1", "MET", "CD44",
    "MUC1", "CEACAM5", "CEACAM6", "CXCL8", "IL6", "MMP9"
}


def normalize_ensembl_gene_id(val):
    """
    Normalize Ensembl gene IDs to canonical form: ENSG + 11 digits.

    Handles examples like:
      ENSG0000000041913   -> ENSG00000000419
      ENSG00000141510.15 -> ENSG00000141510
      ENSG00000141510    -> ENSG00000141510
    """
    if pd.isna(val):
        return None

    s = str(val).strip().upper()

    # Extract canonical ENSG + 11 digits
    m = re.search(r"(ENSG\d{11})", s)
    if m:
        return m.group(1)

    return s


def load_ensembl_mapping(mapping_path: Path) -> pd.DataFrame:
    """
    Load Ensembl-to-symbol mapping from TSV and normalize Ensembl IDs.
    """
    m = pd.read_csv(mapping_path, sep="\t")
    original_cols = list(m.columns)
    m.columns = [str(c).strip().lower() for c in m.columns]

    ensembl_candidates = ["ensembl_gene_id", "ensembl_id", "gene_id", "ensembl"]
    symbol_candidates = ["symbol", "gene_symbol", "hgnc_symbol", "gene_name"]

    ensembl_col = next((c for c in ensembl_candidates if c in m.columns), None)
    symbol_col = next((c for c in symbol_candidates if c in m.columns), None)

    if ensembl_col is None or symbol_col is None:
        raise ValueError(
            f"Could not detect Ensembl/Symbol columns in mapping file. "
            f"Original columns: {original_cols}"
        )

    out = m[[ensembl_col, symbol_col]].copy()
    out.columns = ["gene", "gene_symbol"]

    out["gene"] = out["gene"].apply(normalize_ensembl_gene_id)
    out["gene_symbol"] = out["gene_symbol"].astype(str).str.strip().str.upper()

    out = out.dropna(subset=["gene", "gene_symbol"])
    out = out[out["gene_symbol"] != ""]
    out = out.drop_duplicates(subset=["gene"])

    return out


def build_ml_dataset_curated(project_root, negative_ratio=2.0, random_state=42):
    """
    Build a curated supervised ML dataset.

    Parameters
    ----------
    project_root : str or Path
        Project root path.
    negative_ratio : float
        Number of negatives relative to positives.
    random_state : int
        Random seed.

    Returns
    -------
    pd.DataFrame
        Curated training dataset.
    """
    project_root = Path(project_root)

    input_path = project_root / "data" / "processed" / "classifier_dataset.csv"
    mapping_path = project_root / "data" / "ensembl_to_symbol_cache.tsv"
    output_path = project_root / "data" / "processed" / "ml_training_dataset_curated.csv"
    summary_path = project_root / "data" / "processed" / "ml_dataset_curated_summary.txt"

    print("\nLoading classifier dataset...")
    df = pd.read_csv(input_path)
    print("Classifier dataset shape:", df.shape)

    print("\nLoading Ensembl-to-symbol mapping...")
    mapping_df = load_ensembl_mapping(mapping_path)
    print("Mapping rows:", len(mapping_df))

    # Normalize classifier gene IDs before merge
    df["gene"] = df["gene"].apply(normalize_ensembl_gene_id)

    # Merge mapping
    df = df.merge(mapping_df, on="gene", how="left")

    mapped_count = df["gene_symbol"].notna().sum()
    print("\nMapped gene symbols:", mapped_count)

    print("\nExample mapped rows:")
    print(df[["gene", "gene_symbol"]].dropna().head(10))

    overlap = sorted(set(df["gene_symbol"].dropna().astype(str).str.upper()) & CURATED_CRC_GENES)
    print("\nCurated overlap count:", len(overlap))
    print("Curated overlap symbols:", overlap[:30])

    # Curated labels
    df["is_curated_crc_gene"] = df["gene_symbol"].isin(CURATED_CRC_GENES).astype(int)

    pos = df[df["is_curated_crc_gene"] == 1].copy()
    neg_pool = df[df["is_curated_crc_gene"] == 0].copy()

    print("\nCurated positives:", len(pos))
    print("Negative pool:", len(neg_pool))

    if len(pos) == 0:
        raise ValueError("No curated positive genes found after Ensembl-to-symbol mapping.")

    target_neg = max(int(round(len(pos) * negative_ratio)), 1)

    # Prefer negatives with some evidence so they are not trivial
    evidence_cols = [
        c for c in [
            "n_sources",
            "max_abs_logfc",
            "significance_strength",
            "pathway_total",
            "total_score_sum",
            "max_score",
        ]
        if c in neg_pool.columns
    ]

    if evidence_cols:
        neg_pool["negative_priority_score"] = neg_pool[evidence_cols].fillna(0).sum(axis=1)
        neg_pool = neg_pool.sort_values("negative_priority_score", ascending=False)

    if len(neg_pool) < target_neg:
        neg = neg_pool.copy()
    else:
        hardest_slice_n = min(len(neg_pool), max(target_neg * 5, target_neg))
        candidate_neg = neg_pool.head(hardest_slice_n).copy()
        neg = candidate_neg.sample(n=target_neg, random_state=random_state).copy()

    # Final labels
    pos["crc_label"] = 1
    pos["label_reason"] = "curated_crc_gene"

    neg["crc_label"] = 0
    neg["label_reason"] = "non_curated_gene"

    train_df = pd.concat([pos, neg], ignore_index=True)
    train_df = train_df.sample(frac=1, random_state=random_state).reset_index(drop=True)

    print("\nFinal class balance:")
    print(train_df["crc_label"].value_counts(dropna=False))

    train_df.to_csv(output_path, index=False)

    found_symbols = sorted(pos["gene_symbol"].dropna().astype(str).str.upper().unique())

    summary_text = "\n".join([
        "Curated ML Dataset Summary",
        "==========================",
        f"Classifier dataset rows: {len(df)}",
        f"Mapping rows: {len(mapping_df)}",
        f"Mapped gene symbols: {mapped_count}",
        f"Curated overlap count: {len(overlap)}",
        f"Curated positives: {len(pos)}",
        f"Negative pool: {len(neg_pool)}",
        f"Selected negatives: {len(neg)}",
        f"Negative ratio: {negative_ratio}",
        "",
        "Class balance:",
        str(train_df['crc_label'].value_counts(dropna=False)),
        "",
        "Curated positive gene symbols found:",
        ", ".join(found_symbols) if found_symbols else "None"
    ])
    summary_path.write_text(summary_text, encoding="utf-8")

    print("\nSaved curated ML dataset to:")
    print(output_path)

    print("\nSaved summary to:")
    print(summary_path)

    return train_df


def run_build_ml_dataset_curated(project_root, negative_ratio=2.0):
    return build_ml_dataset_curated(
        project_root=project_root,
        negative_ratio=negative_ratio,
        random_state=42
    )


if __name__ == "__main__":
    run_build_ml_dataset_curated(DEFAULT_PROJECT_ROOT)