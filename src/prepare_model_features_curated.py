"""
prepare_model_features_curated.py

Module-9D
Prepare model-ready features from the curated ML dataset.

Inputs
------
data/processed/ml_training_dataset_curated.csv

Outputs
-------
data/processed/model_features_curated.csv
data/processed/model_labels_curated.csv
data/processed/model_feature_list_curated.txt
"""

from pathlib import Path
import pandas as pd
import numpy as np


DEFAULT_PROJECT_ROOT = Path(r"D:\Ayan\PSU Capstone Project\AI_Capstone")


def prepare_model_features_curated(project_root):

    project_root = Path(project_root)

    input_path = project_root / "data/processed/ml_training_dataset_curated.csv"

    X_out = project_root / "data/processed/model_features_curated.csv"
    y_out = project_root / "data/processed/model_labels_curated.csv"
    feature_list_out = project_root / "data/processed/model_feature_list_curated.txt"

    print("\nLoading curated ML dataset...")
    df = pd.read_csv(input_path)

    print("Input dataset shape:", df.shape)

    if "crc_label" not in df.columns:
        raise ValueError("crc_label column missing.")

    y = df["crc_label"]

    # Columns to explicitly drop
    drop_cols = [
        "gene",
        "gene_symbol",
        "label_reason",
        "is_curated_crc_gene",
        "negative_priority_score",
        "lit_evidence_title",
        "lit_evidence_snippet",
        "path_evidence_top_pathways",
        "omics_evidence_gene_symbol"
    ]

    drop_cols = [c for c in drop_cols if c in df.columns]

    print("\nDropped identifier/text columns:")
    for c in drop_cols:
        print(" -", c)

    X = df.drop(columns=drop_cols)

    print("\nShape after explicit drop:", X.shape)

    # Keep only numeric features
    X = X.select_dtypes(include=[np.number])

    print("Numeric-only feature shape:", X.shape)

    # Drop label column if still present
    if "crc_label" in X.columns:
        X = X.drop(columns=["crc_label"])

    # Remove sparse columns (>40% missing)
    missing_fraction = X.isna().mean()

    sparse_cols = missing_fraction[missing_fraction > 0.4].index.tolist()

    if sparse_cols:
        print("\nDropping sparse columns:")
        for c in sparse_cols:
            print(" -", c)

        X = X.drop(columns=sparse_cols)

    # Impute remaining numeric NA
    X = X.fillna(X.median())

    print("\nFinal feature matrix:", X.shape)

    # Save outputs
    X.to_csv(X_out, index=False)
    y.to_csv(y_out, index=False)

    with open(feature_list_out, "w") as f:
        for col in X.columns:
            f.write(col + "\n")

    print("\nSaved:")
    print(X_out)
    print(y_out)
    print(feature_list_out)

    return X, y


def run_prepare_model_features_curated(project_root):
    return prepare_model_features_curated(project_root)


if __name__ == "__main__":
    run_prepare_model_features_curated(DEFAULT_PROJECT_ROOT)