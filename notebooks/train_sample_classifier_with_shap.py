"""
train_sample_classifier_with_shap.py

Train a sample-level classifier (tumor vs normal) using TCGA-COAD
expression data and identify biomarker genes using SHAP.

This avoids circular modeling on pipeline-derived gene scores.
"""

from __future__ import annotations

from pathlib import Path
import re
import numpy as np
import pandas as pd
import joblib
import shap

from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score, classification_report


DEFAULT_PROJECT_ROOT = Path(r"D:\Ayan\PSU Capstone Project\AI_Capstone")


def normalize_ensembl_gene_id(val):
    """
    Normalize Ensembl gene IDs to canonical ENSG + 11 digits.
    Handles IDs such as:
      ENSG00000000003.15 -> ENSG00000000003
    """
    if pd.isna(val):
        return None
    s = str(val).strip().upper()
    m = re.search(r"(ENSG\d{11})", s)
    return m.group(1) if m else s


def load_ensembl_mapping(mapping_path: Path) -> pd.DataFrame:
    """
    Load Ensembl-to-symbol mapping file.
    Expected columns include:
      ensembl_id, gene_symbol
    """
    m = pd.read_csv(mapping_path, sep="\t")
    m.columns = [str(c).strip().lower() for c in m.columns]

    ensembl_candidates = ["ensembl_gene_id", "ensembl_id", "gene_id", "ensembl"]
    symbol_candidates = ["symbol", "gene_symbol", "hgnc_symbol", "gene_name"]

    ensembl_col = next((c for c in ensembl_candidates if c in m.columns), None)
    symbol_col = next((c for c in symbol_candidates if c in m.columns), None)

    if ensembl_col is None or symbol_col is None:
        raise ValueError(f"Could not detect mapping columns. Found: {list(m.columns)}")

    out = m[[ensembl_col, symbol_col]].copy()
    out.columns = ["Ensembl_ID", "gene_symbol"]

    out["Ensembl_ID"] = out["Ensembl_ID"].apply(normalize_ensembl_gene_id)
    out["gene_symbol"] = out["gene_symbol"].astype(str).str.strip().str.upper()

    out = out.dropna(subset=["Ensembl_ID", "gene_symbol"])
    out = out[out["gene_symbol"] != ""]
    out = out.drop_duplicates(subset=["Ensembl_ID"])

    return out


def get_sample_type_from_tcga_barcode(barcode: str):
    """
    TCGA sample type code is the 4th field in the barcode.
    Example:
      TCGA-XX-XXXX-01A -> tumor
      TCGA-XX-XXXX-11A -> normal
    """
    parts = str(barcode).split("-")
    if len(parts) < 4:
        return None

    sample_code = parts[3][:2]

    if sample_code == "01":
        return "tumor"
    if sample_code == "11":
        return "normal"

    return None


def train_sample_classifier(project_root):
    root = Path(project_root)

    data_path = root / "data" / "raw" / "TCGA-COAD.star_counts.tsv"
    mapping_path = root / "data" / "ensembl_to_symbol_cache.tsv"

    output_dir = root / "outputs"
    output_dir.mkdir(parents=True, exist_ok=True)

    model_dir = output_dir / "models"
    model_dir.mkdir(parents=True, exist_ok=True)

    print("Loading expression dataset...")
    expr = pd.read_csv(data_path, sep="\t")

    print("Raw expression shape:", expr.shape)

    if "Ensembl_ID" not in expr.columns:
        raise ValueError("Expected column 'Ensembl_ID' not found in expression file.")

    print("\nLoading Ensembl mapping...")
    mapping_df = load_ensembl_mapping(mapping_path)
    print("Mapping rows:", len(mapping_df))

    expr["Ensembl_ID"] = expr["Ensembl_ID"].apply(normalize_ensembl_gene_id)

    # Merge gene symbols
    expr = expr.merge(mapping_df, on="Ensembl_ID", how="left")

    print("Mapped gene symbols:", expr["gene_symbol"].notna().sum())

    # Keep genes with symbols
    expr = expr.dropna(subset=["gene_symbol"]).copy()

    # Remove duplicate gene symbols by averaging expression
    sample_cols = [c for c in expr.columns if c not in ["Ensembl_ID", "gene_symbol"]]
    expr = expr.groupby("gene_symbol", as_index=False)[sample_cols].mean()

    print("Expression shape after symbol mapping/grouping:", expr.shape)

    # Set gene symbols as index and transpose to samples x genes
    expr = expr.set_index("gene_symbol")
    X = expr.T.copy()

    # Create sample labels
    X["sample_type"] = [get_sample_type_from_tcga_barcode(idx) for idx in X.index]

    # Keep only tumor and normal
    X = X[X["sample_type"].isin(["tumor", "normal"])].copy()

    y = X["sample_type"].map({"tumor": 1, "normal": 0})
    X = X.drop(columns=["sample_type"])

    print("\nSample-level feature matrix:", X.shape)
    print("Class distribution:")
    print(y.value_counts())

    # Train/test split
    X_train, X_test, y_train, y_test = train_test_split(
        X,
        y,
        test_size=0.2,
        stratify=y,
        random_state=42
    )

    print("\nTrain samples:", len(X_train))
    print("Test samples:", len(X_test))

    # Random Forest with anti-overfitting settings
    rf = RandomForestClassifier(
        n_estimators=400,
        class_weight="balanced",
        min_samples_leaf=2,
        n_jobs=-1,
        random_state=42
    )

    print("\nTraining Random Forest...")
    rf.fit(X_train, y_train)

    preds = rf.predict(X_test)
    probs = rf.predict_proba(X_test)[:, 1]

    auc = roc_auc_score(y_test, probs)

    print("\nROC-AUC:", round(auc, 4))
    print("\nClassification Report\n")
    print(classification_report(y_test, preds))

    # Save model
    model_path = model_dir / "rf_sample_classifier.pkl"
    joblib.dump(rf, model_path)

    print("\nSaved model:")
    print(model_path)

    # SHAP
    print("\nComputing SHAP values...")
    explainer = shap.TreeExplainer(rf)
    shap_values = explainer.shap_values(X_train)

    # Handle different SHAP output formats robustly
    if isinstance(shap_values, list):
        # Older SHAP binary-classification format: [class0_matrix, class1_matrix]
        shap_matrix = shap_values[1]
    else:
        shap_matrix = np.array(shap_values)

        # Newer SHAP may return:
        # (n_samples, n_features) or
        # (n_samples, n_features, n_classes)
        if shap_matrix.ndim == 3:
            # use class 1 = tumor
            shap_matrix = shap_matrix[:, :, 1]
        elif shap_matrix.ndim != 2:
            raise ValueError(f"Unexpected SHAP output shape: {shap_matrix.shape}")

    print("SHAP matrix shape:", shap_matrix.shape)

    shap_importance = np.abs(shap_matrix).mean(axis=0)

    shap_df = pd.DataFrame({
        "gene": X_train.columns.tolist(),
        "shap_importance": shap_importance
    }).sort_values("shap_importance", ascending=False)

    shap_path = output_dir / "shap_gene_importance.csv"
    shap_df.to_csv(shap_path, index=False)

    top_path = output_dir / "top_shap_biomarker_genes.csv"
    shap_df.head(50).to_csv(top_path, index=False)

    print("\nSaved SHAP gene importance:")
    print(shap_path)

    print("\nSaved top SHAP biomarker genes:")
    print(top_path)

    return rf, shap_df


def run_sample_classifier(project_root):
    return train_sample_classifier(project_root)


if __name__ == "__main__":
    run_sample_classifier(DEFAULT_PROJECT_ROOT)