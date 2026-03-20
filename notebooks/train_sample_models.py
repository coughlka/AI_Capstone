"""
train_sample_models.py

Workflow Module 9
Run sample-level model experiments and compare:
- Logistic Regression
- Random Forest
- XGBoost

Experiments
-----------
1. Tumor vs Normal using all genes
2. Tumor vs Normal using top 1000 variable genes
3. Tumor vs Normal using top 500 variable genes

Outputs
-------
outputs/model_comparison_results.csv
"""

from __future__ import annotations

from pathlib import Path
import re
import warnings

import numpy as np
import pandas as pd

from sklearn.model_selection import train_test_split
from sklearn.metrics import (
    roc_auc_score,
    average_precision_score,
    accuracy_score,
    precision_score,
    recall_score,
    f1_score,
)
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier

try:
    from xgboost import XGBClassifier
    XGBOOST_AVAILABLE = True
except Exception:
    XGBOOST_AVAILABLE = False


DEFAULT_PROJECT_ROOT = Path(__file__).resolve().parent.parent


def normalize_ensembl_gene_id(val):
    """
    Normalize Ensembl gene IDs to canonical ENSG + 11 digits.
    Example:
        ENSG00000123456.12 -> ENSG00000123456
    """
    if pd.isna(val):
        return None

    s = str(val).strip().upper()
    m = re.search(r"(ENSG\d{11})", s)
    return m.group(1) if m else s


def load_ensembl_mapping(mapping_path: Path) -> pd.DataFrame:
    """
    Load Ensembl-to-symbol mapping.
    """
    mapping = pd.read_csv(mapping_path, sep="\t")
    mapping.columns = [str(c).strip().lower() for c in mapping.columns]

    ensembl_candidates = ["ensembl_gene_id", "ensembl_id", "gene_id", "ensembl"]
    symbol_candidates = ["symbol", "gene_symbol", "hgnc_symbol", "gene_name"]

    ensembl_col = next((c for c in ensembl_candidates if c in mapping.columns), None)
    symbol_col = next((c for c in symbol_candidates if c in mapping.columns), None)

    if ensembl_col is None or symbol_col is None:
        raise ValueError(
            f"Could not detect mapping columns. Found columns: {list(mapping.columns)}"
        )

    out = mapping[[ensembl_col, symbol_col]].copy()
    out.columns = ["Ensembl_ID", "gene_symbol"]

    out["Ensembl_ID"] = out["Ensembl_ID"].apply(normalize_ensembl_gene_id)
    out["gene_symbol"] = out["gene_symbol"].astype(str).str.strip().str.upper()

    out = out.dropna(subset=["Ensembl_ID", "gene_symbol"])
    out = out[out["gene_symbol"] != ""]
    out = out.drop_duplicates(subset=["Ensembl_ID"])

    return out


def get_sample_type_from_tcga_barcode(barcode: str):
    """
    TCGA sample type code is the 4th field:
      01 = tumor
      11 = normal
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


def load_sample_expression_matrix(project_root) -> tuple[pd.DataFrame, pd.Series]:
    """
    Load TCGA-COAD expression data and convert to sample-level matrix:
      rows = samples
      cols = genes
    """
    root = Path(project_root)
    data_path = root / "data" / "TCGA-COAD.star_counts.tsv"
    mapping_path = root / "data" / "ensembl_to_symbol_cache.tsv"

    print("Loading expression dataset...")
    expr = pd.read_csv(data_path, sep="\t")
    print("Raw expression shape:", expr.shape)

    if "Ensembl_ID" not in expr.columns:
        raise ValueError("Expected column 'Ensembl_ID' not found in expression file.")

    print("\nLoading Ensembl mapping...")
    mapping_df = load_ensembl_mapping(mapping_path)
    print("Mapping rows:", len(mapping_df))

    expr["Ensembl_ID"] = expr["Ensembl_ID"].apply(normalize_ensembl_gene_id)

    expr = expr.merge(mapping_df, on="Ensembl_ID", how="left")
    print("Mapped gene symbols:", expr["gene_symbol"].notna().sum())

    expr = expr.dropna(subset=["gene_symbol"]).copy()

    sample_cols = [c for c in expr.columns if c not in ["Ensembl_ID", "gene_symbol"]]

    # Average duplicate symbols if any
    expr = expr.groupby("gene_symbol", as_index=False)[sample_cols].mean()
    print("Expression shape after mapping/grouping:", expr.shape)

    expr = expr.set_index("gene_symbol")
    X = expr.T.copy()

    X["sample_type"] = [get_sample_type_from_tcga_barcode(idx) for idx in X.index]
    X = X[X["sample_type"].isin(["tumor", "normal"])].copy()

    y = X["sample_type"].map({"tumor": 1, "normal": 0})
    X = X.drop(columns=["sample_type"])

    print("\nSample-level matrix:", X.shape)
    print("Class distribution:")
    print(y.value_counts())

    return X, y


def select_top_variable_genes(X: pd.DataFrame, n_genes: int) -> pd.DataFrame:
    """
    Select top n most variable genes across samples.
    """
    variances = X.var(axis=0)
    top_genes = variances.sort_values(ascending=False).head(n_genes).index.tolist()
    return X[top_genes].copy()


def evaluate_model(model, X_train, X_test, y_train, y_test) -> dict:
    """
    Fit and evaluate a binary classifier.
    """
    model.fit(X_train, y_train)

    probs = model.predict_proba(X_test)[:, 1]
    preds = (probs >= 0.5).astype(int)

    metrics = {
        "roc_auc": roc_auc_score(y_test, probs),
        "pr_auc": average_precision_score(y_test, probs),
        "accuracy": accuracy_score(y_test, preds),
        "precision": precision_score(y_test, preds, zero_division=0),
        "recall": recall_score(y_test, preds, zero_division=0),
        "f1": f1_score(y_test, preds, zero_division=0),
    }
    return metrics


def build_models():
    """
    Create model dictionary.
    """
    models = {}

    models["LogisticRegression"] = Pipeline([
        ("imputer", SimpleImputer(strategy="median")),
        ("scaler", StandardScaler()),
        ("model", LogisticRegression(
            max_iter=3000,
            class_weight="balanced",
            random_state=42
        ))
    ])

    models["RandomForest"] = Pipeline([
        ("imputer", SimpleImputer(strategy="median")),
        ("model", RandomForestClassifier(
            n_estimators=400,
            class_weight="balanced",
            min_samples_leaf=2,
            n_jobs=-1,
            random_state=42
        ))
    ])

    if XGBOOST_AVAILABLE:
        models["XGBoost"] = Pipeline([
            ("imputer", SimpleImputer(strategy="median")),
            ("model", XGBClassifier(
                n_estimators=300,
                max_depth=6,
                learning_rate=0.05,
                subsample=0.8,
                colsample_bytree=0.8,
                objective="binary:logistic",
                eval_metric="logloss",
                random_state=42
            ))
        ])

    return models


def run_experiment(experiment_name, X, y, models) -> list[dict]:
    """
    Train/test split and evaluate all models on a given feature matrix.
    """
    print(f"\n========== {experiment_name} ==========")
    print("Feature matrix:", X.shape)

    X_train, X_test, y_train, y_test = train_test_split(
        X,
        y,
        test_size=0.2,
        stratify=y,
        random_state=42
    )

    print("Train samples:", len(X_train))
    print("Test samples:", len(X_test))

    rows = []

    for model_name, model in models.items():
        print(f"\nTraining {model_name}...")
        metrics = evaluate_model(model, X_train, X_test, y_train, y_test)

        row = {
            "experiment": experiment_name,
            "n_features": X.shape[1],
            "model": model_name,
            **metrics
        }
        rows.append(row)

        print(
            f"{model_name}: "
            f"ROC-AUC={metrics['roc_auc']:.4f}, "
            f"PR-AUC={metrics['pr_auc']:.4f}, "
            f"F1={metrics['f1']:.4f}"
        )

    return rows


def train_sample_models(project_root):
    """
    Main Module 9 runner.
    """
    root = Path(project_root)
    output_path = root / "outputs" / "model_comparison_results.csv"

    X, y = load_sample_expression_matrix(root)
    models = build_models()

    results = []

    # Experiment A: all genes
    results.extend(run_experiment(
        experiment_name="Tumor_vs_Normal_AllGenes",
        X=X,
        y=y,
        models=models
    ))

    # Experiment B: top 1000 variable genes
    X_1000 = select_top_variable_genes(X, 1000)
    results.extend(run_experiment(
        experiment_name="Tumor_vs_Normal_Top1000VarGenes",
        X=X_1000,
        y=y,
        models=models
    ))

    # Experiment C: top 500 variable genes
    X_500 = select_top_variable_genes(X, 500)
    results.extend(run_experiment(
        experiment_name="Tumor_vs_Normal_Top500VarGenes",
        X=X_500,
        y=y,
        models=models
    ))

    results_df = pd.DataFrame(results)
    results_df = results_df.sort_values(
        by=["experiment", "roc_auc", "f1"],
        ascending=[True, False, False]
    ).reset_index(drop=True)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    results_df.to_csv(output_path, index=False)

    print("\nSaved model comparison results to:")
    print(output_path)

    return results_df


def run_train_sample_models(project_root):
    return train_sample_models(project_root)


if __name__ == "__main__":
    run_train_sample_models(DEFAULT_PROJECT_ROOT)