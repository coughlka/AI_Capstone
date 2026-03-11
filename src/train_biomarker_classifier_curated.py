"""
train_biomarker_classifier_curated.py

Module-10E
Train curated CRC biomarker classifiers using cross-validation.

Inputs
------
data/processed/model_features_curated.csv
data/processed/model_labels_curated.csv

Outputs
-------
outputs/curated_cv_metrics.csv
outputs/curated_feature_importance.csv
outputs/curated_oof_predictions.csv
outputs/models/random_forest_curated.pkl
outputs/models/logistic_curated.pkl
"""

from pathlib import Path
import pandas as pd
import numpy as np
import joblib

from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    roc_auc_score,
    average_precision_score,
    accuracy_score,
    precision_score,
    recall_score,
    f1_score,
    classification_report,
)


DEFAULT_PROJECT_ROOT = Path(r"D:\Ayan\PSU Capstone Project\AI_Capstone")


def get_leakage_columns(columns):
    """
    Drop columns that are too close to label construction/ranking logic.
    """
    explicit_drop = {
        "crc_evidence_score",
        "crc_evidence_zscore",
        "evidence_breadth",
        "ranked_candidates_final_score",
        "is_known_crc_gene",
        "has_literature_evidence",
        "has_expression_evidence",
        "has_significance_evidence",
        "has_pathway_evidence",
        "has_network_evidence",
        "has_mutation_evidence",
        "has_methylation_evidence",
        "has_cnv_evidence",
        "has_clinical_perf_evidence",
        "expr_sig_interaction",
        "lit_expr_interaction",
        "network_pathway_interaction",
        "missing_numeric_fraction",
        "data_completeness",
    }

    prefix_drop = [
        "ranked_candidates_",
    ]

    suffix_drop = [
        "_norm",
        "_z",
    ]

    leakage_cols = set()

    for col in columns:
        if col in explicit_drop:
            leakage_cols.add(col)
            continue

        if any(col.startswith(prefix) for prefix in prefix_drop):
            leakage_cols.add(col)
            continue

        if any(col.endswith(suffix) for suffix in suffix_drop):
            leakage_cols.add(col)
            continue

    return sorted(leakage_cols)


def evaluate_predictions(y_true, y_prob, threshold=0.5):
    y_pred = (y_prob >= threshold).astype(int)

    metrics = {
        "roc_auc": roc_auc_score(y_true, y_prob),
        "pr_auc": average_precision_score(y_true, y_prob),
        "accuracy": accuracy_score(y_true, y_pred),
        "precision": precision_score(y_true, y_pred, zero_division=0),
        "recall": recall_score(y_true, y_pred, zero_division=0),
        "f1": f1_score(y_true, y_pred, zero_division=0),
    }
    return metrics, y_pred


def train_biomarker_classifier_curated(project_root):
    project_root = Path(project_root)

    X_path = project_root / "data" / "processed" / "model_features_curated.csv"
    y_path = project_root / "data" / "processed" / "model_labels_curated.csv"

    outputs_dir = project_root / "outputs"
    models_dir = outputs_dir / "models"
    outputs_dir.mkdir(parents=True, exist_ok=True)
    models_dir.mkdir(parents=True, exist_ok=True)

    metrics_path = outputs_dir / "curated_cv_metrics.csv"
    fi_path = outputs_dir / "curated_feature_importance.csv"
    oof_path = outputs_dir / "curated_oof_predictions.csv"
    rf_model_path = models_dir / "random_forest_curated.pkl"
    lr_model_path = models_dir / "logistic_curated.pkl"

    print("\nLoading curated model data...")
    X = pd.read_csv(X_path)
    y = pd.read_csv(y_path).squeeze("columns")

    print("Original feature matrix:", X.shape)
    print("Label distribution:")
    print(y.value_counts())

    # Drop leakage-prone columns
    leakage_cols = get_leakage_columns(X.columns)
    X = X.drop(columns=leakage_cols, errors="ignore")

    print(f"\nDropped leakage-prone columns: {len(leakage_cols)}")
    for c in leakage_cols:
        print(" -", c)

    X = X.select_dtypes(include=["number"]).copy()
    print("\nLeakage-controlled numeric feature matrix:", X.shape)

    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

    rf = Pipeline([
        ("imputer", SimpleImputer(strategy="median")),
        ("model", RandomForestClassifier(
            n_estimators=400,
            random_state=42,
            n_jobs=-1
        ))
    ])

    lr = Pipeline([
        ("imputer", SimpleImputer(strategy="median")),
        ("scaler", StandardScaler()),
        ("model", LogisticRegression(
            max_iter=3000,
            random_state=42
        ))
    ])

    print("\nRunning 5-fold CV for Random Forest...")
    rf_oof_prob = cross_val_predict(
        rf, X, y, cv=cv, method="predict_proba", n_jobs=1
    )[:, 1]
    rf_metrics, rf_oof_pred = evaluate_predictions(y, rf_oof_prob)

    print("Random Forest CV metrics:")
    for k, v in rf_metrics.items():
        print(f"  {k}: {v:.4f}")

    print("\nRunning 5-fold CV for Logistic Regression...")
    lr_oof_prob = cross_val_predict(
        lr, X, y, cv=cv, method="predict_proba", n_jobs=1
    )[:, 1]
    lr_metrics, lr_oof_pred = evaluate_predictions(y, lr_oof_prob)

    print("Logistic Regression CV metrics:")
    for k, v in lr_metrics.items():
        print(f"  {k}: {v:.4f}")

    # Fit final models on full curated dataset
    print("\nFitting final models on full curated dataset...")
    rf.fit(X, y)
    lr.fit(X, y)

    joblib.dump(rf, rf_model_path)
    joblib.dump(lr, lr_model_path)

    # Feature importance from fitted RF
    rf_model = rf.named_steps["model"]
    importance_df = pd.DataFrame({
        "feature": X.columns,
        "importance": rf_model.feature_importances_
    }).sort_values("importance", ascending=False)
    importance_df.to_csv(fi_path, index=False)

    # Save metrics
    metrics_df = pd.DataFrame([
        {"model": "RandomForest", **rf_metrics},
        {"model": "LogisticRegression", **lr_metrics},
    ])
    metrics_df.to_csv(metrics_path, index=False)

    # Save out-of-fold predictions
    oof_df = pd.DataFrame({
        "true_label": y,
        "rf_oof_probability": rf_oof_prob,
        "rf_oof_pred": rf_oof_pred,
        "lr_oof_probability": lr_oof_prob,
        "lr_oof_pred": lr_oof_pred,
    })
    oof_df.to_csv(oof_path, index=False)

    print("\nSaved metrics:")
    print(metrics_path)

    print("\nSaved feature importance:")
    print(fi_path)

    print("\nSaved OOF predictions:")
    print(oof_path)

    print("\nSaved final curated models:")
    print(rf_model_path)
    print(lr_model_path)

    print("\nRandom Forest OOF classification report:\n")
    print(classification_report(y, rf_oof_pred))

    print("\nLogistic Regression OOF classification report:\n")
    print(classification_report(y, lr_oof_pred))

    return rf, lr, metrics_df, importance_df


def run_training_curated(project_root):
    return train_biomarker_classifier_curated(project_root)


if __name__ == "__main__":
    run_training_curated(DEFAULT_PROJECT_ROOT)