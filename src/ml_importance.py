"""ML-based gene importance scoring using Random Forest + SHAP.

Trains a Random Forest classifier to distinguish tumor from normal
samples, then uses SHAP (SHapley Additive exPlanations) values to quantify
each gene's contribution to the classification.  This captures non-linear
effects and gene-gene interactions that univariate t-tests miss.
"""

import os

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score, classification_report, confusion_matrix
from sklearn.model_selection import StratifiedKFold

from src.omics import load_expression_data
from src.utils import load_config, ensure_dirs, write_csv

TAG = "ml_importance"


def _train_and_explain(X: np.ndarray, y: np.ndarray, gene_names: list,
                       n_folds: int = 5, seed: int = 42) -> np.ndarray:
    """Train Random Forest with cross-validation and return mean |SHAP| per gene.

    Uses StratifiedKFold so every fold preserves the tumor/normal ratio.
    Random Forest with 1000 trees uses far more features per ensemble than
    XGBoost, yielding non-zero SHAP importance for hundreds of genes —
    critical for biomarker discovery where sparse models miss candidates.

    Returns:
        1-D array of shape (n_genes,) with mean absolute SHAP values.
    """
    import shap

    skf = StratifiedKFold(n_splits=n_folds, shuffle=True, random_state=seed)
    shap_accum = np.zeros(X.shape[1], dtype=np.float64)
    n_samples_seen = 0
    fold_metrics = []
    all_y_true, all_y_pred, all_y_prob = [], [], []

    for fold_idx, (train_idx, val_idx) in enumerate(skf.split(X, y), 1):
        X_train, X_val = X[train_idx], X[val_idx]
        y_train, y_val = y[train_idx], y[val_idx]

        model = RandomForestClassifier(
            n_estimators=1000,
            max_features='sqrt',
            min_samples_leaf=2,
            class_weight='balanced',
            random_state=seed + fold_idx,
            n_jobs=-1,
        )
        model.fit(X_train, y_train)

        y_pred = model.predict(X_val)
        y_prob = model.predict_proba(X_val)[:, 1]
        acc = (y_pred == y_val).mean()
        auc = roc_auc_score(y_val, y_prob)
        fold_metrics.append({'fold': fold_idx, 'accuracy': acc, 'auc': auc,
                             'n_tumor': y_val.sum(), 'n_normal': (y_val == 0).sum()})
        all_y_true.extend(y_val)
        all_y_pred.extend(y_pred)
        all_y_prob.extend(y_prob)
        print(f"[{TAG}]   Fold {fold_idx}/{n_folds}: accuracy = {acc:.3f}, AUC = {auc:.3f}")

        explainer = shap.TreeExplainer(model)
        shap_vals = explainer.shap_values(X_val)

        # Handle different SHAP output formats for binary classification:
        # - Older SHAP: list of [class_0_array, class_1_array]
        # - Newer SHAP: 3-D array of shape (n_samples, n_genes, 2)
        # Use class_1 (tumor) importances in either case.
        if isinstance(shap_vals, list):
            shap_vals = shap_vals[1]
        elif shap_vals.ndim == 3:
            shap_vals = shap_vals[:, :, 1]

        # shap_vals shape: (n_val_samples, n_genes)
        shap_accum += np.abs(shap_vals).sum(axis=0)
        n_samples_seen += len(val_idx)

    mean_abs_shap = shap_accum / n_samples_seen

    # Print CV summary
    accs = [m['accuracy'] for m in fold_metrics]
    aucs = [m['auc'] for m in fold_metrics]
    print(f"\n[{TAG}]   CV Summary ({n_folds}-fold):")
    print(f"[{TAG}]     Accuracy: {np.mean(accs):.3f} +/- {np.std(accs):.3f}")
    print(f"[{TAG}]     AUC:      {np.mean(aucs):.3f} +/- {np.std(aucs):.3f}")

    all_y_true = np.array(all_y_true)
    all_y_pred = np.array(all_y_pred)
    cm = confusion_matrix(all_y_true, all_y_pred)
    print(f"[{TAG}]     Confusion matrix (pooled OOF):")
    print(f"[{TAG}]       TN={cm[0,0]}  FP={cm[0,1]}")
    print(f"[{TAG}]       FN={cm[1,0]}  TP={cm[1,1]}")
    print(f"[{TAG}]     Classification report (pooled OOF):")
    report = classification_report(all_y_true, all_y_pred,
                                   target_names=['normal', 'tumor'])
    for line in report.split('\n'):
        print(f"[{TAG}]       {line}")

    return mean_abs_shap, fold_metrics


def run_ml_importance(config_path: str) -> str:
    """Run ML importance analysis: Random Forest classifier + SHAP explanations.

    Args:
        config_path: Path to the configuration YAML file.

    Returns:
        Path to the output CSV file.
    """
    print(f"[{TAG}] Loading expression data...")
    counts_filtered, tumor_samples, normal_samples, config = load_expression_data(
        config_path, tag=TAG
    )

    ml_config = config.get('ml_importance', {})
    n_folds = ml_config.get('n_folds', 5)
    top_n_features = ml_config.get('top_n_features', 5000)
    outputs_dir = config['paths']['outputs_dir']
    output_path = os.path.join(outputs_dir, 'ml_importance.csv')

    all_samples = tumor_samples + normal_samples
    X_full = counts_filtered[all_samples].T.values  # (samples, genes)
    y = np.array([1] * len(tumor_samples) + [0] * len(normal_samples))
    gene_names = counts_filtered.index.tolist()

    print(f"[{TAG}] Feature matrix: {X_full.shape[0]} samples x {X_full.shape[1]} genes")

    # Pre-select top genes by variance to keep training tractable
    variances = np.nanvar(X_full, axis=0)
    if X_full.shape[1] > top_n_features:
        top_idx = np.argsort(variances)[-top_n_features:]
        X = X_full[:, top_idx]
        selected_genes = [gene_names[i] for i in top_idx]
        print(f"[{TAG}] Pre-selected top {top_n_features} genes by variance")
    else:
        X = X_full
        selected_genes = gene_names
        top_idx = np.arange(len(gene_names))

    # Replace NaN with 0
    X = np.nan_to_num(X, nan=0.0)

    print(f"[{TAG}] Training {n_folds}-fold Random Forest + SHAP...")
    shap_scores, fold_metrics = _train_and_explain(X, y, selected_genes, n_folds=n_folds)

    # Save CV metrics
    metrics_path = os.path.join(outputs_dir, 'ml_cv_metrics.csv')
    pd.DataFrame(fold_metrics).to_csv(metrics_path, index=False)
    print(f"[{TAG}] Wrote CV metrics to {metrics_path}")

    # Build results DataFrame
    results = pd.DataFrame({
        'gene': selected_genes,
        'shap_importance': shap_scores,
    })

    # Add genes that were filtered out (zero importance)
    excluded_genes = [g for g in gene_names if g not in set(selected_genes)]
    if excluded_genes:
        excluded_df = pd.DataFrame({
            'gene': excluded_genes,
            'shap_importance': 0.0,
        })
        results = pd.concat([results, excluded_df], ignore_index=True)

    results = results.sort_values('shap_importance', ascending=False).reset_index(drop=True)
    results['feature_rank'] = range(1, len(results) + 1)

    write_csv(results, output_path)
    print(f"[{TAG}] Wrote {len(results)} gene importance scores to {output_path}")

    # SHAP distribution summary
    nonzero = results[results['shap_importance'] > 0]
    print(f"\n[{TAG}] SHAP distribution:")
    print(f"[{TAG}]   Non-zero genes: {len(nonzero)} / {len(results)}")
    print(f"[{TAG}]   Max:  {nonzero['shap_importance'].max():.6f}")
    print(f"[{TAG}]   Mean: {nonzero['shap_importance'].mean():.6f}")
    print(f"[{TAG}]   Median: {nonzero['shap_importance'].median():.6f}")

    # Print top 15
    print(f"\n[{TAG}] Top 15 genes by SHAP importance:")
    for _, row in results.head(15).iterrows():
        print(f"        {row['feature_rank']:3.0f}. {row['gene']}: "
              f"SHAP={row['shap_importance']:.6f}")

    print(f"\n[{TAG}] Done.")
    return output_path
