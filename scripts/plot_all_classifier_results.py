from __future__ import annotations

from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, precision_recall_curve, auc


# ============================================================
# Helpers
# ============================================================

def _ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def _safe_read_csv(path: Path) -> Optional[pd.DataFrame]:
    if path.exists():
        return pd.read_csv(path)
    print(f"[WARN] Missing file: {path}")
    return None


def _save_current_figure(path: Path) -> None:
    _ensure_dir(path.parent)
    plt.tight_layout()
    plt.savefig(path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved: {path}")


# ============================================================
# 1) Tumor vs Normal plots
# ============================================================

def plot_tumor_vs_normal_shap(project_root: str, top_n: int = 20) -> None:
    root = Path(project_root)
    output_dir = root / "outputs"
    plot_dir = output_dir / "plots" / "tumor_vs_normal"
    _ensure_dir(plot_dir)

    shap_candidates = [
        output_dir / "shap_gene_importance.csv",
        output_dir / "top_shap_biomarker_genes.csv",
    ]

    shap_path = None
    for p in shap_candidates:
        if p.exists():
            shap_path = p
            break

    if shap_path is None:
        print("[WARN] No tumor-vs-normal SHAP file found.")
        return

    df = pd.read_csv(shap_path)

    gene_col = "gene" if "gene" in df.columns else df.columns[0]
    value_col = "shap_importance" if "shap_importance" in df.columns else df.columns[1]

    top_df = df.sort_values(value_col, ascending=False).head(top_n).copy()
    top_df = top_df.iloc[::-1]

    plt.figure(figsize=(10, 7))
    plt.barh(top_df[gene_col].astype(str), top_df[value_col])
    plt.xlabel("Mean Absolute SHAP Importance")
    plt.ylabel("Gene")
    plt.title(f"Tumor vs Normal: Top {top_n} SHAP Genes")
    _save_current_figure(plot_dir / f"tumor_vs_normal_top_{top_n}_shap_genes.png")


# ============================================================
# 2) CRC vs STAD plots
# ============================================================

def plot_crc_vs_stad_model_comparison(project_root: str) -> None:
    root = Path(project_root)
    output_path = root / "outputs" / "crc_vs_stad_model_comparison.csv"
    plot_dir = root / "outputs" / "plots" / "crc_vs_stad"
    _ensure_dir(plot_dir)

    df = _safe_read_csv(output_path)
    if df is None:
        return

    required_cols = {"experiment", "model", "roc_auc", "pr_auc", "f1"}
    if not required_cols.issubset(df.columns):
        print(f"[WARN] Missing expected columns in {output_path}. Found: {list(df.columns)}")
        return

    experiments = df["experiment"].unique().tolist()

    for experiment in experiments:
        sub = df[df["experiment"] == experiment].copy().reset_index(drop=True)

        # ROC-AUC / PR-AUC grouped bar chart
        x = np.arange(len(sub))
        width = 0.35

        plt.figure(figsize=(10, 6))
        plt.bar(x - width / 2, sub["roc_auc"], width, label="ROC-AUC")
        plt.bar(x + width / 2, sub["pr_auc"], width, label="PR-AUC")
        plt.xticks(x, sub["model"], rotation=20)
        plt.ylabel("Score")
        plt.ylim(0, 1.05)
        plt.title(f"CRC vs STAD Model Comparison\n{experiment}")
        plt.legend()
        safe_name = experiment.replace(" ", "_").replace("/", "_")
        _save_current_figure(plot_dir / f"{safe_name}_roc_pr_comparison.png")

        # F1 chart
        plt.figure(figsize=(8, 5))
        plt.bar(sub["model"], sub["f1"])
        plt.ylabel("F1 Score")
        plt.ylim(0, 1.05)
        plt.title(f"CRC vs STAD F1 Comparison\n{experiment}")
        plt.xticks(rotation=20)
        _save_current_figure(plot_dir / f"{safe_name}_f1_comparison.png")


# ============================================================
# 3) COAD dMMR plots
# ============================================================

def plot_coad_mmr_cv_results(project_root: str) -> None:
    root = Path(project_root)
    cv_path = root / "outputs" / "coad_mmr_classifier" / "cv_results.csv"
    plot_dir = root / "outputs" / "coad_mmr_classifier" / "plots"
    _ensure_dir(plot_dir)

    df = _safe_read_csv(cv_path)
    if df is None:
        return

    # ROC-AUC mean
    if {"model", "roc_auc_mean"}.issubset(df.columns):
        plt.figure(figsize=(8, 5))
        plt.bar(df["model"], df["roc_auc_mean"])
        plt.ylabel("ROC-AUC (CV Mean)")
        plt.ylim(0, 1.0)
        plt.title("COAD dMMR vs pMMR: Cross-Validation ROC-AUC")
        plt.xticks(rotation=20)
        _save_current_figure(plot_dir / "cv_roc_auc_comparison.png")

    # F1 mean
    if {"model", "f1_mean"}.issubset(df.columns):
        plt.figure(figsize=(8, 5))
        plt.bar(df["model"], df["f1_mean"])
        plt.ylabel("F1 Score (CV Mean)")
        plt.ylim(0, 1.0)
        plt.title("COAD dMMR vs pMMR: Cross-Validation F1")
        plt.xticks(rotation=20)
        _save_current_figure(plot_dir / "cv_f1_comparison.png")

    # Balanced accuracy mean
    if {"model", "balanced_accuracy_mean"}.issubset(df.columns):
        plt.figure(figsize=(8, 5))
        plt.bar(df["model"], df["balanced_accuracy_mean"])
        plt.ylabel("Balanced Accuracy (CV Mean)")
        plt.ylim(0, 1.0)
        plt.title("COAD dMMR vs pMMR: Cross-Validation Balanced Accuracy")
        plt.xticks(rotation=20)
        _save_current_figure(plot_dir / "cv_balanced_accuracy_comparison.png")


def plot_coad_mmr_holdout_results(project_root: str) -> None:
    root = Path(project_root)
    holdout_path = root / "outputs" / "coad_mmr_classifier" / "holdout_results.csv"
    plot_dir = root / "outputs" / "coad_mmr_classifier" / "plots"
    _ensure_dir(plot_dir)

    df = _safe_read_csv(holdout_path)
    if df is None:
        return

    metrics = ["roc_auc", "pr_auc", "f1", "balanced_accuracy"]
    available_metrics = [m for m in metrics if m in df.columns]

    if not available_metrics or "model" not in df.columns:
        print("[WARN] Holdout results file missing required columns.")
        return

    x = np.arange(len(df))
    width = 0.18 if len(available_metrics) >= 4 else 0.25

    plt.figure(figsize=(11, 6))
    for i, metric in enumerate(available_metrics):
        plt.bar(x + (i - (len(available_metrics)-1)/2) * width, df[metric], width, label=metric)

    plt.xticks(x, df["model"], rotation=20)
    plt.ylabel("Score")
    plt.ylim(0, 1.05)
    plt.title("COAD dMMR vs pMMR: Holdout Model Performance")
    plt.legend()
    _save_current_figure(plot_dir / "holdout_model_performance.png")

    # Threshold comparison
    if {"model", "default_f1", "tuned_f1"}.issubset(df.columns):
        plt.figure(figsize=(9, 5))
        x = np.arange(len(df))
        width = 0.35
        plt.bar(x - width / 2, df["default_f1"], width, label="Default F1")
        plt.bar(x + width / 2, df["tuned_f1"], width, label="Tuned F1")
        plt.xticks(x, df["model"], rotation=20)
        plt.ylabel("F1 Score")
        plt.ylim(0, 1.05)
        plt.title("COAD dMMR vs pMMR: Effect of Threshold Tuning")
        plt.legend()
        _save_current_figure(plot_dir / "holdout_threshold_tuning_f1.png")


def plot_coad_mmr_roc_pr_curves(project_root: str) -> None:
    root = Path(project_root)
    pred_path = root / "outputs" / "coad_mmr_classifier" / "holdout_predictions.csv"
    plot_dir = root / "outputs" / "coad_mmr_classifier" / "plots"
    _ensure_dir(plot_dir)

    df = _safe_read_csv(pred_path)
    if df is None:
        return

    required_cols = {"model", "y_true", "y_prob"}
    if not required_cols.issubset(df.columns):
        print("[WARN] holdout_predictions.csv missing required columns.")
        return

    # ROC curves
    plt.figure(figsize=(8, 6))
    for model_name, sub in df.groupby("model"):
        fpr, tpr, _ = roc_curve(sub["y_true"], sub["y_prob"])
        roc_val = auc(fpr, tpr)
        plt.plot(fpr, tpr, label=f"{model_name} (AUC={roc_val:.3f})")

    plt.plot([0, 1], [0, 1], linestyle="--")
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title("COAD dMMR vs pMMR: ROC Curves")
    plt.legend()
    _save_current_figure(plot_dir / "holdout_roc_curves.png")

    # PR curves
    plt.figure(figsize=(8, 6))
    for model_name, sub in df.groupby("model"):
        precision, recall, _ = precision_recall_curve(sub["y_true"], sub["y_prob"])
        pr_val = auc(recall, precision)
        plt.plot(recall, precision, label=f"{model_name} (AUC={pr_val:.3f})")

    plt.xlabel("Recall")
    plt.ylabel("Precision")
    plt.title("COAD dMMR vs pMMR: Precision-Recall Curves")
    plt.legend()
    _save_current_figure(plot_dir / "holdout_pr_curves.png")


def plot_coad_mmr_shap(project_root: str, top_n: int = 20) -> None:
    root = Path(project_root)
    shap_path = root / "outputs" / "coad_mmr_classifier" / "shap_feature_importance.csv"
    plot_dir = root / "outputs" / "coad_mmr_classifier" / "plots"
    _ensure_dir(plot_dir)

    df = _safe_read_csv(shap_path)
    if df is None:
        return

    if not {"feature", "mean_abs_shap"}.issubset(df.columns):
        print("[WARN] SHAP file missing expected columns.")
        return

    top_df = df.sort_values("mean_abs_shap", ascending=False).head(top_n).copy()
    top_df = top_df.iloc[::-1]

    plt.figure(figsize=(10, 7))
    plt.barh(top_df["feature"].astype(str), top_df["mean_abs_shap"])
    plt.xlabel("Mean Absolute SHAP Value")
    plt.ylabel("Feature")
    plt.title(f"COAD dMMR vs pMMR: Top {top_n} SHAP Features")
    _save_current_figure(plot_dir / f"top_{top_n}_shap_features.png")


# ============================================================
# 4) SHAP vs biomarker comparison plots
# ============================================================

def plot_shap_biomarker_overlap_summary(project_root: str) -> None:
    root = Path(project_root)
    summary_path = root / "outputs" / "coad_mmr_classifier" / "shap_vs_biomarker_comparison_summary.csv"
    plot_dir = root / "outputs" / "coad_mmr_classifier" / "plots"
    _ensure_dir(plot_dir)

    df = _safe_read_csv(summary_path)
    if df is None:
        return

    if not {"comparison_type", "n_overlap"}.issubset(df.columns):
        print("[WARN] Overlap summary file missing expected columns.")
        return

    plt.figure(figsize=(8, 5))
    plt.bar(df["comparison_type"], df["n_overlap"])
    plt.ylabel("Number of Overlapping Genes")
    plt.title("SHAP vs Biomarker Pipeline Overlap")
    plt.xticks(rotation=20)
    _save_current_figure(plot_dir / "shap_biomarker_overlap_counts.png")

    if "overlap_fraction_of_shap_top" in df.columns:
        plt.figure(figsize=(8, 5))
        plt.bar(df["comparison_type"], df["overlap_fraction_of_shap_top"])
        plt.ylabel("Overlap Fraction of Top SHAP Genes")
        plt.ylim(0, 1.0)
        plt.title("SHAP vs Biomarker Pipeline Overlap Fraction")
        plt.xticks(rotation=20)
        _save_current_figure(plot_dir / "shap_biomarker_overlap_fraction.png")


# ============================================================
# Main orchestration
# ============================================================

def plot_all_classifier_results(project_root: str) -> None:
    print("\n=== Plotting tumor vs normal results ===")
    plot_tumor_vs_normal_shap(project_root)

    print("\n=== Plotting CRC vs STAD results ===")
    plot_crc_vs_stad_model_comparison(project_root)

    print("\n=== Plotting COAD dMMR results ===")
    plot_coad_mmr_cv_results(project_root)
    plot_coad_mmr_holdout_results(project_root)
    plot_coad_mmr_roc_pr_curves(project_root)
    plot_coad_mmr_shap(project_root)

    print("\n=== Plotting SHAP vs biomarker overlap ===")
    plot_shap_biomarker_overlap_summary(project_root)

    print("\nAll available plots generated.")


if __name__ == "__main__":
    # Default to the project root (two levels up from scripts/), overridable via CLI.
    import argparse

    default_root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description="Generate classifier result plots.")
    parser.add_argument("--project-root", type=str, default=str(default_root),
                        help="Project root (default: parent of scripts/).")
    args = parser.parse_args()
    plot_all_classifier_results(args.project_root)
