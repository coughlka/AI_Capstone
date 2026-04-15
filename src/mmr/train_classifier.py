"""Train within-COAD dMMR vs pMMR classifier on TCGA-COAD expression data.

Original author: Ayan Choudhury (train_coad_mmr_classifier.py).
Adapted by Keith Coughlin to:
  - read paths and hyperparameters from config/config.yaml
  - load expression data using the same TSV format the rest of the
    pipeline already consumes (data/TCGA-COAD.star_counts.tsv)
  - drop hardcoded Windows project root and joblib snapshots that
    don't fit the rest of the codebase
  - expose run_mmr_classifier(config_path) so the runner script and
    other modules can call it the same way they call run_omics, etc.

Trains and compares Logistic Regression, Random Forest, and (when
available) XGBoost. Saves CV results, holdout metrics, predictions, and
SHAP feature importance for the best holdout model.
"""

from __future__ import annotations

import json
import os
import warnings
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import joblib
import numpy as np
import pandas as pd
from sklearn.base import clone
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import SelectKBest, VarianceThreshold, f_classif
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    average_precision_score,
    balanced_accuracy_score,
    f1_score,
    precision_score,
    recall_score,
    roc_auc_score,
)
from sklearn.model_selection import StratifiedKFold, cross_validate, train_test_split
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from src.utils import load_config

warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

try:
    from xgboost import XGBClassifier
    XGBOOST_AVAILABLE = True
except ImportError:
    XGBOOST_AVAILABLE = False

try:
    import shap
    SHAP_AVAILABLE = True
except ImportError:
    SHAP_AVAILABLE = False

TAG = "mmr_train"


@dataclass
class MMRConfig:
    """Resolved MMR training configuration.

    Built from config['mmr'] in config/config.yaml so the runner script
    can pass everything in one place.
    """
    expression_path: str
    phenotype_path: str
    output_dir: str
    phenotype_barcode_col: str = "submitter_id.samples"
    phenotype_label_col: str = "loss_expression_of_mismatch_repair_proteins_by_ihc"
    test_size: float = 0.20
    random_state: int = 42
    cv_folds: int = 5
    min_nonzero_fraction: float = 0.10
    variance_threshold: float = 0.0
    top_k_features: int = 1000
    n_jobs: int = -1
    max_samples: Optional[int] = None
    biomarker_candidates_file: Optional[str] = None
    biomarker_gene_col: Optional[str] = None


def _strip_ensembl_version(series: pd.Series) -> pd.Series:
    return series.astype(str).str.replace(r"\.\d+$", "", regex=True)


def _sample_type_from_tcga_barcode(barcode: str) -> Optional[str]:
    if not isinstance(barcode, str):
        return None
    parts = barcode.split("-")
    if len(parts) < 4:
        return None
    sample_field = parts[3]
    if len(sample_field) < 2:
        return None
    return sample_field[:2]


def _participant_id_from_barcode(barcode: str) -> str:
    if not isinstance(barcode, str):
        return ""
    parts = barcode.split("-")
    return "-".join(parts[:3]) if len(parts) >= 3 else barcode


def _sample_id_short(barcode: str) -> str:
    if not isinstance(barcode, str):
        return ""
    parts = barcode.split("-")
    return "-".join(parts[:4]) if len(parts) >= 4 else barcode


def _compute_metrics(y_true: np.ndarray, y_prob: np.ndarray, threshold: float = 0.5) -> Dict[str, float]:
    y_pred = (y_prob >= threshold).astype(int)
    return {
        "roc_auc": float(roc_auc_score(y_true, y_prob)),
        "pr_auc": float(average_precision_score(y_true, y_prob)),
        "f1": float(f1_score(y_true, y_pred, zero_division=0)),
        "balanced_accuracy": float(balanced_accuracy_score(y_true, y_pred)),
        "precision": float(precision_score(y_true, y_pred, zero_division=0)),
        "recall": float(recall_score(y_true, y_pred, zero_division=0)),
    }


def _find_best_threshold(y_true: np.ndarray, y_prob: np.ndarray) -> Tuple[float, float]:
    grid = np.linspace(0.10, 0.90, 81)
    best_threshold, best_f1 = 0.50, -1.0
    for t in grid:
        y_pred = (y_prob >= t).astype(int)
        score = f1_score(y_true, y_pred, zero_division=0)
        if score > best_f1:
            best_f1 = score
            best_threshold = float(t)
    return best_threshold, float(best_f1)


def _load_expression_matrix(expression_path: Path) -> pd.DataFrame:
    """Load TCGA-COAD STAR counts as a genes x samples matrix.

    The project's star_counts.tsv has gene Ensembl IDs in the first
    column and sample TCGA barcodes as the remaining columns. We keep
    Ensembl IDs as features so they line up with the rest of the pipeline.
    """
    expr_df = pd.read_csv(expression_path, sep="\t", low_memory=False)
    cols = expr_df.columns.tolist()
    sample_cols = [c for c in cols if str(c).startswith("TCGA-")]
    if not sample_cols:
        raise ValueError(f"No TCGA sample columns found in {expression_path}")

    gene_id_col = cols[0]
    work_df = expr_df[[gene_id_col] + sample_cols].copy()
    work_df = work_df.rename(columns={gene_id_col: "gene_id"})
    work_df["gene_id"] = _strip_ensembl_version(work_df["gene_id"])

    numeric = work_df[sample_cols].apply(pd.to_numeric, errors="coerce")
    work_df = pd.concat([work_df[["gene_id"]], numeric], axis=1)
    work_df = work_df.groupby("gene_id", as_index=False)[sample_cols].mean()
    return work_df.set_index("gene_id")


def _cpm_log1p_transform(expr_genes_x_samples: pd.DataFrame) -> pd.DataFrame:
    """Convert raw counts to CPM and log1p.

    Note: the project's star_counts.tsv is already log2(count+1) per the
    Xena dataset metadata. Reapplying CPM/log1p on top is a no-op-ish
    monotone transform that preserves rank-based feature selection,
    matching Ayan's original behavior. We do not change the math.
    """
    counts = expr_genes_x_samples.copy()
    lib_sizes = counts.sum(axis=0).replace(0, np.nan)
    cpm = counts.divide(lib_sizes, axis=1) * 1_000_000.0
    cpm = cpm.fillna(0.0)
    return np.log1p(cpm)


def _derive_mmr_labels(phenotype_path: Path, barcode_col: str, label_col: str) -> pd.DataFrame:
    pheno = pd.read_csv(phenotype_path, sep="\t", low_memory=False)
    if barcode_col not in pheno.columns:
        raise ValueError(f"Barcode column '{barcode_col}' not found in phenotype file.")
    if label_col not in pheno.columns:
        raise ValueError(f"Label column '{label_col}' not found in phenotype file.")

    out = pheno[[barcode_col, label_col]].copy()
    out.columns = ["sample_barcode_raw", "label_source"]
    out["sample_barcode_raw"] = out["sample_barcode_raw"].astype(str).str.strip()
    out["label_source"] = out["label_source"].astype(str).str.strip().str.upper()

    label_map = {"YES": 1, "NO": 0, "POSITIVE": 1, "NEGATIVE": 0}
    out = out[out["label_source"].isin(label_map)].copy()
    out["target"] = out["label_source"].map(label_map)
    out["participant_id"] = out["sample_barcode_raw"].map(_participant_id_from_barcode)
    out["sample_id_short"] = out["sample_barcode_raw"].map(_sample_id_short)
    return out.drop_duplicates(subset=["participant_id"], keep="first")


def _build_modeling_table(
    expr_genes_x_samples: pd.DataFrame,
    labels_df: pd.DataFrame,
    max_samples: Optional[int] = None,
) -> Tuple[pd.DataFrame, pd.Series, pd.DataFrame]:
    all_samples = list(expr_genes_x_samples.columns)
    print(f"[{TAG}] Total expression sample columns: {len(all_samples)}")

    initial_meta = pd.DataFrame({"sample_barcode_full": all_samples})
    initial_meta["sample_type_code"] = initial_meta["sample_barcode_full"].map(
        _sample_type_from_tcga_barcode
    )
    initial_meta["participant_id"] = initial_meta["sample_barcode_full"].map(
        _participant_id_from_barcode
    )

    print(f"[{TAG}] Sample type counts:")
    print(initial_meta["sample_type_code"].value_counts(dropna=False).sort_index())

    sample_meta = initial_meta[initial_meta["sample_type_code"] == "01"].copy()
    print(f"[{TAG}] Primary tumor samples: {len(sample_meta)}")

    sample_meta = sample_meta.merge(
        labels_df[["participant_id", "label_source", "target"]],
        on="participant_id",
        how="inner",
    )
    print(f"[{TAG}] After phenotype join: {len(sample_meta)} rows")

    sample_meta = sample_meta.drop_duplicates(subset=["participant_id"], keep="first")
    print(f"[{TAG}] Class counts after alignment:")
    print(sample_meta["target"].value_counts().sort_index())

    if max_samples is not None:
        sample_meta = sample_meta.iloc[:max_samples].copy()

    selected = sample_meta["sample_barcode_full"].tolist()
    X = expr_genes_x_samples[selected].T.copy()
    y = sample_meta.set_index("sample_barcode_full").loc[X.index, "target"].astype(int)

    sample_meta = (
        sample_meta.set_index("sample_barcode_full")
        .loc[X.index]
        .reset_index()
        .rename(columns={"index": "sample_barcode_full"})
    )
    return X, y, sample_meta


def _filter_low_information_genes(
    X_train: pd.DataFrame, X_test: pd.DataFrame, min_nonzero_fraction: float
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    keep_mask = (X_train > 0).mean(axis=0) >= min_nonzero_fraction
    cols = X_train.columns[keep_mask]
    return X_train[cols].copy(), X_test[cols].copy()


def _make_models(config: MMRConfig) -> Dict[str, object]:
    models = {
        "logistic_regression": LogisticRegression(
            class_weight="balanced",
            max_iter=3000,
            solver="liblinear",
            random_state=config.random_state,
        ),
        "random_forest": RandomForestClassifier(
            n_estimators=500,
            max_depth=None,
            min_samples_leaf=2,
            class_weight="balanced",
            random_state=config.random_state,
            n_jobs=config.n_jobs,
        ),
    }
    if XGBOOST_AVAILABLE:
        models["xgboost"] = XGBClassifier(
            n_estimators=400,
            max_depth=4,
            learning_rate=0.05,
            subsample=0.9,
            colsample_bytree=0.8,
            reg_lambda=1.0,
            objective="binary:logistic",
            eval_metric="logloss",
            random_state=config.random_state,
            n_jobs=config.n_jobs,
        )
    return models


def _make_pipeline(model_name: str, model: object, top_k_features: int, variance_threshold: float) -> Pipeline:
    steps = [
        ("imputer", SimpleImputer(strategy="median")),
        ("variance", VarianceThreshold(threshold=variance_threshold)),
        ("select", SelectKBest(score_func=f_classif, k=top_k_features)),
    ]
    if model_name == "logistic_regression":
        steps.append(("scaler", StandardScaler()))
    steps.append(("model", model))
    return Pipeline(steps)


def _get_cv_results(X: pd.DataFrame, y: pd.Series, models: Dict[str, object], config: MMRConfig) -> pd.DataFrame:
    scoring = {
        "roc_auc": "roc_auc",
        "average_precision": "average_precision",
        "f1": "f1",
        "balanced_accuracy": "balanced_accuracy",
    }
    cv = StratifiedKFold(
        n_splits=config.cv_folds, shuffle=True, random_state=config.random_state
    )
    rows = []
    for model_name, model in models.items():
        print(f"[{TAG}] CV: {model_name}")
        pipe = _make_pipeline(
            model_name=model_name,
            model=model,
            top_k_features=min(config.top_k_features, X.shape[1]),
            variance_threshold=config.variance_threshold,
        )
        scores = cross_validate(
            pipe, X, y, cv=cv, scoring=scoring, n_jobs=1, return_train_score=False
        )
        row = {"model": model_name}
        for metric_name, values in scores.items():
            if metric_name.startswith("test_"):
                clean = metric_name.replace("test_", "")
                row[f"{clean}_mean"] = float(np.mean(values))
                row[f"{clean}_std"] = float(np.std(values))
        rows.append(row)
    return pd.DataFrame(rows).sort_values("roc_auc_mean", ascending=False).reset_index(drop=True)


def _train_and_evaluate_holdout(
    X: pd.DataFrame, y: pd.Series, models: Dict[str, object], config: MMRConfig, output_dir: Path
) -> Tuple[pd.DataFrame, Dict[str, Pipeline], pd.DataFrame]:
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=config.test_size, stratify=y, random_state=config.random_state
    )
    print(f"[{TAG}] Holdout train: {X_train.shape}  test: {X_test.shape}")

    X_train, X_test = _filter_low_information_genes(
        X_train, X_test, min_nonzero_fraction=config.min_nonzero_fraction
    )
    print(f"[{TAG}] After low-information filter: train {X_train.shape}  test {X_test.shape}")

    if X_train.shape[1] == 0:
        raise ValueError(
            "No genes remain after low-information filtering. "
            "Lower min_nonzero_fraction in config['mmr']."
        )

    results = []
    prediction_frames = []
    fitted_pipelines: Dict[str, Pipeline] = {}

    for model_name, model in models.items():
        print(f"[{TAG}] Training {model_name}")
        pipe = _make_pipeline(
            model_name=model_name,
            model=model,
            top_k_features=min(config.top_k_features, X_train.shape[1]),
            variance_threshold=config.variance_threshold,
        )
        pipe.fit(X_train, y_train)
        y_prob = pipe.predict_proba(X_test)[:, 1]

        default = _compute_metrics(y_test.values, y_prob, threshold=0.50)
        best_t, _ = _find_best_threshold(y_test.values, y_prob)
        tuned = _compute_metrics(y_test.values, y_prob, threshold=best_t)

        results.append({
            "model": model_name,
            "default_threshold": 0.50,
            "best_threshold": best_t,
            "default_f1": default["f1"],
            "tuned_f1": tuned["f1"],
            "roc_auc": tuned["roc_auc"],
            "pr_auc": tuned["pr_auc"],
            "f1": tuned["f1"],
            "balanced_accuracy": tuned["balanced_accuracy"],
            "precision": tuned["precision"],
            "recall": tuned["recall"],
        })

        prediction_frames.append(pd.DataFrame({
            "sample_id": X_test.index,
            "y_true": y_test.values,
            "y_prob": y_prob,
            "y_pred_default_0_50": (y_prob >= 0.50).astype(int),
            "y_pred_tuned": (y_prob >= best_t).astype(int),
            "model": model_name,
            "best_threshold": best_t,
        }))
        fitted_pipelines[model_name] = pipe

    holdout_df = pd.DataFrame(results).sort_values("roc_auc", ascending=False).reset_index(drop=True)
    predictions_df = pd.concat(prediction_frames, ignore_index=True)

    joblib.dump(
        {"X_train": X_train, "X_test": X_test, "y_train": y_train, "y_test": y_test},
        output_dir / "train_test_split.joblib",
    )
    return holdout_df, fitted_pipelines, predictions_df


def _get_selected_feature_names(fitted_pipeline: Pipeline, original_feature_names: List[str]) -> List[str]:
    feature_names = np.array(original_feature_names)
    feature_names = feature_names[fitted_pipeline.named_steps["variance"].get_support()]
    feature_names = feature_names[fitted_pipeline.named_steps["select"].get_support()]
    return feature_names.tolist()


def _coerce_numeric_array(x):
    arr = np.asarray(x)
    if np.issubdtype(arr.dtype, np.number):
        return arr.astype(float)
    flat = arr.reshape(-1)
    cleaned = [float(str(v).strip().strip("[]")) for v in flat]
    return np.array(cleaned, dtype=float).reshape(arr.shape)


def _run_shap_analysis(
    fitted_pipeline: Pipeline,
    X_reference: pd.DataFrame,
    X_explain: pd.DataFrame,
    output_dir: Path,
) -> Optional[pd.DataFrame]:
    if not SHAP_AVAILABLE:
        print(f"[{TAG}] SHAP not installed; skipping SHAP analysis.")
        return None

    model = fitted_pipeline.named_steps["model"]
    preprocess = Pipeline(fitted_pipeline.steps[:-1])
    X_ref_t = _coerce_numeric_array(preprocess.transform(X_reference))
    X_exp_t = _coerce_numeric_array(preprocess.transform(X_explain))

    selected = _get_selected_feature_names(fitted_pipeline, X_reference.columns.tolist())
    max_explain = min(100, X_exp_t.shape[0])
    X_exp_t_small = X_exp_t[:max_explain]

    try:
        if model.__class__.__name__ in {"RandomForestClassifier", "XGBClassifier"}:
            explainer = shap.TreeExplainer(model)
            shap_values = explainer.shap_values(X_exp_t_small)
            if isinstance(shap_values, list):
                shap_matrix = shap_values[1]
            else:
                shap_matrix = np.asarray(shap_values)
                if shap_matrix.ndim == 3:
                    shap_matrix = shap_matrix[:, :, 1]
        else:
            explainer = shap.Explainer(model, X_ref_t)
            shap_obj = explainer(X_exp_t_small)
            shap_matrix = shap_obj.values
            if shap_matrix.ndim == 3:
                shap_matrix = shap_matrix[:, :, 1]

        shap_matrix = _coerce_numeric_array(shap_matrix)
        mean_abs_shap = np.abs(shap_matrix).mean(axis=0)

        shap_df = pd.DataFrame({
            "feature": selected,
            "mean_abs_shap": mean_abs_shap,
        }).sort_values("mean_abs_shap", ascending=False)

        shap_df.to_csv(output_dir / "shap_feature_importance.csv", index=False)
        shap_df.head(100).to_csv(output_dir / "shap_top100_features.csv", index=False)
        return shap_df
    except Exception as e:
        print(f"[{TAG}] SHAP failed: {e}")
        return None


def _select_best_model(holdout_df: pd.DataFrame, primary: str = "f1", secondary: str = "roc_auc") -> str:
    ranked = holdout_df.sort_values(by=[primary, secondary], ascending=[False, False]).reset_index(drop=True)
    return ranked.iloc[0]["model"]


def _build_mmr_config(config_path: str) -> MMRConfig:
    config = load_config(config_path)
    mmr_block = config.get("mmr", {})
    if "phenotype_path" not in mmr_block:
        raise ValueError(
            "config['mmr']['phenotype_path'] is required. "
            "Add an 'mmr' block to config/config.yaml pointing to the GDC phenotype TSV."
        )
    return MMRConfig(
        expression_path=mmr_block.get("expression_path", config["omics"]["counts_path"]),
        phenotype_path=mmr_block["phenotype_path"],
        output_dir=mmr_block.get("output_dir", "outputs/coad_mmr_classifier"),
        phenotype_barcode_col=mmr_block.get(
            "phenotype_barcode_col", "submitter_id.samples"
        ),
        phenotype_label_col=mmr_block.get(
            "phenotype_label_col", "loss_expression_of_mismatch_repair_proteins_by_ihc"
        ),
        test_size=float(mmr_block.get("test_size", 0.20)),
        random_state=int(mmr_block.get("random_state", 42)),
        cv_folds=int(mmr_block.get("cv_folds", 5)),
        min_nonzero_fraction=float(mmr_block.get("min_nonzero_fraction", 0.10)),
        variance_threshold=float(mmr_block.get("variance_threshold", 0.0)),
        top_k_features=int(mmr_block.get("top_k_features", 1000)),
        n_jobs=int(mmr_block.get("n_jobs", -1)),
        max_samples=mmr_block.get("max_samples"),
        biomarker_candidates_file=mmr_block.get("biomarker_candidates_file", "outputs/ranked_candidates.csv"),
        biomarker_gene_col=mmr_block.get("biomarker_gene_col", "gene"),
    )


def run_mmr_classifier(config_path: str) -> Dict[str, object]:
    """Train and evaluate the within-COAD MMR classifier.

    Args:
        config_path: Path to config/config.yaml.

    Returns:
        A summary dict containing class counts, the best model name, and
        the paths of the artifacts written under outputs/coad_mmr_classifier/.
    """
    cfg = _build_mmr_config(config_path)

    expression_path = Path(cfg.expression_path)
    phenotype_path = Path(cfg.phenotype_path)
    output_dir = Path(cfg.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if not expression_path.exists():
        raise FileNotFoundError(f"Expression file not found: {expression_path}")
    if not phenotype_path.exists():
        raise FileNotFoundError(
            f"Phenotype file not found: {phenotype_path}. "
            "Download TCGA-COAD.GDC_phenotype.tsv from the GDC/Xena hub "
            "and place it in the data/ directory."
        )

    print(f"[{TAG}] Loading expression: {expression_path}")
    expr_raw = _load_expression_matrix(expression_path)
    print(f"[{TAG}] Expression matrix (genes x samples): {expr_raw.shape}")

    print(f"[{TAG}] Applying CPM + log1p transform")
    expr_norm = _cpm_log1p_transform(expr_raw)

    print(f"[{TAG}] Deriving labels from {phenotype_path}")
    labels_df = _derive_mmr_labels(
        phenotype_path=phenotype_path,
        barcode_col=cfg.phenotype_barcode_col,
        label_col=cfg.phenotype_label_col,
    )
    print(f"[{TAG}] Label counts: {labels_df['target'].value_counts().sort_index().to_dict()}")

    X, y, sample_meta = _build_modeling_table(
        expr_genes_x_samples=expr_norm,
        labels_df=labels_df,
        max_samples=cfg.max_samples,
    )
    print(f"[{TAG}] Modeling table (samples x genes): {X.shape}")

    if X.shape[0] == 0:
        raise ValueError("No aligned samples found between expression and phenotype.")
    if y.nunique() < 2:
        raise ValueError("Need at least 2 classes after alignment to train classifier.")

    X.to_csv(output_dir / "X_modeling_table.csv")
    pd.DataFrame({"sample_id": y.index, "target": y.values}).to_csv(
        output_dir / "y_labels.csv", index=False
    )
    sample_meta.to_csv(output_dir / "sample_metadata.csv", index=False)

    models = _make_models(cfg)
    if not XGBOOST_AVAILABLE:
        print(f"[{TAG}] xgboost not installed; running LR + RF only.")

    print(f"[{TAG}] Running cross-validation")
    cv_df = _get_cv_results(X=X, y=y, models=models, config=cfg)
    cv_df.to_csv(output_dir / "cv_results.csv", index=False)
    print(f"[{TAG}] CV results:\n{cv_df}")

    print(f"[{TAG}] Training holdout models")
    holdout_df, fitted_pipelines, predictions_df = _train_and_evaluate_holdout(
        X=X, y=y, models=models, config=cfg, output_dir=output_dir
    )
    holdout_df.to_csv(output_dir / "holdout_results.csv", index=False)
    predictions_df.to_csv(output_dir / "holdout_predictions.csv", index=False)
    print(f"[{TAG}] Holdout results:\n{holdout_df}")

    best_model_name = _select_best_model(holdout_df, primary="f1", secondary="roc_auc")
    best_pipeline = fitted_pipelines[best_model_name]
    print(f"[{TAG}] Best model for explanation: {best_model_name}")
    joblib.dump(best_pipeline, output_dir / f"best_model_{best_model_name}.joblib")

    final_pipeline = _make_pipeline(
        model_name=best_model_name,
        model=clone(models[best_model_name]),
        top_k_features=min(cfg.top_k_features, X.shape[1]),
        variance_threshold=cfg.variance_threshold,
    )
    final_pipeline.fit(X, y)
    joblib.dump(final_pipeline, output_dir / f"final_model_full_data_{best_model_name}.joblib")

    print(f"[{TAG}] Running SHAP on best holdout model")
    split_obj = joblib.load(output_dir / "train_test_split.joblib")
    shap_df = _run_shap_analysis(
        fitted_pipeline=best_pipeline,
        X_reference=split_obj["X_train"],
        X_explain=split_obj["X_test"],
        output_dir=output_dir,
    )

    summary = {
        "config": asdict(cfg),
        "n_samples": int(X.shape[0]),
        "n_genes": int(X.shape[1]),
        "class_counts": {
            str(k): int(v) for k, v in y.value_counts().sort_index().to_dict().items()
        },
        "best_model": best_model_name,
        "cv_results_path": str(output_dir / "cv_results.csv"),
        "holdout_results_path": str(output_dir / "holdout_results.csv"),
        "predictions_path": str(output_dir / "holdout_predictions.csv"),
        "shap_available": bool(shap_df is not None),
    }
    with open(output_dir / "run_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    print(f"[{TAG}] Done. Summary written to {output_dir / 'run_summary.json'}")
    return summary


if __name__ == "__main__":
    import sys
    cfg = sys.argv[1] if len(sys.argv) > 1 else "config/config.yaml"
    run_mmr_classifier(cfg)
