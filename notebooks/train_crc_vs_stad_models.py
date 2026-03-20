"""
train_crc_vs_stad_models.py

Harder sample-level classification:
    COAD tumor vs STAD tumor

Compares:
    - Logistic Regression
    - Random Forest
    - XGBoost (if available)

Inputs
------
data/raw/TCGA-COAD.star_counts.tsv
data/raw/HiSeqV2.tsv   (or TCGA-STAD_xena.tsv / TCGA-STAD.star_counts.tsv)
data/ensembl_to_symbol_cache.tsv

Outputs
-------
outputs/crc_vs_stad_model_comparison.csv
"""

from __future__ import annotations

from pathlib import Path
import re

import pandas as pd

from sklearn.model_selection import StratifiedKFold, cross_val_predict
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


DEFAULT_PROJECT_ROOT = Path(r"D:\Ayan\PSU Capstone Project\AI_Capstone")


def normalize_ensembl_gene_id(val):
    """
    Normalize Ensembl gene IDs to canonical ENSG + 11 digits.

    Examples:
        ENSG00000123456.12 -> ENSG00000123456
        ENSG0000012345612  -> ENSG00000123456
    """
    if pd.isna(val):
        return None

    s = str(val).strip().upper()
    m = re.search(r"(ENSG\d{11})", s)
    return m.group(1) if m else s


def load_ensembl_mapping(mapping_path: Path) -> pd.DataFrame:
    """
    Load Ensembl -> gene symbol mapping.
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
    out.columns = ["gene_id", "gene_symbol"]

    out["gene_id"] = out["gene_id"].apply(normalize_ensembl_gene_id)
    out["gene_symbol"] = out["gene_symbol"].astype(str).str.strip().str.upper()

    out = out.dropna(subset=["gene_id", "gene_symbol"])
    out = out[out["gene_symbol"] != ""]
    out = out.drop_duplicates(subset=["gene_id"])

    return out


def get_tcga_sample_type(barcode: str):
    """
    TCGA sample type:
      01 = primary tumor
      11 = normal
    """
    parts = str(barcode).split("-")
    if len(parts) < 4:
        return None

    code = parts[3][:2]

    if code == "01":
        return "tumor"
    if code == "11":
        return "normal"

    return None


def detect_id_column(df: pd.DataFrame) -> str:
    """
    Detect the identifier column in an expression matrix.
    """
    candidates = [
        "Ensembl_ID",
        "ensembl_id",
        "gene_id",
        "Gene_ID",
        "sample",   # Xena STAD file uses this for gene symbols
        df.columns[0],
    ]

    for c in candidates:
        if c in df.columns:
            return c

    raise ValueError("Could not detect identifier column.")


def read_expression_file(path: Path) -> pd.DataFrame:
    """
    Read a TSV expression matrix.
    """
    return pd.read_csv(path, sep="\t", low_memory=False)


def convert_to_sample_matrix(
    expr_df: pd.DataFrame,
    mapping_df: pd.DataFrame,
    cohort_name: str,
) -> pd.DataFrame:
    """
    Convert gene x sample matrix to sample x gene matrix.
    Keep tumor samples only.

    Handles either:
    - Ensembl IDs in first column
    - gene symbols in first column
    """
    id_col = detect_id_column(expr_df)
    expr_df = expr_df.rename(columns={id_col: "feature_id"}).copy()

    # Keep only TCGA sample columns
    sample_cols = [c for c in expr_df.columns if str(c).startswith("TCGA-")]
    expr_df = expr_df[["feature_id"] + sample_cols].copy()

    # Ensure numeric
    expr_df[sample_cols] = expr_df[sample_cols].apply(pd.to_numeric, errors="coerce")

    # Detect identifier type from first column
    preview = expr_df["feature_id"].astype(str).head(20).tolist()
    ensembl_like = sum(bool(re.search(r"ENSG\d{11}", x.upper())) for x in preview) >= 5

    if ensembl_like:
        expr_df["feature_id"] = expr_df["feature_id"].apply(normalize_ensembl_gene_id)
        expr_df = expr_df.rename(columns={"feature_id": "gene_id"})
        expr_df = expr_df.merge(mapping_df, on="gene_id", how="left")
        expr_df = expr_df.dropna(subset=["gene_symbol"]).copy()
    else:
        # Already gene symbols
        expr_df["gene_symbol"] = (
            expr_df["feature_id"]
            .astype(str)
            .str.strip()
            .str.upper()
        )
        expr_df = expr_df.drop(columns=["feature_id"])

    # Average duplicate symbols if any
    expr_df = expr_df.groupby("gene_symbol", as_index=False)[sample_cols].mean()

    # Transpose to samples x genes
    expr_df = expr_df.set_index("gene_symbol")
    X = expr_df.T.copy()

    # Keep tumor samples only
    X["sample_type"] = [get_tcga_sample_type(idx) for idx in X.index]
    X = X[X["sample_type"] == "tumor"].copy()
    X = X.drop(columns=["sample_type"])

    X["cohort"] = cohort_name
    return X


def resolve_stad_file(root: Path) -> Path:
    """
    Try common STAD file locations/names.
    """
    candidates = [
        root / "data" / "raw" / "TCGA-STAD_xena.tsv",
        root / "data" / "raw" / "TCGA-STAD.star_counts.tsv",
        root / "data" / "raw" / "HiSeqV2.tsv",
        root / "data" / "raw" / "HiSeqV2",
    ]

    for p in candidates:
        if p.exists():
            return p

    raise FileNotFoundError(
        "Could not find STAD expression file. Expected one of:\n"
        + "\n".join(str(p) for p in candidates)
    )


def percentile_rank_transform(X: pd.DataFrame) -> pd.DataFrame:
    """
    Reduce platform/source effects by converting each sample's expression
    values to within-sample percentile ranks.
    Output range: 0..1
    """
    return X.rank(axis=1, pct=True)


def select_top_variable_genes(X: pd.DataFrame, n_genes: int) -> pd.DataFrame:
    """
    Select top n most variable genes.
    """
    variances = X.var(axis=0)
    top_genes = variances.sort_values(ascending=False).head(n_genes).index.tolist()
    return X[top_genes].copy()


def build_models():
    """
    Build comparison models.
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


def evaluate_oof_predictions(y_true, prob):
    """
    Compute metrics from out-of-fold probabilities.
    """
    pred = (prob >= 0.5).astype(int)

    return {
        "roc_auc": roc_auc_score(y_true, prob),
        "pr_auc": average_precision_score(y_true, prob),
        "accuracy": accuracy_score(y_true, pred),
        "precision": precision_score(y_true, pred, zero_division=0),
        "recall": recall_score(y_true, pred, zero_division=0),
        "f1": f1_score(y_true, pred, zero_division=0),
    }


def run_experiment(experiment_name, X, y, models):
    """
    Run 5-fold stratified CV for all models.
    """
    print(f"\n========== {experiment_name} ==========")
    print("Feature matrix:", X.shape)
    print("Class counts:")
    print(y.value_counts())

    if X.shape[1] == 0:
        raise ValueError(f"No features available for experiment: {experiment_name}")

    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

    rows = []

    for model_name, model in models.items():
        print(f"\nRunning {model_name}...")
        prob = cross_val_predict(
            model,
            X,
            y,
            cv=cv,
            method="predict_proba",
            n_jobs=1
        )[:, 1]

        metrics = evaluate_oof_predictions(y, prob)

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


def train_crc_vs_stad_models(project_root):
    """
    Main runner.
    """
    root = Path(project_root)

    coad_path = root / "data" / "raw" / "TCGA-COAD.star_counts.tsv"
    stad_path = resolve_stad_file(root)
    mapping_path = root / "data" / "ensembl_to_symbol_cache.tsv"
    output_path = root / "outputs" / "crc_vs_stad_model_comparison.csv"

    print("Loading gene mapping...")
    mapping_df = load_ensembl_mapping(mapping_path)
    print("Mapping rows:", len(mapping_df))

    print("\nLoading COAD...")
    coad_raw = read_expression_file(coad_path)
    print("COAD raw shape:", coad_raw.shape)

    print("\nLoading STAD...")
    stad_raw = read_expression_file(stad_path)
    print("STAD raw shape:", stad_raw.shape)
    print("Using STAD file:", stad_path)

    coad = convert_to_sample_matrix(coad_raw, mapping_df, "COAD")
    stad = convert_to_sample_matrix(stad_raw, mapping_df, "STAD")

    print("\nCOAD tumor samples:", len(coad))
    print("STAD tumor samples:", len(stad))

    print("\nExample COAD columns:", list(coad.columns[:10]))
    print("Example STAD columns:", list(stad.columns[:10]))

    # Align common genes
    common_genes = sorted((set(coad.columns) & set(stad.columns)) - {"cohort"})
    print("Common genes:", len(common_genes))

    if len(common_genes) == 0:
        raise ValueError(
            "No common genes found between COAD and STAD after preprocessing."
        )

    coad = coad[common_genes].copy()
    stad = stad[common_genes].copy()

    # Combine
    X = pd.concat([coad, stad], axis=0)
    y = pd.Series(
        [1] * len(coad) + [0] * len(stad),
        index=X.index,
        name="crc_label"
    )

    print("\nCombined sample matrix:", X.shape)
    print("Combined class distribution:")
    print(y.value_counts())

    # Within-sample percentile transform to reduce source/platform effects
    print("\nApplying within-sample percentile-rank transform...")
    X_rank = percentile_rank_transform(X)

    models = build_models()
    results = []

    # Experiment A: all common genes
    results.extend(run_experiment(
        experiment_name="COAD_vs_STAD_AllCommonGenes_RankTransformed",
        X=X_rank,
        y=y,
        models=models
    ))

    # Experiment B: top 1000 variable genes
    X_1000 = select_top_variable_genes(X_rank, min(1000, X_rank.shape[1]))
    results.extend(run_experiment(
        experiment_name="COAD_vs_STAD_Top1000VarGenes_RankTransformed",
        X=X_1000,
        y=y,
        models=models
    ))

    # Experiment C: top 500 variable genes
    X_500 = select_top_variable_genes(X_rank, min(500, X_rank.shape[1]))
    results.extend(run_experiment(
        experiment_name="COAD_vs_STAD_Top500VarGenes_RankTransformed",
        X=X_500,
        y=y,
        models=models
    ))

    results_df = pd.DataFrame(results).sort_values(
        by=["experiment", "roc_auc", "f1"],
        ascending=[True, False, False]
    ).reset_index(drop=True)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    results_df.to_csv(output_path, index=False)

    print("\nSaved results to:")
    print(output_path)

    return results_df


def run_train_crc_vs_stad_models(project_root):
    return train_crc_vs_stad_models(project_root)


if __name__ == "__main__":
    run_train_crc_vs_stad_models(DEFAULT_PROJECT_ROOT)