"""Compare MMR-classifier SHAP top features against ranked biomarker candidates.

Original author: Ayan Choudhury (compare_shap_vs_biomarkers.py).
Adapted by Keith Coughlin to run from the project config and write its
outputs alongside the rest of the MMR classifier artifacts.

Loads the SHAP feature importance produced by the MMR classifier and
intersects the top-N features with the project's ranked_candidates.csv
output. Reports overlap against the full candidate pool and against the
top-N best-scoring biomarkers, so we can see whether the MMR-driving
genes are also rising in the omics-driven ranking.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, Tuple

import pandas as pd

from src.utils import load_config

TAG = "mmr_compare"


def _clean_gene_ids(series: pd.Series) -> pd.Series:
    return (
        series.astype(str)
        .str.strip()
        .str.replace(r"\.\d+$", "", regex=True)
    )


def _validate_columns(df: pd.DataFrame, required_cols: Iterable[str], df_name: str) -> None:
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise ValueError(
            f"Missing required column(s) {missing} in {df_name}. "
            f"Available columns: {list(df.columns)}"
        )


def _compute_overlap(
    shap_df: pd.DataFrame,
    biomarker_df: pd.DataFrame,
    biomarker_gene_col: str,
    top_n_shap: int,
    label: str,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    _validate_columns(shap_df, ["feature"], "shap_df")
    _validate_columns(biomarker_df, [biomarker_gene_col], "biomarker_df")

    shap_top = shap_df.head(top_n_shap).copy()
    shap_top["feature_clean"] = _clean_gene_ids(shap_top["feature"])

    biomarker_work = biomarker_df.copy()
    biomarker_work["biomarker_gene"] = biomarker_work[biomarker_gene_col].astype(str).str.strip()
    biomarker_work["biomarker_gene_clean"] = _clean_gene_ids(biomarker_work["biomarker_gene"])

    overlap_df = shap_top.merge(
        biomarker_work,
        left_on="feature_clean",
        right_on="biomarker_gene_clean",
        how="inner",
    )

    summary_df = pd.DataFrame({
        "comparison_type": [label],
        "top_n_shap": [top_n_shap],
        "n_shap_top_unique": [shap_top["feature_clean"].nunique()],
        "n_biomarkers_input": [biomarker_work["biomarker_gene_clean"].nunique()],
        "n_overlap": [overlap_df["feature_clean"].nunique()],
        "overlap_fraction_of_shap_top": [
            overlap_df["feature_clean"].nunique() / max(1, shap_top["feature_clean"].nunique())
        ],
    })
    return overlap_df, summary_df


def run_mmr_shap_vs_biomarkers(config_path: str) -> str:
    """Compare MMR-classifier SHAP features against ranked biomarker candidates.

    Args:
        config_path: Path to config/config.yaml.

    Returns:
        Path to the combined comparison summary CSV.
    """
    config = load_config(config_path)
    mmr_block = config.get("mmr", {})

    output_dir = Path(mmr_block.get("output_dir", "outputs/coad_mmr_classifier"))
    shap_path = output_dir / "shap_feature_importance.csv"

    biomarker_file = mmr_block.get("biomarker_candidates_file", "outputs/ranked_candidates.csv")
    biomarker_gene_col = mmr_block.get("biomarker_gene_col", "gene")
    biomarker_score_col = mmr_block.get("biomarker_score_col", "final_score")
    top_n_shap = int(mmr_block.get("compare_top_n_shap", 100))
    top_n_biomarkers = int(mmr_block.get("compare_top_n_biomarkers", 500))

    if not shap_path.exists():
        raise FileNotFoundError(
            f"SHAP file not found: {shap_path}. Run the MMR classifier first."
        )
    if not Path(biomarker_file).exists():
        raise FileNotFoundError(f"Biomarker candidates file not found: {biomarker_file}")

    print(f"[{TAG}] Loading SHAP features: {shap_path}")
    shap_df = pd.read_csv(shap_path)

    print(f"[{TAG}] Loading biomarker candidates: {biomarker_file}")
    biomarker_df = pd.read_csv(biomarker_file)
    _validate_columns(shap_df, ["feature"], "SHAP file")
    _validate_columns(biomarker_df, [biomarker_gene_col, biomarker_score_col], "biomarker file")

    overlap_all_df, summary_all_df = _compute_overlap(
        shap_df=shap_df,
        biomarker_df=biomarker_df,
        biomarker_gene_col=biomarker_gene_col,
        top_n_shap=top_n_shap,
        label="all_biomarker_candidates",
    )

    biomarker_top_df = (
        biomarker_df.sort_values(biomarker_score_col, ascending=False)
        .head(top_n_biomarkers)
        .copy()
    )
    overlap_top_df, summary_top_df = _compute_overlap(
        shap_df=shap_df,
        biomarker_df=biomarker_top_df,
        biomarker_gene_col=biomarker_gene_col,
        top_n_shap=top_n_shap,
        label=f"top_{top_n_biomarkers}_biomarkers",
    )

    combined_summary_df = pd.concat([summary_all_df, summary_top_df], ignore_index=True)

    overlap_all_path = output_dir / "shap_vs_all_biomarker_candidates_overlap.csv"
    overlap_top_path = output_dir / f"shap_vs_top_{top_n_biomarkers}_biomarkers_overlap.csv"
    summary_path = output_dir / "shap_vs_biomarker_comparison_summary.csv"

    overlap_all_df.to_csv(overlap_all_path, index=False)
    overlap_top_df.to_csv(overlap_top_path, index=False)
    combined_summary_df.to_csv(summary_path, index=False)

    print(f"[{TAG}] Wrote {overlap_all_path}")
    print(f"[{TAG}] Wrote {overlap_top_path}")
    print(f"[{TAG}] Wrote {summary_path}")
    print(f"[{TAG}] Combined summary:\n{combined_summary_df}")

    return str(summary_path)


if __name__ == "__main__":
    import sys
    cfg = sys.argv[1] if len(sys.argv) > 1 else "config/config.yaml"
    run_mmr_shap_vs_biomarkers(cfg)
