"""
build_classifier_dataset.py

Module-7D
Assemble pipeline evidence outputs into a per-gene feature matrix for ML training.

Joins omics, literature, pathway, LLM, and scoring outputs on gene (Ensembl ID)
to produce a single tabular dataset consumed by build_ml_dataset_curated.py.

Inputs
------
outputs/omics_evidence.csv
outputs/ranked_candidates.csv
outputs/lit_evidence.csv       (optional)
outputs/pathway_evidence.csv   (optional)
outputs/llm_scores.csv         (optional)

Outputs
-------
data/processed/classifier_dataset.csv
data/processed/classifier_dataset_summary.txt
"""

from pathlib import Path

import numpy as np
import pandas as pd


DEFAULT_PROJECT_ROOT = Path(__file__).resolve().parent.parent


def build_classifier_dataset(project_root):
    """Assemble pipeline evidence outputs into a per-gene feature matrix.

    Parameters
    ----------
    project_root : str or Path
        Project root directory.

    Returns
    -------
    pd.DataFrame
        Per-gene feature matrix saved to data/processed/classifier_dataset.csv.
    """
    project_root = Path(project_root)
    outputs_dir = project_root / "outputs"
    processed_dir = project_root / "data" / "processed"
    processed_dir.mkdir(parents=True, exist_ok=True)

    output_path = processed_dir / "classifier_dataset.csv"
    summary_path = processed_dir / "classifier_dataset_summary.txt"

    # ------------------------------------------------------------------ #
    # 1. Omics evidence — base gene universe
    # ------------------------------------------------------------------ #
    omics_path = outputs_dir / "omics_evidence.csv"
    if not omics_path.exists():
        raise FileNotFoundError(
            f"omics_evidence.csv not found at {omics_path}. "
            "Run the omics pipeline step first."
        )

    print("\nLoading omics evidence...")
    omics = pd.read_csv(omics_path)
    print(f"  Omics genes: {len(omics)}")

    # Normalize versioned Ensembl IDs (ENSG00000141510.15 -> ENSG00000141510)
    omics["gene"] = omics["gene"].str.split(".").str[0]

    # Derived omics features
    epsilon = 1e-300
    omics["max_abs_logfc"] = omics["log2fc"].abs()
    omics["significance_strength"] = -np.log10(omics["fdr"].clip(lower=epsilon))

    # ------------------------------------------------------------------ #
    # 2. Ranked candidates — composite scores
    # ------------------------------------------------------------------ #
    ranked_path = outputs_dir / "ranked_candidates.csv"
    if not ranked_path.exists():
        raise FileNotFoundError(
            f"ranked_candidates.csv not found at {ranked_path}. "
            "Run the scoring pipeline step first."
        )

    print("Loading ranked candidates...")
    ranked = pd.read_csv(ranked_path)
    ranked["gene"] = ranked["gene"].str.split(".").str[0]
    ranked = ranked.rename(columns={
        "final_score": "total_score_sum",
        "omics_score": "omics_score",
        "literature_score": "literature_score",
        "pathway_score": "pathway_score",
    })
    ranked["max_score"] = ranked[["omics_score", "literature_score", "pathway_score"]].max(axis=1)
    print(f"  Ranked genes: {len(ranked)}")

    # ------------------------------------------------------------------ #
    # 3. Pathway evidence — pathway count per gene
    # ------------------------------------------------------------------ #
    pathway_path = outputs_dir / "pathway_evidence.csv"
    pathway_total_map: dict = {}
    top_pathways_map: dict = {}
    if pathway_path.exists():
        print("Loading pathway evidence...")
        pathway = pd.read_csv(pathway_path)
        pathway["gene"] = pathway["gene"].str.split(".").str[0]
        pathway_total_map = dict(zip(pathway["gene"], pathway["pathway_count"]))
        if "top_pathways" in pathway.columns:
            top_pathways_map = dict(zip(pathway["gene"], pathway["top_pathways"]))
        print(f"  Pathway genes: {len(pathway)}")
    else:
        print("  pathway_evidence.csv not found — pathway features will be 0")

    # ------------------------------------------------------------------ #
    # 4. Literature evidence — paper count per gene
    # ------------------------------------------------------------------ #
    lit_path = outputs_dir / "lit_evidence.csv"
    n_sources_map: dict = {}
    if lit_path.exists():
        print("Loading literature evidence...")
        lit = pd.read_csv(lit_path)
        lit["gene"] = lit["gene"].str.split(".").str[0]
        n_sources_map = lit.groupby("gene").size().to_dict()
        print(f"  Literature genes with papers: {len(n_sources_map)}")
    else:
        print("  lit_evidence.csv not found — n_sources will be 0")

    # ------------------------------------------------------------------ #
    # 5. LLM scores
    # ------------------------------------------------------------------ #
    llm_path = outputs_dir / "llm_scores.csv"
    llm_score_map: dict = {}
    llm_rationale_map: dict = {}
    if llm_path.exists():
        print("Loading LLM scores...")
        llm = pd.read_csv(llm_path)
        llm["gene"] = llm["gene"].str.split(".").str[0]
        llm_score_map = dict(zip(llm["gene"], llm["llm_score"]))
        if "rationale" in llm.columns:
            llm_rationale_map = dict(zip(llm["gene"], llm["rationale"]))
        print(f"  LLM scored genes: {len(llm_score_map)}")
    else:
        print("  llm_scores.csv not found — llm_score will be 0")

    # ------------------------------------------------------------------ #
    # 6. Assemble feature matrix
    # ------------------------------------------------------------------ #
    print("\nAssembling classifier dataset...")

    # Start from omics (full gene universe).
    # Omit gene_symbol here — build_ml_dataset_curated.py adds it via its own
    # ensembl_to_symbol_cache merge, and a pre-existing gene_symbol column
    # would cause a collision that breaks that merge.
    df = omics[["gene", "log2fc", "p_value", "fdr",
                 "direction", "tumor_mean", "normal_mean",
                 "max_abs_logfc", "significance_strength"]].copy()

    # Join ranked scores
    ranked_cols = ["gene", "omics_score", "literature_score", "pathway_score",
                   "total_score_sum", "max_score"]
    df = df.merge(ranked[ranked_cols], on="gene", how="left")

    # Map evidence counts from optional sources
    df["pathway_total"] = df["gene"].map(pathway_total_map).fillna(0).astype(int)
    df["path_evidence_top_pathways"] = df["gene"].map(top_pathways_map)

    df["n_sources"] = df["gene"].map(n_sources_map).fillna(0).astype(int)

    df["llm_score"] = df["gene"].map(llm_score_map).fillna(0)
    df["lit_evidence_snippet"] = df["gene"].map(llm_rationale_map)

    # Fill any remaining NaN in numeric columns
    numeric_cols = ["omics_score", "literature_score", "pathway_score",
                    "total_score_sum", "max_score"]
    df[numeric_cols] = df[numeric_cols].fillna(0)

    print(f"Final classifier dataset shape: {df.shape}")
    print("\nClass-ready evidence column stats:")
    for col in ["n_sources", "max_abs_logfc", "significance_strength",
                "pathway_total", "total_score_sum", "max_score"]:
        if col in df.columns:
            nonzero = (df[col] > 0).sum()
            print(f"  {col}: nonzero={nonzero}, mean={df[col].mean():.3f}")

    # ------------------------------------------------------------------ #
    # 7. Save outputs
    # ------------------------------------------------------------------ #
    df.to_csv(output_path, index=False)
    print(f"\nSaved classifier dataset to: {output_path}")

    summary_lines = [
        "Classifier Dataset Summary",
        "==========================",
        f"Total genes: {len(df)}",
        f"Genes with literature (n_sources > 0): {(df['n_sources'] > 0).sum()}",
        f"Genes with pathway annotations: {(df['pathway_total'] > 0).sum()}",
        f"Genes with LLM scores: {(df['llm_score'] > 0).sum()}",
        f"Columns: {list(df.columns)}",
    ]
    summary_path.write_text("\n".join(summary_lines), encoding="utf-8")
    print(f"Saved summary to: {summary_path}")

    return df


def run_build_classifier_dataset(project_root):
    return build_classifier_dataset(project_root)


if __name__ == "__main__":
    run_build_classifier_dataset(DEFAULT_PROJECT_ROOT)
