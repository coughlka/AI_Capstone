"""Scoring module for combining evidence and ranking candidates."""

import os

import numpy as np
import pandas as pd

from src.utils import load_config, ensure_dirs, write_csv, read_csv


def _min_max_normalize(series: pd.Series, target_min: float = 0, target_max: float = 100) -> pd.Series:
    """Normalize a series to a target range using min-max scaling.

    Handles constant values (min == max) by returning target_min.

    Args:
        series: Pandas Series to normalize.
        target_min: Minimum value of target range.
        target_max: Maximum value of target range.

    Returns:
        Normalized Pandas Series.
    """
    min_val = series.min()
    max_val = series.max()

    if min_val == max_val:
        # All values are the same; return target_min to avoid division by zero
        return pd.Series([target_min] * len(series), index=series.index)

    return (series - min_val) / (max_val - min_val) * (target_max - target_min) + target_min


def run_scoring(config_path: str) -> str:
    """Run the scoring and ranking pipeline step.

    Combines evidence from omics, literature, and pathway modules to compute
    a weighted final score for each gene.

    Omics scoring is based on differential expression signal:
    omics_signal = |log2FC| * -log10(FDR + epsilon)

    Args:
        config_path: Path to the configuration YAML file.

    Returns:
        Path to the output CSV file.

    Raises:
        FileNotFoundError: If any required upstream output file is missing.
    """
    print("[scoring] Loading configuration...")
    config = load_config(config_path)
    ensure_dirs(config)

    outputs_dir = config['paths']['outputs_dir']
    output_path = os.path.join(outputs_dir, 'ranked_candidates.csv')

    # Define required input files
    omics_path = os.path.join(outputs_dir, 'omics_evidence.csv')
    lit_path = os.path.join(outputs_dir, 'lit_evidence.csv')
    pathway_path = os.path.join(outputs_dir, 'pathway_evidence.csv')

    # Check for missing upstream outputs
    missing = []
    if not os.path.exists(omics_path):
        missing.append("omics_evidence.csv (run omics step first)")
    if not os.path.exists(lit_path):
        missing.append("lit_evidence.csv (run pubmed step first)")
    if not os.path.exists(pathway_path):
        missing.append("pathway_evidence.csv (run pathway step first)")

    if missing:
        raise FileNotFoundError(
            f"Missing required upstream outputs in {outputs_dir}/:\n  - " +
            "\n  - ".join(missing)
        )

    # Load upstream outputs
    print(f"[scoring] Reading omics evidence from: {omics_path}")
    omics_df = read_csv(omics_path)

    print(f"[scoring] Reading literature evidence from: {lit_path}")
    lit_df = read_csv(lit_path)

    print(f"[scoring] Reading pathway evidence from: {pathway_path}")
    pathway_df = read_csv(pathway_path)

    # Get weights from config
    weights = config.get('scoring', {}).get('weights', {})
    w_omics = weights.get('omics', 0.45)
    w_lit = weights.get('literature', 0.35)
    w_path = weights.get('pathway', 0.20)

    print(f"[scoring] Using weights: omics={w_omics}, literature={w_lit}, pathway={w_path}")

    # Start with omics data as the base (gene list)
    if omics_df.empty:
        print("[scoring] Warning: omics evidence is empty. Creating empty output.")
        ranked = pd.DataFrame(columns=[
            'gene', 'final_score', 'omics_score', 'literature_score', 'pathway_score'
        ])
        write_csv(ranked, output_path)
        return output_path

    # Compute omics_score based on differential expression signal
    # omics_signal = |log2FC| * -log10(FDR + epsilon)
    print("[scoring] Computing omics scores from DE signal...")

    # Check if we have the new DE columns or old mean_expr column
    if 'log2fc' in omics_df.columns and 'fdr' in omics_df.columns:
        # New DE-based scoring
        epsilon = 1e-300  # Prevent log10(0)
        omics_df['omics_signal'] = (
            omics_df['log2fc'].abs() *
            (-np.log10(omics_df['fdr'] + epsilon))
        )
        scored = omics_df[['gene']].copy()
        scored['omics_score'] = _min_max_normalize(omics_df['omics_signal'])
        print("[scoring] Using DE-based scoring: |log2FC| * -log10(FDR)")
    elif 'mean_expr' in omics_df.columns:
        # Fallback to old mean expression scoring
        scored = omics_df[['gene']].copy()
        scored['omics_score'] = _min_max_normalize(omics_df['mean_expr'])
        print("[scoring] Warning: Using legacy mean expression scoring")
    else:
        # No valid scoring columns found
        print("[scoring] Warning: No valid omics scoring columns found. Using zeros.")
        scored = omics_df[['gene']].copy()
        scored['omics_score'] = 0.0

    # Compute literature_score: count rows per gene in lit_evidence, normalize to 0-100
    print("[scoring] Computing literature scores...")
    if lit_df.empty:
        scored['literature_score'] = 0.0
    else:
        lit_counts = lit_df.groupby('gene').size().reset_index(name='lit_count')
        scored = scored.merge(lit_counts, on='gene', how='left')
        scored['lit_count'] = scored['lit_count'].fillna(0)
        scored['literature_score'] = _min_max_normalize(scored['lit_count'])
        scored = scored.drop(columns=['lit_count'])

    # Compute pathway_score: use pathway_count, normalize to 0-100
    print("[scoring] Computing pathway scores...")
    if pathway_df.empty:
        scored['pathway_score'] = 0.0
    else:
        pathway_counts = pathway_df[['gene', 'pathway_count']].copy()
        scored = scored.merge(pathway_counts, on='gene', how='left')
        scored['pathway_count'] = scored['pathway_count'].fillna(0)
        scored['pathway_score'] = _min_max_normalize(scored['pathway_count'])
        scored = scored.drop(columns=['pathway_count'])

    # Compute final weighted score
    print("[scoring] Computing final scores...")
    scored['final_score'] = (
        w_omics * scored['omics_score'] +
        w_lit * scored['literature_score'] +
        w_path * scored['pathway_score']
    )

    # Sort by final_score descending
    ranked = scored.sort_values('final_score', ascending=False)

    # Reorder columns
    ranked = ranked[['gene', 'final_score', 'omics_score', 'literature_score', 'pathway_score']]

    print(f"[scoring] Writing output to: {output_path}")
    write_csv(ranked, output_path)

    # Print top 10 ranked genes
    print(f"\n[scoring] Top 10 ranked genes:")
    for i, (_, row) in enumerate(ranked.head(10).iterrows()):
        print(f"        {i+1}. {row['gene']}: final={row['final_score']:.2f}, "
              f"omics={row['omics_score']:.2f}")

    print(f"\n[scoring] Done. {len(ranked)} genes ranked.")
    return output_path
