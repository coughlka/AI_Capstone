"""Omics evidence module for gene expression analysis."""

import os

import pandas as pd

from src.utils import load_config, ensure_dirs, write_csv


def run_omics(config_path: str) -> str:
    """Run the omics analysis pipeline step.

    Loads gene expression counts, computes per-gene mean and variance,
    and outputs a summary CSV.

    Args:
        config_path: Path to the configuration YAML file.

    Returns:
        Path to the output CSV file.

    Raises:
        FileNotFoundError: If the counts TSV file is not found.
    """
    print("[omics] Loading configuration...")
    config = load_config(config_path)
    ensure_dirs(config)

    counts_path = config['omics']['counts_path']
    dataset_label = config['omics']['dataset_label']
    outputs_dir = config['paths']['outputs_dir']
    output_path = os.path.join(outputs_dir, 'omics_evidence.csv')

    print(f"[omics] Reading counts from: {counts_path}")
    if not os.path.exists(counts_path):
        raise FileNotFoundError(
            f"Omics counts file not found: {counts_path}. "
            "Please ensure the data file exists in the data/ directory."
        )

    # Read TSV: first column is gene id, remaining columns are samples
    counts_df = pd.read_csv(counts_path, sep='\t', index_col=0)

    print(f"[omics] Loaded {counts_df.shape[0]} genes x {counts_df.shape[1]} samples")

    # Compute per-gene statistics across samples
    print("[omics] Computing per-gene mean and variance...")
    gene_stats = pd.DataFrame({
        'gene': counts_df.index,
        'mean_expr': counts_df.mean(axis=1).values,
        'var_expr': counts_df.var(axis=1).values,
        'dataset': dataset_label
    })

    print(f"[omics] Writing output to: {output_path}")
    write_csv(gene_stats, output_path)

    print(f"[omics] Done. {len(gene_stats)} genes processed.")
    return output_path
