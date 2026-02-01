"""Omics evidence module for tumor vs normal differential expression analysis."""

import os

import numpy as np
import pandas as pd
from scipy import stats

from src.utils import load_config, ensure_dirs, write_csv
from src.gene_mapping import map_ensembl_to_symbols, strip_ensembl_version


def parse_tcga_sample_labels(sample_ids: list) -> pd.DataFrame:
    """Parse TCGA barcodes to extract sample type labels.

    TCGA barcode format: TCGA-XX-XXXX-SSA where SS is the sample type code
    at positions 14-15 (1-indexed) of the barcode.

    Sample type codes:
    - 01-09: Tumor samples (01=primary, 02=recurrent, 06=metastatic, etc.)
    - 10-19: Normal samples (10=blood normal, 11=solid tissue normal, etc.)

    Args:
        sample_ids: List of TCGA sample barcode strings.

    Returns:
        DataFrame with columns: sample_id, sample_type_code, sample_type_group
        where sample_type_group is 'tumor', 'normal', or 'other'.
    """
    records = []
    for sample_id in sample_ids:
        # Extract sample type code from positions 14-15 (0-indexed: 13-14)
        if len(sample_id) >= 15 and sample_id.startswith('TCGA-'):
            sample_type_code = sample_id[13:15]
            try:
                code_int = int(sample_type_code)
                if 1 <= code_int <= 9:
                    group = 'tumor'
                elif 10 <= code_int <= 19:
                    group = 'normal'
                else:
                    group = 'other'
            except ValueError:
                group = 'other'
                sample_type_code = 'XX'
        else:
            sample_type_code = 'XX'
            group = 'other'

        records.append({
            'sample_id': sample_id,
            'sample_type_code': sample_type_code,
            'sample_type_group': group
        })

    return pd.DataFrame(records)


def benjamini_hochberg(p_values: np.ndarray) -> np.ndarray:
    """Apply Benjamini-Hochberg FDR correction to p-values.

    Args:
        p_values: Array of p-values.

    Returns:
        Array of FDR-adjusted p-values (q-values).
    """
    n = len(p_values)
    if n == 0:
        return np.array([])

    # Handle NaN values
    valid_mask = ~np.isnan(p_values)
    fdr = np.full(n, np.nan)

    if not valid_mask.any():
        return fdr

    valid_pvals = p_values[valid_mask]
    n_valid = len(valid_pvals)

    # Sort p-values and track original indices
    sorted_indices = np.argsort(valid_pvals)
    sorted_pvals = valid_pvals[sorted_indices]

    # Compute BH adjusted p-values
    # q_i = p_i * n / rank
    ranks = np.arange(1, n_valid + 1)
    adjusted = sorted_pvals * n_valid / ranks

    # Ensure monotonicity: q_i = min(q_i, q_{i+1}, ..., q_n)
    adjusted_monotonic = np.minimum.accumulate(adjusted[::-1])[::-1]

    # Cap at 1.0
    adjusted_monotonic = np.minimum(adjusted_monotonic, 1.0)

    # Map back to original order
    unsorted_adjusted = np.empty(n_valid)
    unsorted_adjusted[sorted_indices] = adjusted_monotonic

    # Place back into full array
    fdr[valid_mask] = unsorted_adjusted

    return fdr


def run_omics(config_path: str) -> str:
    """Run tumor vs normal differential expression analysis.

    Loads gene expression counts, classifies samples as tumor or normal
    based on TCGA barcodes, performs differential expression analysis
    using Welch's t-test, and applies FDR correction.

    Args:
        config_path: Path to the configuration YAML file.

    Returns:
        Path to the output CSV file.

    Raises:
        FileNotFoundError: If the counts TSV file is not found.
        ValueError: If insufficient tumor or normal samples are found.
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
    n_genes_raw = counts_df.shape[0]
    n_samples = counts_df.shape[1]

    print(f"[omics] Loaded {n_genes_raw} genes x {n_samples} samples")

    # Parse sample labels from TCGA barcodes
    print("[omics] Parsing TCGA sample barcodes...")
    sample_labels = parse_tcga_sample_labels(counts_df.columns.tolist())

    tumor_samples = sample_labels[sample_labels['sample_type_group'] == 'tumor']['sample_id'].tolist()
    normal_samples = sample_labels[sample_labels['sample_type_group'] == 'normal']['sample_id'].tolist()

    n_tumor = len(tumor_samples)
    n_normal = len(normal_samples)

    print(f"[omics] Sample classification:")
    print(f"        - Tumor samples: {n_tumor}")
    print(f"        - Normal samples: {n_normal}")

    # Check we have enough samples for DE analysis
    if n_tumor < 3:
        raise ValueError(f"Insufficient tumor samples for DE analysis: {n_tumor} (need at least 3)")
    if n_normal < 3:
        raise ValueError(f"Insufficient normal samples for DE analysis: {n_normal} (need at least 3)")

    # Print sample type breakdown
    type_counts = sample_labels.groupby('sample_type_code').size().sort_values(ascending=False)
    print("[omics] Sample type breakdown:")
    for code, count in type_counts.items():
        print(f"        - {code}: {count}")

    # Ensure numeric values
    counts_df = counts_df.apply(pd.to_numeric, errors='coerce')

    # Filter genes to reduce noise and runtime
    print("[omics] Filtering genes...")
    gene_means = counts_df.mean(axis=1)
    nonzero_fraction = (counts_df > 0).sum(axis=1) / n_samples

    # Keep genes where mean >= 1.0 OR nonzero in >= 20% of samples
    keep_mask = (gene_means >= 1.0) | (nonzero_fraction >= 0.20)
    counts_filtered = counts_df[keep_mask]
    n_genes_filtered = counts_filtered.shape[0]

    print(f"[omics] Gene filtering: {n_genes_raw} -> {n_genes_filtered} genes")
    print(f"        (kept genes with mean >= 1.0 or nonzero in >= 20% samples)")

    # Split into tumor and normal expression matrices
    tumor_expr = counts_filtered[tumor_samples]
    normal_expr = counts_filtered[normal_samples]

    # Compute differential expression statistics
    print("[omics] Computing differential expression statistics...")

    results = []
    genes = counts_filtered.index.tolist()

    for gene in genes:
        tumor_vals = tumor_expr.loc[gene].values
        normal_vals = normal_expr.loc[gene].values

        tumor_mean = np.nanmean(tumor_vals)
        normal_mean = np.nanmean(normal_vals)

        # log2 fold change (already on log2 scale, so just subtract)
        log2fc = tumor_mean - normal_mean

        # Direction
        direction = 'up' if log2fc > 0 else 'down'

        # Welch's t-test
        # Handle cases where variance is zero
        tumor_var = np.nanvar(tumor_vals, ddof=1)
        normal_var = np.nanvar(normal_vals, ddof=1)

        if tumor_var == 0 and normal_var == 0:
            # No variance in either group
            p_value = 1.0
        else:
            try:
                _, p_value = stats.ttest_ind(tumor_vals, normal_vals, equal_var=False, nan_policy='omit')
                if np.isnan(p_value):
                    p_value = 1.0
            except Exception:
                p_value = 1.0

        results.append({
            'gene': gene,
            'log2fc': log2fc,
            'p_value': p_value,
            'direction': direction,
            'tumor_mean': tumor_mean,
            'normal_mean': normal_mean
        })

    results_df = pd.DataFrame(results)

    # Apply Benjamini-Hochberg FDR correction
    print("[omics] Applying Benjamini-Hochberg FDR correction...")
    results_df['fdr'] = benjamini_hochberg(results_df['p_value'].values)

    # Add dataset label
    results_df['dataset'] = dataset_label

    # Add gene symbol mapping
    print("[omics] Mapping Ensembl IDs to gene symbols...")
    gene_mapping_config = config.get('gene_mapping', {})
    cache_path = gene_mapping_config.get('cache_path', 'data/ensembl_to_symbol_cache.tsv')
    use_api = gene_mapping_config.get('use_api', True)

    if use_api:
        results_df = map_ensembl_to_symbols(results_df, ensembl_col='gene', cache_path=cache_path)
    else:
        # Just add empty column if API disabled
        results_df['gene_symbol'] = ''
        print("[omics] Gene symbol mapping disabled (use_api=false)")

    # Reorder columns to match contract (gene_symbol after gene)
    results_df = results_df[['gene', 'gene_symbol', 'log2fc', 'p_value', 'fdr', 'direction',
                              'tumor_mean', 'normal_mean', 'dataset']]

    # Sort by FDR for output
    results_df = results_df.sort_values('fdr')

    print(f"[omics] Writing output to: {output_path}")
    write_csv(results_df, output_path)

    # Generate candidate list for downstream modules
    candidates_config = config.get('candidates', {})
    top_n = candidates_config.get('top_n', 500)
    fdr_threshold = candidates_config.get('fdr_threshold', 0.05)
    candidates_path = candidates_config.get('output_path', os.path.join(outputs_dir, 'candidates.csv'))

    print(f"\n[omics] Generating candidate list (top {top_n}, FDR < {fdr_threshold})...")

    # Filter significant genes and take top N by DE signal
    sig_genes = results_df[results_df['fdr'] < fdr_threshold].copy()
    sig_genes['de_signal'] = sig_genes['log2fc'].abs() * (-np.log10(sig_genes['fdr'] + 1e-300))
    candidates = sig_genes.nlargest(min(top_n, len(sig_genes)), 'de_signal')

    # Create candidate list with essential columns
    candidates_out = candidates[['gene', 'gene_symbol', 'log2fc', 'fdr', 'direction']].copy()
    candidates_out['rank'] = range(1, len(candidates_out) + 1)

    write_csv(candidates_out, candidates_path)
    print(f"[omics] Wrote {len(candidates_out)} candidates to {candidates_path}")

    # Print summary statistics
    sig_005 = (results_df['fdr'] < 0.05).sum()
    sig_001 = (results_df['fdr'] < 0.01).sum()
    up_reg = ((results_df['fdr'] < 0.05) & (results_df['direction'] == 'up')).sum()
    down_reg = ((results_df['fdr'] < 0.05) & (results_df['direction'] == 'down')).sum()

    print(f"\n[omics] Summary statistics:")
    print(f"        - Total genes tested: {n_genes_filtered}")
    print(f"        - Significant at FDR < 0.05: {sig_005}")
    print(f"        - Significant at FDR < 0.01: {sig_001}")
    print(f"        - Upregulated (FDR < 0.05): {up_reg}")
    print(f"        - Downregulated (FDR < 0.05): {down_reg}")

    # Print top 10 genes by absolute log2fc
    print(f"\n[omics] Top 10 genes by |log2FC| (FDR < 0.05):")
    top_sig = results_df[results_df['fdr'] < 0.05].copy()
    if len(top_sig) > 0:
        top_sig['abs_log2fc'] = top_sig['log2fc'].abs()
        top_by_fc = top_sig.nlargest(10, 'abs_log2fc')
        for _, row in top_by_fc.iterrows():
            symbol = row['gene_symbol'] if row['gene_symbol'] else '(unmapped)'
            print(f"        {symbol} ({row['gene']}): log2FC={row['log2fc']:.3f}, FDR={row['fdr']:.2e}, {row['direction']}")
    else:
        print("        (no significant genes at FDR < 0.05)")

    # Validation check for known CRC genes
    known_crc_genes = ['APC', 'KRAS', 'TP53', 'SMAD4', 'PIK3CA', 'BRAF']
    print(f"\n[omics] Known CRC gene check:")
    for gene_symbol in known_crc_genes:
        match = results_df[results_df['gene_symbol'] == gene_symbol]
        if len(match) > 0:
            row = match.iloc[0]
            status = "SIG" if row['fdr'] < 0.05 else "n.s."
            print(f"        {gene_symbol}: log2FC={row['log2fc']:.3f}, FDR={row['fdr']:.2e}, {row['direction']} [{status}]")
        else:
            print(f"        {gene_symbol}: not found in dataset")

    print(f"\n[omics] Done. {len(results_df)} genes processed.")
    return output_path
