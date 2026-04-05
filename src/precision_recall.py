"""Precision/recall evaluation at rank thresholds.

Computes precision, recall, and F1 for known CRC biomarkers at various
rank cutoffs. Breaks down results by biomarker tier (FDA, COSMIC, literature).
"""

import os

import pandas as pd

from src.validation import KNOWN_CRC_BIOMARKERS, BIOMARKER_TIERS
from src.utils import load_config, read_csv


def _pr_at_k(ranked_symbols, known_set, k):
    """Compute precision, recall, F1 at rank threshold k."""
    top_k = set(ranked_symbols[:k])
    tp = len(top_k & known_set)
    precision = tp / k if k > 0 else 0
    recall = tp / len(known_set) if known_set else 0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0
    return tp, precision, recall, f1


def run_precision_recall(config_path: str) -> str:
    """Run precision/recall evaluation at multiple rank thresholds.

    Args:
        config_path: Path to the configuration YAML file.

    Returns:
        Path to the output CSV file.
    """
    config = load_config(config_path)
    outputs_dir = config['paths']['outputs_dir']

    ranked_path = os.path.join(outputs_dir, 'ranked_candidates.csv')
    omics_path = os.path.join(outputs_dir, 'omics_evidence.csv')
    output_path = os.path.join(outputs_dir, 'precision_recall.csv')

    print("[pr_eval] Loading data...")
    ranked = read_csv(ranked_path)
    omics = read_csv(omics_path)

    # Sort by final score and map symbols
    ranked = ranked.sort_values('final_score', ascending=False).reset_index(drop=True)
    symbol_lookup = omics.set_index('gene')['gene_symbol'].to_dict()
    ranked_symbols = ranked['gene'].map(symbol_lookup).tolist()

    # Build known biomarker sets by tier
    all_known = set(KNOWN_CRC_BIOMARKERS.keys())
    tier1 = {g for g, t in BIOMARKER_TIERS.items() if t == 1}
    tier2 = {g for g, t in BIOMARKER_TIERS.items() if t == 2}
    tier3 = {g for g, t in BIOMARKER_TIERS.items() if t == 3}

    # Filter to genes actually in the dataset
    dataset_symbols = set(symbol_lookup.values())
    all_found = all_known & dataset_symbols
    tier1_found = tier1 & dataset_symbols
    tier2_found = tier2 & dataset_symbols
    tier3_found = tier3 & dataset_symbols

    thresholds = [10, 25, 50, 100, 250, 500, 1000, 2500, 5000, 10000]

    results = []
    for k in thresholds:
        tp, p, r, f1 = _pr_at_k(ranked_symbols, all_found, k)
        tp1, p1, r1, f1_1 = _pr_at_k(ranked_symbols, tier1_found, k)
        tp2, p2, r2, f1_2 = _pr_at_k(ranked_symbols, tier2_found, k)
        tp3, p3, r3, f1_3 = _pr_at_k(ranked_symbols, tier3_found, k)

        results.append({
            'k': k,
            'tp_all': tp, 'precision_all': round(p, 4), 'recall_all': round(r, 4), 'f1_all': round(f1, 4),
            'tp_tier1': tp1, 'precision_tier1': round(p1, 4), 'recall_tier1': round(r1, 4), 'f1_tier1': round(f1_1, 4),
            'tp_tier2': tp2, 'precision_tier2': round(p2, 4), 'recall_tier2': round(r2, 4), 'f1_tier2': round(f1_2, 4),
            'tp_tier3': tp3, 'precision_tier3': round(p3, 4), 'recall_tier3': round(r3, 4), 'f1_tier3': round(f1_3, 4),
        })

    results_df = pd.DataFrame(results)
    results_df.to_csv(output_path, index=False)

    # Print results
    print(f"\n[pr_eval] Known biomarkers in dataset: {len(all_found)} "
          f"(Tier 1: {len(tier1_found)}, Tier 2: {len(tier2_found)}, Tier 3: {len(tier3_found)})")

    print(f"\n[pr_eval] All biomarkers ({len(all_found)}):")
    print(f"  {'K':>7} | {'TP':>4} {'Prec':>7} {'Recall':>7} {'F1':>7}")
    print(f"  {'-'*7} | {'-'*4} {'-'*7} {'-'*7} {'-'*7}")
    for _, row in results_df.iterrows():
        print(f"  {row['k']:>7,} | {row['tp_all']:>4} {row['precision_all']:>7.4f} {row['recall_all']:>7.4f} {row['f1_all']:>7.4f}")

    print(f"\n[pr_eval] Tier 1 - FDA/guideline ({len(tier1_found)}):")
    print(f"  {'K':>7} | {'TP':>4} {'Prec':>7} {'Recall':>7} {'F1':>7}")
    print(f"  {'-'*7} | {'-'*4} {'-'*7} {'-'*7} {'-'*7}")
    for _, row in results_df.iterrows():
        print(f"  {row['k']:>7,} | {row['tp_tier1']:>4} {row['precision_tier1']:>7.4f} {row['recall_tier1']:>7.4f} {row['f1_tier1']:>7.4f}")

    print(f"\n[pr_eval] Tier 2 - COSMIC CGC ({len(tier2_found)}):")
    print(f"  {'K':>7} | {'TP':>4} {'Prec':>7} {'Recall':>7} {'F1':>7}")
    print(f"  {'-'*7} | {'-'*4} {'-'*7} {'-'*7} {'-'*7}")
    for _, row in results_df.iterrows():
        print(f"  {row['k']:>7,} | {row['tp_tier2']:>4} {row['precision_tier2']:>7.4f} {row['recall_tier2']:>7.4f} {row['f1_tier2']:>7.4f}")

    print(f"\n[pr_eval] Tier 3 - Validated literature ({len(tier3_found)}):")
    print(f"  {'K':>7} | {'TP':>4} {'Prec':>7} {'Recall':>7} {'F1':>7}")
    print(f"  {'-'*7} | {'-'*4} {'-'*7} {'-'*7} {'-'*7}")
    for _, row in results_df.iterrows():
        print(f"  {row['k']:>7,} | {row['tp_tier3']:>4} {row['precision_tier3']:>7.4f} {row['recall_tier3']:>7.4f} {row['f1_tier3']:>7.4f}")

    print(f"\n[pr_eval] Results saved to {output_path}")
    print(f"[pr_eval] Done.")
    return output_path
