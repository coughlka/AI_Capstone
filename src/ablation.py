"""Scoring weight ablation study.

Sweeps over weight combinations for the four evidence streams and measures
enrichment ratio against known CRC biomarkers at each configuration.
Reads pre-computed component scores from ranked_candidates.csv so no
upstream pipeline steps need to re-run.
"""

import os
from itertools import product

import numpy as np
import pandas as pd

from src.validation import KNOWN_CRC_BIOMARKERS
from src.utils import load_config, read_csv


def _enrichment_at_weights(scored_df, symbol_lookup, w_omics, w_ml, w_lit, w_path):
    """Recompute final scores with given weights and return validation metrics."""
    df = scored_df.copy()
    df['final_score'] = (
        w_omics * df['omics_score'] +
        w_ml * df['ml_importance_score'] +
        w_lit * df['literature_score'] +
        w_path * df['pathway_score']
    )
    df = df.sort_values('final_score', ascending=False).reset_index(drop=True)
    df['rank'] = df.index + 1
    total = len(df)

    # Map symbols and find known biomarkers
    df['gene_symbol'] = df['gene'].map(symbol_lookup)
    known_symbols = set(KNOWN_CRC_BIOMARKERS.keys())
    found = df[df['gene_symbol'].isin(known_symbols)]

    if found.empty:
        return None

    ranks = found['rank']
    median_rank = float(ranks.median())
    expected_median = total / 2
    enrichment = expected_median / median_rank if median_rank > 0 else 0

    # BMP3 rank
    bmp3 = found[found['gene_symbol'] == 'BMP3']
    bmp3_rank = int(bmp3['rank'].iloc[0]) if not bmp3.empty else None

    return {
        'w_omics': w_omics,
        'w_ml': w_ml,
        'w_literature': w_lit,
        'w_pathway': w_path,
        'enrichment_ratio': round(enrichment, 3),
        'median_rank': int(median_rank),
        'in_top_500': int((ranks <= 500).sum()),
        'in_top_1000': int((ranks <= 1000).sum()),
        'in_top_5000': int((ranks <= 5000).sum()),
        'bmp3_rank': bmp3_rank,
        'found': len(found),
    }


def run_ablation(config_path: str) -> str:
    """Run weight ablation study over a grid of weight combinations.

    Args:
        config_path: Path to the configuration YAML file.

    Returns:
        Path to the ablation results CSV.
    """
    config = load_config(config_path)
    outputs_dir = config['paths']['outputs_dir']

    ranked_path = os.path.join(outputs_dir, 'ranked_candidates.csv')
    omics_path = os.path.join(outputs_dir, 'omics_evidence.csv')
    output_path = os.path.join(outputs_dir, 'ablation_results.csv')

    print("[ablation] Loading data...")
    scored = read_csv(ranked_path)
    omics = read_csv(omics_path)
    symbol_lookup = omics.set_index('gene')['gene_symbol'].to_dict()

    # Current weights from config
    weights = config.get('scoring', {}).get('weights', {})
    cur_omics = weights.get('omics', 0.40)
    cur_ml = weights.get('ml_importance', 0.10)
    cur_lit = weights.get('literature', 0.30)
    cur_path = weights.get('pathway', 0.20)

    # Generate weight grid: step 0.10, each in [0.0, 0.70], must sum to 1.0
    step = 0.10
    values = np.arange(0.0, 0.71, step)
    values = np.round(values, 2)

    combos = []
    for wo, wm, wl, wp in product(values, repeat=4):
        if abs(wo + wm + wl + wp - 1.0) < 1e-6:
            combos.append((round(wo, 2), round(wm, 2), round(wl, 2), round(wp, 2)))

    print(f"[ablation] Testing {len(combos)} weight combinations...")

    results = []
    for i, (wo, wm, wl, wp) in enumerate(combos):
        r = _enrichment_at_weights(scored, symbol_lookup, wo, wm, wl, wp)
        if r:
            results.append(r)
        if (i + 1) % 50 == 0:
            print(f"[ablation]   {i+1}/{len(combos)} complete...")

    results_df = pd.DataFrame(results)
    results_df = results_df.sort_values('enrichment_ratio', ascending=False).reset_index(drop=True)
    results_df.to_csv(output_path, index=False)

    # Print results
    print(f"\n[ablation] Tested {len(results)} valid combinations")
    print(f"[ablation] Results saved to {output_path}")

    # Current config performance
    cur = _enrichment_at_weights(scored, symbol_lookup, cur_omics, cur_ml, cur_lit, cur_path)
    print(f"\n[ablation] Current weights: omics={cur_omics}, ml={cur_ml}, lit={cur_lit}, path={cur_path}")
    print(f"[ablation]   Enrichment: {cur['enrichment_ratio']}x, "
          f"BMP3 rank: {cur['bmp3_rank']}, "
          f"Top 500: {cur['in_top_500']}, Top 1000: {cur['in_top_1000']}")

    print(f"\n[ablation] Top 10 weight combinations by enrichment ratio:")
    print(f"  {'omics':>6} {'ml':>6} {'lit':>6} {'path':>6} | {'enrich':>7} {'med_rank':>9} {'top500':>6} {'top1k':>6} {'top5k':>6} {'BMP3':>6}")
    print(f"  {'-'*6} {'-'*6} {'-'*6} {'-'*6} | {'-'*7} {'-'*9} {'-'*6} {'-'*6} {'-'*6} {'-'*6}")
    for _, row in results_df.head(10).iterrows():
        print(f"  {row['w_omics']:6.2f} {row['w_ml']:6.2f} {row['w_literature']:6.2f} {row['w_pathway']:6.2f} | "
              f"{row['enrichment_ratio']:7.3f} {row['median_rank']:9,} {row['in_top_500']:6} {row['in_top_1000']:6} {row['in_top_5000']:6} {row['bmp3_rank']:6}")

    print(f"\n[ablation] Done.")
    return output_path
