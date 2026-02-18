"""Validation module: checks ranked results against known CRC biomarkers."""

import json
import os

import pandas as pd

from src.utils import load_config, ensure_dirs, read_csv


# Curated list of known CRC biomarkers with functional categories
KNOWN_CRC_BIOMARKERS = {
    # Mutation-driven oncogenes/tumor suppressors
    "APC": "mutation-driven",
    "KRAS": "mutation-driven",
    "TP53": "mutation-driven",
    "BRAF": "mutation-driven",
    "PIK3CA": "mutation-driven",
    "SMAD4": "mutation-driven",
    "FBXW7": "mutation-driven",
    "PTEN": "mutation-driven",
    "NRAS": "mutation-driven",
    "CTNNB1": "mutation-driven",
    "RNF43": "mutation-driven",
    "POLE": "mutation-driven",
    # Expression/methylation biomarkers
    "BMP3": "expression/methylation",
    "NDRG4": "expression/methylation",
    "VIM": "expression/methylation",
    "CDH1": "expression/methylation",
    "CEACAM5": "expression/methylation",
    "MKI67": "expression/methylation",
    "CDX2": "expression/methylation",
    # Signaling/pathway genes
    "EGFR": "signaling/pathway",
    "ERBB2": "signaling/pathway",
    "VEGFA": "signaling/pathway",
    "TGFBR2": "signaling/pathway",
    "MYC": "signaling/pathway",
    "MLH1": "signaling/pathway",
    "MSH2": "signaling/pathway",
    "MSH6": "signaling/pathway",
    "PMS2": "signaling/pathway",
    "TCF7L2": "signaling/pathway",
    "AXIN2": "signaling/pathway",
    "SOX9": "signaling/pathway",
    "DCC": "signaling/pathway",
    "MSH3": "signaling/pathway",
}


def run_validation(config_path: str) -> str:
    """Validate ranked results against known CRC biomarkers.

    Args:
        config_path: Path to the configuration YAML file.

    Returns:
        Path to the validation report CSV.
    """
    print("[validation] Loading configuration...")
    config = load_config(config_path)
    ensure_dirs(config)

    outputs_dir = config['paths']['outputs_dir']
    ranked_path = os.path.join(outputs_dir, 'ranked_candidates.csv')
    omics_path = os.path.join(outputs_dir, 'omics_evidence.csv')
    report_path = os.path.join(outputs_dir, 'validation_report.csv')
    summary_path = os.path.join(outputs_dir, 'validation_summary.json')

    if not os.path.exists(ranked_path):
        raise FileNotFoundError(
            f"Missing {ranked_path} (run scoring step first)"
        )
    if not os.path.exists(omics_path):
        raise FileNotFoundError(
            f"Missing {omics_path} (run omics step first)"
        )

    # Load data
    print(f"[validation] Reading ranked candidates from: {ranked_path}")
    ranked = read_csv(ranked_path)

    print(f"[validation] Reading omics evidence from: {omics_path}")
    omics = read_csv(omics_path)

    # Build symbol -> ensembl mapping from omics
    symbol_to_ensembl = omics.set_index('gene_symbol')['gene'].to_dict()

    # Add rank column (1-indexed)
    ranked = ranked.sort_values('final_score', ascending=False).reset_index(drop=True)
    ranked['rank'] = ranked.index + 1
    total_genes = len(ranked)

    # Build ensembl -> symbol mapping
    ensembl_to_symbol = omics.set_index('gene')['gene_symbol'].to_dict()
    ranked['gene_symbol'] = ranked['gene'].map(ensembl_to_symbol)

    # Validate each known biomarker
    print(f"[validation] Checking {len(KNOWN_CRC_BIOMARKERS)} known CRC biomarkers...")
    results = []

    for symbol, category in KNOWN_CRC_BIOMARKERS.items():
        ensembl_id = symbol_to_ensembl.get(symbol)
        if ensembl_id is None:
            results.append({
                'gene_symbol': symbol,
                'gene': 'NOT_FOUND',
                'category': category,
                'rank': None,
                'percentile': None,
                'final_score': None,
                'omics_score': None,
                'literature_score': None,
                'pathway_score': None,
                'in_top_500': False,
                'in_top_1000': False,
                'in_top_5000': False,
            })
            continue

        match = ranked[ranked['gene'] == ensembl_id]
        if match.empty:
            results.append({
                'gene_symbol': symbol,
                'gene': ensembl_id,
                'category': category,
                'rank': None,
                'percentile': None,
                'final_score': None,
                'omics_score': None,
                'literature_score': None,
                'pathway_score': None,
                'in_top_500': False,
                'in_top_1000': False,
                'in_top_5000': False,
            })
            continue

        row = match.iloc[0]
        rank = int(row['rank'])
        percentile = round((1 - rank / total_genes) * 100, 1)

        results.append({
            'gene_symbol': symbol,
            'gene': ensembl_id,
            'category': category,
            'rank': rank,
            'percentile': percentile,
            'final_score': round(row['final_score'], 2),
            'omics_score': round(row['omics_score'], 2),
            'literature_score': round(row['literature_score'], 2),
            'pathway_score': round(row['pathway_score'], 2),
            'in_top_500': rank <= 500,
            'in_top_1000': rank <= 1000,
            'in_top_5000': rank <= 5000,
        })

    report_df = pd.DataFrame(results)
    report_df = report_df.sort_values('rank', na_position='last')

    # Write per-gene report
    print(f"[validation] Writing validation report to: {report_path}")
    report_df.to_csv(report_path, index=False)

    # Compute summary metrics
    found = report_df[report_df['rank'].notna()]
    not_found = report_df[report_df['rank'].isna()]
    ranks = found['rank'].astype(int)

    expected_median = total_genes / 2
    actual_median = float(ranks.median()) if not ranks.empty else None

    in_top_500 = int(found['in_top_500'].sum())
    in_top_1000 = int(found['in_top_1000'].sum())
    in_top_5000 = int(found['in_top_5000'].sum())

    # Category breakdown
    category_stats = {}
    for cat in ['mutation-driven', 'expression/methylation', 'signaling/pathway']:
        cat_df = found[found['category'] == cat]
        if not cat_df.empty:
            category_stats[cat] = {
                'count': len(cat_df),
                'median_rank': int(cat_df['rank'].median()),
                'mean_score': round(float(cat_df['final_score'].mean()), 2),
                'best_gene': cat_df.iloc[0]['gene_symbol'],
                'best_rank': int(cat_df.iloc[0]['rank']),
            }

    # Top known biomarkers
    top_known = []
    for _, row in found.head(10).iterrows():
        top_known.append({
            'gene_symbol': row['gene_symbol'],
            'rank': int(row['rank']),
            'final_score': row['final_score'],
            'category': row['category'],
        })

    summary = {
        'total_known_biomarkers': len(KNOWN_CRC_BIOMARKERS),
        'found_in_dataset': len(found),
        'not_found': len(not_found),
        'not_found_genes': not_found['gene_symbol'].tolist(),
        'total_ranked_genes': total_genes,
        'median_rank': int(actual_median) if actual_median else None,
        'expected_median_rank': int(expected_median),
        'mean_rank': int(ranks.mean()) if not ranks.empty else None,
        'enrichment_ratio': round(expected_median / actual_median, 2) if actual_median else None,
        'in_top_500': in_top_500,
        'in_top_1000': in_top_1000,
        'in_top_5000': in_top_5000,
        'category_breakdown': category_stats,
        'top_known_biomarkers': top_known,
    }

    print(f"[validation] Writing summary to: {summary_path}")
    with open(summary_path, 'w') as f:
        json.dump(summary, f, indent=2)

    # Print summary
    print(f"\n[validation] === Known CRC Biomarker Validation ===")
    print(f"  Known biomarkers checked:  {len(KNOWN_CRC_BIOMARKERS)}")
    print(f"  Found in dataset:          {len(found)}/{len(KNOWN_CRC_BIOMARKERS)}")
    print(f"  In top 500:                {in_top_500}")
    print(f"  In top 1,000:              {in_top_1000}")
    print(f"  In top 5,000:              {in_top_5000}")
    print(f"  Median rank:               {int(actual_median):,} (expected random: {int(expected_median):,})")
    if actual_median:
        print(f"  Enrichment ratio:          {expected_median / actual_median:.1f}x above random")
    print(f"\n  Top 10 known biomarkers by rank:")
    for _, row in found.head(10).iterrows():
        print(f"    {row['gene_symbol']:12s}  rank={int(row['rank']):>5,}  "
              f"score={row['final_score']:.2f}  ({row['category']})")

    print(f"\n[validation] Done.")
    return report_path


if __name__ == '__main__':
    import sys
    config = sys.argv[1] if len(sys.argv) > 1 else 'config/config.yaml'
    run_validation(config)
