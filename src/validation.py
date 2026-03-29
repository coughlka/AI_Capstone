"""Validation module: checks ranked results against known CRC biomarkers.

Uses a 3-tier biomarker set:
  Tier 1: FDA-approved / guideline-recommended (highest confidence)
  Tier 2: COSMIC Cancer Gene Census (CRC-specific entries)
  Tier 3: Validated clinical literature (prognostic, immune, epigenetic)
"""

import json
import os

import pandas as pd

from src.utils import load_config, ensure_dirs, read_csv


# --- 3-tier known CRC biomarker definitions ---

TIER1_FDA_GUIDELINE = {
    "KRAS", "NRAS", "BRAF", "EGFR", "ERBB2", "MLH1", "MSH2", "MSH6",
    "PMS2", "EPCAM", "VEGFA", "BMP3", "NDRG4", "DPYD",
}

TIER2_COSMIC_CGC = {
    "APC", "TP53", "PIK3CA", "SMAD4", "FBXW7", "PTEN", "CTNNB1",
    "RNF43", "POLE", "TGFBR2", "AXIN2", "SOX9", "TCF7L2",
    "AMER1", "ARID1A", "SMAD2", "ACVR2A", "ATM", "ZNRF3",
    "BCL9L", "RBM10", "PCBP1", "RPL22", "PTPRT", "FAM123B",
}

TIER3_VALIDATED = {
    "VIM", "CDH1", "CEACAM5", "MKI67", "CDX2", "DCC", "MYC", "MSH3",
    "ERBB3", "MET", "VEGFB", "FLT1", "KDR", "PDGFRA", "IGF2",
    "MAP2K1", "TGFB1", "SMAD3", "NOTCH1", "JAG1", "DLL4",
    "GSK3B", "CSNK1A1", "DVL2", "MYB", "ETV4", "FOXQ1",
    "CD274", "PDCD1", "CTLA4", "LAG3",
    "DNMT1", "DNMT3B", "TET2", "KDM6A",
    "MUTYH", "CHEK2", "BRCA1", "BRCA2", "PALB2", "PIK3R1",
    "CTNNB1",
}

# Build tier lookup
BIOMARKER_TIERS = {}
for g in TIER1_FDA_GUIDELINE:
    BIOMARKER_TIERS[g] = 1
for g in TIER2_COSMIC_CGC:
    BIOMARKER_TIERS.setdefault(g, 2)
for g in TIER3_VALIDATED:
    BIOMARKER_TIERS.setdefault(g, 3)

# All known biomarkers (union)
ALL_BIOMARKERS = TIER1_FDA_GUIDELINE | TIER2_COSMIC_CGC | TIER3_VALIDATED

# Functional category mapping
BIOMARKER_CATEGORIES = {}

for g in ["APC", "TP53", "PIK3CA", "SMAD4", "FBXW7", "PTEN", "CTNNB1",
          "RNF43", "POLE", "AMER1", "ARID1A", "SMAD2", "ACVR2A", "ATM",
          "ZNRF3", "BCL9L", "RBM10", "PCBP1", "RPL22", "PTPRT", "FAM123B",
          "MUTYH", "CHEK2", "BRCA1", "BRCA2", "PALB2", "PIK3R1"]:
    BIOMARKER_CATEGORIES[g] = "mutation/tumor-suppressor"

for g in ["KRAS", "NRAS", "BRAF", "EGFR", "ERBB2", "ERBB3", "MET",
          "VEGFA", "VEGFB", "FLT1", "KDR", "PDGFRA", "IGF2",
          "MAP2K1", "TGFBR2", "TGFB1", "SMAD3",
          "NOTCH1", "JAG1", "DLL4", "TCF7L2", "AXIN2", "SOX9",
          "GSK3B", "CSNK1A1", "DVL2", "MYC", "MYB", "ETV4", "FOXQ1"]:
    BIOMARKER_CATEGORIES[g] = "signaling/pathway"

for g in ["MLH1", "MSH2", "MSH6", "PMS2", "EPCAM", "MSH3"]:
    BIOMARKER_CATEGORIES[g] = "MMR/MSI"

for g in ["BMP3", "NDRG4", "VIM", "CDH1", "CEACAM5", "MKI67", "CDX2", "DCC"]:
    BIOMARKER_CATEGORIES[g] = "expression/methylation"

for g in ["CD274", "PDCD1", "CTLA4", "LAG3"]:
    BIOMARKER_CATEGORIES[g] = "immune-checkpoint"

for g in ["DNMT1", "DNMT3B", "TET2", "KDM6A"]:
    BIOMARKER_CATEGORIES[g] = "epigenetic"

for g in ["DPYD"]:
    BIOMARKER_CATEGORIES[g] = "pharmacogenomic"

# Final dict: gene -> {tier, category}
KNOWN_CRC_BIOMARKERS = {}
for gene in ALL_BIOMARKERS:
    KNOWN_CRC_BIOMARKERS[gene] = {
        "tier": BIOMARKER_TIERS.get(gene, 3),
        "category": BIOMARKER_CATEGORIES.get(gene, "other"),
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

    def _empty_result(symbol, info, ensembl_id='NOT_FOUND'):
        return {
            'gene_symbol': symbol,
            'gene': ensembl_id,
            'tier': info['tier'],
            'category': info['category'],
            'rank': None,
            'percentile': None,
            'final_score': None,
            'omics_score': None,
            'ml_importance_score': None,
            'literature_score': None,
            'pathway_score': None,
            'in_top_500': False,
            'in_top_1000': False,
            'in_top_5000': False,
        }

    for symbol, info in KNOWN_CRC_BIOMARKERS.items():
        ensembl_id = symbol_to_ensembl.get(symbol)
        if ensembl_id is None:
            results.append(_empty_result(symbol, info))
            continue

        match = ranked[ranked['gene'] == ensembl_id]
        if match.empty:
            results.append(_empty_result(symbol, info, ensembl_id))
            continue

        row = match.iloc[0]
        rank = int(row['rank'])
        percentile = round((1 - rank / total_genes) * 100, 1)

        results.append({
            'gene_symbol': symbol,
            'gene': ensembl_id,
            'tier': info['tier'],
            'category': info['category'],
            'rank': rank,
            'percentile': percentile,
            'final_score': round(row['final_score'], 2),
            'omics_score': round(row['omics_score'], 2),
            'ml_importance_score': round(row.get('ml_importance_score', 0), 2),
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

    # Tier breakdown
    tier_stats = {}
    for tier_num, tier_label in [(1, "FDA/guideline"), (2, "COSMIC CGC"), (3, "validated literature")]:
        tier_df = found[found['tier'] == tier_num]
        tier_total = report_df[report_df['tier'] == tier_num]
        if not tier_df.empty:
            tier_stats[f"tier_{tier_num}"] = {
                'label': tier_label,
                'found': len(tier_df),
                'total': len(tier_total),
                'median_rank': int(tier_df['rank'].median()),
                'mean_score': round(float(tier_df['final_score'].mean()), 2),
                'in_top_500': int(tier_df['in_top_500'].sum()),
                'best_gene': tier_df.iloc[0]['gene_symbol'],
                'best_rank': int(tier_df.iloc[0]['rank']),
            }

    # Category breakdown
    category_stats = {}
    for cat in sorted(found['category'].unique()):
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
            'tier': int(row['tier']),
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
        'tier_breakdown': tier_stats,
        'category_breakdown': category_stats,
        'top_known_biomarkers': top_known,
    }

    print(f"[validation] Writing summary to: {summary_path}")
    with open(summary_path, 'w') as f:
        json.dump(summary, f, indent=2)

    # Print summary
    print(f"\n[validation] === Known CRC Biomarker Validation ({len(KNOWN_CRC_BIOMARKERS)} markers) ===")
    print(f"  Found in dataset:          {len(found)}/{len(KNOWN_CRC_BIOMARKERS)}")
    print(f"  In top 500:                {in_top_500}")
    print(f"  In top 1,000:              {in_top_1000}")
    print(f"  In top 5,000:              {in_top_5000}")
    print(f"  Median rank:               {int(actual_median):,} (expected random: {int(expected_median):,})")
    if actual_median:
        print(f"  Enrichment ratio:          {expected_median / actual_median:.1f}x above random")

    print(f"\n  By tier:")
    for key in ['tier_1', 'tier_2', 'tier_3']:
        if key in tier_stats:
            t = tier_stats[key]
            print(f"    {key.replace('_', ' ').title()} ({t['label']}): "
                  f"{t['found']}/{t['total']} found, "
                  f"median rank={t['median_rank']:,}, "
                  f"top 500={t['in_top_500']}")

    print(f"\n  Top 10 known biomarkers by rank:")
    for _, row in found.head(10).iterrows():
        print(f"    {row['gene_symbol']:12s}  rank={int(row['rank']):>5,}  "
              f"score={row['final_score']:.2f}  tier={int(row['tier'])}  ({row['category']})")

    print(f"\n[validation] Done.")
    return report_path


if __name__ == '__main__':
    import sys
    config = sys.argv[1] if len(sys.argv) > 1 else 'config/config.yaml'
    run_validation(config)
