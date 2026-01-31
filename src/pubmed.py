"""PubMed literature evidence module."""

import os

import pandas as pd

from src.utils import load_config, ensure_dirs, write_csv


def run_pubmed(config_path: str) -> str:
    """Run the PubMed literature search pipeline step.

    Currently a stub that creates the output schema without calling external APIs.
    Future implementation will query PubMed for gene-disease associations.

    Args:
        config_path: Path to the configuration YAML file.

    Returns:
        Path to the output CSV file.
    """
    print("[pubmed] Loading configuration...")
    config = load_config(config_path)
    ensure_dirs(config)

    outputs_dir = config['paths']['outputs_dir']
    output_path = os.path.join(outputs_dir, 'lit_evidence.csv')
    candidate_list_path = config.get('candidates', {}).get('candidate_list_path')

    # Try to read candidate list if specified
    candidates = []
    if candidate_list_path and os.path.exists(candidate_list_path):
        print(f"[pubmed] Reading candidate list from: {candidate_list_path}")
        candidates_df = pd.read_csv(candidate_list_path)
        if 'gene' in candidates_df.columns:
            candidates = candidates_df['gene'].tolist()
        print(f"[pubmed] Found {len(candidates)} candidate genes")
    else:
        print("[pubmed] No candidate list found, creating empty output schema")

    # Create output DataFrame with correct schema (stub - no API calls)
    # Schema: gene, pmid, year, study_type, role, sample_type, directionality, snippet
    lit_evidence = pd.DataFrame(columns=[
        'gene', 'pmid', 'year', 'study_type', 'role',
        'sample_type', 'directionality', 'snippet'
    ])

    print(f"[pubmed] Writing output to: {output_path}")
    write_csv(lit_evidence, output_path)

    print("[pubmed] Done. (Stub - no external API calls made)")
    return output_path
