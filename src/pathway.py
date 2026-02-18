"""Pathway enrichment evidence module using g:Profiler."""

import os
import time
from typing import Optional

import httpx
import pandas as pd

from src.utils import load_config, ensure_dirs, write_csv

GPROFILER_URL = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"
MAX_RETRIES = 3
RETRY_DELAYS = [1, 3, 10]  # seconds


def _request_with_retry(
    method: str,
    url: str,
    timeout: float = 60.0,
    **kwargs
) -> Optional[httpx.Response]:
    """Make HTTP request with exponential backoff retry."""
    for attempt in range(MAX_RETRIES):
        try:
            if method.upper() == "POST":
                response = httpx.post(url, timeout=timeout, **kwargs)
            else:
                response = httpx.get(url, timeout=timeout, **kwargs)
            response.raise_for_status()
            return response
        except (httpx.HTTPError, httpx.TimeoutException) as e:
            if attempt < MAX_RETRIES - 1:
                delay = RETRY_DELAYS[attempt]
                print(f"[pathway] Request failed (attempt {attempt + 1}/{MAX_RETRIES}), retrying in {delay}s... ({e})")
                time.sleep(delay)
            else:
                print(f"[pathway] All {MAX_RETRIES} attempts failed: {e}")
                return None
    return None


def run_pathway(config_path: str) -> str:
    """Run pathway enrichment using g:Profiler.

    Queries g:Profiler to find enriched pathways across GO, KEGG,
    Reactome, and WikiPathways for the candidate genes.

    Args:
        config_path: Path to the configuration YAML file.

    Returns:
        Path to the output CSV file.
    """
    print("[pathway] Loading configuration...")
    config = load_config(config_path)
    ensure_dirs(config)

    outputs_dir = config['paths']['outputs_dir']
    output_path = os.path.join(outputs_dir, 'pathway_evidence.csv')

    # Get candidate list
    candidates_config = config.get('candidates', {})
    candidate_list_path = candidates_config.get('output_path', os.path.join(outputs_dir, 'candidates.csv'))

    print(f"[pathway] Reading candidate list from: {candidate_list_path}")
    if not os.path.exists(candidate_list_path):
        print(f"[pathway] Warning: Candidate list not found at {candidate_list_path}. Creating empty output.")
        pd.DataFrame(columns=['gene', 'pathway_count', 'top_pathways']).to_csv(output_path, index=False)
        return output_path

    candidates_df = pd.read_csv(candidate_list_path)

    # Build gene list - prefer symbols, fallback to Ensembl IDs
    # g:Profiler handles both, but symbols tend to have better coverage
    if 'gene_symbol' in candidates_df.columns:
        genes_to_query = candidates_df['gene_symbol'].dropna().unique().tolist()
        print(f"[pathway] Using {len(genes_to_query)} gene symbols for enrichment")
    else:
        genes_to_query = candidates_df['gene'].dropna().unique().tolist()
        print(f"[pathway] Using {len(genes_to_query)} gene IDs for enrichment")

    if not genes_to_query:
        print("[pathway] No genes to map. Exiting.")
        pd.DataFrame(columns=['gene', 'pathway_count', 'top_pathways']).to_csv(output_path, index=False)
        return output_path

    # Get config
    pathway_config = config.get('pathway', {})
    fdr_threshold = pathway_config.get('fdr_threshold', 0.05)
    sources = pathway_config.get('sources', ["GO:BP", "KEGG", "REAC", "WP"])

    print(f"[pathway] Querying g:Profiler (sources: {sources}, FDR < {fdr_threshold})...")

    # Build g:Profiler request
    payload = {
        "organism": "hsapiens",
        "query": genes_to_query,
        "sources": sources,
        "user_threshold": fdr_threshold,
        "significance_threshold_method": "g_SCS",
        "no_evidences": False,
    }

    response = _request_with_retry("POST", GPROFILER_URL, json=payload)

    if response is None:
        print("[pathway] Failed to query g:Profiler API after retries")
        pd.DataFrame(columns=['gene', 'pathway_count', 'top_pathways']).to_csv(output_path, index=False)
        return output_path

    data = response.json()
    results = data.get("result", [])
    print(f"[pathway] Received {len(results)} significant terms")

    # Get the gene alignment from metadata
    # intersections[i] aligns to ensgs[i]; use mapping to get symbols
    meta = data.get("meta", {})
    genes_meta = meta.get("genes_metadata", {}).get("query", {}).get("query_1", {})
    ensgs = genes_meta.get("ensgs", [])
    symbol_mapping = genes_meta.get("mapping", {})

    # Build reverse map: ensg -> symbol
    ensg_to_symbol = {}
    for sym, ensg_list in symbol_mapping.items():
        for eid in ensg_list:
            ensg_to_symbol[eid] = sym

    print(f"[pathway] {len(ensgs)} genes recognized by g:Profiler")

    # Map each gene to its pathway hits
    gene_counts = {}
    gene_top_paths = {}

    for term in results:
        pathway_name = term.get("name", "Unknown")
        source = term.get("source", "")
        label = f"{pathway_name} ({source})"

        # intersections is aligned to ensgs list
        # each element is a list of evidence codes; non-empty = gene is in this term
        intersections = term.get("intersections", [])

        for i, evidence in enumerate(intersections):
            if not evidence or i >= len(ensgs):
                continue
            # This gene is in this term
            symbol = ensg_to_symbol.get(ensgs[i], ensgs[i])

            gene_counts[symbol] = gene_counts.get(symbol, 0) + 1

            if symbol not in gene_top_paths:
                gene_top_paths[symbol] = []
            if len(gene_top_paths[symbol]) < 5:
                gene_top_paths[symbol].append(label)

    genes_with_hits = len(gene_counts)
    print(f"[pathway] {genes_with_hits}/{len(genes_to_query)} genes mapped to at least one pathway")

    # Build output, preserving original candidate order
    final_rows = []

    if 'gene_symbol' in candidates_df.columns:
        original_symbols = candidates_df['gene_symbol'].tolist()
        original_ids = candidates_df['gene'].tolist()
    else:
        original_ids = candidates_df['gene'].tolist()
        original_symbols = original_ids

    for i, gene_id in enumerate(original_ids):
        symbol = str(original_symbols[i])

        if 'gene_symbol' in candidates_df.columns:
            lookup_key = symbol
        else:
            lookup_key = str(gene_id)

        count = gene_counts.get(lookup_key, 0)
        tops = gene_top_paths.get(lookup_key, [])

        final_rows.append({
            'gene': gene_id,
            'pathway_count': count,
            'top_pathways': "; ".join(tops)
        })

    pathway_evidence = pd.DataFrame(final_rows)

    print(f"[pathway] Writing output to: {output_path}")
    write_csv(pathway_evidence, output_path)

    print("[pathway] Done.")
    return output_path
