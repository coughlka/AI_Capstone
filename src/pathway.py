"""Pathway enrichment evidence module."""

import io
import os
import csv
import time
from typing import Optional

import httpx
import pandas as pd

from src.utils import load_config, ensure_dirs, write_csv

REACTOME_API_URL = "https://reactome.org/AnalysisService"
MAX_RETRIES = 3
RETRY_DELAYS = [1, 3, 10]  # seconds


def _request_with_retry(
    method: str,
    url: str,
    timeout: float = 30.0,
    **kwargs
) -> Optional[httpx.Response]:
    """Make HTTP request with exponential backoff retry.

    Args:
        method: HTTP method ('GET' or 'POST')
        url: Request URL
        timeout: Request timeout in seconds
        **kwargs: Additional arguments to pass to httpx

    Returns:
        Response object if successful, None if all retries failed
    """
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
    """Run the pathway enrichment pipeline step.

    Queries Reactome Analysis Service to find enriched pathways for the candidate genes.

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
    
    # Get candidate list path from config (matching omics.py output)
    candidates_config = config.get('candidates', {})
    candidate_list_path = candidates_config.get('output_path', os.path.join(outputs_dir, 'candidates.csv'))

    print(f"[pathway] Reading candidate list from: {candidate_list_path}")
    if not os.path.exists(candidate_list_path):
        print(f"[pathway] Warning: Candidate list not found at {candidate_list_path}. creating empty output.")
        pd.DataFrame(columns=['gene', 'pathway_count', 'top_pathways']).to_csv(output_path, index=False)
        return output_path

    candidates_df = pd.read_csv(candidate_list_path)
    
    # Use gene symbols if available, otherwise Ensembl IDs
    if 'gene_symbol' in candidates_df.columns:
        # Filter out empty symbols
        genes_to_map = candidates_df['gene_symbol'].dropna().unique().tolist()
        print(f"[pathway] Using {len(genes_to_map)} gene symbols for enrichment")
    else:
        genes_to_map = candidates_df['gene'].dropna().unique().tolist()
        print(f"[pathway] Using {len(genes_to_map)} gene IDs for enrichment")

    if not genes_to_map:
        print("[pathway] No genes to map. Exiting.")
        pd.DataFrame(columns=['gene', 'pathway_count', 'top_pathways']).to_csv(output_path, index=False)
        return output_path

    # Step 2: Submit to Reactome API
    print("[pathway] Submitting genes to Reactome Analysis Service...")

    # Join genes with newlines
    payload = "\n".join(genes_to_map)

    # Explicitly set Content-Type to text/plain as required by Reactome
    headers = {"Content-Type": "text/plain"}

    # We use projection to get the token (with retry logic)
    response = _request_with_retry(
        "POST",
        f"{REACTOME_API_URL}/identifiers/projection",
        timeout=30.0,
        content=payload,
        headers=headers,
        params={"pageSize": 1, "page": 1}
    )

    if response is None:
        print("[pathway] Failed to query Reactome API after retries")
        pd.DataFrame(columns=['gene', 'pathway_count', 'top_pathways']).to_csv(output_path, index=False)
        return output_path

    analysis_data = response.json()
    print(f"[pathway] Analysis submitted successfully")

    token = analysis_data['summary']['token']
    print(f"[pathway] Analysis token received: {token}")

    # Step 3: Fetch pathways for each gene (mapping)
    print("[pathway] Retrieving gene-to-pathway mappings...")

    # Download results as CSV (with retry logic)
    mapping_response = _request_with_retry(
        "GET",
        f"{REACTOME_API_URL}/download/{token}/pathways/TOTAL/result.csv",
        timeout=60.0
    )

    if mapping_response is None:
        print("[pathway] Failed to retrieve results after retries")
        pd.DataFrame(columns=['gene', 'pathway_count', 'top_pathways']).to_csv(output_path, index=False)
        return output_path

    # Parse the response (CSV)
    # Columns expected:
    # Pathway identifier, Pathway name, #Entities found, ..., Entities FDR, ..., Submitted entities found
    lines = mapping_response.text.strip().split('\n')
    
    if not lines:
        print("[pathway] No significant pathways found.")
        pd.DataFrame(columns=['gene', 'pathway_count', 'top_pathways']).to_csv(output_path, index=False)
        return output_path

    # We can use csv.reader to handle quoted fields (e.g. "KRAS;EGFR")
    reader = csv.reader(io.StringIO(mapping_response.text))
    header = next(reader, None)

    print(f"[pathway] Header: {header}")
    
    if not header:
        print("[pathway] Empty result CSV.")
        pd.DataFrame(columns=['gene', 'pathway_count', 'top_pathways']).to_csv(output_path, index=False)
        return output_path

    # Find columns
    try:
        idx_pathway_name = header.index("Pathway name")
        idx_fdr = header.index("Entities FDR")
        idx_submitted = header.index("Submitted entities found")
    except ValueError:
        # Fallback to known indices if headers change slightly
        # Based on verification: 
        # 1: Pathway name, 6: Entities FDR, 12: Submitted entities found
        idx_pathway_name = 1
        idx_fdr = 6
        idx_submitted = 12

    gene_counts = {}
    gene_top_paths = {}
    count_sig_pathways = 0
    
    # Get FDR threshold from config
    pathway_config = config.get('pathway', {})
    fdr_threshold = pathway_config.get('fdr_threshold', 0.05)
    print(f"[pathway] Using FDR threshold: {fdr_threshold}")

    for row in reader:
        if len(row) <= idx_submitted:
            continue
            
        pathway_name = row[idx_pathway_name]
        try:
            fdr = float(row[idx_fdr])
        except ValueError:
            continue
            
        if fdr < fdr_threshold:
            count_sig_pathways += 1
            # Genes are semicolon separated
            entities_str = row[idx_submitted]
            matched_genes = entities_str.split(';')
            
            for gene in matched_genes:
                gene = gene.strip()
                if not gene: continue
                
                gene_counts[gene] = gene_counts.get(gene, 0) + 1
                
                if gene not in gene_top_paths:
                    gene_top_paths[gene] = []
                
                # Limit to top 5 pathways per gene
                if len(gene_top_paths[gene]) < 5:
                    gene_top_paths[gene].append(pathway_name)
    
    print(f"[pathway] Found {count_sig_pathways} significant pathways (FDR < {fdr_threshold})")

    # Prepare output
    final_rows = []
    
    # Iterate through original candidates to preserve order
    if 'gene_symbol' in candidates_df.columns:
        original_symbols = candidates_df['gene_symbol'].tolist()
        original_ids = candidates_df['gene'].tolist()
    else:
        original_ids = candidates_df['gene'].tolist()
        original_symbols = original_ids

    for i, gene_id in enumerate(original_ids):
        symbol = str(original_symbols[i])
        
        # The API returns whatever we sent it (symbols or IDs) in "Submitted entities found"
        # We sent 'genes_to_map'.
        # If we sent symbols, lookup key is symbol.
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
