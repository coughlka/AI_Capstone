"""Pathway enrichment evidence module."""

import io
import os
import csv
import time
import httpx
import pandas as pd

from src.utils import load_config, ensure_dirs, write_csv

REACTOME_API_URL = "https://reactome.org/AnalysisService"


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
    
    try:
        # Join genes with newlines.
        payload = "\n".join(genes_to_map)
        
        # Explicitly set Content-Type to text/plain as required by Reactome
        headers = {"Content-Type": "text/plain"}
        
        # We use projection to get the token
        response = httpx.post(
            f"{REACTOME_API_URL}/identifiers/projection",
            content=payload,
            headers=headers,
            params={"pageSize": 1, "page": 1}, 
            timeout=30.0
        )
        response.raise_for_status()
        
        analysis_data = response.json()

        print(f"[pathway] Analysis data:\n{analysis_data}")

        token = analysis_data['summary']['token']
        print(f"[pathway] Analysis token received: {token}")
        
    except httpx.HTTPError as e:
        print(f"[pathway] Error querying Reactome API: {e}")
        # Build empty df so pipeline doesn't crash completely
        pd.DataFrame(columns=['gene', 'pathway_count', 'top_pathways']).to_csv(output_path, index=False)
        return output_path

    # Step 3: Fetch pathways for each gene (mapping)
    print("[pathway] Retrieving gene-to-pathway mappings...")
    
    try:
        # Download results as CSV to get the mapping
        # Endpoint: /download/{token}/pathways/TOTAL/result.csv
        
        mapping_response = httpx.get(
            f"{REACTOME_API_URL}/download/{token}/pathways/TOTAL/result.csv",
            timeout=60.0
        )
        mapping_response.raise_for_status()
        
    except httpx.HTTPError as e:
        print(f"[pathway] Error retrieving results: {e}")
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
