"""Gene ID mapping utilities for Ensembl to gene symbol conversion."""

import json
import os
import urllib.request
import urllib.error
from typing import Dict, Optional

import pandas as pd


def strip_ensembl_version(ensembl_id: str) -> str:
    """Remove version suffix from Ensembl ID.

    Args:
        ensembl_id: Ensembl ID with or without version (e.g., 'ENSG00000141510.18')

    Returns:
        Ensembl ID without version (e.g., 'ENSG00000141510')
    """
    if '.' in ensembl_id:
        return ensembl_id.split('.')[0]
    return ensembl_id


def load_mapping_file(mapping_path: str) -> Dict[str, str]:
    """Load a local Ensembl to gene symbol mapping file.

    Expected format: TSV with columns 'ensembl_id' and 'gene_symbol'

    Args:
        mapping_path: Path to the mapping TSV file.

    Returns:
        Dictionary mapping Ensembl IDs (without version) to gene symbols.
    """
    if not os.path.exists(mapping_path):
        return {}

    df = pd.read_csv(mapping_path, sep='\t')

    if 'ensembl_id' not in df.columns or 'gene_symbol' not in df.columns:
        print(f"[gene_mapping] Warning: mapping file missing required columns")
        return {}

    mapping = {}
    for _, row in df.iterrows():
        ens_id = strip_ensembl_version(str(row['ensembl_id']))
        symbol = str(row['gene_symbol'])
        if ens_id and symbol and symbol != 'nan':
            mapping[ens_id] = symbol

    return mapping


def fetch_symbols_from_mygene(ensembl_ids: list, batch_size: int = 1000) -> Dict[str, str]:
    """Fetch gene symbols from mygene.info API.

    Args:
        ensembl_ids: List of Ensembl IDs (with or without versions).
        batch_size: Number of IDs to query per API call.

    Returns:
        Dictionary mapping Ensembl IDs (without version) to gene symbols.
    """
    # Strip versions for querying
    id_map = {strip_ensembl_version(eid): eid for eid in ensembl_ids}
    unique_ids = list(id_map.keys())

    mapping = {}
    total = len(unique_ids)

    print(f"[gene_mapping] Fetching symbols for {total} unique Ensembl IDs from mygene.info...")

    for i in range(0, total, batch_size):
        batch = unique_ids[i:i + batch_size]
        batch_num = i // batch_size + 1
        total_batches = (total + batch_size - 1) // batch_size

        print(f"[gene_mapping] Batch {batch_num}/{total_batches} ({len(batch)} IDs)...")

        try:
            # mygene.info POST query
            url = "https://mygene.info/v3/query"
            data = {
                "q": ",".join(batch),
                "scopes": "ensembl.gene",
                "fields": "symbol",
                "species": "human"
            }

            # URL encode the data
            encoded_data = urllib.parse.urlencode(data).encode('utf-8')

            req = urllib.request.Request(url, data=encoded_data, method='POST')
            req.add_header('Content-Type', 'application/x-www-form-urlencoded')

            with urllib.request.urlopen(req, timeout=60) as response:
                results = json.loads(response.read().decode('utf-8'))

            for result in results:
                if isinstance(result, dict) and 'symbol' in result and 'query' in result:
                    ens_id = result['query']
                    symbol = result['symbol']
                    if symbol:
                        mapping[ens_id] = symbol

        except urllib.error.URLError as e:
            print(f"[gene_mapping] Warning: API request failed: {e}")
            continue
        except json.JSONDecodeError as e:
            print(f"[gene_mapping] Warning: Failed to parse API response: {e}")
            continue

    print(f"[gene_mapping] Retrieved symbols for {len(mapping)}/{total} genes")
    return mapping


def save_mapping_file(mapping: Dict[str, str], output_path: str) -> None:
    """Save mapping dictionary to TSV file for future use.

    Args:
        mapping: Dictionary mapping Ensembl IDs to gene symbols.
        output_path: Path to save the TSV file.
    """
    df = pd.DataFrame([
        {'ensembl_id': ens_id, 'gene_symbol': symbol}
        for ens_id, symbol in mapping.items()
    ])
    df.to_csv(output_path, sep='\t', index=False)
    print(f"[gene_mapping] Saved mapping to {output_path}")


def get_gene_symbols(
    ensembl_ids: list,
    cache_path: Optional[str] = None,
    use_api: bool = True
) -> Dict[str, str]:
    """Get gene symbols for Ensembl IDs, using cache if available.

    Args:
        ensembl_ids: List of Ensembl IDs (with or without versions).
        cache_path: Optional path to cache file. If provided, will load from
                   cache first and save new mappings to cache.
        use_api: Whether to fetch missing mappings from mygene.info API.

    Returns:
        Dictionary mapping Ensembl IDs (without version) to gene symbols.
    """
    mapping = {}

    # Try to load from cache first
    if cache_path and os.path.exists(cache_path):
        print(f"[gene_mapping] Loading cached mapping from {cache_path}")
        mapping = load_mapping_file(cache_path)
        print(f"[gene_mapping] Loaded {len(mapping)} cached mappings")

    # Find IDs that still need mapping
    stripped_ids = [strip_ensembl_version(eid) for eid in ensembl_ids]
    missing_ids = [eid for eid in stripped_ids if eid not in mapping]

    if missing_ids and use_api:
        print(f"[gene_mapping] {len(missing_ids)} IDs need mapping from API")
        api_mapping = fetch_symbols_from_mygene(missing_ids)
        mapping.update(api_mapping)

        # Update cache if we got new mappings
        if cache_path and api_mapping:
            save_mapping_file(mapping, cache_path)

    return mapping


def map_ensembl_to_symbols(df: pd.DataFrame, ensembl_col: str = 'gene',
                           cache_path: Optional[str] = None) -> pd.DataFrame:
    """Add gene_symbol column to DataFrame with Ensembl IDs.

    Args:
        df: DataFrame containing Ensembl IDs.
        ensembl_col: Name of column containing Ensembl IDs.
        cache_path: Optional path to cache file.

    Returns:
        DataFrame with added 'gene_symbol' column.
    """
    ensembl_ids = df[ensembl_col].tolist()
    mapping = get_gene_symbols(ensembl_ids, cache_path=cache_path)

    # Map IDs (strip version for lookup)
    df = df.copy()
    df['gene_symbol'] = df[ensembl_col].apply(
        lambda x: mapping.get(strip_ensembl_version(x), '')
    )

    mapped_count = (df['gene_symbol'] != '').sum()
    print(f"[gene_mapping] Mapped {mapped_count}/{len(df)} genes to symbols")

    return df
