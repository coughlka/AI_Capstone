"""Gene ID mapping utilities for Ensembl to gene symbol conversion."""

import json
import os
import urllib.error
import urllib.parse
import urllib.request
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


def load_mapping_file(mapping_path: str) -> Dict[str, dict]:
    """Load a local Ensembl to gene symbol mapping file.

    Expected format: TSV with columns 'ensembl_id', 'gene_symbol', and
    optionally 'aliases' (semicolon-separated) and 'type_of_gene'.

    Args:
        mapping_path: Path to the mapping TSV file.

    Returns:
        Dictionary mapping Ensembl IDs (without version) to dicts with
        keys 'gene_symbol' and 'type_of_gene'.
    """
    if not os.path.exists(mapping_path):
        return {}

    df = pd.read_csv(mapping_path, sep='\t')

    if 'ensembl_id' not in df.columns or 'gene_symbol' not in df.columns:
        print(f"[gene_mapping] Warning: mapping file missing required columns")
        return {}

    has_type = 'type_of_gene' in df.columns
    mapping = {}
    for _, row in df.iterrows():
        ens_id = strip_ensembl_version(str(row['ensembl_id']))
        symbol = str(row['gene_symbol'])
        if ens_id and symbol and symbol != 'nan':
            mapping[ens_id] = {
                'gene_symbol': symbol,
                'type_of_gene': str(row['type_of_gene']) if has_type else '',
            }

    return mapping


def load_extended_mapping(mapping_path: str) -> Dict[str, dict]:
    """Load extended mapping including aliases and gene type.

    Handles both 3-column (legacy: ensembl_id, gene_symbol, type_of_gene) and
    4-column (ensembl_id, gene_symbol, aliases, type_of_gene) cache formats.

    Returns:
        Dictionary mapping Ensembl IDs to dicts with keys:
        'gene_symbol', 'aliases' (list), 'type_of_gene' (str).
    """
    if not os.path.exists(mapping_path):
        return {}

    df = pd.read_csv(mapping_path, sep='\t')
    if 'ensembl_id' not in df.columns or 'gene_symbol' not in df.columns:
        return {}

    has_aliases = 'aliases' in df.columns
    has_type = 'type_of_gene' in df.columns

    result = {}
    for _, row in df.iterrows():
        ens_id = strip_ensembl_version(str(row['ensembl_id']))
        symbol = str(row['gene_symbol'])
        if not ens_id or symbol == 'nan':
            continue

        if has_aliases:
            aliases_str = str(row.get('aliases', ''))
            aliases = [a.strip() for a in aliases_str.split(';') if a.strip() and a.strip() != 'nan']
        else:
            aliases = []

        gene_type = str(row['type_of_gene']) if has_type else ''
        if gene_type == 'nan':
            gene_type = ''

        result[ens_id] = {
            'gene_symbol': symbol,
            'aliases': aliases,
            'type_of_gene': gene_type,
        }
    return result


def fetch_symbols_from_mygene(ensembl_ids: list, batch_size: int = 1000) -> Dict[str, dict]:
    """Fetch gene symbols and types from mygene.info API.

    Args:
        ensembl_ids: List of Ensembl IDs (with or without versions).
        batch_size: Number of IDs to query per API call.

    Returns:
        Dictionary mapping Ensembl IDs (without version) to dicts with
        keys 'gene_symbol' and 'type_of_gene'.
    """
    extended = fetch_extended_from_mygene(ensembl_ids, batch_size=batch_size)
    return {eid: {'gene_symbol': info['gene_symbol'], 'type_of_gene': info.get('type_of_gene', '')}
            for eid, info in extended.items()}


def fetch_extended_from_mygene(ensembl_ids: list, batch_size: int = 1000) -> Dict[str, dict]:
    """Fetch gene symbols, aliases, and gene type from mygene.info API.

    Args:
        ensembl_ids: List of Ensembl IDs (with or without versions).
        batch_size: Number of IDs to query per API call.

    Returns:
        Dictionary mapping Ensembl IDs (without version) to dicts with keys:
        'gene_symbol' (str), 'aliases' (list[str]), 'type_of_gene' (str).
    """
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
            url = "https://mygene.info/v3/query"
            data = {
                "q": ",".join(batch),
                "scopes": "ensembl.gene",
                "fields": "symbol,alias,type_of_gene",
                "species": "human"
            }

            encoded_data = urllib.parse.urlencode(data).encode('utf-8')

            req = urllib.request.Request(url, data=encoded_data, method='POST')
            req.add_header('Content-Type', 'application/x-www-form-urlencoded')

            with urllib.request.urlopen(req, timeout=60) as response:
                results = json.loads(response.read().decode('utf-8'))

            for result in results:
                if isinstance(result, dict) and 'symbol' in result and 'query' in result:
                    ens_id = result['query']
                    symbol = result['symbol']
                    if not symbol:
                        continue

                    # Parse aliases: can be a string or list
                    raw_alias = result.get('alias', [])
                    if isinstance(raw_alias, str):
                        aliases = [raw_alias]
                    elif isinstance(raw_alias, list):
                        aliases = [a for a in raw_alias if isinstance(a, str)]
                    else:
                        aliases = []

                    gene_type = result.get('type_of_gene', '')

                    mapping[ens_id] = {
                        'gene_symbol': symbol,
                        'aliases': aliases,
                        'type_of_gene': gene_type or '',
                    }

        except urllib.error.URLError as e:
            print(f"[gene_mapping] Warning: API request failed: {e}")
            continue
        except json.JSONDecodeError as e:
            print(f"[gene_mapping] Warning: Failed to parse API response: {e}")
            continue

    print(f"[gene_mapping] Retrieved symbols for {len(mapping)}/{total} genes")
    return mapping


def fetch_aliases_from_ensembl(ensembl_ids: list) -> Dict[str, dict]:
    """Fetch gene synonyms from the Ensembl REST API (xrefs endpoint).

    This is a fallback for genes where mygene.info returned no aliases,
    which commonly happens for lncRNAs and other non-coding genes.
    The Ensembl xrefs endpoint returns HGNC synonyms (e.g., NIHCOLE -> LINC02163).

    Args:
        ensembl_ids: List of Ensembl IDs (without versions).

    Returns:
        Dictionary mapping Ensembl IDs to dicts with 'aliases' key.
    """
    results = {}
    total = len(ensembl_ids)
    if not total:
        return results

    print(f"[gene_mapping] Fetching aliases from Ensembl REST API for {total} genes...")

    for i, ens_id in enumerate(ensembl_ids):
        try:
            url = f"https://rest.ensembl.org/xrefs/id/{ens_id}?content-type=application/json"
            req = urllib.request.Request(url)
            with urllib.request.urlopen(req, timeout=30) as response:
                xrefs = json.loads(response.read().decode('utf-8'))

            aliases = []
            for xref in xrefs:
                if not isinstance(xref, dict):
                    continue
                if xref.get('dbname') == 'HGNC':
                    synonyms = xref.get('synonyms', [])
                    aliases.extend(s for s in synonyms if isinstance(s, str))

            if aliases:
                results[ens_id] = {'aliases': aliases}
                print(f"[gene_mapping]   {ens_id}: found aliases {aliases}")

        except (urllib.error.URLError, json.JSONDecodeError) as e:
            print(f"[gene_mapping]   {ens_id}: Ensembl lookup failed ({e})")
            continue

    print(f"[gene_mapping] Ensembl fallback: found aliases for {len(results)}/{total} genes")
    return results


def save_mapping_file(mapping: Dict[str, dict], output_path: str) -> None:
    """Save mapping dictionary to TSV file for future use.

    Accepts extended {id: {gene_symbol, aliases, type_of_gene}} format.

    Args:
        mapping: Dictionary mapping Ensembl IDs to extended info dicts.
        output_path: Path to save the TSV file.
    """
    rows = []
    for ens_id, value in mapping.items():
        if isinstance(value, dict):
            rows.append({
                'ensembl_id': ens_id,
                'gene_symbol': value.get('gene_symbol', ''),
                'aliases': ';'.join(value.get('aliases', [])),
                'type_of_gene': value.get('type_of_gene', ''),
            })
        else:
            rows.append({
                'ensembl_id': ens_id,
                'gene_symbol': value,
                'aliases': '',
                'type_of_gene': '',
            })
    df = pd.DataFrame(rows)
    df.to_csv(output_path, sep='\t', index=False)
    print(f"[gene_mapping] Saved mapping to {output_path}")


def get_gene_symbols(
    ensembl_ids: list,
    cache_path: Optional[str] = None,
    use_api: bool = True
) -> Dict[str, dict]:
    """Get gene symbols and types for Ensembl IDs, using cache if available.

    Args:
        ensembl_ids: List of Ensembl IDs (with or without versions).
        cache_path: Optional path to cache file. If provided, will load from
                   cache first and save new mappings to cache.
        use_api: Whether to fetch missing mappings from mygene.info API.

    Returns:
        Dictionary mapping Ensembl IDs (without version) to dicts with
        keys 'gene_symbol' and 'type_of_gene'.
    """
    extended = get_extended_gene_info(ensembl_ids, cache_path=cache_path, use_api=use_api)
    return {eid: {'gene_symbol': info['gene_symbol'], 'type_of_gene': info.get('type_of_gene', '')}
            for eid, info in extended.items()}


def get_extended_gene_info(
    ensembl_ids: list,
    cache_path: Optional[str] = None,
    use_api: bool = True
) -> Dict[str, dict]:
    """Get gene symbols, aliases, and type for Ensembl IDs.

    Args:
        ensembl_ids: List of Ensembl IDs (with or without versions).
        cache_path: Optional path to cache file.
        use_api: Whether to fetch missing mappings from mygene.info API.

    Returns:
        Dictionary mapping Ensembl IDs to dicts with keys:
        'gene_symbol', 'aliases' (list), 'type_of_gene' (str).
    """
    mapping = {}

    # Try to load from cache first
    if cache_path and os.path.exists(cache_path):
        print(f"[gene_mapping] Loading cached mapping from {cache_path}")
        extended = load_extended_mapping(cache_path)
        if extended:
            mapping = extended
        else:
            # Fall back to simple mapping for legacy cache files
            simple = load_mapping_file(cache_path)
            mapping = {eid: {'gene_symbol': info['gene_symbol'], 'aliases': [],
                             'type_of_gene': info.get('type_of_gene', '')}
                       for eid, info in simple.items()}
        print(f"[gene_mapping] Loaded {len(mapping)} cached mappings")

    # Find IDs that still need mapping
    stripped_ids = [strip_ensembl_version(eid) for eid in ensembl_ids]
    missing_ids = [eid for eid in stripped_ids if eid not in mapping]

    if missing_ids and use_api:
        print(f"[gene_mapping] {len(missing_ids)} IDs need mapping from API")
        api_mapping = fetch_extended_from_mygene(missing_ids)
        mapping.update(api_mapping)

        # Update cache if we got new mappings
        if cache_path and api_mapping:
            save_mapping_file(mapping, cache_path)

    return mapping


def map_ensembl_to_symbols(df: pd.DataFrame, ensembl_col: str = 'gene',
                           cache_path: Optional[str] = None,
                           include_type: bool = False) -> pd.DataFrame:
    """Add gene_symbol column to DataFrame with Ensembl IDs.

    Args:
        df: DataFrame containing Ensembl IDs.
        ensembl_col: Name of column containing Ensembl IDs.
        cache_path: Optional path to cache file.
        include_type: If True, also add a 'type_of_gene' column
                      (e.g. 'protein-coding', 'ncRNA', 'pseudo').

    Returns:
        DataFrame with added 'gene_symbol' column, and optionally
        'type_of_gene'.
    """
    ensembl_ids = df[ensembl_col].tolist()
    mapping = get_extended_gene_info(ensembl_ids, cache_path=cache_path)

    df = df.copy()
    df['gene_symbol'] = df[ensembl_col].apply(
        lambda x: mapping.get(strip_ensembl_version(x), {}).get('gene_symbol', '')
    )

    if include_type:
        df['type_of_gene'] = df[ensembl_col].apply(
            lambda x: mapping.get(strip_ensembl_version(x), {}).get('type_of_gene', '')
        )

    mapped_count = (df['gene_symbol'] != '').sum()
    print(f"[gene_mapping] Mapped {mapped_count}/{len(df)} genes to symbols")

    return df
