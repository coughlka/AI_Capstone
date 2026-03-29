"""PubMed literature evidence module.

Queries NCBI E-Utilities to retrieve abstracts for candidate genes,
producing per-gene literature evidence for the scoring pipeline.
"""

import os
import re
import time
import xml.etree.ElementTree as ET
from datetime import datetime, timezone
from typing import Dict, List, Optional

import pandas as pd
import requests
from dotenv import load_dotenv

from src.gene_mapping import load_extended_mapping, fetch_aliases_from_ensembl, save_mapping_file
from src.utils import load_config, ensure_dirs, write_csv

# Load .env from project root
_PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
load_dotenv(os.path.join(_PROJECT_ROOT, ".env"))

EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


class _RateLimiter:
    def __init__(self, max_per_second: float):
        self.min_interval = 1.0 / max_per_second if max_per_second > 0 else 0.0
        self._last = 0.0

    def wait(self):
        now = time.time()
        elapsed = now - self._last
        if elapsed < self.min_interval:
            time.sleep(self.min_interval - elapsed)
        self._last = time.time()


class PubMedClient:
    """Thin wrapper around NCBI E-Utilities with rate limiting."""

    def __init__(self, email: str, tool: str, api_key: Optional[str] = None, timeout: int = 30):
        self.email = email
        self.tool = tool
        self.api_key = api_key or os.getenv("NCBI_API_KEY")
        self.timeout = timeout
        self.session = requests.Session()
        max_rps = 3.0 if self.api_key else 1.0
        self.limiter = _RateLimiter(max_per_second=max_rps)

    def _params(self, extra: dict) -> dict:
        params = {"tool": self.tool, "email": self.email}
        if self.api_key:
            params["api_key"] = self.api_key
        params.update(extra)
        return params

    def esearch(self, query: str, retmax: int = 5) -> List[str]:
        self.limiter.wait()
        params = self._params({
            "db": "pubmed", "term": query,
            "retmode": "json", "retmax": str(retmax),
        })
        r = self.session.get(f"{EUTILS_BASE}/esearch.fcgi", params=params, timeout=self.timeout)
        r.raise_for_status()
        return (r.json().get("esearchresult", {}) or {}).get("idlist", []) or []

    def efetch_xml(self, pmids: List[str]) -> str:
        if not pmids:
            return ""
        self.limiter.wait()
        params = self._params({
            "db": "pubmed", "id": ",".join(pmids),
            "retmode": "xml", "rettype": "abstract",
        })
        r = self.session.get(f"{EUTILS_BASE}/efetch.fcgi", params=params, timeout=self.timeout)
        r.raise_for_status()
        return r.text


def _strip(text: Optional[str]) -> str:
    return re.sub(r"\s+", " ", (text or "")).strip()


def _parse_pubmed_xml(xml_text: str) -> Dict[str, dict]:
    """Parse PubMed XML into dict of pmid -> record."""
    if not xml_text.strip():
        return {}

    root = ET.fromstring(xml_text)
    retrieved_at = datetime.now(timezone.utc).isoformat()
    out = {}

    for article in root.findall(".//PubmedArticle"):
        pmid_el = article.find(".//MedlineCitation/PMID")
        if pmid_el is None or not pmid_el.text:
            continue
        pmid = pmid_el.text.strip()

        title = _strip(article.findtext(".//Article/ArticleTitle") or "")

        # Build abstract
        abstract_parts = []
        for abs_el in article.findall(".//Article/Abstract/AbstractText"):
            label = abs_el.attrib.get("Label") or abs_el.attrib.get("NlmCategory") or ""
            chunk = _strip("".join(abs_el.itertext()))
            if chunk:
                abstract_parts.append(f"{label}: {chunk}" if label else chunk)
        abstract = _strip(" ".join(abstract_parts))

        # Extract year
        year = ""
        pubdate_el = article.find(".//Article/Journal/JournalIssue/PubDate")
        if pubdate_el is not None:
            y = _strip(pubdate_el.findtext("Year") or "")
            medline = _strip(pubdate_el.findtext("MedlineDate") or "")
            if y:
                year = y
            elif medline:
                # Try to extract year from MedlineDate (e.g., "2020 Jan-Feb")
                m = re.search(r"(\d{4})", medline)
                if m:
                    year = m.group(1)

        journal = _strip(article.findtext(".//Article/Journal/Title") or "")

        # MeSH terms
        mesh_terms = []
        for mh in article.findall(".//MeshHeadingList/MeshHeading/DescriptorName"):
            t = _strip("".join(mh.itertext()))
            if t:
                mesh_terms.append(t)

        out[pmid] = {
            "pmid": pmid,
            "title": title,
            "abstract": abstract,
            "year": year,
            "journal": journal,
            "mesh_terms": mesh_terms,
        }

    return out


def _build_query(gene_symbol: str, templates: List[str], aliases: Optional[List[str]] = None) -> str:
    """Build a PubMed query from gene symbol, aliases, and templates.

    When aliases are provided, the query searches for any of the names,
    ensuring genes like NIHCOLE (alias: LINC02163) are found under all
    known identifiers.
    """
    all_names = [gene_symbol]
    if aliases:
        for alias in aliases:
            if alias.upper() != gene_symbol.upper() and alias not in all_names:
                all_names.append(alias)

    if templates:
        parts = []
        for name in all_names:
            parts.extend(f'({t.format(gene=name)})' for t in templates)
        return " OR ".join(parts)

    name_parts = [f'("{name}"[Title/Abstract] OR "{name}"[All Fields])' for name in all_names]
    return " OR ".join(name_parts)


def _snippet(abstract: str, max_len: int = 300) -> str:
    """First max_len chars of abstract as a snippet."""
    if len(abstract) <= max_len:
        return abstract
    return abstract[:max_len].rsplit(" ", 1)[0] + "..."


def run_pubmed(config_path: str) -> str:
    """Run the PubMed literature search pipeline step.

    Queries PubMed for each candidate gene and retrieves abstracts.
    Uses gene aliases from the mapping cache to broaden search coverage,
    and retries zero-result genes via Ensembl REST API alias lookup.

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

    # Get candidate list
    candidates_config = config.get('candidates', {})
    candidate_list_path = candidates_config.get('output_path', os.path.join(outputs_dir, 'candidates.csv'))

    print(f"[pubmed] Reading candidate list from: {candidate_list_path}")
    if not os.path.exists(candidate_list_path):
        print(f"[pubmed] Warning: Candidate list not found at {candidate_list_path}. Creating empty output.")
        pd.DataFrame(columns=['gene', 'pmid', 'year', 'title', 'abstract', 'snippet']).to_csv(output_path, index=False)
        return output_path

    candidates_df = pd.read_csv(candidate_list_path)

    # Get gene symbols
    if 'gene_symbol' in candidates_df.columns:
        gene_ids = candidates_df['gene'].tolist()
        gene_symbols = candidates_df['gene_symbol'].tolist()
    else:
        gene_ids = candidates_df['gene'].tolist()
        gene_symbols = gene_ids

    # PubMed config
    pubmed_config = config.get('pubmed', {})
    email = pubmed_config.get('email', 'user@example.com')
    tool = pubmed_config.get('tool', 'ai-capstone')
    api_key = pubmed_config.get('api_key') or os.getenv("NCBI_API_KEY")
    max_abstracts = pubmed_config.get('max_abstracts_per_gene', 5)
    query_templates = pubmed_config.get('query_templates', [])

    print(f"[pubmed] Email: {email}, API key: {'yes' if api_key else 'no'}")
    print(f"[pubmed] Max abstracts per gene: {max_abstracts}")
    print(f"[pubmed] Query templates: {query_templates}")

    client = PubMedClient(email=email, tool=tool, api_key=api_key)

    # Load gene aliases from mapping cache for improved search coverage
    alias_map = {}  # gene_id -> list of aliases
    gene_mapping_config = config.get('gene_mapping', {})
    mapping_cache = gene_mapping_config.get('cache_path', 'data/ensembl_to_symbol_cache.tsv')
    if os.path.exists(mapping_cache):
        extended = load_extended_mapping(mapping_cache)
        for ens_id, info in extended.items():
            alias_map[ens_id] = info.get('aliases', [])
        print(f"[pubmed] Loaded aliases for {len(alias_map)} genes")

    rows = []
    total = len(gene_symbols)
    zero_result_genes = []  # Track genes with no PubMed results for retry

    for i, (gene_id, symbol) in enumerate(zip(gene_ids, gene_symbols)):
        symbol_str = str(symbol).strip()
        if not symbol_str or symbol_str.lower() in ('nan', 'none'):
            continue

        # Get aliases for this gene to broaden PubMed search
        gene_id_stripped = gene_id.split('.')[0] if '.' in str(gene_id) else str(gene_id)
        aliases = alias_map.get(gene_id_stripped, [])

        query = _build_query(symbol_str, query_templates, aliases=aliases)

        try:
            pmids = client.esearch(query=query, retmax=max_abstracts)
        except Exception as e:
            print(f"[pubmed] [{i+1}/{total}] {symbol_str}: search failed ({e})")
            continue

        if not pmids:
            if (i + 1) % 50 == 0 or i == 0:
                print(f"[pubmed] [{i+1}/{total}] {symbol_str}: 0 results")
            # Track for alias-based retry if no aliases were available
            if not aliases:
                zero_result_genes.append((gene_id, gene_id_stripped, symbol_str))
            continue

        try:
            xml = client.efetch_xml(pmids)
            parsed = _parse_pubmed_xml(xml)
        except Exception as e:
            print(f"[pubmed] [{i+1}/{total}] {symbol_str}: fetch failed ({e})")
            continue

        for pmid in pmids:
            if pmid not in parsed:
                continue
            rec = parsed[pmid]
            rows.append({
                'gene': gene_id,
                'pmid': rec['pmid'],
                'year': rec['year'],
                'title': rec['title'],
                'abstract': rec['abstract'],
                'snippet': _snippet(rec['abstract']),
            })

        if (i + 1) % 50 == 0 or i == 0:
            print(f"[pubmed] [{i+1}/{total}] {symbol_str}: {len(pmids)} abstracts")

    # Retry genes with zero results by fetching aliases from Ensembl REST API.
    # This catches lncRNAs and other genes where mygene.info had no alias data
    # (e.g., NIHCOLE -> LINC02163).
    if zero_result_genes:
        print(f"\n[pubmed] {len(zero_result_genes)} genes had 0 results and no aliases — "
              f"trying Ensembl alias lookup...")
        retry_ids = [eid for _, eid, _ in zero_result_genes]
        ensembl_aliases = fetch_aliases_from_ensembl(retry_ids)

        # Update the mapping cache with newly discovered aliases
        if ensembl_aliases and os.path.exists(mapping_cache):
            extended = load_extended_mapping(mapping_cache)
            for eid, info in ensembl_aliases.items():
                if eid in extended and info.get('aliases'):
                    extended[eid]['aliases'] = info['aliases']
            save_mapping_file(extended, mapping_cache)
            print(f"[pubmed] Updated alias cache with {len(ensembl_aliases)} new entries")

        retried = 0
        for gene_id, gene_id_stripped, symbol_str in zero_result_genes:
            new_aliases = ensembl_aliases.get(gene_id_stripped, {}).get('aliases', [])
            if not new_aliases:
                continue

            query = _build_query(symbol_str, query_templates, aliases=new_aliases)
            try:
                pmids = client.esearch(query=query, retmax=max_abstracts)
            except Exception as e:
                print(f"[pubmed] [retry] {symbol_str}: search failed ({e})")
                continue

            if not pmids:
                continue

            retried += 1
            print(f"[pubmed] [retry] {symbol_str} (aliases: {new_aliases}): "
                  f"{len(pmids)} abstracts found!")

            try:
                xml = client.efetch_xml(pmids)
                parsed = _parse_pubmed_xml(xml)
            except Exception as e:
                print(f"[pubmed] [retry] {symbol_str}: fetch failed ({e})")
                continue

            for pmid in pmids:
                if pmid not in parsed:
                    continue
                rec = parsed[pmid]
                rows.append({
                    'gene': gene_id,
                    'pmid': rec['pmid'],
                    'year': rec['year'],
                    'title': rec['title'],
                    'abstract': rec['abstract'],
                    'snippet': _snippet(rec['abstract']),
                })

        print(f"[pubmed] Alias retry recovered abstracts for {retried} genes")

    lit_evidence = pd.DataFrame(rows) if rows else pd.DataFrame(columns=['gene', 'pmid', 'year', 'title', 'abstract', 'snippet'])

    unique_genes = lit_evidence['gene'].nunique() if not lit_evidence.empty else 0
    print(f"[pubmed] Retrieved {len(lit_evidence)} abstracts for {unique_genes} genes")
    print(f"[pubmed] Writing output to: {output_path}")
    write_csv(lit_evidence, output_path)

    print("[pubmed] Done.")
    return output_path


def run_pubmed_retry(config_path: str) -> str:
    """Retry PubMed searches only for genes with zero abstracts in existing output.

    Reads the existing lit_evidence.csv, identifies genes with no results,
    fetches aliases from the Ensembl REST API, and retries PubMed searches.
    Appends new results to the existing output file.

    Args:
        config_path: Path to the configuration YAML file.

    Returns:
        Path to the updated output CSV file.
    """
    print("[pubmed-retry] Loading configuration...")
    config = load_config(config_path)
    ensure_dirs(config)

    outputs_dir = config['paths']['outputs_dir']
    output_path = os.path.join(outputs_dir, 'lit_evidence.csv')

    if not os.path.exists(output_path):
        print(f"[pubmed-retry] No existing output at {output_path}. Run full pubmed first.")
        return output_path

    existing = pd.read_csv(output_path)
    genes_with_results = set(existing['gene'].unique()) if not existing.empty else set()

    # Get candidate list
    candidates_config = config.get('candidates', {})
    candidate_list_path = candidates_config.get(
        'output_path', os.path.join(outputs_dir, 'candidates.csv'))
    candidates_df = pd.read_csv(candidate_list_path)

    if 'gene_symbol' in candidates_df.columns:
        gene_ids = candidates_df['gene'].tolist()
        gene_symbols = candidates_df['gene_symbol'].tolist()
    else:
        gene_ids = candidates_df['gene'].tolist()
        gene_symbols = gene_ids

    # Find genes with zero results
    zero_genes = []
    for gene_id, symbol in zip(gene_ids, gene_symbols):
        symbol_str = str(symbol).strip()
        if not symbol_str or symbol_str.lower() in ('nan', 'none'):
            continue
        if gene_id not in genes_with_results:
            gene_id_stripped = gene_id.split('.')[0] if '.' in str(gene_id) else str(gene_id)
            zero_genes.append((gene_id, gene_id_stripped, symbol_str))

    if not zero_genes:
        print("[pubmed-retry] All genes already have PubMed results. Nothing to retry.")
        return output_path

    print(f"[pubmed-retry] {len(zero_genes)} genes have zero abstracts")

    # Load existing aliases
    gene_mapping_config = config.get('gene_mapping', {})
    mapping_cache = gene_mapping_config.get('cache_path', 'data/ensembl_to_symbol_cache.tsv')
    alias_map = {}
    if os.path.exists(mapping_cache):
        extended = load_extended_mapping(mapping_cache)
        for ens_id, info in extended.items():
            alias_map[ens_id] = info.get('aliases', [])

    # Separate: genes that already have aliases vs genes that need Ensembl lookup
    need_ensembl = [g for g in zero_genes if not alias_map.get(g[1], [])]
    have_aliases = [g for g in zero_genes if alias_map.get(g[1], [])]

    print(f"[pubmed-retry] {len(have_aliases)} already have aliases, "
          f"{len(need_ensembl)} need Ensembl lookup")

    # Fetch aliases from Ensembl for genes that have none
    if need_ensembl:
        ensembl_ids = [eid for _, eid, _ in need_ensembl]
        ensembl_aliases = fetch_aliases_from_ensembl(ensembl_ids)

        # Update cache
        if ensembl_aliases and os.path.exists(mapping_cache):
            ext = load_extended_mapping(mapping_cache)
            for eid, info in ensembl_aliases.items():
                if eid in ext and info.get('aliases'):
                    ext[eid]['aliases'] = info['aliases']
                    alias_map[eid] = info['aliases']
            save_mapping_file(ext, mapping_cache)
            print(f"[pubmed-retry] Updated alias cache with {len(ensembl_aliases)} entries")

    # PubMed config
    pubmed_config = config.get('pubmed', {})
    email = pubmed_config.get('email', 'user@example.com')
    tool = pubmed_config.get('tool', 'ai-capstone')
    api_key = pubmed_config.get('api_key') or os.getenv("NCBI_API_KEY")
    max_abstracts = pubmed_config.get('max_abstracts_per_gene', 5)
    query_templates = pubmed_config.get('query_templates', [])

    client = PubMedClient(email=email, tool=tool, api_key=api_key)

    new_rows = []
    recovered = 0
    for gene_id, gene_id_stripped, symbol_str in zero_genes:
        aliases = alias_map.get(gene_id_stripped, [])
        query = _build_query(symbol_str, query_templates, aliases=aliases)

        try:
            pmids = client.esearch(query=query, retmax=max_abstracts)
        except Exception as e:
            print(f"[pubmed-retry] {symbol_str}: search failed ({e})")
            continue

        if not pmids:
            continue

        recovered += 1
        print(f"[pubmed-retry] {symbol_str} (aliases: {aliases}): {len(pmids)} abstracts")

        try:
            xml = client.efetch_xml(pmids)
            parsed = _parse_pubmed_xml(xml)
        except Exception as e:
            print(f"[pubmed-retry] {symbol_str}: fetch failed ({e})")
            continue

        for pmid in pmids:
            if pmid not in parsed:
                continue
            rec = parsed[pmid]
            new_rows.append({
                'gene': gene_id,
                'pmid': rec['pmid'],
                'year': rec['year'],
                'title': rec['title'],
                'abstract': rec['abstract'],
                'snippet': _snippet(rec['abstract']),
            })

    if new_rows:
        new_df = pd.DataFrame(new_rows)
        combined = pd.concat([existing, new_df], ignore_index=True)
        write_csv(combined, output_path)
        print(f"[pubmed-retry] Added {len(new_rows)} abstracts for {recovered} genes")
    else:
        print("[pubmed-retry] No new abstracts found")

    print("[pubmed-retry] Done.")
    return output_path
