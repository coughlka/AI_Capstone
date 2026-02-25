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

from src.utils import load_config, ensure_dirs, write_csv

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


def _build_query(gene_symbol: str, templates: List[str]) -> str:
    """Build a PubMed query from gene symbol and templates."""
    if templates:
        parts = [f'({t.format(gene=gene_symbol)})' for t in templates]
        return " OR ".join(parts)
    return f'("{gene_symbol}"[Title/Abstract] OR "{gene_symbol}"[All Fields])'


def run_pubmed(config_path: str) -> str:
    """Run the PubMed literature search pipeline step.

    Queries PubMed for each candidate gene and retrieves abstracts.

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
        pd.DataFrame(columns=['gene', 'pmid', 'year', 'title', 'abstract']).to_csv(output_path, index=False)
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
    api_key = pubmed_config.get('api_key')
    max_abstracts = pubmed_config.get('max_abstracts_per_gene', 5)
    query_templates = pubmed_config.get('query_templates', [])

    print(f"[pubmed] Email: {email}, API key: {'yes' if api_key else 'no'}")
    print(f"[pubmed] Max abstracts per gene: {max_abstracts}")
    print(f"[pubmed] Query templates: {query_templates}")

    client = PubMedClient(email=email, tool=tool, api_key=api_key)

    rows = []
    total = len(gene_symbols)

    for i, (gene_id, symbol) in enumerate(zip(gene_ids, gene_symbols)):
        symbol_str = str(symbol).strip()
        if not symbol_str or symbol_str.lower() in ('nan', 'none'):
            continue

        query = _build_query(symbol_str, query_templates)

        try:
            pmids = client.esearch(query=query, retmax=max_abstracts)
        except Exception as e:
            print(f"[pubmed] [{i+1}/{total}] {symbol_str}: search failed ({e})")
            continue

        if not pmids:
            if (i + 1) % 50 == 0 or i == 0:
                print(f"[pubmed] [{i+1}/{total}] {symbol_str}: 0 results")
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
            })

        if (i + 1) % 50 == 0 or i == 0:
            print(f"[pubmed] [{i+1}/{total}] {symbol_str}: {len(pmids)} abstracts")

    lit_evidence = pd.DataFrame(rows) if rows else pd.DataFrame(columns=['gene', 'pmid', 'year', 'title', 'abstract'])

    unique_genes = lit_evidence['gene'].nunique() if not lit_evidence.empty else 0
    print(f"[pubmed] Retrieved {len(lit_evidence)} abstracts for {unique_genes} genes")
    print(f"[pubmed] Writing output to: {output_path}")
    write_csv(lit_evidence, output_path)

    print("[pubmed] Done.")
    return output_path
