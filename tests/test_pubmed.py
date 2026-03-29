"""Tests for pubmed module — alias-aware queries, XML parsing, and retry logic."""

import os
import json
import tempfile
from unittest.mock import patch, MagicMock

import pandas as pd
import pytest

from src.pubmed import (
    _build_query,
    _snippet,
    _strip,
    _parse_pubmed_xml,
    run_pubmed,
    run_pubmed_retry,
)


# --- _strip ---

def test_strip_normalizes_whitespace():
    assert _strip("  hello   world  ") == "hello world"


def test_strip_handles_none():
    assert _strip(None) == ""


# --- _snippet ---

def test_snippet_short_text():
    assert _snippet("Short text") == "Short text"


def test_snippet_truncates():
    text = "word " * 100  # 500 chars
    result = _snippet(text, max_len=50)
    assert len(result) <= 55  # 50 + room for "..."
    assert result.endswith("...")


# --- _build_query ---

def test_build_query_no_aliases():
    result = _build_query("TP53", ["colorectal cancer {gene}"])
    assert result == "(colorectal cancer TP53)"


def test_build_query_with_templates_no_aliases():
    result = _build_query("BMP3", ["colorectal cancer {gene}", "colon cancer {gene}"])
    assert result == "(colorectal cancer BMP3) OR (colon cancer BMP3)"


def test_build_query_with_aliases():
    result = _build_query("NIHCOLE", ["colorectal cancer {gene}"], aliases=["LINC02163"])
    assert "(colorectal cancer NIHCOLE)" in result
    assert "(colorectal cancer LINC02163)" in result


def test_build_query_aliases_deduped():
    """Alias matching primary symbol (case-insensitive) should be skipped."""
    result = _build_query("TP53", ["cancer {gene}"], aliases=["tp53", "LFS1"])
    assert result.count("TP53") == 1  # not duplicated
    assert "LFS1" in result


def test_build_query_no_templates_with_aliases():
    result = _build_query("GENE1", [], aliases=["ALT1"])
    assert '"GENE1"' in result
    assert '"ALT1"' in result


def test_build_query_no_templates_no_aliases():
    result = _build_query("GENE1", [])
    assert '"GENE1"[Title/Abstract]' in result


# --- _parse_pubmed_xml ---

SAMPLE_XML = """<?xml version="1.0"?>
<PubmedArticleSet>
  <PubmedArticle>
    <MedlineCitation>
      <PMID>12345678</PMID>
      <Article>
        <ArticleTitle>TP53 in colorectal cancer</ArticleTitle>
        <Abstract>
          <AbstractText>This study examines TP53 mutations in CRC patients.</AbstractText>
        </Abstract>
        <Journal>
          <Title>Cancer Research</Title>
          <JournalIssue>
            <PubDate><Year>2023</Year></PubDate>
          </JournalIssue>
        </Journal>
      </Article>
      <MeshHeadingList>
        <MeshHeading>
          <DescriptorName>Colorectal Neoplasms</DescriptorName>
        </MeshHeading>
      </MeshHeadingList>
    </MedlineCitation>
  </PubmedArticle>
</PubmedArticleSet>"""


def test_parse_xml_basic():
    parsed = _parse_pubmed_xml(SAMPLE_XML)
    assert "12345678" in parsed
    rec = parsed["12345678"]
    assert rec["title"] == "TP53 in colorectal cancer"
    assert "TP53 mutations" in rec["abstract"]
    assert rec["year"] == "2023"
    assert rec["journal"] == "Cancer Research"
    assert "Colorectal Neoplasms" in rec["mesh_terms"]


STRUCTURED_ABSTRACT_XML = """<?xml version="1.0"?>
<PubmedArticleSet>
  <PubmedArticle>
    <MedlineCitation>
      <PMID>99999999</PMID>
      <Article>
        <ArticleTitle>Structured abstract test</ArticleTitle>
        <Abstract>
          <AbstractText Label="BACKGROUND">Background text here.</AbstractText>
          <AbstractText Label="METHODS">Methods text here.</AbstractText>
          <AbstractText Label="RESULTS">Results text here.</AbstractText>
        </Abstract>
        <Journal>
          <Title>Test Journal</Title>
          <JournalIssue>
            <PubDate><Year>2024</Year></PubDate>
          </JournalIssue>
        </Journal>
      </Article>
    </MedlineCitation>
  </PubmedArticle>
</PubmedArticleSet>"""


def test_parse_xml_structured_abstract():
    parsed = _parse_pubmed_xml(STRUCTURED_ABSTRACT_XML)
    rec = parsed["99999999"]
    assert "BACKGROUND: Background text here." in rec["abstract"]
    assert "METHODS: Methods text here." in rec["abstract"]
    assert "RESULTS: Results text here." in rec["abstract"]


def test_parse_xml_empty():
    assert _parse_pubmed_xml("") == {}
    assert _parse_pubmed_xml("   ") == {}


MEDLINE_DATE_XML = """<?xml version="1.0"?>
<PubmedArticleSet>
  <PubmedArticle>
    <MedlineCitation>
      <PMID>11111111</PMID>
      <Article>
        <ArticleTitle>MedlineDate test</ArticleTitle>
        <Abstract><AbstractText>Abstract.</AbstractText></Abstract>
        <Journal>
          <Title>J Test</Title>
          <JournalIssue>
            <PubDate><MedlineDate>2020 Jan-Feb</MedlineDate></PubDate>
          </JournalIssue>
        </Journal>
      </Article>
    </MedlineCitation>
  </PubmedArticle>
</PubmedArticleSet>"""


def test_parse_xml_medline_date():
    parsed = _parse_pubmed_xml(MEDLINE_DATE_XML)
    assert parsed["11111111"]["year"] == "2020"


# --- run_pubmed integration ---

def test_run_pubmed_no_candidates(tmp_dir):
    """run_pubmed with missing candidate file creates empty output."""
    config_path = os.path.join(tmp_dir, "config.yaml")
    with open(config_path, "w") as f:
        f.write(f"""
paths:
  data_dir: {tmp_dir}
  outputs_dir: {tmp_dir}
candidates:
  output_path: {os.path.join(tmp_dir, 'candidates.csv')}
pubmed:
  email: test@test.com
  tool: test
  max_abstracts_per_gene: 5
  query_templates:
    - "colorectal cancer {{gene}}"
gene_mapping:
  cache_path: {os.path.join(tmp_dir, 'cache.tsv')}
""")
    result = run_pubmed(config_path)
    assert os.path.exists(result)
    df = pd.read_csv(result)
    assert list(df.columns) == ['gene', 'pmid', 'year', 'title', 'abstract', 'snippet']
    assert len(df) == 0


def test_run_pubmed_with_mocked_api(tmp_dir):
    """run_pubmed with mocked PubMed API returns expected results."""
    # Create candidates file
    candidates_path = os.path.join(tmp_dir, "candidates.csv")
    pd.DataFrame({
        "gene": ["ENSG00000141510"],
        "gene_symbol": ["TP53"],
    }).to_csv(candidates_path, index=False)

    # Create config
    config_path = os.path.join(tmp_dir, "config.yaml")
    with open(config_path, "w") as f:
        f.write(f"""
paths:
  data_dir: {tmp_dir}
  outputs_dir: {tmp_dir}
candidates:
  output_path: {candidates_path}
pubmed:
  email: test@test.com
  tool: test
  max_abstracts_per_gene: 5
  query_templates:
    - "colorectal cancer {{gene}}"
gene_mapping:
  cache_path: {os.path.join(tmp_dir, 'cache.tsv')}
""")

    mock_search_resp = MagicMock()
    mock_search_resp.status_code = 200
    mock_search_resp.json.return_value = {"esearchresult": {"idlist": ["12345678"]}}
    mock_search_resp.raise_for_status = MagicMock()

    mock_fetch_resp = MagicMock()
    mock_fetch_resp.status_code = 200
    mock_fetch_resp.text = SAMPLE_XML
    mock_fetch_resp.raise_for_status = MagicMock()

    with patch("src.pubmed.requests.Session") as mock_session_cls:
        mock_session = MagicMock()
        mock_session.get.side_effect = [mock_search_resp, mock_fetch_resp]
        mock_session_cls.return_value = mock_session

        result = run_pubmed(config_path)

    df = pd.read_csv(result)
    assert len(df) == 1
    assert df.iloc[0]["gene"] == "ENSG00000141510"
    assert df.iloc[0]["title"] == "TP53 in colorectal cancer"
    assert "abstract" in df.columns
    assert "snippet" in df.columns
    assert "TP53 mutations" in df.iloc[0]["abstract"]


def test_run_pubmed_with_aliases(tmp_dir):
    """Aliases from cache should be used in queries."""
    # Create cache with aliases
    cache_path = os.path.join(tmp_dir, "cache.tsv")
    pd.DataFrame([{
        "ensembl_id": "ENSG00000257151",
        "gene_symbol": "LINC02163",
        "aliases": "NIHCOLE",
        "type_of_gene": "ncRNA",
    }]).to_csv(cache_path, sep="\t", index=False)

    candidates_path = os.path.join(tmp_dir, "candidates.csv")
    pd.DataFrame({
        "gene": ["ENSG00000257151"],
        "gene_symbol": ["LINC02163"],
    }).to_csv(candidates_path, index=False)

    config_path = os.path.join(tmp_dir, "config.yaml")
    with open(config_path, "w") as f:
        f.write(f"""
paths:
  data_dir: {tmp_dir}
  outputs_dir: {tmp_dir}
candidates:
  output_path: {candidates_path}
pubmed:
  email: test@test.com
  tool: test
  max_abstracts_per_gene: 5
  query_templates:
    - "colorectal cancer {{gene}}"
gene_mapping:
  cache_path: {cache_path}
""")

    queries_captured = []

    mock_search_resp = MagicMock()
    mock_search_resp.status_code = 200
    mock_search_resp.json.return_value = {"esearchresult": {"idlist": []}}
    mock_search_resp.raise_for_status = MagicMock()

    def capture_get(url, params=None, timeout=None):
        if params and "term" in params:
            queries_captured.append(params["term"])
        return mock_search_resp

    with patch("src.pubmed.requests.Session") as mock_session_cls:
        mock_session = MagicMock()
        mock_session.get.side_effect = capture_get
        mock_session_cls.return_value = mock_session

        run_pubmed(config_path)

    # The query should include both LINC02163 and its alias NIHCOLE
    assert len(queries_captured) >= 1
    q = queries_captured[0]
    assert "LINC02163" in q
    assert "NIHCOLE" in q


# --- run_pubmed_retry ---

def test_run_pubmed_retry_no_existing_output(tmp_dir):
    """Retry should exit gracefully when no existing output."""
    config_path = os.path.join(tmp_dir, "config.yaml")
    with open(config_path, "w") as f:
        f.write(f"""
paths:
  data_dir: {tmp_dir}
  outputs_dir: {tmp_dir}
candidates:
  output_path: {os.path.join(tmp_dir, 'candidates.csv')}
gene_mapping:
  cache_path: {os.path.join(tmp_dir, 'cache.tsv')}
""")
    result = run_pubmed_retry(config_path)
    assert result.endswith("lit_evidence.csv")


def test_run_pubmed_retry_all_genes_have_results(tmp_dir):
    """Retry should do nothing when all genes have results."""
    # Create existing output
    lit_path = os.path.join(tmp_dir, "lit_evidence.csv")
    pd.DataFrame({
        "gene": ["ENSG00000141510"],
        "pmid": ["12345678"],
        "year": ["2023"],
        "title": ["Test"],
        "abstract": ["Abstract text"],
        "snippet": ["Abstract text"],
    }).to_csv(lit_path, index=False)

    # Create candidates (same gene)
    candidates_path = os.path.join(tmp_dir, "candidates.csv")
    pd.DataFrame({
        "gene": ["ENSG00000141510"],
        "gene_symbol": ["TP53"],
    }).to_csv(candidates_path, index=False)

    config_path = os.path.join(tmp_dir, "config.yaml")
    with open(config_path, "w") as f:
        f.write(f"""
paths:
  data_dir: {tmp_dir}
  outputs_dir: {tmp_dir}
candidates:
  output_path: {candidates_path}
gene_mapping:
  cache_path: {os.path.join(tmp_dir, 'cache.tsv')}
pubmed:
  email: test@test.com
  tool: test
  max_abstracts_per_gene: 5
  query_templates:
    - "cancer {{gene}}"
""")
    result = run_pubmed_retry(config_path)
    df = pd.read_csv(result)
    assert len(df) == 1  # unchanged
