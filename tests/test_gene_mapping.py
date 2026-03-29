"""Tests for gene_mapping module — alias support and backward compatibility."""

import os
import json
from unittest.mock import patch, MagicMock

import pandas as pd
import pytest

from src.gene_mapping import (
    strip_ensembl_version,
    load_mapping_file,
    load_extended_mapping,
    fetch_extended_from_mygene,
    fetch_aliases_from_ensembl,
    save_mapping_file,
    get_gene_symbols,
    get_extended_gene_info,
    map_ensembl_to_symbols,
)


# --- strip_ensembl_version ---

def test_strip_version_with_version():
    assert strip_ensembl_version("ENSG00000141510.18") == "ENSG00000141510"


def test_strip_version_without_version():
    assert strip_ensembl_version("ENSG00000141510") == "ENSG00000141510"


# --- load_mapping_file (3-column legacy) ---

def test_load_mapping_file_3col(sample_cache_3col):
    mapping = load_mapping_file(sample_cache_3col)
    assert len(mapping) == 3
    assert mapping["ENSG00000141510"]["gene_symbol"] == "TP53"
    assert mapping["ENSG00000141510"]["type_of_gene"] == "protein-coding"


def test_load_mapping_file_missing():
    mapping = load_mapping_file("/nonexistent/path.tsv")
    assert mapping == {}


# --- load_extended_mapping ---

def test_load_extended_4col(sample_cache_4col):
    ext = load_extended_mapping(sample_cache_4col)
    assert len(ext) == 4
    tp53 = ext["ENSG00000141510"]
    assert tp53["gene_symbol"] == "TP53"
    assert tp53["aliases"] == ["LFS1", "TRP53", "p53"]
    assert tp53["type_of_gene"] == "protein-coding"


def test_load_extended_4col_empty_aliases(sample_cache_4col):
    ext = load_extended_mapping(sample_cache_4col)
    assert ext["ENSG00000099999"]["aliases"] == []


def test_load_extended_3col_backward_compat(sample_cache_3col):
    """3-column cache should load with empty aliases."""
    ext = load_extended_mapping(sample_cache_3col)
    assert len(ext) == 3
    assert ext["ENSG00000141510"]["gene_symbol"] == "TP53"
    assert ext["ENSG00000141510"]["aliases"] == []


def test_load_extended_missing():
    assert load_extended_mapping("/nonexistent.tsv") == {}


# --- save_mapping_file ---

def test_save_extended_format(tmp_dir):
    mapping = {
        "ENSG00000141510": {
            "gene_symbol": "TP53",
            "aliases": ["LFS1", "p53"],
            "type_of_gene": "protein-coding",
        },
        "ENSG00000012048": {
            "gene_symbol": "BRCA1",
            "aliases": [],
            "type_of_gene": "protein-coding",
        },
    }
    out_path = os.path.join(tmp_dir, "saved.tsv")
    save_mapping_file(mapping, out_path)

    df = pd.read_csv(out_path, sep='\t')
    assert len(df) == 2
    assert list(df.columns) == ["ensembl_id", "gene_symbol", "aliases", "type_of_gene"]
    tp53_row = df[df["ensembl_id"] == "ENSG00000141510"].iloc[0]
    assert tp53_row["gene_symbol"] == "TP53"
    assert tp53_row["aliases"] == "LFS1;p53"


def test_save_roundtrip(tmp_dir):
    """Save then load should preserve data."""
    original = {
        "ENSG00000141510": {
            "gene_symbol": "TP53",
            "aliases": ["LFS1", "TRP53"],
            "type_of_gene": "protein-coding",
        },
    }
    path = os.path.join(tmp_dir, "roundtrip.tsv")
    save_mapping_file(original, path)
    loaded = load_extended_mapping(path)
    assert loaded["ENSG00000141510"]["gene_symbol"] == "TP53"
    assert loaded["ENSG00000141510"]["aliases"] == ["LFS1", "TRP53"]
    assert loaded["ENSG00000141510"]["type_of_gene"] == "protein-coding"


# --- fetch_extended_from_mygene ---

def test_fetch_extended_from_mygene_parses_aliases():
    """Test that aliases are correctly parsed from mygene.info response."""
    mock_response = json.dumps([
        {
            "query": "ENSG00000141510",
            "symbol": "TP53",
            "alias": ["LFS1", "TRP53", "p53"],
            "type_of_gene": "protein-coding",
        },
        {
            "query": "ENSG00000012048",
            "symbol": "BRCA1",
            "alias": "RNF53",  # single string, not list
            "type_of_gene": "protein-coding",
        },
    ])

    mock_resp = MagicMock()
    mock_resp.read.return_value = mock_response.encode('utf-8')
    mock_resp.__enter__ = lambda s: s
    mock_resp.__exit__ = MagicMock(return_value=False)

    with patch("src.gene_mapping.urllib.request.urlopen", return_value=mock_resp):
        result = fetch_extended_from_mygene(["ENSG00000141510", "ENSG00000012048"])

    assert result["ENSG00000141510"]["gene_symbol"] == "TP53"
    assert result["ENSG00000141510"]["aliases"] == ["LFS1", "TRP53", "p53"]
    assert result["ENSG00000012048"]["aliases"] == ["RNF53"]


def test_fetch_extended_handles_no_alias():
    """Genes with no aliases should get empty list."""
    mock_response = json.dumps([
        {
            "query": "ENSG00000141510",
            "symbol": "TP53",
            "type_of_gene": "protein-coding",
        },
    ])

    mock_resp = MagicMock()
    mock_resp.read.return_value = mock_response.encode('utf-8')
    mock_resp.__enter__ = lambda s: s
    mock_resp.__exit__ = MagicMock(return_value=False)

    with patch("src.gene_mapping.urllib.request.urlopen", return_value=mock_resp):
        result = fetch_extended_from_mygene(["ENSG00000141510"])

    assert result["ENSG00000141510"]["aliases"] == []


# --- fetch_aliases_from_ensembl ---

def test_fetch_aliases_from_ensembl():
    """Test Ensembl xrefs parsing for HGNC synonyms."""
    mock_xrefs = json.dumps([
        {"dbname": "HGNC", "synonyms": ["NIHCOLE"]},
        {"dbname": "EntrezGene", "synonyms": []},
    ])

    mock_resp = MagicMock()
    mock_resp.read.return_value = mock_xrefs.encode('utf-8')
    mock_resp.__enter__ = lambda s: s
    mock_resp.__exit__ = MagicMock(return_value=False)

    with patch("src.gene_mapping.urllib.request.urlopen", return_value=mock_resp):
        result = fetch_aliases_from_ensembl(["ENSG00000257151"])

    assert "ENSG00000257151" in result
    assert result["ENSG00000257151"]["aliases"] == ["NIHCOLE"]


def test_fetch_aliases_no_hgnc():
    """Gene with no HGNC xrefs should not appear in results."""
    mock_xrefs = json.dumps([
        {"dbname": "EntrezGene", "synonyms": []},
    ])

    mock_resp = MagicMock()
    mock_resp.read.return_value = mock_xrefs.encode('utf-8')
    mock_resp.__enter__ = lambda s: s
    mock_resp.__exit__ = MagicMock(return_value=False)

    with patch("src.gene_mapping.urllib.request.urlopen", return_value=mock_resp):
        result = fetch_aliases_from_ensembl(["ENSG00000257151"])

    assert len(result) == 0


# --- get_gene_symbols (backward compat) ---

def test_get_gene_symbols_returns_dict_of_dicts(sample_cache_4col):
    """get_gene_symbols should return {id: {gene_symbol, type_of_gene}}."""
    result = get_gene_symbols(
        ["ENSG00000141510.18"], cache_path=sample_cache_4col, use_api=False
    )
    assert "ENSG00000141510" in result
    assert result["ENSG00000141510"]["gene_symbol"] == "TP53"
    assert result["ENSG00000141510"]["type_of_gene"] == "protein-coding"


# --- get_extended_gene_info ---

def test_get_extended_gene_info(sample_cache_4col):
    result = get_extended_gene_info(
        ["ENSG00000141510.18"], cache_path=sample_cache_4col, use_api=False
    )
    info = result["ENSG00000141510"]
    assert info["gene_symbol"] == "TP53"
    assert info["aliases"] == ["LFS1", "TRP53", "p53"]
    assert info["type_of_gene"] == "protein-coding"


def test_get_extended_gene_info_legacy_cache(sample_cache_3col):
    """3-column cache should work with empty aliases."""
    result = get_extended_gene_info(
        ["ENSG00000141510"], cache_path=sample_cache_3col, use_api=False
    )
    info = result["ENSG00000141510"]
    assert info["gene_symbol"] == "TP53"
    assert info["aliases"] == []


# --- map_ensembl_to_symbols ---

def test_map_ensembl_to_symbols(sample_cache_4col):
    df = pd.DataFrame({"gene": ["ENSG00000141510.18", "ENSG00000012048.23"]})
    result = map_ensembl_to_symbols(df, cache_path=sample_cache_4col, include_type=True)
    assert result["gene_symbol"].tolist() == ["TP53", "BRCA1"]
    assert result["type_of_gene"].tolist() == ["protein-coding", "protein-coding"]


def test_map_ensembl_to_symbols_without_type(sample_cache_4col):
    df = pd.DataFrame({"gene": ["ENSG00000141510.18"]})
    result = map_ensembl_to_symbols(df, cache_path=sample_cache_4col, include_type=False)
    assert "gene_symbol" in result.columns
    assert "type_of_gene" not in result.columns
