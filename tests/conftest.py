"""Shared fixtures for the AI Capstone test suite.

Fixtures provide realistic test data that mirrors actual pipeline outputs,
including known CRC genes so we can validate ranking behavior.
"""

import os
import tempfile

import numpy as np
import pandas as pd
import pytest
import yaml


@pytest.fixture
def tmp_dir(tmp_path):
    """Provide a temporary directory for test outputs."""
    return str(tmp_path)


@pytest.fixture
def sample_config(tmp_path):
    """Minimal valid pipeline config dict and corresponding YAML file."""
    outputs_dir = str(tmp_path / "outputs")
    config = {
        "paths": {
            "data_dir": str(tmp_path / "data"),
            "outputs_dir": outputs_dir,
        },
        "omics": {
            "counts_path": str(tmp_path / "data" / "counts.tsv"),
            "cohort": "TCGA-COAD",
            "dataset_label": "test-dataset",
        },
        "gene_mapping": {
            "cache_path": str(tmp_path / "data" / "cache.tsv"),
            "use_api": False,
        },
        "candidates": {
            "top_n": 500,
            "fdr_threshold": 0.05,
            "output_path": str(tmp_path / "outputs" / "candidates.csv"),
        },
        "pubmed": {
            "email": "test@example.com",
            "tool": "test-tool",
            "api_key": None,
            "max_abstracts_per_gene": 5,
            "query_templates": [
                "colorectal cancer {gene}",
                "colon cancer {gene}",
            ],
        },
        "pathway": {
            "source": "gprofiler",
            "fdr_threshold": 0.05,
            "sources": ["GO:BP", "KEGG", "REAC", "WP"],
        },
        "scoring": {
            "weights": {
                "omics": 0.45,
                "literature": 0.35,
                "pathway": 0.20,
            }
        },
    }

    config_path = str(tmp_path / "config.yaml")
    with open(config_path, "w") as f:
        yaml.dump(config, f)

    return config, config_path


@pytest.fixture
def sample_omics_df():
    """10 genes with realistic DE values. Includes KRAS and APC (known CRC genes)
    which should have strong signals, and several non-CRC genes with weak signals.
    """
    return pd.DataFrame([
        # Known CRC genes — strong signal
        {"gene": "ENSG00000133703", "gene_symbol": "KRAS", "log2fc": 2.5, "p_value": 1e-20, "fdr": 1e-18, "direction": "up", "tumor_mean": 8.5, "normal_mean": 6.0, "dataset": "test"},
        {"gene": "ENSG00000134982", "gene_symbol": "APC", "log2fc": -3.1, "p_value": 1e-25, "fdr": 1e-23, "direction": "down", "tumor_mean": 3.2, "normal_mean": 6.3, "dataset": "test"},
        {"gene": "ENSG00000141510", "gene_symbol": "TP53", "log2fc": 1.8, "p_value": 1e-15, "fdr": 1e-13, "direction": "up", "tumor_mean": 7.8, "normal_mean": 6.0, "dataset": "test"},
        # Moderate signal
        {"gene": "ENSG00000112715", "gene_symbol": "VEGFA", "log2fc": 1.2, "p_value": 1e-8, "fdr": 1e-6, "direction": "up", "tumor_mean": 6.2, "normal_mean": 5.0, "dataset": "test"},
        {"gene": "ENSG00000148773", "gene_symbol": "MKI67", "log2fc": -0.9, "p_value": 1e-5, "fdr": 1e-3, "direction": "down", "tumor_mean": 4.1, "normal_mean": 5.0, "dataset": "test"},
        # Weak / non-significant
        {"gene": "ENSG00000111640", "gene_symbol": "GAPDH", "log2fc": 0.1, "p_value": 0.5, "fdr": 0.8, "direction": "up", "tumor_mean": 5.1, "normal_mean": 5.0, "dataset": "test"},
        {"gene": "ENSG00000012048", "gene_symbol": "BRCA1", "log2fc": -0.2, "p_value": 0.3, "fdr": 0.7, "direction": "down", "tumor_mean": 4.8, "normal_mean": 5.0, "dataset": "test"},
        {"gene": "ENSG00000099999", "gene_symbol": "GENE_A", "log2fc": 0.5, "p_value": 0.01, "fdr": 0.1, "direction": "up", "tumor_mean": 5.5, "normal_mean": 5.0, "dataset": "test"},
        {"gene": "ENSG00000088888", "gene_symbol": "GENE_B", "log2fc": -0.3, "p_value": 0.2, "fdr": 0.6, "direction": "down", "tumor_mean": 4.7, "normal_mean": 5.0, "dataset": "test"},
        {"gene": "ENSG00000077777", "gene_symbol": "GENE_C", "log2fc": 0.05, "p_value": 0.9, "fdr": 0.95, "direction": "up", "tumor_mean": 5.05, "normal_mean": 5.0, "dataset": "test"},
    ])


@pytest.fixture
def sample_lit_df():
    """Literature evidence for 5 genes with realistic abstracts."""
    rows = []
    genes_with_lit = [
        ("ENSG00000133703", "KRAS mutations in colorectal cancer", "KRAS mutations are found in 40% of CRC cases and drive resistance to anti-EGFR therapy."),
        ("ENSG00000134982", "APC and Wnt signaling in CRC", "Truncating APC mutations activate Wnt/beta-catenin signaling in colorectal cancer initiation."),
        ("ENSG00000141510", "p53 pathway in colon cancer", "TP53 loss of function promotes CRC progression through impaired DNA damage response."),
        ("ENSG00000112715", "VEGFA in tumor angiogenesis", "VEGFA promotes angiogenesis in various solid tumors including colorectal cancer."),
        ("ENSG00000148773", "Ki-67 proliferation index", "MKI67 is a general proliferation marker used across multiple cancer types."),
    ]
    for gene_id, title, abstract in genes_with_lit:
        rows.append({"gene": gene_id, "pmid": f"PM{hash(gene_id) % 100000}", "year": "2023", "title": title, "abstract": abstract})
    return pd.DataFrame(rows)


@pytest.fixture
def sample_pathway_df():
    """Pathway evidence for 5 genes."""
    return pd.DataFrame([
        {"gene": "ENSG00000133703", "pathway_count": 15, "top_pathways": "MAPK signaling (KEGG); RAS signaling (REAC)"},
        {"gene": "ENSG00000134982", "pathway_count": 22, "top_pathways": "Wnt signaling (KEGG); beta-catenin (GO:BP)"},
        {"gene": "ENSG00000141510", "pathway_count": 30, "top_pathways": "p53 signaling (KEGG); Apoptosis (REAC)"},
        {"gene": "ENSG00000112715", "pathway_count": 10, "top_pathways": "VEGF signaling (KEGG)"},
        {"gene": "ENSG00000148773", "pathway_count": 3, "top_pathways": "Cell cycle (GO:BP)"},
    ])


@pytest.fixture
def sample_llm_scores_df():
    """LLM scores for 5 genes with realistic distributions."""
    return pd.DataFrame([
        {"gene": "ENSG00000133703", "gene_symbol": "KRAS", "llm_score": 92, "rationale": "Strong CRC evidence with EGFR resistance."},
        {"gene": "ENSG00000134982", "gene_symbol": "APC", "llm_score": 95, "rationale": "Key Wnt pathway driver in CRC initiation."},
        {"gene": "ENSG00000141510", "gene_symbol": "TP53", "llm_score": 88, "rationale": "Major tumor suppressor in CRC progression."},
        {"gene": "ENSG00000112715", "gene_symbol": "VEGFA", "llm_score": 55, "rationale": "Angiogenesis role, not CRC-specific."},
        {"gene": "ENSG00000148773", "gene_symbol": "MKI67", "llm_score": 35, "rationale": "General proliferation marker."},
    ])


@pytest.fixture
def sample_ranked_df():
    """Pre-computed ranked candidates (10 genes, sorted by final_score desc)."""
    return pd.DataFrame([
        {"gene": "ENSG00000134982", "final_score": 85.2, "omics_score": 100.0, "literature_score": 95.0, "pathway_score": 22.0},
        {"gene": "ENSG00000133703", "final_score": 80.1, "omics_score": 90.0, "literature_score": 92.0, "pathway_score": 15.0},
        {"gene": "ENSG00000141510", "final_score": 72.5, "omics_score": 75.0, "literature_score": 88.0, "pathway_score": 30.0},
        {"gene": "ENSG00000112715", "final_score": 45.3, "omics_score": 40.0, "literature_score": 55.0, "pathway_score": 10.0},
        {"gene": "ENSG00000148773", "final_score": 30.2, "omics_score": 25.0, "literature_score": 35.0, "pathway_score": 3.0},
        {"gene": "ENSG00000099999", "final_score": 20.1, "omics_score": 15.0, "literature_score": 0.0, "pathway_score": 0.0},
        {"gene": "ENSG00000111640", "final_score": 10.5, "omics_score": 5.0, "literature_score": 0.0, "pathway_score": 0.0},
        {"gene": "ENSG00000012048", "final_score": 8.2, "omics_score": 3.0, "literature_score": 0.0, "pathway_score": 0.0},
        {"gene": "ENSG00000088888", "final_score": 5.1, "omics_score": 2.0, "literature_score": 0.0, "pathway_score": 0.0},
        {"gene": "ENSG00000077777", "final_score": 1.0, "omics_score": 0.5, "literature_score": 0.0, "pathway_score": 0.0},
    ])
