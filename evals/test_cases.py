"""Gold standard test cases for LLM-as-a-judge evaluation.

Real test cases load abstracts from the pipeline's lit_evidence.csv so that
evals reflect the actual evidence the scorer sees in production. Synthetic
edge cases (empty input, conflicting evidence, fictional gene) remain
hardcoded because they can't come from real data.

Gene selection rationale:
  - True positives:  LY6G6D, SCARA5, DPYD — high CRC relevance in pipeline
  - Moderate:        GSTM5, HDAC9 — mixed/indirect CRC evidence
  - True negatives:  CNTN2, UCN3 — abstracts retrieved but no CRC relevance
  - Edge cases:      EMPTY_GENE, CONFLICTED_GENE, FICTIONAL_GENE_X — synthetic
"""

import os
from pathlib import Path

import pandas as pd


# ---------------------------------------------------------------------------
# Load real abstracts from pipeline output
# ---------------------------------------------------------------------------

_LIT_PATH = Path(__file__).resolve().parent.parent / "outputs" / "lit_evidence.csv"
_OMICS_PATH = Path(__file__).resolve().parent.parent / "outputs" / "omics_evidence.csv"


def _load_abstracts_for_gene(symbol: str) -> list[dict]:
    """Load real abstracts from lit_evidence.csv for a gene symbol."""
    if not _LIT_PATH.exists() or not _OMICS_PATH.exists():
        raise FileNotFoundError(
            f"Pipeline outputs required: {_LIT_PATH} and {_OMICS_PATH}. "
            "Run the pipeline first (python run_pipeline.py)."
        )

    omics = pd.read_csv(_OMICS_PATH)
    sym_to_ens = omics.set_index("gene_symbol")["gene"].to_dict()
    ensembl_id = sym_to_ens.get(symbol)
    if ensembl_id is None:
        return []

    lit = pd.read_csv(_LIT_PATH)
    rows = lit[lit["gene"] == ensembl_id]

    abstracts = []
    for _, row in rows.iterrows():
        title = str(row.get("title", ""))
        abstract = str(row.get("abstract", ""))
        if abstract and abstract != "nan":
            abstracts.append({"title": title, "abstract": abstract})
    return abstracts


# ---------------------------------------------------------------------------
# True Positives — strong CRC biomarker evidence (expect 70-100)
# ---------------------------------------------------------------------------

TRUE_POSITIVES = [
    {
        "gene_symbol": "LY6G6D",
        "category": "true_positive",
        "expected_score_range": (70, 100),
        "rationale_should_mention": ["colorectal", "antigen", "immune"],
        "abstracts": _load_abstracts_for_gene("LY6G6D"),
    },
    {
        "gene_symbol": "SCARA5",
        "category": "true_positive",
        "expected_score_range": (70, 100),
        "rationale_should_mention": ["colorectal", "tumor", "suppressor"],
        "abstracts": _load_abstracts_for_gene("SCARA5"),
    },
    {
        "gene_symbol": "DPYD",
        "category": "true_positive",
        "expected_score_range": (70, 100),
        "rationale_should_mention": ["fluorouracil", "toxicity", "colorectal"],
        "abstracts": _load_abstracts_for_gene("DPYD"),
    },
]

# ---------------------------------------------------------------------------
# Moderate Cases — mixed or indirect CRC evidence (expect 30-80)
# ---------------------------------------------------------------------------

MODERATE_CASES = [
    {
        "gene_symbol": "GSTM5",
        "category": "moderate",
        "expected_score_range": (30, 80),
        "rationale_should_mention": ["colorectal"],
        "abstracts": _load_abstracts_for_gene("GSTM5"),
    },
    {
        "gene_symbol": "HDAC9",
        "category": "moderate",
        "expected_score_range": (30, 80),
        "rationale_should_mention": ["epigenetic", "colorectal"],
        "abstracts": _load_abstracts_for_gene("HDAC9"),
    },
]

# ---------------------------------------------------------------------------
# True Negatives — abstracts exist but no CRC relevance (expect 0-25)
# ---------------------------------------------------------------------------

TRUE_NEGATIVES = [
    {
        "gene_symbol": "CNTN2",
        "category": "true_negative",
        "expected_score_range": (0, 25),
        "rationale_should_mention": [],
        "abstracts": _load_abstracts_for_gene("CNTN2"),
    },
    {
        "gene_symbol": "UCN3",
        "category": "true_negative",
        "expected_score_range": (0, 25),
        "rationale_should_mention": [],
        "abstracts": _load_abstracts_for_gene("UCN3"),
    },
    {
        "gene_symbol": "FICTIONAL_GENE_X",
        "category": "true_negative",
        "expected_score_range": (0, 15),
        "rationale_should_mention": [],
        "abstracts": [
            {
                "title": "Characterization of zinc finger protein interactions in Drosophila wing development",
                "abstract": "We report the identification of novel zinc finger protein complexes regulating wing disc patterning in Drosophila melanogaster. These transcription factors coordinate hedgehog and decapentaplegic signaling during larval imaginal disc morphogenesis, with no orthologs detected in mammalian genomes.",
            },
        ],
    },
]

# ---------------------------------------------------------------------------
# Edge Cases — test boundary conditions (synthetic, hardcoded)
# ---------------------------------------------------------------------------

EDGE_CASES = [
    {
        "gene_symbol": "EMPTY_GENE",
        "category": "edge_no_abstracts",
        "expected_score_range": (0, 10),
        "rationale_should_mention": [],
        "abstracts": [],
    },
    {
        "gene_symbol": "CONFLICTED_GENE",
        "category": "edge_conflicting",
        "expected_score_range": (20, 60),
        "rationale_should_mention": ["conflict", "inconsistent"],
        "abstracts": [
            {
                "title": "CONFLICTED_GENE as a novel CRC biomarker: a pilot study",
                "abstract": "We report that CONFLICTED_GENE is significantly overexpressed in colorectal tumor tissue compared to adjacent normal mucosa (p<0.001). High expression was associated with advanced stage and poor prognosis, suggesting a role as a diagnostic and prognostic CRC biomarker.",
            },
            {
                "title": "Re-evaluation of CONFLICTED_GENE in colorectal cancer: no independent prognostic value",
                "abstract": "In our large multicenter validation study (n=2,450), we could not replicate the previously reported association between CONFLICTED_GENE expression and CRC outcomes. After adjustment for established prognostic factors, CONFLICTED_GENE showed no independent predictive value (HR=1.02, p=0.89).",
            },
        ],
    },
]


def get_all_test_cases():
    """Return all test cases as a flat list."""
    return TRUE_POSITIVES + MODERATE_CASES + TRUE_NEGATIVES + EDGE_CASES


def get_test_cases_by_category():
    """Return test cases organized by category."""
    return {
        "true_positive": TRUE_POSITIVES,
        "moderate": MODERATE_CASES,
        "true_negative": TRUE_NEGATIVES,
        "edge_cases": EDGE_CASES,
    }
