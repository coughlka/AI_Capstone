# AI Capstone: Colorectal Cancer Biomarker Evidence Synthesis

An explainable AI pipeline for synthesizing and ranking colorectal cancer biomarker candidates from multiple evidence sources, featuring LLM-powered literature analysis.

## Overview

This project integrates four evidence streams to identify and rank potential CRC biomarkers:

- **Omics evidence** — Tumor vs normal differential expression from TCGA-COAD RNA-seq (473 tumor, 41 normal samples)
- **Literature evidence** — PubMed abstract retrieval via NCBI E-Utilities with CRC-specific query templates
- **Pathway evidence** — Multi-database pathway enrichment via g:Profiler (GO:BP, KEGG, Reactome, WikiPathways)
- **LLM scoring** — Claude analyzes each gene's PubMed abstracts and scores biomarker relevance (0-100) with rationale

The pipeline produces a ranked list of 500 biomarker candidates with transparent, auditable scoring and an interactive web UI for exploration.

## Team

| Member | Module | Responsibility |
|--------|--------|----------------|
| Keith | `src/omics.py`, `src/scoring.py`, `src/llm_scoring.py`, `src/validation.py`, `evals/` | Omics analysis, LLM scoring, evaluation framework, pipeline orchestration |
| Ayan | `src/pubmed.py` | PubMed literature retrieval |
| Gabriel | `src/pathway.py` | Pathway enrichment analysis |

## Quick Start

### 1. Clone and Setup

```bash
git clone git@github.com:coughlka/AI_Capstone.git
cd AI_Capstone

# Create virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt
```

### 2. Environment Variables

Create a `web-app/.env` file with your Anthropic API key (required for LLM scoring and chat):

```
ANTHROPIC_API_KEY=sk-ant-...
```

### 3. Download Data

Download TCGA-COAD STAR counts from [UCSC Xena Browser](https://xenabrowser.net/datapages/?dataset=TCGA-COAD.star_counts.tsv&host=https%3A%2F%2Fgdc.xenahubs.net):

1. Go to the link above
2. Click "Download" to get `TCGA-COAD.star_counts.tsv`
3. Place the file in `data/TCGA-COAD.star_counts.tsv`

### 4. Run the Pipeline

```bash
python run_pipeline.py --config config/config.yaml
```

The pipeline runs 6 steps in sequence:
1. **Omics** — Differential expression analysis (Welch's t-test, BH FDR correction)
2. **PubMed** — Literature retrieval for top 500 candidate genes (~5 min)
3. **Pathway** — g:Profiler enrichment analysis
4. **LLM Scoring** — Claude scores literature relevance for each gene with `temperature=0` for deterministic results (~15 min)
5. **Scoring** — Weighted combination: omics (45%), literature (35%), pathway (20%)
6. **Validation** — Ranks 33 known CRC biomarkers against pipeline output to measure enrichment

### 5. Launch the Web UI

```bash
cd web-app
python -m uvicorn main:app --host 0.0.0.0 --port 8000
```

Open http://localhost:8000 to browse results.

## Project Structure

```
AI_Capstone/
├── run_pipeline.py              # CLI entrypoint (6-step pipeline)
├── config/
│   └── config.yaml              # Central configuration
├── src/
│   ├── utils.py                 # Shared utilities (config, CSV I/O)
│   ├── gene_mapping.py          # Ensembl -> gene symbol mapping
│   ├── omics.py                 # Differential expression analysis
│   ├── pubmed.py                # PubMed literature retrieval (NCBI E-Utilities)
│   ├── pathway.py               # Pathway enrichment (g:Profiler API)
│   ├── llm_scoring.py           # LLM-based literature relevance scoring (temperature=0)
│   ├── scoring.py               # Evidence combination and ranking
│   └── validation.py            # Known CRC biomarker validation (33 biomarkers)
├── evals/
│   ├── test_cases.py            # Gold standard test cases (real + synthetic)
│   ├── judge.py                 # LLM-as-a-judge scorer and evaluator
│   ├── consistency.py           # Multi-run consistency testing
│   ├── robustness.py            # Perturbation robustness testing
│   └── run_evals.py             # Eval suite orchestrator
├── tests/                       # Unit tests (81 tests)
├── web-app/
│   ├── main.py                  # FastAPI application
│   ├── .env                     # API keys (not committed)
│   ├── api/
│   │   ├── biomarkers.py        # Gene data API endpoints
│   │   └── chat.py              # LLM chat endpoint (Claude)
│   └── static/                  # Frontend (HTML/CSS/JS)
├── notebooks/
│   └── eda_omics.ipynb          # Exploratory data analysis
├── docs/
│   └── module_contracts.md      # Output schema documentation
├── data/                        # Input data (gitignored)
└── outputs/                     # Generated outputs (gitignored)
```

## Output Files

| File | Description |
|------|-------------|
| `outputs/omics_evidence.csv` | Differential expression results (log2FC, FDR, direction) for all 36,519 genes |
| `outputs/candidates.csv` | Top 500 candidate genes selected for downstream analysis |
| `outputs/lit_evidence.csv` | PubMed abstracts (~1,500 abstracts for ~395 genes) |
| `outputs/pathway_evidence.csv` | Pathway enrichment (gene, pathway_count, top_pathways) |
| `outputs/llm_scores.csv` | LLM relevance scores (0-100) with rationale per gene |
| `outputs/ranked_candidates.csv` | Final ranked list with combined weighted scores |
| `outputs/validation_summary.json` | Known biomarker enrichment analysis (33 CRC biomarkers) |

## Web UI Features

- **Evidence Browser** — Paginated, sortable, filterable table of all ranked genes
- **Gene Detail Modal** — Click "View" to see scores, omics evidence, pathway data, and AI literature analysis with rationale
- **LLM Chat** — Ask questions about the data using the chat panel (powered by Claude)
- **Search & Filter** — Filter by gene symbol, direction (up/down), minimum score

## Module Contracts

See [docs/module_contracts.md](docs/module_contracts.md) for detailed output schemas.

## Configuration

All pipeline parameters are in `config/config.yaml`:

- **Candidates**: `top_n` (default 500), `fdr_threshold` (0.05)
- **PubMed**: API credentials, `max_abstracts_per_gene` (5), query templates
- **Pathway**: g:Profiler sources (GO:BP, KEGG, REAC, WP), FDR threshold
- **Scoring weights**: omics (0.45), literature (0.35), pathway (0.20)

## Evaluation Framework

The `evals/` directory contains an LLM-as-a-judge evaluation suite that measures scoring quality:

```bash
python -m evals.run_evals
```

**What it tests:**
- **Score calibration** — Are numeric scores appropriate for the evidence?
- **Rubric adherence** — Does the rationale address all 4 weighted criteria?
- **Groundedness** — Are claims supported by the provided abstracts only?
- **Consistency** — Same input scored 3x produces identical results (temperature=0)
- **Robustness** — Score stability under perturbations (shuffle, remove, duplicate abstracts)

**Test cases** load real abstracts from `lit_evidence.csv` so evals reflect actual production behavior. Synthetic edge cases (empty input, conflicting evidence, fictional gene) test boundary conditions.

**Latest results:** 10/10 cases pass, mean judge score 8.4/10, consistency std=0.0, robustness 92%.

## Key Results

- **473 tumor** vs **41 normal** samples from TCGA-COAD
- **25,800 significant genes** (FDR < 0.05)
- **500 candidate genes** selected by strongest differential expression signal
- **~395 genes** with PubMed literature evidence (~1,500 full abstracts)
- **~382 genes** with pathway annotations (79% coverage via g:Profiler)
- **LLM scores** range 0-95 with mean ~49 (well-differentiated, deterministic with temperature=0)
- **Biomarker validation**: 33/33 known CRC biomarkers found, 2.1x enrichment ratio (median rank 8,812 vs expected 18,259), 9 in top 5,000
- **Top validated biomarkers**: BMP3 (rank 23), APC (rank 2,872), PTEN (rank 3,524), SMAD4 (rank 3,681), KRAS (rank 3,736)
