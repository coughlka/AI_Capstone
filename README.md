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
| Keith Coughlin | `src/omics.py`, `src/scoring.py`, `src/llm_scoring.py`, `src/ml_importance.py`, `src/validation.py`, `run_pipeline.py` | Omics, LLM scoring, ML classifier, pipeline orchestration, MMR integration |
| Ayan Choudhury | `src/pubmed.py`, `src/mmr/` | PubMed literature retrieval, MMR deficiency classifier |
| Gabriel Saenz | `src/pathway.py`, `web-app/`, `terraform/` | Pathway enrichment, FastAPI web app, Azure / Terraform deployment, CI/CD |

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

### 4. Run the Main Pipeline

```bash
python run_pipeline.py --config config/config.yaml
```

The pipeline runs 6 steps in sequence:
1. **Omics** — Differential expression analysis (Welch's t-test, BH FDR correction)
2. **PubMed** — Literature retrieval for top 500 candidate genes
3. **Pathway** — g:Profiler enrichment analysis
4. **LLM Scoring** — Claude scores literature relevance for each gene (content-hash cached)
5. **Scoring** — Weighted combination: omics (40%), literature (30%), pathway (20%), ML importance (10%)
6. **Validation** — Enrichment ratio vs 80 known CRC biomarkers (FDA/guideline + COSMIC + literature tiers)

### 5. Run the MMR Deficiency Classifier (optional, separate pipeline)

Requires `data/TCGA-COAD.GDC_phenotype.tsv` in addition to the expression matrix.

```bash
python run_mmr_classifier.py --config config/config.yaml
```

Trains LR / RF / XGBoost classifiers on IHC-labeled MMR deficiency status, runs SHAP on the best model, and compares SHAP top features against the ranked biomarker list.

### 6. Launch the Web UI

```bash
cd web-app
python -m uvicorn main:app --host 0.0.0.0 --port 8000
```

Open http://localhost:8000 to browse results.

## Project Structure

```
AI_Capstone/
├── run_pipeline.py              # CLI entrypoint (6-step main pipeline)
├── run_mmr_classifier.py        # CLI entrypoint (3-step MMR pipeline)
├── config/
│   └── config.yaml              # Central configuration
├── src/
│   ├── utils.py                 # Shared utilities (config, CSV I/O)
│   ├── gene_mapping.py          # Ensembl -> gene symbol mapping
│   ├── omics.py                 # Differential expression analysis
│   ├── pubmed.py                # PubMed literature retrieval (NCBI E-Utilities)
│   ├── pathway.py               # Pathway enrichment (g:Profiler API)
│   ├── llm_scoring.py           # LLM-based literature relevance scoring
│   ├── ml_importance.py         # Supervised ML classifier, feature importance stream
│   ├── scoring.py               # Evidence combination and ranking
│   ├── validation.py            # Enrichment vs 80 known CRC biomarkers
│   └── mmr/
│       ├── validate_labels.py   # Phenotype label QA
│       ├── train_classifier.py  # LR/RF/XGB + SHAP on MMR-deficiency labels
│       └── compare_shap.py      # SHAP top features vs ranked biomarker overlap
├── web-app/
│   ├── main.py                  # FastAPI application
│   ├── .env                     # API keys (not committed)
│   ├── api/
│   │   ├── biomarkers.py        # Gene data API endpoints
│   │   └── chat.py              # LLM chat endpoint (rate limited + validated)
│   └── static/                  # Frontend (HTML/CSS/JS)
├── terraform/                   # Azure infrastructure as code
├── .github/workflows/           # CI/CD (Docker build + ACR push on merge to main)
├── notebooks/
│   └── eda_omics.ipynb          # Exploratory data analysis
├── docs/                        # Report drafts, module contracts, supplementary notes
├── deliverables/                # Final submission bundle (README.txt, staged docs)
├── data/                        # Input data (gitignored)
└── outputs/                     # Generated outputs (tracked for canonical result CSVs)
```

## Output Files

| File | Description |
|------|-------------|
| `outputs/omics_evidence.csv` | Differential expression results (log2FC, FDR, direction) for all 36,519 genes |
| `outputs/candidates.csv` | Top 500 candidate genes selected for downstream analysis |
| `outputs/lit_evidence.csv` | PubMed abstracts per gene (query + retrieval log) |
| `outputs/pathway_evidence.csv` | Pathway enrichment (gene, pathway_count, top_pathways) |
| `outputs/llm_scores.csv` | LLM relevance scores (0-100) with rationale per gene |
| `outputs/ml_importance.csv` | Supervised RF/LR feature importance per gene |
| `outputs/ranked_candidates.csv` | Final ranked list with combined weighted scores |
| `outputs/validation_report.csv` | Per-gene validation result against 80-biomarker panel |
| `outputs/validation_summary.json` | Enrichment ratio, tier breakdowns, top-N recall |
| `outputs/coad_mmr_classifier/` | MMR classifier outputs (CV metrics, holdout, SHAP, overlap) |

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
- **Scoring weights**: omics (0.40), literature (0.30), pathway (0.20), ml_importance (0.10)
- **MMR classifier**: expression path, phenotype path, label column, CV folds, top-k features

## Key Results

- **473 tumor** vs **41 normal** samples from TCGA-COAD
- **500 candidate genes** selected by strongest differential expression signal
- **Validation enrichment: 1.46x** above random against an 80-gene known-biomarker panel
- **79 / 80** known biomarkers recovered from the expression matrix
- **4 known biomarkers** in the top 500
- Best Tier 1 FDA / guideline biomarker: **BMP3 at rank #14**
- Best Tier 3 literature biomarker: **FOXQ1 at rank #53**
- **MMR classifier**: 340 labeled samples (284 pMMR / 56 dMMR), logistic regression holdout ROC-AUC 0.737, tuned F1 0.56
