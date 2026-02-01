# AI Capstone: Colorectal Cancer Biomarker Evidence Synthesis

An explainable AI pipeline for synthesizing and ranking colorectal cancer biomarker candidates from multiple evidence sources.

## Overview

This project integrates:
- **Omics evidence** — Tumor vs normal differential expression from TCGA-COAD RNA-seq
- **Literature evidence** — PubMed abstract mining for gene-disease associations
- **Pathway evidence** — Pathway enrichment from Reactome/KEGG

The pipeline produces a ranked list of biomarker candidates with transparent, auditable scoring.

## Team

| Member | Module | Responsibility |
|--------|--------|----------------|
| Keith | `src/omics.py`, `src/scoring.py` | Omics analysis, pipeline orchestration |
| Ayan | `src/pubmed.py` | Literature evidence extraction |
| Gabriel | `src/pathway.py` | Pathway mapping and enrichment |

## Quick Start

### 1. Clone and Setup

```bash
git clone git@github.com:coughlka/AI_Capstone.git
cd AI_Capstone

# Create virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install pyyaml pandas numpy scipy
```

### 2. Download Data

Download TCGA-COAD STAR counts from [UCSC Xena Browser](https://xenabrowser.net/datapages/?dataset=TCGA-COAD.star_counts.tsv&host=https%3A%2F%2Fgdc.xenahubs.net):

1. Go to the link above
2. Click "Download" to get `TCGA-COAD.star_counts.tsv`
3. Place the file in `data/TCGA-COAD.star_counts.tsv`

### 3. Run the Pipeline

```bash
python run_pipeline.py --config config/config.yaml
```

## Project Structure

```
AI_Capstone/
├── run_pipeline.py          # CLI entrypoint
├── config/
│   └── config.yaml          # Central configuration
├── src/
│   ├── utils.py             # Shared utilities
│   ├── gene_mapping.py      # Ensembl → gene symbol mapping
│   ├── omics.py             # Differential expression analysis
│   ├── pubmed.py            # Literature mining (stub)
│   ├── pathway.py           # Pathway enrichment (stub)
│   └── scoring.py           # Evidence combination and ranking
├── docs/
│   └── module_contracts.md  # Output schema documentation
├── data/                    # Input data (gitignored)
└── outputs/                 # Generated outputs (gitignored)
```

## Output Files

| File | Description |
|------|-------------|
| `outputs/omics_evidence.csv` | Differential expression results (log2FC, FDR, direction) |
| `outputs/candidates.csv` | Top 500 candidate genes for downstream analysis |
| `outputs/lit_evidence.csv` | Literature evidence (Ayan) |
| `outputs/pathway_evidence.csv` | Pathway evidence (Gabriel) |
| `outputs/ranked_candidates.csv` | Final ranked list with combined scores |

## Module Contracts

See [docs/module_contracts.md](docs/module_contracts.md) for detailed output schemas.

### For Ayan (pubmed.py)

Your module should:
1. Read candidate genes from `outputs/candidates.csv`
2. Query PubMed for colorectal cancer associations
3. Output `outputs/lit_evidence.csv` with columns:
   - `gene`, `pmid`, `year`, `study_type`, `role`, `sample_type`, `directionality`, `snippet`

### For Gabriel (pathway.py)

Your module should:
1. Read candidate genes from `outputs/candidates.csv`
2. Query Reactome/KEGG for pathway membership
3. Output `outputs/pathway_evidence.csv` with columns:
   - `gene`, `pathway_count`, `top_pathways`

## Current Results (Omics)

- **473 tumor** vs **41 normal** samples
- **25,800 significant genes** (FDR < 0.05)
- Known CRC genes validated: APC, KRAS, SMAD4, PIK3CA, BRAF

## Dependencies

- Python 3.9+
- pandas
- numpy
- scipy
- pyyaml
