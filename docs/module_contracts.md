# Module Contracts

This document describes the output schemas for each pipeline module.

## Directory Structure

- `data/` - Input data files (gitignored, do not commit)
- `outputs/` - Pipeline output files (gitignored, do not commit)

## Pipeline Flow

```
omics.py -> pubmed.py -> pathway.py -> llm_scoring.py -> scoring.py
```

Each step reads from `config/config.yaml` and writes to `outputs/`.

## Output Schemas

### outputs/omics_evidence.csv

Tumor vs normal differential expression evidence from TCGA-COAD RNA-seq data.

| Column | Type | Description |
|--------|------|-------------|
| gene | string | Gene identifier (Ensembl ID with version) |
| gene_symbol | string | HGNC gene symbol (mapped via Ensembl API) |
| log2fc | float | Log2 fold change (tumor - normal on log2 scale) |
| p_value | float | Welch's t-test p-value |
| fdr | float | Benjamini-Hochberg FDR-adjusted p-value |
| direction | string | Direction of change: 'up' or 'down' |
| tumor_mean | float | Mean log2(count+1) expression in tumor samples |
| normal_mean | float | Mean log2(count+1) expression in normal samples |
| dataset | string | Dataset label from config |

Rows are sorted by FDR (ascending). Contains all 36,519 genes.

### outputs/candidates.csv

Top 500 candidate genes selected by strongest differential expression signal (`|log2FC| * -log10(FDR)`).

| Column | Type | Description |
|--------|------|-------------|
| gene | string | Gene identifier (Ensembl ID) |
| gene_symbol | string | HGNC gene symbol |
| log2fc | float | Log2 fold change |
| fdr | float | FDR-adjusted p-value |
| direction | string | 'up' or 'down' |

### outputs/lit_evidence.csv

PubMed abstracts retrieved via NCBI E-Utilities (esearch + efetch).

| Column | Type | Description |
|--------|------|-------------|
| gene | string | Gene identifier (Ensembl ID) |
| pmid | string | PubMed ID |
| year | string | Publication year |
| title | string | Article title |
| snippet | string | First 300 characters of abstract |

Up to `max_abstracts_per_gene` (default 5) abstracts per gene, queried using CRC-specific templates (e.g., "colorectal cancer {gene}").

### outputs/pathway_evidence.csv

Pathway enrichment results from g:Profiler, querying GO:BP, KEGG, Reactome, and WikiPathways.

| Column | Type | Description |
|--------|------|-------------|
| gene | string | Gene identifier (Ensembl ID) |
| pathway_count | int | Number of significant pathways the gene appears in |
| top_pathways | string | Pipe-separated list of top pathway names |

### outputs/llm_scores.csv

LLM-based literature relevance scores. Claude analyzes each gene's PubMed abstracts and scores biomarker relevance on a 0-100 scale.

| Column | Type | Description |
|--------|------|-------------|
| gene | string | Gene identifier (Ensembl ID) |
| gene_symbol | string | HGNC gene symbol |
| llm_score | int | Biomarker relevance score (0-100) |
| rationale | string | 1-2 sentence explanation of the score |

Scoring criteria: direct CRC evidence (40%), biomarker potential (30%), mechanistic insight (20%), evidence quality (10%).

### outputs/ranked_candidates.csv

Final ranked list combining all evidence streams.

| Column | Type | Description |
|--------|------|-------------|
| gene | string | Gene identifier (Ensembl ID) |
| final_score | float | Weighted composite score (0-100) |
| omics_score | float | Normalized omics evidence score (0-100) |
| literature_score | float | LLM-based literature score (0-100), falls back to count-based if unavailable |
| pathway_score | float | Normalized pathway evidence score (0-100) |

Rows are sorted by `final_score` in descending order. Default weights: omics (45%), literature (35%), pathway (20%).

## Notes

- The `data/` and `outputs/` directories are gitignored and should not be committed
- Place input data files in `data/` before running the pipeline
- All output files are regenerated on each pipeline run
- LLM scoring requires an `ANTHROPIC_API_KEY` in `web-app/.env`
- PubMed queries are rate-limited (3 req/sec with API key, 1 req/sec without)
