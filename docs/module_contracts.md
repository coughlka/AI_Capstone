# Module Contracts

This document describes the required output schemas for each pipeline module.

## Directory Structure

- `data/` - Input data files (gitignored, do not commit)
- `outputs/` - Pipeline output files (gitignored, do not commit)

## Output Schemas

### outputs/omics_evidence.csv

Gene expression statistics from RNA-seq data.

| Column | Type | Description |
|--------|------|-------------|
| gene | string | Gene identifier |
| mean_expr | float | Mean expression across samples |
| var_expr | float | Variance of expression across samples |
| dataset | string | Dataset label from config |

### outputs/lit_evidence.csv

Literature evidence from PubMed searches.

| Column | Type | Description |
|--------|------|-------------|
| gene | string | Gene symbol |
| pmid | string | PubMed ID |
| year | int | Publication year |
| study_type | string | Type of study (e.g., clinical, in vitro) |
| role | string | Gene's role described in the paper |
| sample_type | string | Sample type used in the study |
| directionality | string | Direction of effect (up/down/mixed) |
| snippet | string | Relevant text excerpt |

### outputs/pathway_evidence.csv

Pathway enrichment data from pathway databases.

| Column | Type | Description |
|--------|------|-------------|
| gene | string | Gene symbol |
| pathway_count | int | Number of pathways the gene belongs to |
| top_pathways | string | Comma-separated list of top pathway names |

### outputs/ranked_candidates.csv

Final ranked list of candidate genes.

| Column | Type | Description |
|--------|------|-------------|
| gene | string | Gene identifier |
| final_score | float | Weighted composite score (0-100) |
| omics_score | float | Normalized omics evidence score (0-100) |
| literature_score | float | Normalized literature evidence score (0-100) |
| pathway_score | float | Normalized pathway evidence score (0-100) |

Rows are sorted by `final_score` in descending order.

## Notes

- The `data/` and `outputs/` directories are gitignored and should not be committed
- Place input data files in `data/` before running the pipeline
- All output files are regenerated on each pipeline run
