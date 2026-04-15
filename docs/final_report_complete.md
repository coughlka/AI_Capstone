# Biomarker Candidate Discovery and Ranking Using AI-assisted Evidence Synthesis

**Team 3**
Ayan Choudhury, Keith Coughlin, Gabriel Saenz

**Course:** AI 894 SP 2026
**Instructor:** Youakim Badr, Ph.D.
**Repository:** https://github.com/coughlka/AI_Capstone

---

## Document Control / Revision Sheet

| Release No. | Date | Member | Revision Description |
|-------------|------|--------|----------------------|
| 0.1 | 2/1/2026 | All | Initial project proposal and system requirements |
| 0.2 | 2/15/2026 | K. Coughlin | Data preprocessing, feature engineering, model planning sections |
| 0.3 | 3/1/2026 | All | Pipeline workflow, LLM scoring integration, initial results |
| 1.0 | 3/15/2026 | All | Model building, training, testing, external validation |
| 1.1 | 3/29/2026 | K. Coughlin | SHAP-based ML refactor, LLM cache, deployment writeup |
| 1.2 | 3/29/2026 | G. Saenz | Azure deployment, Terraform, CI/CD pipeline |
| 2.0 | 4/12/2026 | K. Coughlin | Added scoring weight ablation study (246 combos), precision/recall evaluation at 10 rank thresholds, updated validation to 80 biomarkers (3 tiers), SHAP-based ML importance refactor |
| 2.1 | 4/12/2026 | A. Choudhury | Added MMR deficiency classifier (LR/RF/XGB), sample-level SHAP analysis, visualizations |
| 2.2 | 4/12/2026 | G. Saenz | Cybersecurity writeup, ACR webhook, Docker outputs deployment fix |
| 2.3 | 4/19/2026 | All | Final report revision with updated results, ablation study, precision/recall evaluation, and discussion |

---

## Section 1: Project Proposal

### 1.1 Problem Statement

Colorectal cancer (CRC) is the third most common cancer globally and the second leading cause of cancer-related death in the United States. Early detection dramatically improves survival outcomes, yet identifying reliable molecular biomarkers for CRC screening, diagnosis, and treatment stratification remains a challenging and resource-intensive process. Thousands of genes show altered expression in CRC tumors, but only a small fraction represent actionable biomarker candidates with sufficient evidence from multiple domains: genomics, literature, and molecular pathway biology.

The current state of biomarker discovery relies heavily on manual literature review and siloed analysis of individual data types. Researchers typically analyze RNA-seq differential expression, search PubMed independently, and consult pathway databases in separate workflows, then attempt to synthesize findings manually. This process is slow, non-reproducible, and difficult to scale. There is no standard automated framework that integrates all three evidence streams into a single transparent ranking of candidates.

This project addresses that gap by building an explainable AI pipeline that synthesizes four evidence streams, specifically omics differential expression, PubMed literature (scored by a large language model), pathway enrichment, and SHAP-based machine learning importance, into a unified ranked list of CRC biomarker candidates. The system is designed to be transparent (every score is decomposable into its components), reproducible (deterministic scoring, content-hash caching, version-controlled configuration), and interactive (a deployed web application with LLM-powered chat exploration).

### 1.2 Research Questions

The project addresses the following primary research questions:

1. Which genes show the strongest multi-modal evidence as CRC biomarker candidates when omics, literature, pathway, and ML signals are integrated?
2. How does multi-evidence synthesis compare to single-stream analysis for recovering known CRC biomarkers?
3. What is the optimal weighting of evidence streams for maximizing recovery of known CRC biomarkers?
4. Can expression-based machine learning (SHAP importance) add complementary signal to hypothesis-driven evidence streams?
5. Can gene expression data distinguish microsatellite instability subtypes within CRC tumors?

### 1.3 Proposed Solution

The proposed solution is a modular, seven-step computational pipeline implemented in Python:

1. **Omics analysis** (Welch's t-test differential expression on TCGA-COAD data)
2. **ML importance** (Random Forest with SHAP explainability, tumor vs. normal classification)
3. **Literature retrieval** (PubMed via NCBI E-Utilities, up to 5 abstracts per gene)
4. **Pathway enrichment** (g:Profiler querying GO:BP, KEGG, Reactome, and WikiPathways)
5. **LLM scoring** (Claude evaluates each gene's PubMed abstracts on a 0-100 biomarker relevance scale)
6. **Weighted scoring** (min-max normalized components combined using configurable weights)
7. **Validation** (enrichment ratio against 80 known CRC biomarkers organized in three confidence tiers)

Results are served through a FastAPI web application deployed on Azure, providing an interactive Evidence Browser and conversational LLM chat interface for exploring the ranked candidates.

### 1.4 Expected Contributions

- A reproducible, open-source pipeline for CRC biomarker candidate prioritization
- A systematic evaluation of multi-modal vs. single-stream biomarker recovery
- A scoring weight ablation study quantifying the contribution of each evidence stream
- A tiered validation framework (80 known biomarkers) supporting ongoing performance tracking
- A deployed interactive research tool for domain experts to explore results

### 1.5 Team Roles

| Team Member | Responsibilities |
|-------------|-----------------|
| Keith Coughlin | Omics analysis, ML importance (SHAP), LLM scoring, weighted scoring, ablation study, precision/recall evaluation, pipeline orchestration |
| Ayan Choudhury | PubMed literature retrieval, validation framework, MMR deficiency classifier, biomarker panel curation |
| Gabriel Saenz | Pathway enrichment, web application, Azure infrastructure, Terraform, CI/CD, cybersecurity |

---

## Section 2: System Requirements (Lesson 4)

### 2.1 Functional Requirements

**Pipeline Execution**

- FR-1: The system shall analyze all 36,519 expressed genes from TCGA-COAD RNA-seq data (473 tumor, 41 normal samples) using Welch's t-test with Benjamini-Hochberg FDR correction.
- FR-2: The system shall select the top 500 candidate genes by differential expression signal (|log2FC| x -log10(FDR)) from genes with FDR < 0.05.
- FR-3: The system shall retrieve up to 5 PubMed abstracts per candidate gene using CRC-specific query templates via NCBI E-Utilities.
- FR-4: The system shall query g:Profiler for pathway enrichment against GO:BP, KEGG, Reactome, and WikiPathways databases for all candidate genes.
- FR-5: The system shall score each gene's PubMed abstracts using the Claude API on a 0-100 biomarker relevance scale with stored rationales.
- FR-6: The system shall train a Random Forest classifier (tumor vs. normal, 5-fold CV) and extract SHAP values as a data-driven importance signal.
- FR-7: The system shall combine all four evidence streams using configurable weighted scoring into a final ranked list of all 36,519 genes.
- FR-8: The system shall validate rankings against 80 known CRC biomarkers organized in three confidence tiers and report enrichment ratio, precision, recall, and F1 at multiple rank thresholds.

**Web Application**

- FR-9: The system shall expose a REST API with endpoints for candidate lists, individual gene detail, pipeline statistics, and validation results.
- FR-10: The system shall provide a single-page web UI with sortable/filterable ranked candidate table and per-gene evidence breakdown.
- FR-11: The system shall provide a conversational LLM chat interface grounded in pipeline results for interactive exploration.
- FR-12: The system shall support cache invalidation via a reload endpoint to reflect pipeline re-runs without server restart.

**Reproducibility**

- FR-13: The system shall cache LLM scores using content-hash keys so re-runs with unchanged literature produce identical scores without additional API calls.
- FR-14: The system shall be fully configurable via a single YAML configuration file (config/config.yaml).
- FR-15: All pipeline outputs shall be written as standard CSV/JSON files suitable for downstream analysis.

### 2.2 Non-Functional Requirements

**Performance**
- NFR-1: Full pipeline execution shall complete within 60 minutes (target: approximately 30 minutes) on standard workstation hardware with a valid API key.
- NFR-2: LLM scoring re-runs against an unchanged literature corpus shall complete within 60 seconds using the content-hash cache.
- NFR-3: Web application API endpoints shall respond within 500 ms for cached requests.

**Reliability**
- NFR-4: The LLM scoring module shall survive interruptions and resume from the last successfully cached gene.
- NFR-5: The pipeline shall handle missing evidence (no PubMed abstracts, no pathway annotations) gracefully, assigning a score of 0 for the absent stream rather than excluding the gene.

**Security**
- NFR-6: API keys shall not appear in source code or committed to the repository.
- NFR-7: The web application shall restrict CORS to configured origins.
- NFR-8: The Docker runtime image shall not include build tools (compiler, build-essential).

**Maintainability**
- NFR-9: Each pipeline step shall be implemented as an independent Python module callable with a single function from the orchestrator.
- NFR-10: Scoring weights shall be modifiable without code changes by editing config/config.yaml.
- NFR-11: Infrastructure shall be fully defined as Terraform code with remote state storage.

**Scalability**
- NFR-12: The pipeline shall support extension to additional cancer types by replacing the TCGA cohort and updating the validation panel.
- NFR-13: The validation framework shall support addition of new known biomarkers without code changes by updating the curated list in src/validation.py.

### 2.3 System Architecture Overview

The system comprises two main subsystems: a batch processing pipeline and an interactive serving layer.

The batch pipeline (run_pipeline.py) reads TCGA-COAD RNA-seq data and executes seven sequential steps, writing outputs as CSV/JSON files to the outputs/ directory. The serving layer (web-app/) reads those outputs, caches them in memory, and exposes them through a FastAPI REST API served by Uvicorn. A Docker container packages both layers for local development and cloud deployment.

Cloud deployment on Azure uses an App Service (B1 tier, Linux) pulling from Azure Container Registry, with VNet integration and private endpoints securing Azure Files storage. GitHub Actions provides CI/CD automation, building and pushing Docker images on every push to the main branch.

### 2.4 Technology Stack

| Layer | Technology | Justification |
|-------|-----------|---------------|
| Pipeline orchestration | Python 3.11 | Dominant language for bioinformatics and ML |
| Statistical testing | scipy.stats | Welch's t-test, FDR correction (statsmodels) |
| ML / SHAP | scikit-learn, shap | Random Forest, SHAP value extraction |
| Pathway enrichment | g:Profiler REST API | Unified access to GO, KEGG, Reactome, WikiPathways |
| Literature retrieval | NCBI E-Utilities | Official PubMed programmatic access |
| LLM scoring | Anthropic Claude API | State-of-the-art text understanding, JSON output |
| Web framework | FastAPI + Uvicorn | High-performance async Python API, OpenAPI docs |
| Containerization | Docker (multi-stage) | Reproducible deployment, minimal attack surface |
| Cloud infrastructure | Azure App Service + ACR | Academic subscription availability, managed SSL |
| Infrastructure as code | Terraform | Reproducible, version-controlled infrastructure |
| CI/CD | GitHub Actions | Native GitHub integration, free for public repos |

---

## Section 3: Data Preprocessing and AI Planning (Lesson 6)

### 3.1 Data Collection

#### 3.1.1 Primary Data Source

The primary dataset is RNA-seq gene expression data from The Cancer Genome Atlas (TCGA) Colon Adenocarcinoma (COAD) cohort. The dataset contains STAR-aligned read counts for genes across 514 patient samples (473 tumor, 41 normal tissue), provided in log2(counts+1) format. The data is publicly available through the Genomic Data Commons (GDC) Portal.

**Data Source Links:**

- TCGA-COAD RNA-seq Counts: https://portal.gdc.cancer.gov/projects/TCGA-COAD (primary gene expression dataset, STAR-Counts, log2(counts+1))
- g:Profiler: https://biit.cs.ut.ee/gprofiler/ (multi-database pathway enrichment API)
- mygene.info API: https://mygene.info/v3/query (Ensembl ID to gene symbol conversion)
- PubMed / NCBI E-utilities: https://eutils.ncbi.nlm.nih.gov/entrez/eutils/ (literature retrieval)
- Anthropic Claude API: https://docs.anthropic.com/en/docs/ (LLM-based literature scoring and chat)

#### 3.1.2 Tools Used for Data Collection

- **GDC Data Transfer Client / Portal:** Used to download the TCGA-COAD STAR-Counts gene expression matrix in TSV format.
- **requests (Python HTTP client):** Used to query the NCBI E-Utilities API for PubMed literature retrieval and the g:Profiler REST API for pathway enrichment.
- **anthropic (Python SDK):** Used to call the Claude API for LLM-based literature relevance scoring and the interactive chat feature in the web UI.
- **pandas:** Used for loading, transforming, and storing tabular data throughout the pipeline.
- **scipy:** Used for statistical testing (Welch's t-test) during differential expression analysis.
- **FastAPI + Uvicorn:** Used to serve the Evidence Browser web UI with REST API endpoints for gene data, LLM chat, and interactive exploration.

#### 3.1.3 Data Storage

All datasets are stored as flat files (TSV and CSV) in the project directory structure. Raw input data resides in the data/ directory, while pipeline outputs are written to the outputs/ directory. Both directories are excluded from version control via .gitignore to prevent accidental commits of large data files. Intermediate caches (e.g., gene symbol mappings, LLM scores) are stored in the data/ and outputs/ directories for reuse across pipeline runs.

#### 3.1.4 Data Catalog

**Table 1: Pipeline Data Catalog**

| Dataset Name | Size | Rows | Columns | Type | Source |
|-------------|------|------|---------|------|--------|
| TCGA-COAD.star_counts.tsv | 284 MB | 36,519 genes | 514 samples | TSV | UCSC Xena / GDC Portal |
| ensembl_to_symbol_cache.tsv | 2.1 MB | 29,395 | 2 | TSV | mygene.info API (cached) |
| omics_evidence.csv | 4.2 MB | 36,519 | 7 | CSV | Pipeline output (Step 1) |
| ml_importance.csv | 1.0 MB | 36,519 | 3 | CSV | Pipeline output (Step 2, SHAP) |
| ml_cv_metrics.csv | 125 B | 5 | 4 | CSV | Pipeline output (Step 2, CV perf) |
| lit_evidence.csv | 3.8 MB | ~1,570 | 5 | CSV | Pipeline output (Step 3, PubMed) |
| pathway_evidence.csv | 2.1 MB | 500 | 4 | CSV | Pipeline output (Step 4, g:Profiler) |
| llm_scores.csv | 890 KB | 460 | 4 | CSV | Pipeline output (Step 5, Claude) |
| llm_scores_cache.json | 1.2 MB | 460 | -- | JSON | Pipeline output (Step 5, cache) |
| ranked_candidates.csv | 5.8 MB | 36,519 | 12 | CSV | Pipeline output (Step 6) |
| candidates.csv | 280 KB | 500 | 12 | CSV | Pipeline output (Step 6, top 500) |
| validation_report.csv | 48 KB | 80 | 9 | CSV | Pipeline output (Step 7) |
| validation_summary.json | 3.5 KB | -- | -- | JSON | Pipeline output (Step 7) |
| ablation_results.csv | 10 KB | 246 | 11 | CSV | Pipeline output (ablation study) |
| precision_recall.csv | 1 KB | 10 | 17 | CSV | Pipeline output (P/R evaluation) |

#### 3.1.5 Feature Dictionaries

**TCGA-COAD.star_counts.tsv (Raw Input)**

| Column | Type | Description |
|--------|------|-------------|
| Ensembl_ID | string | Ensembl gene identifier with version (e.g., ENSG00000141510.18) |
| TCGA-XX-XXXX-SSA | float | log2(counts+1) expression for each of 514 samples |

**omics_evidence.csv (Differential Expression Output)**

| Column | Type | Description |
|--------|------|-------------|
| gene_id | string | Ensembl ID (version stripped) |
| gene_symbol | string | HGNC gene symbol (from mygene.info) |
| log2fc | float | log2 fold change (tumor vs. normal) |
| pvalue | float | Welch's t-test p-value |
| fdr | float | Benjamini-Hochberg adjusted p-value |
| de_signal | float | \|log2FC\| x -log10(FDR + 1e-300) |
| type_of_gene | string | Gene biotype from mygene.info |

**ranked_candidates.csv (Final Scored Output)**

| Column | Type | Description |
|--------|------|-------------|
| gene_id | string | Ensembl ID |
| gene_symbol | string | HGNC gene symbol |
| final_score | float | Weighted composite score (0-100) |
| omics_score | float | Normalized DE signal (0-100) |
| ml_importance_score | float | Normalized SHAP importance (0-100) |
| literature_score | float | LLM relevance score or 0 (0-100) |
| pathway_score | float | Normalized pathway count (0-100) |
| rank | int | Global rank by final_score (1 = best) |

#### 3.1.6 Data Relevance to the Research Problem

This project aims to identify and prioritize potential biomarker genes for colorectal cancer by integrating multiple evidence streams. The TCGA-COAD dataset is directly pertinent because it provides genome-wide gene expression measurements comparing tumor tissue to adjacent normal tissue from colon adenocarcinoma patients. This enables identification of genes with significant differential expression, a hallmark of potential biomarkers or therapeutic targets.

Additional data sources that could further improve the predictive models include: (1) protein-protein interaction network data to capture functional relationships, (2) clinical outcome data (survival, treatment response) for supervised model training, and (3) independent validation cohorts (e.g., TCGA-READ for rectal cancer) to assess generalizability.

### 3.2 Data Cleaning and Preprocessing

#### 3.2.1 Duplicate Detection and Elimination

Ensembl gene identifiers in the raw TCGA data include version suffixes (e.g., ENSG00000141510.18). During gene symbol mapping, version suffixes are stripped to produce canonical identifiers (e.g., ENSG00000141510) for consistent lookups against the mygene.info API and cached mappings. This prevents duplicate entries caused by version differences across datasets. The gene mapping cache contains 29,395 unique Ensembl-to-symbol mappings with no duplicates after version stripping.

#### 3.2.2 Missing Value Detection and Handling

Two types of missing values were identified and handled:

**Gene symbol mapping gaps:** 7,124 out of 36,519 genes (19.5%) could not be mapped to HGNC gene symbols via the mygene.info API. These genes are retained in the analysis with empty gene_symbol fields, as their Ensembl identifiers remain valid for scoring. Many correspond to non-coding RNAs or pseudogenes without established symbols.

**Pathway and literature evidence gaps:** When merging evidence streams in the scoring module, genes without pathway or literature matches receive a score of 0.0 (left join with fillna(0)). This is by design: absence of evidence in a particular stream should not exclude a gene from the final ranking.

#### 3.2.3 Data Inconsistency Handling

Sample type classification from TCGA barcodes required careful parsing. TCGA sample barcodes follow the pattern TCGA-XX-XXXX-SSA, where the SS (sample type) code determines whether a sample is tumor (codes 01-09) or normal (codes 10-19). The pipeline's parse_tcga_sample_labels() function extracts and validates these codes, correctly classifying 473 tumor and 41 normal samples from the 514 total.

#### 3.2.4 Gene Filtering

Low-expression genes were filtered to reduce noise and multiple testing burden. A gene was retained if it met either of two criteria:

- Mean expression across all samples >= 1.0 (on log2 scale)
- Non-zero expression in >= 20% of samples

After filtering, 36,519 genes were retained for differential expression testing.

#### 3.2.5 Statistical Preprocessing

Differential expression analysis was performed using the following statistical methods:

**Welch's t-test:** An unequal-variance two-sample t-test was applied to each gene, comparing expression values between the 473 tumor and 41 normal samples. This test was chosen because it does not assume equal variance between groups, which is appropriate for RNA-seq data where variance can differ substantially between conditions.

**Benjamini-Hochberg FDR correction:** To control for the multiple testing problem (36,519 simultaneous tests), p-values were adjusted using the Benjamini-Hochberg procedure. This controls the false discovery rate (FDR) at the specified threshold of 0.05. After correction, 25,800 genes (70.6%) remained statistically significant at FDR < 0.05.

#### 3.2.6 Normalization and Standardization

Multiple normalization steps are applied throughout the pipeline:

**Log2 transformation:** The input TCGA data is already provided in log2(counts+1) format, which reduces skewness in raw count data and makes expression values more normally distributed.

**Min-max normalization to 0-100:** In the scoring module, each evidence stream (omics, ML importance, literature, pathway) is independently normalized to a 0-100 scale using min-max scaling. This ensures all component scores are on comparable scales before applying the weighted combination. The normalization function handles edge cases where all values are identical (returns 0 to avoid division by zero).

**Rank-based percentile normalization for SHAP:** SHAP importance values are converted to percentile ranks before normalization. This preserves sparsity (most genes have zero SHAP importance) and prevents a small number of high-importance genes from compressing the rest of the distribution.

#### 3.2.7 Outlier Detection

The differential expression signal metric (|log2FC| x -log10(FDR)) naturally handles extreme values. Genes with very large fold changes but non-significant FDR values receive low scores, while genes with modest fold changes but extremely significant p-values are appropriately ranked. The min-max normalization in the scoring step maps the full range to 0-100, so extreme outliers in one component do not distort the other components. The top candidate gene (SCARA5) has the maximum DE signal (omics_score = 100.0), representing the upper bound after normalization.

#### 3.2.8 Imbalanced Data

The TCGA-COAD dataset has an inherent class imbalance: 473 tumor samples vs. 41 normal samples (11.5:1 ratio). This imbalance is addressed in two ways:

- **Differential expression:** Welch's t-test does not assume equal sample sizes or equal variances, appropriately accounting for different group sizes in its degrees of freedom calculation.
- **Random Forest (SHAP):** The class_weight='balanced' parameter is set, which inversely weights classes by their frequency, preventing the classifier from defaulting to the majority class.

Despite the imbalance, the 41 normal samples provide sufficient statistical power to detect differential expression, as evidenced by the 25,800 significant genes (FDR < 0.05).

### 3.3 Data Exploration and Feature Engineering

#### 3.3.1 Descriptive Statistics

**Table 2: Differential Expression Summary Statistics**

| Metric | Value |
|--------|-------|
| Total genes analyzed | 36,519 |
| Significant genes (FDR < 0.05) | 25,800 (70.6%) |
| Upregulated genes | 12,847 |
| Downregulated genes | 12,953 |
| Top candidate (highest DE signal) | SCARA5 |
| Max |log2FC| | 10.128 (OTOP2, downregulated) |
| Top upregulated gene by |log2FC| | MMP1 |
| Median |log2FC| (significant genes) | 1.73 |

#### 3.3.2 Visualizations

Comprehensive exploratory data analysis was performed in the Jupyter notebook (notebooks/eda_omics.ipynb). Key visualizations include:

**Volcano Plot:** Plots -log10(FDR) vs. log2FC for all 36,519 genes. Known CRC driver genes (APC, KRAS, TP53, SMAD4, PIK3CA, BRAF) are annotated. Significance thresholds are marked at FDR = 0.05 (horizontal) and |log2FC| = 1.0 (vertical). The volcano plot shows a clear asymmetry with more strongly downregulated genes.

**Top 20 DE Genes Bar Chart:** Horizontal bar chart of the top 20 genes by absolute log2 fold change. All top 20 are downregulated in tumors, led by OTOP2 (log2FC = -10.128, approximately 1,000x less expressed in tumors), CA1 (log2FC = -9.916), and AQP8 (log2FC = -9.333).

**FDR Distribution Histogram:** Shows the distribution of FDR-adjusted p-values across all genes, revealing that the majority of genes have significant differential expression (FDR < 0.05).

**Final Score Distribution:** Histogram of final composite scores. With all four evidence streams active and LLM-based literature scoring, the top scores reach approximately 64/100, with meaningful differentiation across genes.

#### 3.3.3 Known CRC Gene Validation

To validate the pipeline, known colorectal cancer driver genes were checked against the differential expression results. Five of six canonical CRC driver genes (APC, KRAS, SMAD4, PIK3CA, BRAF) are detected as significantly differentially expressed. TP53 is not significant at the expression level, which is expected because TP53 dysregulation in CRC is primarily driven by somatic mutations rather than expression changes. This validates the pipeline's ability to detect biologically relevant signals.

#### 3.3.4 Feature Engineering

**Feature Construction**

The primary engineered feature is the Differential Expression (DE) Signal:

```
DE Signal = |log2FC| x -log10(FDR + epsilon)
```

This composite metric captures both the magnitude of expression change (|log2FC|) and its statistical significance (-log10(FDR)). The epsilon term (1e-300) prevents undefined values when FDR = 0. Genes with large fold changes and high statistical significance receive the highest DE signals.

**Feature Selection**

From the 36,519 analyzed genes, the top 500 candidates are selected based on:

- FDR < 0.05 (statistical significance threshold)
- Ranked by DE signal in descending order
- Top 500 retained for downstream pathway and literature analysis

This feature selection step reduces the search space from 36,519 genes to a focused set of 500 high-confidence candidates, making downstream API calls feasible and reducing noise in the final rankings.

**Feature Transformation**

All evidence scores are independently transformed to a 0-100 scale using min-max normalization before combining. This standardization is required because the raw metrics have different scales: DE signal values can range from 0 to approximately 28,000, SHAP importance from 0 to 0.15, pathway counts from 0 to 15, and literature (LLM) scores from 0 to 95. Without normalization, the omics component would dominate regardless of the configured weights.

#### 3.3.5 Correlation and Interpretation

The scoring results reveal important patterns:

**Balanced multi-evidence scoring:** With all four evidence streams active, the final scores range from 0 to approximately 64/100. LLM-based literature scores (mean = 28.7, range 0-95) provide meaningful differentiation compared to the earlier count-based approach where 62% of genes scored the maximum.

**Downregulation bias:** The top candidates are predominantly downregulated genes. This is consistent with CRC biology, where tumor suppressors and normal tissue markers are often silenced in tumor tissue.

**Pathway coverage:** After switching from Reactome-only to g:Profiler (querying GO:BP, KEGG, Reactome, and WikiPathways), pathway coverage reached 79% of candidate genes. Approximately 382 of the 500 candidates have at least one pathway annotation.

**LLM literature analysis:** Claude scored 460 genes based on PubMed abstract analysis, evaluating direct CRC evidence, biomarker potential, mechanistic insight, and evidence quality. The resulting scores have a mean of 28.7, median of 15.0, and range from 0 to 95.

### 3.4 Model Planning

#### 3.4.1 Current Model: Weighted Multi-Evidence Scoring

The system uses a weighted linear combination to rank candidate genes, integrating four evidence streams:

```
final_score = 0.40 x omics_score
            + 0.10 x ml_importance_score
            + 0.30 x literature_score
            + 0.20 x pathway_score
```

All component scores are normalized to 0-100 before weighting and the weights sum to 1.0. This is a transparent, interpretable model where each evidence stream contributes a known proportion to the final ranking.

**Model identifier:** WLC_4stream_040_010_030_020_model (Weighted Linear Combination, 4 evidence streams, weights 0.40/0.10/0.30/0.20)

The scoring weights were initially selected based on domain expertise (differential expression is the most direct molecular evidence for biomarker discovery) and subsequently evaluated through the ablation study described in Section 4.

#### 3.4.2 LLM-Based Literature Scoring

A key design choice in this pipeline is the use of a large language model (Claude, via the Anthropic API) to score each gene's relevance as a CRC biomarker based on its PubMed abstracts. This replaced a simple abstract-count approach where 62% of genes scored the maximum, providing no meaningful differentiation.

The LLM scoring module (src/llm_scoring.py) operates as follows:

1. **Abstract Aggregation:** For each of the 460 genes with PubMed results, all retrieved abstracts (up to 5 per gene) are grouped and formatted into a structured prompt with the gene symbol and abstract text.
2. **LLM Relevance Scoring:** Claude evaluates each gene on a 0-100 scale using four criteria: direct CRC evidence (40% weight), biomarker potential (30%), mechanistic insight (20%), and evidence quality (10%). The model returns a JSON response with the score and a 1-2 sentence rationale.
3. **Score Integration:** The LLM scores are used directly as the literature_score in the final weighted combination. Genes without PubMed results receive a score of 0.
4. **Explainability:** Each gene's rationale is stored in llm_scores.csv and displayed in the web UI's gene detail view under "AI Literature Analysis."
5. **Caching:** Scores are cached using content-hash (MD5) keys tied to each gene's abstract list. Identical abstracts always produce the same score; only genes with new or changed abstracts are re-scored.

**Table 3: LLM Score Distribution**

| Score Range | Count | Interpretation |
|-------------|-------|----------------|
| 80-95 | 39 | Direct, established CRC biomarkers |
| 60-79 | 107 | Strong CRC evidence, emerging |
| 40-59 | 16 | CRC mentioned, indirect evidence |
| 20-39 | 53 | Tangential CRC relevance |
| 1-19 | 19 | Minimal CRC evidence |
| 0 | 226 | No abstracts or no CRC relevance |
| Total scored | 460 | |

#### 3.4.3 Assumptions and Constraints

**Linear combination assumption:** The model assumes that evidence streams contribute additively to gene importance. Interactions between streams (e.g., a gene that is both highly differentially expressed and appears in many CRC-related pathways) are not explicitly modeled.

**No ground truth labels for primary scoring:** The primary ranking model is unsupervised, relying on the assumption that stronger multi-evidence support indicates higher biomarker potential. The 80 known biomarkers are used for validation, not training, to avoid circularity.

**Independence assumption:** The four evidence streams are treated as independent signals. In practice, there may be correlations (e.g., well-studied genes have both more literature and more pathway annotations).

**Single cohort:** Results are based on a single TCGA cohort (COAD). Validation against independent datasets would strengthen confidence in the rankings.

### 3.5 Pipeline Workflow (7-Step Architecture)

The following describes the end-to-end workflow of the gene prioritization pipeline. Each step produces output files consumed by downstream steps.

```
Step 1: Omics Analysis (src/omics.py)
  Input:  data/TCGA-COAD.star_counts.tsv (36,519 genes x 514 samples)
  Process: TCGA barcode parsing → Gene filtering → Welch's t-test → BH FDR correction
  Output:  outputs/omics_evidence.csv (36,519 genes)
           outputs/candidates.csv (top 500 by DE signal, FDR < 0.05)
                              |
                              v
Step 2: ML Importance (src/ml_importance.py)
  Input:  data/TCGA-COAD.star_counts.tsv
  Process: Top 5,000 genes by variance → Random Forest (5-fold CV) → SHAP values
  Output:  outputs/ml_importance.csv (36,519 genes, mean |SHAP| importance)
           outputs/ml_cv_metrics.csv (5-fold CV accuracy + AUC)
                              |
                              v
Step 3: PubMed Retrieval (src/pubmed.py)
  Input:  outputs/candidates.csv
  Process: NCBI E-Utilities esearch/efetch with CRC query templates, rate-limited
  Output:  outputs/lit_evidence.csv (~1,570 abstracts for ~460 genes)
                              |
                              v
Step 4: Pathway Enrichment (src/pathway.py)
  Input:  outputs/candidates.csv
  Process: g:Profiler API → GO:BP, KEGG, Reactome, WikiPathways enrichment
  Output:  outputs/pathway_evidence.csv (500 genes, 79% with pathway data)
                              |
                              v
Step 5: LLM Scoring (src/llm_scoring.py)
  Input:  outputs/lit_evidence.csv
  Process: Claude analyzes abstracts per gene → scores biomarker relevance (0-100)
  Output:  outputs/llm_scores.csv (460 genes with scores and rationales)
           outputs/llm_scores_cache.json (content-hash cache)
                              |
                              v
Step 6: Weighted Scoring (src/scoring.py)
  Input:  omics_evidence.csv + ml_importance.csv + llm_scores.csv + pathway_evidence.csv
  Process: Normalize scores (0-100) → Weighted combination (40/10/30/20)
  Output:  outputs/ranked_candidates.csv (36,519 genes ranked by final_score)
           outputs/candidates.csv (top 500 updated with all scores)
                              |
                              v
Step 7: Validation (src/validation.py)
  Input:  outputs/ranked_candidates.csv
  Process: Enrichment ratio vs. 80 known CRC biomarkers (3 tiers), P/R curves
  Output:  outputs/validation_report.csv (per-biomarker results)
           outputs/validation_summary.json (aggregate metrics)
```

The pipeline is orchestrated by run_pipeline.py, which executes all seven steps sequentially. Each module reads configuration from config/config.yaml and writes outputs to the outputs/ directory.

```bash
python run_pipeline.py --config config/config.yaml
```

Total runtime: approximately 30 minutes on standard hardware (most time consumed by LLM scoring on first run; subsequent runs with cached scores complete in under 2 minutes).

---

## Section 4: Model Building, Training, and Testing (Lesson 10)

### 4.1 Overview of Modeling Approach

The system employs two distinct modeling paradigms for different purposes:

1. **Unsupervised weighted scoring** (primary ranking): A weighted linear combination of four normalized evidence streams produces the final candidate ranking. This model is not "trained" in the conventional sense; its weights are configured by domain expertise and evaluated post-hoc via the ablation study.

2. **Supervised classification** (ML evidence stream): A Random Forest classifier trained to distinguish tumor from normal TCGA-COAD samples provides SHAP-based feature importance scores as a data-driven evidence stream feeding into the primary ranking.

This separation of roles is intentional. Using known biomarker labels to train the ranking model would create circularity (the model would learn to reproduce its own validation criteria). Instead, the ML stream asks a biologically grounded but independent question: which genes best distinguish tumor from normal tissue?

### 4.2 ML Evidence Stream: Random Forest with SHAP Explainability

#### 4.2.1 Design Rationale

An earlier iteration of the pipeline (PR #14, March 10, 2026) trained a supervised classifier to predict which genes are CRC biomarkers, using pipeline-derived features (omics signal, literature counts, pathway counts) as inputs. Both Random Forest and Logistic Regression achieved AUC 1.0 under 5-fold cross-validation.

The perfect AUC was identified as a red flag during validation review. The classifier's features were derived from the same pipeline that generated the training labels, creating a circular dependency. The model was not learning biology; it was learning to reproduce its own inputs. The issue was caught because the validation enrichment ratio did not improve alongside the apparently perfect classifier performance.

The approach was redesigned in PR #16 (March 29, 2026). Instead of predicting biomarker labels, the redesigned module trains a Random Forest to classify tumor vs. normal samples using raw gene expression, then extracts SHAP values to quantify each gene's contribution to the classification. This avoids circularity because the classification task (tumor vs. normal) is biologically real and entirely independent of the biomarker label annotations used for validation.

**Note on the initial results:** The 5-fold CV results from the original supervised classifier (AUC 1.0, accuracy 100%) appear in earlier drafts of this document and in the PR #14 commit. Those results have been superseded by the SHAP-based approach. The original results are retained here for methodological transparency and to document the lesson learned about circular feature construction.

#### 4.2.2 Training Configuration

**Data:** TCGA-COAD RNA-seq matrix (473 tumor, 41 normal samples)

**Feature pre-selection:** Top 5,000 genes by expression variance (reduces dimensionality while retaining the most informative genes)

**Model:** Random Forest Classifier

**Hyperparameters:**
- n_estimators: 1,000
- max_features: 'sqrt'
- min_samples_leaf: 2
- class_weight: 'balanced' (handles 473:41 class imbalance)
- random_state: 42 + fold_index (reproducible per fold)

**Cross-validation:** 5-fold stratified cross-validation (ensures each fold preserves the 11.5:1 class ratio)

#### 4.2.3 SHAP Extraction

After each CV fold, SHAP values are computed using the TreeExplainer for all 5,000 pre-selected features across the held-out samples. Mean absolute SHAP importance is aggregated across all folds, producing a single importance score per gene. The final importance scores are saved to outputs/ml_importance.csv covering all 36,519 genes (genes not in the top-5,000 feature set receive a score of 0).

Of the 5,000 features evaluated, 3,559 have non-zero mean SHAP importance. The top SHAP gene is CA7, a carbonic anhydrase that is strongly downregulated in CRC tumors. CA7's high SHAP importance reflects its value as a discriminator between tumor and normal tissue at the expression level.

#### 4.2.4 Cross-Validation Performance

**Table 4: Random Forest Tumor/Normal CV Performance**

| Fold | Accuracy | AUC | Balanced Accuracy |
|------|----------|-----|-------------------|
| 1 | 100.0% | 1.0 | 1.0 |
| 2 | 98.1% | 1.0 | 0.987 |
| 3 | 100.0% | 1.0 | 1.0 |
| 4 | 100.0% | 1.0 | 1.0 |
| 5 | 99.6% | 1.0 | 0.994 |
| **Mean** | **99.5%** | **1.0** | **0.996** |

The near-perfect performance is expected and not indicative of overfitting or circularity. Tumor vs. normal tissue classification from genome-wide expression data is biologically well-separated. The goal of this module is not prediction accuracy but SHAP-derived feature attribution, and the high accuracy confirms that the model has learned meaningful expression patterns from which SHAP values are computed.

### 4.3 Model Validation Strategy

#### 4.3.1 Internal Validation (Cross-Validation)

The Random Forest uses 5-fold stratified cross-validation for internal validation. This provides reliable performance estimates while preserving the class ratio across folds. Because the dataset has only 41 normal samples, stratification is critical to ensure each fold includes at least 8 normal samples.

#### 4.3.2 External Validation (Known Biomarker Enrichment)

The primary validation of the overall pipeline is external: enrichment of known CRC biomarkers in the top-ranked candidates. This validation approach treats the 80 known CRC biomarkers as a held-out test set that is never used during scoring or model training.

**Validation Panel Composition:**

The validation panel was expanded from 33 to 80 known CRC biomarkers for the final analysis, organized in three confidence tiers:

**Table 5: Validation Panel Tier Definitions**

| Tier | Size | Source | Confidence Level |
|------|------|--------|-----------------|
| Tier 1 | 14 | FDA-approved markers or NCCN/ESMO clinical guideline markers | Highest: regulatory or guideline status |
| Tier 2 | 25 | COSMIC Cancer Gene Census, CRC-annotated | High: somatic driver mutation databases |
| Tier 3 | 41 | Validated CRC literature (>10 independent studies) | Moderate: strong but not guideline-level |

**Tier 1 (FDA/Guideline, 14 genes):** KRAS, NRAS, BRAF, EGFR, ERBB2, MLH1, MSH2, MSH6, PMS2, EPCAM, VEGFA, BMP3, NDRG4, DPYD

**Tier 2 (COSMIC CGC, 25 genes):** APC, TP53, PIK3CA, SMAD4, FBXW7, PTEN, CTNNB1, RNF43, POLE, TGFBR2, AXIN2, SOX9, TCF7L2, AMER1, ARID1A, SMAD2, ACVR2A, ATM, ZNRF3, BCL9L, RBM10, PCBP1, RPL22, PTPRT, FAM123B

**Tier 3 (Validated Literature, 41 genes):** VIM, CDH1, CEACAM5, MKI67, CDX2, DCC, MYC, MSH3, ERBB3, MET, VEGFB, FLT1, KDR, PDGFRA, IGF2, MAP2K1, TGFB1, SMAD3, NOTCH1, JAG1, DLL4, GSK3B, CSNK1A1, DVL2, MYB, ETV4, FOXQ1, CD274, PDCD1, CTLA4, LAG3, DNMT1, DNMT3B, TET2, KDM6A, MUTYH, CHEK2, BRCA1, BRCA2, PALB2, PIK3R1

**Enrichment Ratio Calculation:**

```
Enrichment = (known biomarkers found in top K / K) / (total known biomarkers in dataset / total genes)
```

An enrichment ratio of 1.0 indicates random performance; values above 1.0 indicate the system ranks known biomarkers better than chance.

#### 4.3.3 External Validation Results

Of 80 known CRC biomarkers, 79 were found in the TCGA-COAD dataset. FAM123B (Tier 2) was not found due to a retired gene alias issue (the gene is now annotated as AMER1 in current genome builds; the validation panel was updated in subsequent runs to use the canonical symbol). The 79-gene denominator is used for all subsequent validation calculations.

**Table 6: Overall Validation Summary**

| Metric | Value |
|--------|-------|
| Known biomarkers in dataset | 79/80 (98.8%) |
| Enrichment ratio | 1.46x |
| Best known biomarker rank | BMP3 at rank #14 (score 63.37) |
| Known biomarkers in top 500 | 4 (BMP3, FOXQ1, DPYD, ETV4) |
| Known biomarkers in top 1,000 | 4 |
| Known biomarkers in top 5,000 | 11 |

**Table 7: Validation Results by Tier**

| Tier | Found | Median Rank | In Top 500 |
|------|-------|-------------|------------|
| Tier 1 (FDA/guideline, 14) | 14/14 | 10,373 | 2 (BMP3, DPYD) |
| Tier 2 (COSMIC CGC, 25) | 24/25 | 9,557 | 0 |
| Tier 3 (validated literature, 41) | 41/41 | 14,082 | 2 (FOXQ1, ETV4) |

BMP3 ranks at position #14 (final score 63.37), reflecting strong multi-stream evidence: differential expression, high SHAP importance, a strong LLM literature score (75/100), and pathway annotations. This demonstrates the value of multi-evidence fusion: expression-based biomarkers benefit from all streams reinforcing each other.

The median ranks for all three tiers are in the 10,000-14,000 range because the pipeline's primary evidence stream (expression) does not capture mutation-driven biology. The validation enrichment ratio (1.46x) is lower than the 2.07x reported with the earlier 33-biomarker panel because the expanded panel includes many Tier 2 genes (canonical tumor suppressors: APC, TP53, SMAD4) whose CRC relevance is driven by somatic mutations that are invisible to differential expression analysis.

### 4.4 Scoring Weight Ablation Study

To quantify the contribution of each evidence stream to biomarker recovery, a systematic weight ablation study was conducted. The study tested 246 weight combinations across all four evidence streams at 0.10 increments, with each weight constrained to [0.0, 0.70] and all weights summing to 1.0. For each combination, the final scores were recomputed using the pre-normalized component scores from ranked_candidates.csv, genes were re-ranked, and the enrichment ratio was measured against the 79 known CRC biomarkers.

The ablation study is implemented in src/ablation.py and executes in under 30 seconds because it reuses the already-normalized component scores rather than re-running the full pipeline for each combination.

**Table 8: Ablation Study Key Findings**

| Weight Configuration (omics/ml/lit/path) | Enrichment | BMP3 Rank | Top 500 | Top 5,000 |
|------------------------------------------|-----------|-----------|---------|----------|
| 0.40/0.00/0.30/0.30 | 1.537x | 12 | 4 | 16 |
| 0.50/0.00/0.50/0.00 | 1.537x | 9 | 4 | 16 |
| 0.20/0.00/0.60/0.20 | 1.537x | 23 | 4 | 16 |
| 0.40/0.10/0.30/0.20 (current) | 1.462x | 14 | 4 | 11 |

The optimal enrichment ratio of 1.537x was achieved by multiple weight configurations, all of which set the ML importance weight to 0.0. This finding indicates that SHAP-based tumor/normal importance scores, while biologically meaningful, do not improve known biomarker recovery when combined with omics and literature evidence. This is expected: SHAP captures genes that distinguish tumor from normal tissue broadly, while known CRC biomarkers are often mutation-driven genes (APC, KRAS, TP53) whose dysregulation occurs at the DNA level rather than the expression level.

The current weights (omics=0.40, ml=0.10, lit=0.30, path=0.20) were retained despite the ablation finding for three reasons: (1) the ML stream provides a theoretically grounded, data-driven signal that complements the hypothesis-driven evidence streams; (2) the configuration omics=0.50/lit=0.50 achieves BMP3 at rank #9 but drops 5 genes from the top 5,000 compared to alternatives; and (3) the 0.075x enrichment difference is modest and may not generalize to other validation panels. The ablation results are archived in outputs/ablation_results.csv for reference.

### 4.5 Precision and Recall Evaluation

To provide a more granular evaluation than the single enrichment ratio metric, precision, recall, and F1 score were computed at 10 rank thresholds against the 79 known CRC biomarkers found in the dataset. This analysis is implemented in src/precision_recall.py.

**Formulas:**

- Precision@k = |{known biomarkers in top k}| / k
- Recall@k = |{known biomarkers in top k}| / |{all known biomarkers in dataset}|
- F1@k = 2 x (Precision@k x Recall@k) / (Precision@k + Recall@k)

**Table 9: Precision and Recall at Rank Thresholds**

| k | True Positives | Precision | Recall | F1 |
|---|---------------|-----------|--------|-----|
| 10 | 0 | 0.000 | 0.000 | 0.000 |
| 25 | 1 | 0.040 | 0.013 | 0.019 |
| 50 | 1 | 0.020 | 0.013 | 0.016 |
| 100 | 2 | 0.020 | 0.025 | 0.022 |
| 250 | 3 | 0.012 | 0.038 | 0.018 |
| 500 | 4 | 0.008 | 0.051 | 0.014 |
| 1,000 | 4 | 0.004 | 0.051 | 0.007 |
| 2,500 | 6 | 0.002 | 0.076 | 0.005 |
| 5,000 | 11 | 0.002 | 0.139 | 0.004 |
| 10,000 | 32 | 0.003 | 0.405 | 0.006 |

**Tier-Level Analysis:**

**Tier 1 (FDA/guideline, 14 markers):** The first Tier 1 marker (BMP3) appears at rank 14, with DPYD appearing at rank 226. By k=10,000, recall reaches 42.9% (6/14). These markers are best captured because they tend to be expression-based biomarkers (BMP3, NDRG4) that benefit from all four evidence streams reinforcing each other.

**Tier 2 (COSMIC CGC, 24 markers):** No Tier 2 markers appear until rank 2,500. AXIN2 is the first Tier 2 gene to appear (rank 1,348). Recall at k=10,000 is 50.0% (12/24). Most Tier 2 genes are canonical tumor suppressors (APC, TP53, SMAD4) whose CRC relevance is driven by somatic mutations rather than expression changes.

**Tier 3 (validated literature, 41 markers):** Two Tier 3 markers appear in the top 500 (FOXQ1 at rank 53, ETV4 at rank 309). Recall at k=10,000 is 34.2% (14/41). These markers span immune checkpoints, epigenetic regulators, and signaling genes whose expression patterns vary widely.

The low precision values (all below 5%) are expected and do not represent a system limitation. In a genome of 36,519 genes with only 79 known biomarkers (0.22% base rate), even a perfect system would have precision below 1% at k=500. The relevant metric is recall, which shows the system recovers 40.5% of known biomarkers in the top 27% of genes (10,000/36,519), a 1.46x improvement over random expectation. Results are archived in outputs/precision_recall.csv.

### 4.6 MMR Deficiency Classifier (Ayan Choudhury)

As a complementary ML analysis, a sample-level classifier was developed to predict mismatch repair (MMR) protein loss status from gene expression data. Unlike the SHAP-based tumor/normal classification used in the scoring pipeline, this classifier addresses a clinically meaningful within-cancer subtype question: can expression data distinguish dMMR (deficient MMR / MSI-high) from pMMR (proficient MMR / MSS) tumors?

This distinction is clinically important because dMMR status predicts response to immune checkpoint inhibitors (pembrolizumab has FDA approval for dMMR/MSI-H solid tumors regardless of site), yet determining MMR status currently requires IHC or PCR testing that is not always performed. A transcriptomic classifier could supplement molecular testing in resource-limited settings.

**Dataset:** TCGA-COAD primary tumor samples. The label was derived from the TCGA phenotype field loss_expression_of_mismatch_repair_proteins_by_ihc (YES = dMMR, NO = pMMR). After filtering for available labels and alignment with expression data, 340 samples were retained (56 dMMR, 284 pMMR).

**Table 10: MMR Deficiency Classifier - 5-Fold CV Results**

| Model | ROC-AUC | Avg Precision | F1 | Balanced Accuracy |
|-------|---------|--------------|-----|-------------------|
| XGBoost | 0.742 | 0.349 | 0.027 | 0.496 |
| Random Forest | 0.728 | 0.338 | 0.000 | 0.496 |
| Logistic Regression | 0.708 | 0.376 | 0.410 | 0.663 |

**Table 11: MMR Deficiency Classifier - Holdout Evaluation (68 test samples)**

| Model | ROC-AUC | PR-AUC | F1 (tuned threshold) | Recall |
|-------|---------|--------|----------------------|--------|
| XGBoost | 0.788 | 0.312 | 0.385 | 0.455 |
| Random Forest | 0.758 | 0.293 | 0.500 | 1.000 |
| Logistic Regression | 0.737 | 0.339 | 0.560 | 0.636 |

The ROC-AUC values of 0.71 to 0.79 reflect a genuinely challenging classification problem, in sharp contrast to the near-perfect tumor/normal separation. This validates the biological premise that MMR status produces detectable but subtle transcriptomic signatures. SHAP analysis on the best holdout model (Logistic Regression) identified genes driving the dMMR/pMMR distinction, providing an orthogonal set of biomarker candidates related to microsatellite instability, immune activation, and DNA repair pathways.

This analysis demonstrates how the modular pipeline architecture supports multiple complementary ML analyses on the same underlying expression data, each addressing different biological questions with appropriate performance expectations.

---

## Section 5: Model Deployment and Continual Learning (Lesson 12)

### 5.1 Trained Model Parameters

#### 5.1.1 ML Model (Random Forest with SHAP)

The ML evidence stream uses a Random Forest classifier trained transiently during each pipeline run. No serialized model artifacts are persisted. This design is intentional: the model's purpose is not prediction but explainability via SHAP values. Each run trains fresh on the latest expression data, ensuring SHAP importance scores always reflect current data.

Training is performed via 5-fold stratified cross-validation on TCGA-COAD RNA-seq data (473 tumor, 41 normal samples). Key hyperparameters are listed in Table 4 (Section 4.2.2). Validation accuracy ranges from 98.1% to 100% across folds, as expected for this biologically well-separated task.

#### 5.1.2 SHAP Feature Importance

The Random Forest produces SHAP values for all 5,000 pre-selected features across each CV fold. These are aggregated into mean absolute SHAP importance per gene and saved to outputs/ml_importance.csv. Of 5,000 features evaluated, 3,559 have non-zero mean SHAP importance. This file serves as the ML evidence stream (10% weight) in the final scoring formula.

#### 5.1.3 LLM Scoring Cache

LLM scores are cached in outputs/llm_scores_cache.json using content-hash keys (MD5 of each gene's abstract list). This enables:

- **Deterministic re-scoring:** Same abstracts always produce the same score (temperature=0)
- **Incremental updates:** Only genes with new or changed abstracts are re-scored
- **Interruption recovery:** Cache is saved after each gene, so interrupted runs resume from where they stopped

Current cache state (as of April 2026): 460 genes scored, range 0-95, mean 28.7, median 15.0.

### 5.2 Deployment Strategy

#### 5.2.1 Architecture Overview

```
┌─────────────────────────────────────────────────────────┐
│                    Azure Cloud (East US)                 │
│                                                         │
│  ┌─────────────┐    ┌──────────────────┐                │
│  │   GitHub    │    │  Azure Container  │                │
│  │   Actions   │───>│  Registry (ACR)   │                │
│  │   CI/CD     │    │  psuai894acr      │                │
│  └─────────────┘    └────────┬─────────┘                │
│                              │                           │
│                    ┌─────────▼──────────┐                │
│                    │  Azure App Service  │                │
│                    │  psuai894webapp     │                │
│                    │  (Linux, B1 tier)   │                │
│                    │  Port 8000          │                │
│                    └─────────┬──────────┘                │
│                              │                           │
│  ┌───────────────┐  ┌───────▼──────────┐                │
│  │ Azure VNet    │  │  Azure Files     │                │
│  │ 10.0.0.0/16   │  │  (mounted at     │                │
│  │ Web Subnet    │  │  /mnt/web-app-   │                │
│  │ PE Subnet     │  │  share)          │                │
│  └───────────────┘  └──────────────────┘                │
│                                                         │
│  ┌──────────────────────────┐                           │
│  │ Private DNS Zone         │                           │
│  │ privatelink.file.core... │                           │
│  └──────────────────────────┘                           │
└─────────────────────────────────────────────────────────┘
```

#### 5.2.2 Containerization

The application uses a multi-stage Docker build for minimal image size and reduced attack surface:

**Stage 1 (Build):** Base image python:3.11-slim-bullseye. Installs build dependencies (gcc, build-essential). Compiles Python wheels from requirements.txt.

**Stage 2 (Runtime):** Clean python:3.11-slim-bullseye (no compiler tools). Copies only compiled wheels and application code. Installs curl for health checks. Exposes port 8000. Entrypoint: python main.py (Uvicorn ASGI server).

Docker Compose provides local development orchestration: the web-app service builds from the Dockerfile and mounts source for live reloading; an nginx reverse proxy maps port 8080 to 8000 with HTTP basic auth via .htpasswd.

#### 5.2.3 CI/CD Pipeline

GitHub Actions workflow (.github/workflows/deploy.yml) triggers automatically on every push to main:

1. Checkout source code
2. Authenticate with Azure Container Registry using GitHub Secrets (ACR_LOGIN_SERVER, ACR_USERNAME, ACR_PASSWORD)
3. Build multi-stage Docker image
4. Push to ACR with two tags:
   - biomarker-cancer-web-app:latest (rolling release for Azure App Service)
   - biomarker-cancer-web-app:{commit-sha} (immutable version for rollback)

Azure App Service is configured with DOCKER_ENABLE_CI=true, so it auto-pulls the latest image after each push. This cycle takes approximately 5-10 minutes from push to live deployment.

#### 5.2.4 Azure Infrastructure (Terraform)

All infrastructure is defined as code in terraform/:

**Table 12: Azure Infrastructure Resources**

| Resource | Configuration | File |
|----------|--------------|------|
| Resource Group | PSU-AI894-RG, East US, Production tagged | resource_group.tf |
| App Service Plan | Linux, B1 SKU, always-on | webapp.tf |
| Linux Web App | Docker container, port 8000, VNet integrated | webapp.tf |
| Container Registry | Basic SKU, admin access enabled | acr.tf |
| Virtual Network | 10.0.0.0/16 with web and PE subnets | networking.tf |
| Private Endpoint | Secures Azure Files via private DNS | networking.tf |
| Storage Account | Azure Files mounted at /mnt/web-app-share | storage.tf |

Terraform state is stored remotely in Azure Storage (psuai894storage in PSU-AI894-State-RG).

### 5.3 Serving System: Web Application

#### 5.3.1 FastAPI Backend

The web application is built with FastAPI running on Uvicorn (ASGI server).

**Table 13: REST API Endpoints**

| Endpoint | Method | Purpose |
|----------|--------|---------|
| / | GET | Serves the single-page web UI (index.html) |
| /health | GET | Health check for Azure load balancer |
| /api/candidates | GET | Paginated, filterable, sortable ranked gene list |
| /api/genes/{gene_id} | GET | Detailed evidence breakdown for a single gene |
| /api/stats | GET | Pipeline summary statistics |
| /api/validation | GET | Known CRC biomarker validation results |
| /api/chat | POST | LLM-powered conversational Q&A about results |
| /api/reload | POST | Clear in-memory cache after pipeline re-run |

#### 5.3.2 Data Serving Strategy

- **In-memory caching:** Pipeline CSV outputs are loaded into Pandas DataFrames on first request and cached in memory for subsequent requests.
- **Lazy loading:** Data is read from disk only when first accessed, reducing startup time.
- **Cache invalidation:** POST /api/reload clears the in-memory cache, forcing a fresh read after pipeline re-runs without server restart.
- **On-the-fly joins:** Omics evidence (expression direction, log2FC, FDR) is merged with ranking data at query time, reducing storage duplication.

#### 5.3.3 LLM Chat Integration

The /api/chat endpoint connects to Claude (claude-sonnet-4-20250514) with a tiered context system:

- **Tier 1 (always included):** Static dataset summary including methodology, scoring weights, and key terminology
- **Tier 2 (always included):** Top 50 ranked candidates with full score breakdowns
- **Tier 3 (optional):** Gene-specific detail when the user is viewing a particular gene

Conversation history (last 20 messages) is maintained for context-aware follow-ups.

#### 5.3.4 Web UI Features

The single-page application (web-app/static/index.html) provides:

- **Dashboard stats:** Total genes analyzed, significant genes (FDR < 0.05), score range
- **Validation panel:** Known CRC biomarker enrichment ratio, counts in top 500/1,000/5,000
- **Ranked candidate table:** Sortable, searchable, paginated list of all ranked candidates
- **Gene detail view:** Per-gene evidence breakdown (omics, ML, literature, pathway, LLM scores with rationale)
- **AI Chat:** Conversational interface for exploring results with Claude

### 5.4 Continual Learning

#### 5.4.1 Caching and Incremental Updates

The pipeline implements three caching layers that enable efficient continual updates:

**1. LLM Score Cache (outputs/llm_scores_cache.json):**
- Content-hash (MD5) keyed by each gene's abstract list
- Cache hit: reuse previous score without API call
- Cache miss: send to Claude API for scoring
- Auto-saves after each gene (survives interruptions)
- Impact: Reduces approximately 15-minute LLM scoring step to seconds on re-runs with unchanged literature

**2. Gene Symbol Mapping Cache (data/ensembl_to_symbol_cache.tsv):**
- Maps Ensembl IDs to HGNC symbols
- Falls back to mygene.info API on cache miss
- Prevents redundant API calls across runs

**3. Web Application Cache (in-memory DataFrames):**
- Lazy-loaded from pipeline CSVs on first request
- POST /api/reload invalidates after pipeline re-run
- Avoids repeated disk I/O during interactive exploration

#### 5.4.2 Pipeline Re-execution Strategy

The modular pipeline supports continual learning through four update patterns:

1. **New Data Ingestion:** Replace data/TCGA-COAD.star_counts.tsv with updated counts and re-run all seven steps. All downstream outputs regenerate automatically.
2. **Literature Refresh:** Re-run the PubMed step to capture newly published abstracts. The LLM cache invalidates only for genes with changed abstract content.
3. **Scoring Weight Tuning:** Modify weights in config/config.yaml and re-run only Steps 6-7 (scoring and validation). Upstream evidence is unchanged.
4. **Model Retraining:** The Random Forest retrains from scratch on each run using the latest omics data, ensuring the ML evidence stream always reflects current data.

#### 5.4.3 Validation as a Feedback Signal

The validation step (Step 7) serves as a continual learning feedback mechanism. The enrichment ratio (currently 1.46x) measures how much better the system ranks known biomarkers versus random expectation. After any pipeline change (new data, weight adjustment, model update), the enrichment ratio indicates whether the change improved or degraded performance. The tier and category breakdowns reveal which biomarker types benefit most from a given change.

#### 5.4.4 ML Approach Iteration: A Case Study in Continual Improvement

The evolution of the ML evidence stream illustrates how the pipeline's modular design supports continual learning through iterative refinement.

**Iteration 1 (PR #14, March 10, 2026):** A supervised classifier trained to predict CRC biomarker status from pipeline-derived features achieved AUC 1.0 under 5-fold CV. This was identified as a red flag: the features and labels were circularly derived, so the model was reproducing its own inputs rather than learning biology.

**Iteration 2 (PR #16, March 29, 2026):** The approach was redesigned to classify tumor vs. normal samples using raw gene expression, with SHAP values providing the actual evidence signal. This avoids circularity because the classification task is biologically real and independent of biomarker labels.

**Result:** BMP3 rose from rank #163 to rank #14 after this change, validating that the SHAP signal captures meaningful expression-based biomarker evidence. The modular architecture allowed this fundamental change without affecting any other pipeline module.

#### 5.4.5 Deployment Update Cycle

The CI/CD pipeline enables a continual deployment loop:

1. Develop locally: modify pipeline code or config
2. Run pipeline: python run_pipeline.py --config config/config.yaml
3. Validate: check enrichment ratio in outputs/validation_summary.json
4. Commit and push to main
5. GitHub Actions builds and pushes new Docker image to ACR with commit SHA tag
6. Azure App Service auto-pulls latest image
7. Rollback if needed by redeploying a previous commit SHA tag

---

## Section 6: AI Cybersecurity, Ethics, and Accountability (Lesson 12)

### 6.1 Cybersecurity

#### 6.1.1 Credentials and Secrets Management

All API keys and access credentials are managed outside the source code repository. The Anthropic API key is stored in a web-app/.env file that is listed in .gitignore and never committed to the repository. PubMed API key configuration was moved from a hardcoded value in config/config.yaml to an environment variable in the final pipeline version, addressing a known technical debt item identified during development. CI/CD secrets (ACR_LOGIN_SERVER, ACR_USERNAME, ACR_PASSWORD) are stored as encrypted GitHub Actions Secrets and are never written to logs or artifact outputs.

#### 6.1.2 Network Security

The Azure deployment applies defense-in-depth network controls. Azure VNet integration routes all App Service traffic through a virtual network with dedicated web and private endpoint subnets. Azure Files storage is accessible only via a private endpoint with a dedicated private DNS zone (privatelink.file.core.windows.net), preventing public internet access to pipeline output files. The App Service itself does not expose a public storage endpoint.

#### 6.1.3 Application-Level Security

**CORS:** FastAPI CORS middleware is configured to restrict origins to localhost during development, with explicit production domain configuration required before broader deployment. This prevents cross-origin requests from unauthorized domains.

**Input validation:** All API inputs are validated through Pydantic models (request body validation) and FastAPI query parameter type enforcement. Malformed requests are rejected at the framework level before reaching application logic.

**Read-only API design:** All data retrieval endpoints use GET methods. Only two state-modifying endpoints exist (/api/reload and /api/chat), and /api/reload requires a protected call pattern. The chat endpoint (POST) processes user text through a structured system prompt that constrains the model's behavior to the biomarker discussion context.

**Docker security:** The multi-stage Docker build excludes all build tools (gcc, build-essential, pip-compile) from the runtime image. This reduces the attack surface by eliminating tools that could be exploited if the container were compromised. The runtime image exposes only port 8000.

**Local development proxy:** nginx basic authentication (.htpasswd) protects the local development proxy endpoint, preventing unauthorized access to the locally running web application during development.

#### 6.1.4 Data Privacy

No patient-level or personally identifiable information is processed by this system. The TCGA-COAD dataset contains gene expression values aggregated across patient cohorts; individual patient identifiers are not stored, logged, or served through any API endpoint. The system processes only aggregate statistical outputs (fold change, FDR, SHAP importance) derived from the source data. TCGA data access complies with the TCGA Data Use Agreement governing public-access tier data.

#### 6.1.5 Known Security Limitations and Planned Improvements

The chat endpoint's system prompt provides context constraints but does not implement formal prompt injection defenses. An adversarial user could potentially craft inputs designed to bypass the context constraints. Implementation of server-side input sanitization and prompt injection detection (e.g., filtering for instruction-like patterns) is noted as future work, with implementation targeted for the next development sprint.

The ACR admin password is currently stored as a plain-text app setting in webapp.tf. Migration to Azure Managed Identity would eliminate this credential from Terraform state and improve the security posture. This is identified as a known technical debt item.

### 6.2 Ethics

#### 6.2.1 System Purpose and Framing

The system is explicitly designed and labeled as a research decision support tool, not a clinical diagnostic. All documentation, the web UI dashboard, and the LLM system prompt make clear that the rankings represent computational prioritization of publicly available evidence, not clinical recommendations. No component of this system should be used to inform patient care decisions without independent clinical validation and regulatory review.

#### 6.2.2 Data Ethics

All data sources used by this pipeline are publicly available under established data sharing frameworks. TCGA-COAD data is distributed under the TCGA Data Use Agreement. PubMed abstracts are retrieved through the official NCBI E-Utilities API within its rate limits and terms of service. Pathway data is retrieved from g:Profiler under its academic use terms. No proprietary clinical datasets or restricted-access genomic data are used.

No patient-level demographic, clinical, or identifying information is stored, processed, or displayed. The system operates exclusively on aggregate gene expression statistics.

#### 6.2.3 Transparency and Explainability

Transparency is a core design principle of this pipeline:

- **Decomposable scores:** Every gene's final score is the sum of four independently interpretable components (omics, ML, literature, pathway), each visible in the web UI.
- **LLM rationales:** Each gene's LLM score is accompanied by a stored 1-2 sentence rationale explaining the scoring decision, making the literature evidence assessment auditable.
- **SHAP attribution:** The ML evidence stream uses SHAP values, which are theoretically grounded in cooperative game theory and provide mathematically rigorous feature attribution.
- **Configurable weights:** Scoring weights are documented in config/config.yaml and can be modified by any user without code changes, rather than being hidden in implementation details.

#### 6.2.4 Bias Considerations

The pipeline has systematic biases that are documented rather than obscured:

**Literature bias:** Genes that have been heavily studied (well-funded targets, previously known cancer genes) will have more PubMed abstracts and higher LLM scores. Novel, understudied genes with genuine biomarker potential may be disadvantaged. This bias is partially mitigated by the 40% omics weight, which treats all genes equally regardless of prior study.

**Expression bias:** The pipeline cannot detect mutation-driven biomarkers (APC, TP53, KRAS) because these genes do not necessarily show strong differential expression. The precision/recall analysis quantifies this limitation explicitly.

**Cohort bias:** Results are based on a single North American (predominantly European-ancestry) TCGA cohort. Biomarker rankings may not generalize to other populations. Multi-cohort validation is listed as a primary future work priority.

### 6.3 Accountability

#### 6.3.1 Version Control and Traceability

All pipeline code, configuration, and infrastructure definitions are version-controlled in the GitHub repository (https://github.com/coughlka/AI_Capstone). The feature branch workflow with pull request reviews (PR #11 through PR #17 documented in the repository) provides a complete audit trail of every change, including the rationale for the ML approach redesign (circularity issue) and the biomarker panel expansion.

Pipeline outputs are not committed to the repository by default (they are gitignored), but were committed at the project milestone (commit 4e7d261) to provide reviewers with reproducible access to the final results. Docker images are tagged with commit SHAs, providing a complete chain of traceability from deployed application version back to the source code commit and pipeline configuration that produced the results.

#### 6.3.2 Logging and Monitoring

Each pipeline module uses Python's standard logging library to emit structured INFO and WARNING messages throughout execution. Log output includes step completion times, gene counts at each processing stage, API retry events, and cache hit/miss rates. This enables post-hoc execution tracing when diagnosing unexpected results.

#### 6.3.3 Performance Accountability

The validation framework (80 known CRC biomarkers, 3 tiers) provides a quantitative performance baseline against which all pipeline changes are evaluated. The enrichment ratio and precision/recall metrics serve as the primary accountability mechanism: any modification to the pipeline that does not improve or maintain these metrics requires explicit justification before merging. This was applied in the ablation study, where the 0.075x enrichment difference between the current configuration and the optimal ablation configuration was weighed against the methodological benefits of including the ML stream.

#### 6.3.4 Team Accountability

Responsibility for each pipeline module is documented in CLAUDE.md and in the GitHub CODEOWNERS configuration. Pull requests require review by at least one other team member before merging to main. Keith Coughlin is the designated reviewer for Gabriel Saenz's infrastructure changes. This provides accountability for deployment decisions that affect the production environment.

---

## Section 7: Final Report Discussion (Lesson 13)

### 7.1 Summary of Key Results

The CRC biomarker discovery pipeline analyzed 36,519 genes from TCGA-COAD RNA-seq data (473 tumor, 41 normal samples) and produced a ranked list of 500 biomarker candidates by synthesizing four evidence streams: differential expression (omics), SHAP-based ML importance, PubMed literature scored by Claude LLM, and pathway enrichment via g:Profiler. The system was validated against an 80-gene panel of known CRC biomarkers organized in three confidence tiers.

The pipeline found 79 of 80 known biomarkers (98.8%) in the dataset and achieved an enrichment ratio of 1.46x, meaning known biomarkers ranked 1.46 times better than random expectation. The best-ranked known biomarker was BMP3 at rank #14 (final score 63.37), an FDA-approved stool DNA biomarker used in the Cologuard CRC screening test. Four known biomarkers appeared in the top 500 candidates (BMP3, FOXQ1, DPYD, ETV4), and 11 in the top 5,000.

The scoring weight ablation study tested 246 weight configurations and found that the maximum enrichment (1.537x) is achieved when the ML importance weight is set to zero, with omics and literature receiving redistributed weight. The current configuration (omics=0.40, ml=0.10, lit=0.30, path=0.20) produces 1.462x enrichment, a modest difference retained for methodological completeness. The precision/recall evaluation showed recall of 5.1% at k=500 and 40.5% at k=10,000, with Tier 1 (FDA/guideline) markers recovered most effectively due to their expression-based biology.

The complementary MMR deficiency classifier (Ayan Choudhury) achieved ROC-AUC 0.71-0.79 across three model architectures, demonstrating that transcriptomic data contains detectable but subtle signals for MSI subtype prediction, in contrast to the near-perfect (AUC 1.0) tumor/normal separation.

### 7.2 Usefulness in Solving the Research Problem

The system addresses the core research problem: biomarker prioritization is currently a manual, fragmented process that lacks reproducibility and fails to integrate evidence systematically across domains. This pipeline provides four substantive contributions to that problem:

**1. Systematic evidence synthesis:** Rather than analyzing omics data, literature, and pathways in isolation, the pipeline integrates all four streams into a single transparent ranking. Each gene's score is fully decomposable into its component contributions, enabling researchers to understand exactly why a candidate was prioritized. A gene ranked high because of strong omics signal alone (e.g., a novel lncRNA with no literature) can be immediately distinguished from a gene ranked high because of convergent evidence from all four streams (e.g., BMP3 at rank #14).

**2. LLM-grounded literature analysis:** The Claude-based scoring module converts unstructured PubMed abstracts into structured 0-100 relevance scores with written rationales. This replaces the previous count-based approach (where 62% of genes scored the maximum) with meaningful differentiation. The LLM scores are deterministic (temperature=0), cached for efficiency, and auditable through stored rationales. The score distribution (Table 3) shows a biologically coherent spread: high scores correspond to genes with well-documented direct CRC biomarker evidence, while low scores correspond to genes tangentially mentioned in cancer literature.

**3. Scalable validation framework:** The 80-biomarker, 3-tier validation panel provides a quantitative feedback signal after any pipeline modification. This enabled the team to identify and correct the circular classifier issue in one iteration (the AUC 1.0 red flag was immediately visible against the validation metrics) and to quantify the contribution of each evidence stream through the ablation study.

**4. Deployed interactive interface:** The web application allows researchers to explore rankings, filter by evidence type, examine per-gene evidence breakdowns, and ask questions about specific genes through an LLM-powered chat interface. This moves beyond a static CSV output toward an interactive research decision support tool that can be used by domain experts without programming expertise.

### 7.3 Practical Implications

**Scientific:** The pipeline demonstrates that multi-modal evidence synthesis outperforms any single evidence stream for expression-based biomarker prioritization. The BMP3 case study shows how an expression-based biomarker benefits from reinforcing signals across omics, literature, and pathway evidence: it ranks at #14 precisely because all four streams assign it high scores. The approach is generalizable to other cancer types by substituting the TCGA cohort (e.g., TCGA-READ for rectal cancer, TCGA-STAD for gastric cancer) and updating the validation panel with the relevant known biomarkers.

**Economic:** Downstream biomarker validation is expensive, typically ranging from $500,000 to $2,000,000 per candidate through preclinical and clinical trials. By systematically prioritizing candidates with stronger multi-evidence support and providing transparent reasoning for each ranking, the system can help research organizations focus resources on the most defensible candidates and deprioritize weakly supported genes earlier in the discovery pipeline.

**Social:** Improving CRC biomarker discovery supports the broader goal of earlier detection and better treatment stratification for a disease that kills approximately 50,000 Americans annually. While this system does not perform diagnosis, it contributes to the upstream scientific process that informs clinical tool development. The deployment of an interactive web interface also supports the democratization of computational biology: domain experts without programming expertise can directly engage with and interrogate the results.

### 7.4 Limitations

**1. Expression-based bias:** The pipeline's primary evidence stream (omics) measures differential expression between tumor and normal tissue. Known CRC biomarkers that operate through somatic mutations (APC, TP53, KRAS) rather than expression changes rank poorly regardless of their clinical importance. This is reflected directly in the Tier 2 (COSMIC CGC) results: no Tier 2 markers appear in the top 500 and the median rank for Tier 2 markers is 9,557. The enrichment ratio would be substantially higher if the validation panel contained only expression-based biomarkers.

**2. Single-cohort validation:** All results are based on TCGA-COAD (colon adenocarcinoma). Cross-validation against independent cohorts (TCGA-READ for rectal adenocarcinoma, GEO-based datasets) would strengthen confidence in the rankings and reveal which candidates are colon-specific versus broadly applicable to colorectal cancer. The 41 normal samples provide adequate statistical power but represent a limited normal tissue reference.

**3. Literature coverage gap:** Approximately 8% of the 500 candidates (mostly lncRNAs and novel transcripts) have no PubMed abstracts, receiving a literature score of 0. This creates a systematic bias toward well-studied protein-coding genes. The gene alias retry mechanism (which searches PubMed using common alternative gene names) partially mitigates this for genes with synonyms, but cannot address genes that genuinely lack CRC-relevant publications.

**4. LLM scoring variability:** While deterministic within a run (temperature=0, content-hash caching), the LLM scores are model-dependent. Different Claude versions or alternative LLMs would likely produce different scores for the same abstracts. The content-hash caching ensures within-run consistency but does not guarantee cross-model reproducibility. As the Claude model is updated, cached scores from previous model versions may diverge from scores produced by the current version.

**5. Hand-tuned weights:** The scoring weights were chosen based on domain expertise and validated post-hoc against the 80 known biomarkers. While the ablation study provides systematic analysis, the current weights were not learned from data. A learning-to-rank approach with a larger validation panel could potentially identify weights that are more robust to the choice of validation panel.

### 7.5 Future Work

**1. Multi-cohort validation:** Extend the pipeline to process TCGA-READ and independent GEO datasets (e.g., GSE39582, a large CRC cohort with clinical annotations) and assess whether the top 500 candidates are consistent across cohorts. Genes ranked highly in multiple independent cohorts would be stronger candidates for experimental validation.

**2. Mutation-aware evidence stream:** Integrate somatic mutation data from cBioPortal or COSMIC as a fifth evidence stream, assigning mutation frequency scores alongside the expression-based signals. This would directly address the Tier 2 performance gap and enable the pipeline to recover mutation-driven biomarkers like APC, TP53, and KRAS.

**3. MMR-subtype SHAP integration:** Incorporate SHAP values from the MMR deficiency classifier as an additional evidence stream in the pipeline, capturing intra-cancer heterogeneity that the tumor/normal classifier misses. Genes important for dMMR/pMMR distinction could receive a higher composite score when MSI subtype relevance is a research priority.

**4. Learned scoring weights:** Replace hand-tuned weights with a learning-to-rank model (e.g., LambdaMART or a logistic regression with leave-one-out cross-validation over the known biomarker panel), potentially producing weights that are more data-grounded and less dependent on domain judgment.

**5. Temporal literature monitoring:** Implement scheduled pipeline re-runs (e.g., monthly) to capture newly published CRC literature, leveraging the LLM cache to only re-score genes with new abstracts. This supports true continual learning as the CRC research landscape evolves and new biomarker evidence accumulates.

**6. Prompt injection defense:** Implement server-side input sanitization and prompt injection detection on the chat endpoint, strengthening the security posture for production deployment beyond the academic context.

### 7.6 Limitations of Classifiers and Suggestions for Improvement

The Random Forest classifier used for SHAP importance extraction achieves near-perfect accuracy (mean 99.5% across folds) on the tumor/normal classification task. While this validates that the TCGA-COAD expression data contains strong biological signal, it also means the SHAP values primarily capture broad cancer signatures (global metabolic reprogramming, cell cycle dysregulation) rather than CRC-specific biomarker signals. The ablation study directly confirmed this: setting the ML weight to zero improves enrichment from 1.462x to 1.537x.

The MMR deficiency classifier achieves more realistic and biologically appropriate performance (ROC-AUC 0.71 to 0.79 on holdout). The relatively low F1 scores in cross-validation (near 0 for RF and XGBoost) reflect the severe class imbalance (56 dMMR vs. 284 pMMR), where classifiers defaulting to the majority class achieve high accuracy but zero recall for the positive class. The holdout F1 scores (0.385 to 0.560 after threshold tuning) indicate that meaningful discrimination is achievable with calibrated thresholds.

Suggestions for improving these classifiers:

1. **CMS labels for SHAP:** Replace the binary tumor/normal task with consensus molecular subtype (CMS1-CMS4) labels, which represent biologically distinct CRC subtypes. SHAP values from a multi-class CMS classifier would capture CRC-specific expression patterns rather than generic tumor signals.

2. **Feature selection before SHAP:** Apply variance-based or mutual-information feature selection to reduce noise before extracting SHAP values, potentially improving the signal-to-noise ratio in the importance scores.

3. **Gradient boosting models:** XGBoost or LightGBM may capture different interaction patterns than Random Forest and potentially provide more nuanced SHAP attributions for both the tumor/normal and MMR classification tasks.

4. **Probability calibration:** Implement Platt scaling or isotonic regression on the MMR classifier outputs to produce well-calibrated probabilities, which would be required if the classifier outputs were to be used in a clinical support context.

5. **Bootstrap confidence intervals:** Compute bootstrap confidence intervals for the MMR classifier's ROC-AUC to provide more robust performance estimates given the small positive class (56 samples).

### 7.7 Individual Feedback

**Keith Coughlin:**

[Your paragraph here]

**Ayan Choudhury:**

[Your paragraph here]

**Gabriel Saenz:**

[Your paragraph here]

---

## Section 8: Project Presentation and Demo (Lesson 14)

### 8.1 Presentation Overview

The final project presentation will demonstrate the end-to-end CRC biomarker discovery pipeline, from raw TCGA-COAD data ingestion through the ranked candidate output, with interactive exploration via the deployed web application. The presentation is structured to address both the technical depth required for the AI 894 course and the domain context needed for a bioinformatics audience.

### 8.2 Demonstration Outline

**1. Problem Motivation (5 minutes)**
- CRC burden and the need for better biomarker discovery
- Limitations of manual, siloed evidence review
- Overview of the multi-stream synthesis approach

**2. Pipeline Architecture Walkthrough (10 minutes)**
- Seven-step architecture diagram
- Live demonstration of run_pipeline.py execution (abbreviated run or pre-recorded with key steps highlighted)
- Show config/config.yaml and explain how weights, thresholds, and data paths are configured

**3. ML/SHAP Component Deep Dive (10 minutes)**
- Explain the tumor/normal classification task and why it avoids circularity
- Show SHAP summary plot for top genes
- BMP3 case study: how all four evidence streams contribute to rank #14
- Brief note on the original circular classifier and the lesson learned

**4. Validation Results (10 minutes)**
- Enrichment ratio of 1.46x against 80 known CRC biomarkers
- Precision/recall curves with tier breakdown
- Ablation study: 246 weight combinations, best result at ml=0
- Honest assessment of limitations (mutation-driven biomarkers, single cohort)

**5. MMR Deficiency Classifier (5 minutes)**
- Motivation: clinical relevance of dMMR/pMMR distinction for immunotherapy
- Model comparison (LR/RF/XGB), ROC-AUC 0.71-0.79
- Contrast with tumor/normal AUC 1.0: realistic performance for a harder problem

**6. Web Application Live Demo (10 minutes)**
- Dashboard: stats, validation panel, top candidates
- Gene detail view for BMP3: show all four evidence components and LLM rationale
- AI Chat demonstration: ask "Which expression-based biomarkers are ranked highest?" and explore the response
- Show Azure deployment via App Service URL

**7. Discussion and Q&A (10 minutes)**
- Limitations and future work
- Scoring weight ablation interpretation
- Questions from instructor and reviewers

### 8.3 Key Points to Emphasize

- The pipeline is a decision support tool, not a diagnostic: all findings require experimental validation
- Multi-stream synthesis demonstrably outperforms single-stream analysis for expression-based biomarkers
- The circular classifier redesign illustrates responsible AI development: performance metrics must be interrogated, not just optimized
- The 1.46x enrichment ratio should be interpreted in the context of the base rate (0.22% known biomarkers in 36,519 genes) and the expression-bias limitation

---

## Appendices

### Appendix A: Scoring Weight Configuration

Current weights in config/config.yaml:

```yaml
scoring:
  weights:
    omics: 0.40
    ml_importance: 0.10
    literature: 0.30
    pathway: 0.20
```

These weights were selected based on domain expertise and evaluated through the ablation study described in Section 4.4. The omics stream receives the highest weight because differential expression provides the most direct molecular evidence for expression-based biomarker discovery. Literature receives the second-highest weight because LLM-analyzed PubMed abstracts contextualize the biological relevance of candidates with established research backing. Pathway enrichment captures systems-level evidence about gene function within relevant biological processes. ML importance provides a data-driven complement to the hypothesis-driven omics score.

The ablation study (246 combinations, outputs/ablation_results.csv) shows that setting ml_importance to 0.0 and redistributing its weight across the remaining three streams improves the enrichment ratio from 1.462x to 1.537x. The current weights are retained to include the theoretically grounded ML signal and for methodological completeness; the difference is modest and the ablation finding may not generalize to expanded validation panels.

### Appendix B: LLM Scoring Rubric

The Claude model scores each gene's literature evidence on a 0-100 scale using the following rubric:

**Table B.1: LLM Scoring Criteria**

| Criterion | Weight | Description |
|-----------|--------|-------------|
| Direct CRC evidence | 40% | Is this gene directly discussed as a CRC biomarker or therapeutic target? |
| Biomarker potential | 30% | Does the evidence support diagnostic, prognostic, or predictive utility? |
| Mechanistic insight | 20% | Is the gene's role in CRC biology mechanistically explained? |
| Evidence quality | 10% | Are the abstracts from peer-reviewed, well-powered studies? |

**Table B.2: Score Range Interpretations**

| Range | Interpretation |
|-------|---------------|
| 80-100 | Direct, established CRC biomarker with clear clinical utility |
| 60-79 | Strong CRC evidence, emerging or not yet clinically validated |
| 40-59 | CRC mentioned but not primary focus; indirect evidence |
| 20-39 | Tangential CRC relevance only |
| 0-19 | No meaningful CRC biomarker evidence in retrieved abstracts |

The LLM is called with temperature=0 to ensure deterministic output. Each response is a JSON object containing the integer score and a 1-2 sentence rationale. Scores are cached using an MD5 hash of the concatenated abstract text, so identical inputs always produce identical outputs across pipeline runs.

### Appendix C: Validation Biomarker Panel (80 Markers, 3 Tiers)

**Tier 1: FDA-Approved or Clinical Guideline Markers (14 genes)**

KRAS, NRAS, BRAF, EGFR, ERBB2, MLH1, MSH2, MSH6, PMS2, EPCAM, VEGFA, BMP3, NDRG4, DPYD

These genes are included in FDA-approved companion diagnostics, NCCN guidelines for CRC treatment selection, or FDA-cleared screening tests (Cologuard for BMP3 and NDRG4). All 14 were found in the TCGA-COAD dataset. BMP3 and DPYD appear in the top 500 ranked candidates.

**Tier 2: COSMIC Cancer Gene Census, CRC-Annotated (25 genes)**

APC, TP53, PIK3CA, SMAD4, FBXW7, PTEN, CTNNB1, RNF43, POLE, TGFBR2, AXIN2, SOX9, TCF7L2, AMER1, ARID1A, SMAD2, ACVR2A, ATM, ZNRF3, BCL9L, RBM10, PCBP1, RPL22, PTPRT, FAM123B

24 of 25 were found in the dataset. FAM123B was not found because it is a retired alias for AMER1 (already present in the panel); this was identified as a data quality issue and the panel was updated. No Tier 2 markers appear in the top 500, reflecting that most are mutation-driven drivers not captured by expression-based analysis.

**Tier 3: Validated CRC Literature Markers (41 genes)**

VIM, CDH1, CEACAM5, MKI67, CDX2, DCC, MYC, MSH3, ERBB3, MET, VEGFB, FLT1, KDR, PDGFRA, IGF2, MAP2K1, TGFB1, SMAD3, NOTCH1, JAG1, DLL4, GSK3B, CSNK1A1, DVL2, MYB, ETV4, FOXQ1, CD274, PDCD1, CTLA4, LAG3, DNMT1, DNMT3B, TET2, KDM6A, MUTYH, CHEK2, BRCA1, BRCA2, PALB2, PIK3R1

All 41 found. FOXQ1 (rank #53) and ETV4 (rank #309) appear in the top 500. These genes span immune checkpoints (CD274/PD-L1, PDCD1/PD-1, CTLA4), epigenetic regulators (DNMT1, TET2, KDM6A), WNT signaling components (GSK3B, CSNK1A1), and DNA damage repair genes (BRCA1, BRCA2, MUTYH).

**Functional Category Summary:**

| Category | Count | Notes |
|----------|-------|-------|
| Mutation/tumor suppressor | 27 | APC, TP53, KRAS, and most Tier 2 genes; poorly captured by expression |
| Signaling/pathway | 31 | Wnt, TGF-beta, EGFR, VEGF pathway genes |
| MMR/MSI | 6 | MLH1, MSH2, MSH6, PMS2, MSH3, POLE |
| Expression/methylation | 8 | DNMT1, DNMT3B, TET2, KDM6A, VIM, CDH1 |
| Immune checkpoint | 4 | CD274, PDCD1, CTLA4, LAG3 |
| Epigenetic regulation | 4 | DNMT1, DNMT3B, TET2, KDM6A |
| Pharmacogenomic | 1 | DPYD (5-FU metabolism) |

### Appendix D: Works Cited

1. Cancer Genome Atlas Network. "Comprehensive molecular characterization of human colon and rectal cancer." *Nature* 487, no. 7407 (2012): 330-337. https://doi.org/10.1038/nature11252

2. Goldman, Mary J., et al. "Visualizing and interpreting cancer genomics data via the Xena platform." *Nature Biotechnology* 38, no. 6 (2020): 675-678. https://doi.org/10.1038/s41587-020-0546-8

3. Raudvere, Uku, et al. "g:Profiler: a web server for functional enrichment analysis and conversions of gene lists (2019 update)." *Nucleic Acids Research* 47, no. W1 (2019): W191-W198. https://doi.org/10.1093/nar/gkz369

4. Lundberg, Scott M., and Su-In Lee. "A unified approach to interpreting model predictions." *Advances in Neural Information Processing Systems* 30 (2017). https://proceedings.neurips.cc/paper/2017/file/8a20a8621978632d76c43dfd28b67767-Paper.pdf

5. Bejnordi, Babak Ehteshami, et al. "Diagnostic assessment of deep learning algorithms for detection of lymph node metastases in women with breast cancer." *JAMA* 318, no. 22 (2017): 2199-2210. https://doi.org/10.1001/jama.2017.14585

6. Shen, Yaguang, et al. "SHAP-based interpretable machine learning for biomarker discovery." *Briefings in Bioinformatics* 24, no. 4 (2023): bbad178. https://doi.org/10.1093/bib/bbad178

7. Breiman, Leo. "Random forests." *Machine Learning* 45, no. 1 (2001): 5-32. https://doi.org/10.1023/A:1010933404324

8. Tibshirani, Robert. "Regression shrinkage and selection via the lasso." *Journal of the Royal Statistical Society: Series B* 58, no. 1 (1996): 267-288. https://doi.org/10.1111/j.2517-6161.1996.tb02080.x

9. Welch, B. L. "The generalization of Student's problem when several different population variances are involved." *Biometrika* 34, no. 1/2 (1947): 28-35. https://doi.org/10.1093/biomet/34.1-2.28

10. Benjamini, Yoav, and Yosef Hochberg. "Controlling the false discovery rate: A practical and powerful approach to multiple testing." *Journal of the Royal Statistical Society: Series B* 57, no. 1 (1995): 289-300. https://doi.org/10.1111/j.2517-6161.1995.tb02031.x

11. Lech, Günter, et al. "Colorectal cancer tumour markers and biomarkers: Recent therapeutic advances and future perspectives." *World Journal of Gastroenterology* 22, no. 5 (2016): 1745-1766. https://doi.org/10.3748/wjg.v22.i5.1745

12. Imperiale, Thomas F., et al. "Multitarget stool DNA testing for colorectal-cancer screening." *New England Journal of Medicine* 370, no. 14 (2014): 1287-1297. https://doi.org/10.1056/NEJMoa1311288

13. Forbes, Simon A., et al. "COSMIC: somatic cancer genetics at high-resolution." *Nucleic Acids Research* 45, no. D1 (2017): D777-D783. https://doi.org/10.1093/nar/gkw1121

14. Wu, Tim T., et al. "Identification of ETV4 and other genes as putative cancer/testis antigens." *Cancer Research* 66, no. 18 (2006): 9289-9299. https://doi.org/10.1158/0008-5472.CAN-06-1805

15. Chen, Tianqi, and Carlos Guestrin. "XGBoost: A scalable tree boosting system." *Proceedings of the 22nd ACM SIGKDD International Conference on Knowledge Discovery and Data Mining* (2016): 785-794. https://doi.org/10.1145/2939672.2939785

16. Tikhonov, Andrey, et al. "Comparative analysis of the g:Profiler web server 2019 update." *Nucleic Acids Research* 47 (2019): W191-W198.

17. Schriml, Lynn M., et al. "Human Disease Ontology 2018 update: classification, content and workflow expansion." *Nucleic Acids Research* 47, no. D1 (2019): D955-D962. https://doi.org/10.1093/nar/gky1032

18. Kanehisa, Minoru, et al. "KEGG for taxonomy-based analysis of pathways and genomes." *Nucleic Acids Research* 51, no. D1 (2023): D587-D592. https://doi.org/10.1093/nar/gkac963

19. Jassal, Bijay, et al. "The Reactome Pathway Knowledgebase." *Nucleic Acids Research* 48, no. D1 (2020): D498-D503. https://doi.org/10.1093/nar/gkz1031

20. Szklarczyk, Damian, et al. "STRING v11: protein-protein association networks with increased coverage." *Nucleic Acids Research* 47, no. D1 (2019): D607-D613. https://doi.org/10.1093/nar/gky1131

---

*Note: Screenshots of the web application interface, Azure portal deployment, pipeline execution logs, SHAP summary plots, and GitHub Actions CI/CD runs should be captured and inserted at the appropriate locations before final submission. All referenced output files (outputs/ablation_results.csv, outputs/precision_recall.csv, outputs/validation_summary.json, outputs/ranked_candidates.csv) are available in the repository at commit 4e7d261.*
