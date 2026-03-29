# AI 894 Biweekly Project Writeup - Week 9
## Model Deployment and Continual Learning

**Team Members:** Keith Coughlin, Ayan Chakraborty, Gabriel Saenz
**Date:** March 29, 2026
**Repository:** https://github.com/coughlka/AI_Capstone

---

## 1. Purpose

This report details the model deployment strategy and continual learning mechanisms for our CRC (colorectal cancer) biomarker discovery pipeline. The system synthesizes four evidence streams — omics differential expression, PubMed literature, pathway enrichment, and LLM-based scoring — to produce a ranked list of 500 biomarker candidates, served through an interactive web application deployed on Azure.

---

## 2. Trained Model Parameters

### 2.1 ML Model (Random Forest + SHAP)

The ML evidence stream uses a Random Forest classifier trained transiently during each pipeline run — no serialized model artifacts are persisted. This is intentional: the model's purpose is not prediction but **explainability** via SHAP values. Each run trains fresh on the latest expression data, ensuring SHAP importance scores always reflect current data.

Training is performed via 5-fold stratified cross-validation on TCGA-COAD RNA-seq data (473 tumor, 41 normal samples). Key hyperparameters:

**Random Forest:**
- `n_estimators`: 1,000
- `max_features`: 'sqrt'
- `min_samples_leaf`: 2
- `class_weight`: 'balanced' (handles 473:41 class imbalance)
- `random_state`: 42 + fold_index (reproducible per fold)
- Feature pre-selection: top 5,000 genes by variance

**Validation accuracy:** 98.1–100% across folds (expected — tumor vs. normal tissue is biologically well-separated). The high accuracy is not the goal; the SHAP values that explain *which genes drive the classification* are the output that feeds into the scoring pipeline.

**Note on prior approach:** An earlier iteration used a supervised classifier (Random Forest + Logistic Regression) trained to predict biomarker labels from pipeline scores. Both achieved AUC 1.0 — a red flag indicating the features and labels were circularly derived. The current SHAP-based approach avoids this by asking a fundamentally different question: not "which genes are biomarkers?" but "which genes best distinguish tumor from normal tissue, and by how much?" This reframing produces defensible, non-circular importance scores grounded in Shapley value theory.

### 2.2 SHAP Feature Importance

The Random Forest model produces SHAP (SHapley Additive exPlanations) values for all 5,000 pre-selected features across each CV fold. These are aggregated into mean absolute SHAP importance per gene and saved to `outputs/ml_importance.csv`. This serves as the ML evidence stream (10% weight) in the final scoring formula.

### 2.3 LLM Scoring Cache

LLM scores are cached in `outputs/llm_scores_cache.json` using content-hash keys (MD5 of each gene's abstract list). This enables:
- **Deterministic re-scoring**: Same abstracts always produce the same score
- **Incremental updates**: Only genes with new/changed abstracts are re-scored
- **Interruption recovery**: Cache is saved after each gene, so interrupted runs resume from where they stopped

---

## 3. Deployment Strategy

### 3.1 Architecture Overview

```
┌─────────────────────────────────────────────────────────┐
│                    Azure Cloud (East US)                 │
│                                                         │
│  ┌─────────────┐    ┌──────────────────┐                │
│  │   GitHub     │    │  Azure Container  │               │
│  │   Actions    │───>│  Registry (ACR)   │               │
│  │   CI/CD      │    │  psuai894acr      │               │
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

### 3.2 Containerization

The application uses a **multi-stage Docker build** for minimal image size and reduced attack surface:

**Stage 1 — Build:**
- Base image: `python:3.11-slim-bullseye`
- Installs build dependencies (`gcc`, `build-essential`)
- Compiles Python wheels from `requirements.txt`

**Stage 2 — Runtime:**
- Clean `python:3.11-slim-bullseye` (no compiler tools)
- Copies only compiled wheels and application code
- Installs `curl` for health checks
- Exposes port 8000
- Entrypoint: `python main.py` (Uvicorn ASGI server)

**Docker Compose** (local development):
- `web-app` service: Builds from Dockerfile, mounts source for live reloading
- `nginx` reverse proxy: Port 8080 -> 8000, HTTP basic auth via `.htpasswd`

**Files:** `Dockerfile`, `docker-compose.yml`

### 3.3 CI/CD Pipeline

GitHub Actions workflow (`.github/workflows/deploy.yml`) triggers automatically on every push to `main`:

1. **Checkout** source code
2. **Authenticate** with Azure Container Registry using GitHub Secrets (`ACR_LOGIN_SERVER`, `ACR_USERNAME`, `ACR_PASSWORD`)
3. **Build** multi-stage Docker image
4. **Push** to ACR with two tags:
   - `biomarker-cancer-web-app:latest` — rolling release for Azure App Service
   - `biomarker-cancer-web-app:<commit-sha>` — immutable version for rollback

Azure App Service is configured with `DOCKER_ENABLE_CI=true`, so it auto-pulls the latest image after each push.

### 3.4 Azure Infrastructure (Terraform)

All infrastructure is defined as code in `terraform/`:

| Resource | Configuration | File |
|----------|--------------|------|
| **Resource Group** | `PSU-AI894-RG`, East US, Production tagged | `resource_group.tf` |
| **App Service Plan** | Linux, B1 SKU, always-on | `webapp.tf` |
| **Linux Web App** | Docker container, port 8000, VNet integrated | `webapp.tf` |
| **Container Registry** | Basic SKU, admin access enabled | `acr.tf` |
| **Virtual Network** | 10.0.0.0/16 with web and PE subnets | `networking.tf` |
| **Private Endpoint** | Secures Azure Files via private DNS | `networking.tf` |
| **Storage Account** | Azure Files mounted at `/mnt/web-app-share` | `storage.tf` |

Terraform state is stored remotely in Azure Storage (`psuai894storage` in `PSU-AI894-State-RG`).

### 3.5 Security Considerations

- **VNet Integration**: App Service routes all traffic through Azure Virtual Network
- **Private Endpoints**: Azure Files storage accessible only via private DNS zone, not public internet
- **Secrets Management**: API keys stored as GitHub Secrets (CI/CD) and Azure App Settings (runtime)
- **CORS Policy**: Restricted to localhost origins (development); configurable for production domains
- **Nginx Basic Auth**: `.htpasswd` protects local development proxy
- **Minimal Docker Image**: Multi-stage build excludes build tools from runtime

---

## 4. Serving System — Web Application

### 4.1 FastAPI Backend

The web application is built with **FastAPI** running on **Uvicorn** (ASGI server), providing:

| Endpoint | Method | Purpose |
|----------|--------|---------|
| `/` | GET | Serves the single-page web UI (`index.html`) |
| `/health` | GET | Health check for Azure load balancer |
| `/api/candidates` | GET | Paginated, filterable, sortable ranked gene list |
| `/api/genes/{gene_id}` | GET | Detailed evidence breakdown for a single gene |
| `/api/stats` | GET | Pipeline summary statistics |
| `/api/validation` | GET | Known CRC biomarker validation results |
| `/api/chat` | POST | LLM-powered conversational Q&A about results |
| `/api/reload` | POST | Clear in-memory cache after pipeline re-run |

**Files:** `web-app/main.py`, `web-app/api/biomarkers.py`, `web-app/api/chat.py`

### 4.2 Data Serving Strategy

- **In-memory caching**: Pipeline CSV outputs are loaded into Pandas DataFrames on first request and cached in memory
- **Lazy loading**: Data is only read from disk when first accessed, reducing startup time
- **Cache invalidation**: `POST /api/reload` clears the cache, forcing a fresh read after pipeline re-runs
- **On-the-fly joins**: Omics evidence (direction, log2FC, FDR) is merged with ranking data at query time

### 4.3 LLM Chat Integration

The `/api/chat` endpoint connects to **Claude** (claude-sonnet-4-20250514) with a tiered context system:

- **Tier 1 (always included)**: Static dataset summary — methodology, scoring weights, key terminology
- **Tier 2 (always included)**: Top 50 ranked candidates with full score breakdowns
- **Tier 3 (optional)**: Gene-specific detail when user is viewing a particular gene

Conversation history (last 20 messages) is maintained for context-aware follow-ups.

### 4.4 Web UI Features

The single-page application (`web-app/static/index.html`) provides:

- **Dashboard stats**: Total genes, significant genes (FDR < 0.05), score range
- **Validation panel**: Known CRC biomarker enrichment ratio, counts in top 500/1000/5000
- **Ranked candidate table**: Sortable, searchable, paginated list of all ranked genes
- **Gene detail view**: Per-gene evidence breakdown (omics, literature, pathway, LLM scores)
- **AI Chat**: Conversational interface for exploring results with Claude

*[Screenshots of the web UI should be captured from http://localhost:8000 and inserted here]*

---

## 5. Software Architecture

### 5.1 Pipeline Architecture

The pipeline (`run_pipeline.py`) executes 7 sequential steps, each implemented as an independent module:

```
Step 1: Omics Analysis          → outputs/omics_evidence.csv
        (src/omics.py)            36,519 genes, Welch's t-test, BH FDR correction

Step 2: ML Importance           → outputs/ml_importance.csv
        (src/ml_importance.py)    Random Forest + SHAP, 5-fold CV, top 5,000 features

Step 3: PubMed Retrieval        → outputs/lit_evidence.csv
        (src/pubmed.py)           NCBI E-Utilities, up to 5 abstracts/gene

Step 4: Pathway Enrichment      → outputs/pathway_evidence.csv
        (src/pathway.py)          g:Profiler: GO:BP, KEGG, Reactome, WikiPathways

Step 5: LLM Scoring             → outputs/llm_scores.csv + llm_scores_cache.json
        (src/llm_scoring.py)      Claude scores 0-100, content-hash caching

Step 6: Weighted Scoring        → outputs/ranked_candidates.csv
        (src/scoring.py)          Multi-evidence combination, min-max normalization

Step 7: Validation              → outputs/validation_report.csv + validation_summary.json
        (src/validation.py)       80 known CRC biomarkers (3 tiers), enrichment ratio
```

### 5.2 Scoring Formula

```
final_score = 0.40 × omics_score
            + 0.10 × ml_importance_score
            + 0.30 × literature_score
            + 0.20 × pathway_score
```

Each component is normalized to [0, 100] before weighting:
- **Omics**: `|log2FC| × -log10(FDR + 1e-300)`, min-max normalized
- **ML Importance**: SHAP values, rank-based percentile normalization (preserves sparsity)
- **Literature**: LLM score (0-100) from Claude analysis of PubMed abstracts
- **Pathway**: Pathway count from g:Profiler, min-max normalized

Weights are configurable in `config/config.yaml`.

### 5.3 Configuration

All pipeline parameters are centralized in `config/config.yaml`:
- Data paths and output directories
- Omics parameters (FDR threshold, candidate count)
- PubMed query templates and rate limits
- Pathway enrichment sources (GO:BP, KEGG, Reactome, WikiPathways)
- ML parameters (CV folds, feature count)
- Scoring weights
- LLM model selection and caching

---

## 6. Continual Learning

### 6.1 Caching and Incremental Updates

The pipeline implements three caching layers that enable efficient continual updates:

**1. LLM Score Cache** (`outputs/llm_scores_cache.json`):
- Content-hash (MD5) keyed by gene's abstract list
- Cache hit: reuse previous score immediately
- Cache miss: send to Claude API for re-scoring
- Auto-saves after each gene (survives interruptions)
- **Impact**: Reduces ~15-minute LLM scoring step to seconds on re-runs with unchanged literature
- **Current state** (as of March 29, 2026): 460 genes scored, range 0-95, mean 28.7, median 15.0

  | Score Range | Count | Interpretation |
  |-------------|-------|----------------|
  | 80-95 | 39 | Direct established CRC biomarkers |
  | 60-79 | 107 | Strong CRC evidence, emerging |
  | 40-59 | 16 | CRC mentioned, indirect evidence |
  | 20-39 | 53 | Tangential CRC relevance |
  | 1-19 | 19 | Minimal CRC evidence |
  | 0 | 226 | No abstracts or no CRC relevance |

**2. Gene Symbol Mapping Cache** (`data/ensembl_to_symbol_cache.tsv`):
- Maps Ensembl IDs to HGNC symbols
- Falls back to mygene.info API on cache miss
- Prevents redundant API calls across runs

**3. Web Application Cache** (in-memory):
- Lazy-loaded DataFrames from pipeline CSVs
- `POST /api/reload` invalidates after pipeline re-run
- Avoids repeated disk I/O during interactive exploration

### 6.2 Pipeline Re-execution Strategy

The modular pipeline design supports continual learning through:

1. **New Data Ingestion**: Replace `data/TCGA-COAD.star_counts.tsv` with updated counts and re-run — all downstream outputs regenerate automatically
2. **Literature Refresh**: Re-run PubMed step to capture newly published abstracts; LLM cache invalidates only for genes with changed abstract content
3. **Scoring Weight Tuning**: Modify weights in `config/config.yaml` and re-run only Steps 6-7 (scoring + validation) — upstream evidence is unchanged
4. **Model Retraining**: The Random Forest and Logistic Regression models retrain from scratch on each run using the latest omics data, ensuring the ML evidence stream reflects current data

### 6.3 Validation as a Feedback Signal

The validation step (Step 7) serves as a continual learning feedback mechanism:

- **80 known CRC biomarkers** curated across 3 confidence tiers and 7 functional categories
- **3-tier hierarchy**: Tier 1 (14 FDA/guideline markers), Tier 2 (25 COSMIC Cancer Gene Census), Tier 3 (41 validated literature)
- **Enrichment ratio** (currently 1.46x) measures how much better the system ranks known biomarkers vs. random
- After any pipeline change (new data, weight adjustment, model update), the enrichment ratio indicates whether the change improved or degraded performance
- **Tier and category breakdowns** reveal which biomarker types benefit most from changes

Current validation results:

| Metric | Value |
|--------|-------|
| Known biomarkers found | 79/80 (98.8%) |
| Enrichment ratio | 1.46x |
| Best known biomarker | BMP3 (rank #14, score 63.37) |
| In top 500 | 4 (BMP3, FOXQ1, DPYD, ETV4) |
| In top 1,000 | 4 |
| In top 5,000 | 11 |

| Tier | Found | Median Rank | In Top 500 |
|------|-------|-------------|------------|
| Tier 1 (FDA/guideline) | 14/14 | 10,373 | 2 |
| Tier 2 (COSMIC CGC) | 24/25 | 9,557 | 0 |
| Tier 3 (validated literature) | 41/41 | 14,082 | 2 |

BMP3 — an FDA-approved stool DNA biomarker used in Cologuard screening — ranks highly because it exhibits strong differential expression *and* high SHAP importance *and* a high LLM literature score (75/100). This demonstrates the value of multi-evidence fusion: expression-based biomarkers benefit from all streams reinforcing each other. The enrichment ratio is lower than the earlier 33-biomarker validation (1.84x) because the expanded set includes many tier 3 genes (immune checkpoints, epigenetic regulators) whose biology is not captured by expression-based evidence alone.

### 6.4 ML Approach Iteration — A Case Study in Continual Improvement

The evolution of the ML evidence stream illustrates how the pipeline's modular design supports continual learning through iterative refinement:

**Iteration 1 — Supervised Classifier (PR #14, March 10):**
An initial approach trained Random Forest and Logistic Regression classifiers to predict which genes are CRC biomarkers, using pipeline-derived features (omics signal, literature counts, pathway counts) as inputs. Both models achieved AUC 1.0 under 5-fold cross-validation.

**Problem identified:** Perfect AUC was a red flag. The classifier's features were derived from the same pipeline that generated the training labels, creating a circular dependency. The model was not learning biology — it was learning to reproduce its own inputs. This was caught during validation review.

**Iteration 2 — SHAP Explainability (PR #16, March 29):**
The approach was redesigned to ask a fundamentally different question. Instead of predicting biomarker labels, the new module:
1. Trains a Random Forest to classify **tumor vs. normal samples** using raw gene expression
2. Extracts SHAP values to quantify each gene's contribution to the classification
3. Uses mean absolute SHAP importance as an evidence stream in scoring

This avoids circularity because the classification task (tumor vs. normal) is biologically real and independent of the biomarker labels. SHAP values are grounded in Shapley value theory, providing theoretically justified feature attribution.

**Result:** BMP3, an FDA-approved stool DNA biomarker, rose from rank #163 to rank #14 after this change — validating that the SHAP signal captures meaningful expression-based biomarker evidence. The enrichment ratio remained stable at 1.84x because most known biomarkers are mutation-driven (not detectable via expression), confirming that the new stream adds complementary signal without degrading existing performance.

**Key takeaway:** The modular pipeline architecture allowed this fundamental change to the ML stream without affecting any other module. Only `src/ml_importance.py`, `src/scoring.py`, and `config/config.yaml` were modified; omics, PubMed, pathway, and LLM scoring were untouched.

### 6.5 Deployment Update Cycle

The CI/CD pipeline enables a continual deployment loop:

1. **Develop** locally — modify pipeline code or config
2. **Run pipeline** — `python run_pipeline.py --config config/config.yaml`
3. **Validate** — check enrichment ratio in `outputs/validation_summary.json`
4. **Commit & push** to `main`
5. **GitHub Actions** builds new Docker image, pushes to ACR with commit SHA tag
6. **Azure App Service** auto-pulls latest image
7. **Rollback** if needed by redeploying a previous commit SHA tag

This cycle takes approximately 5-10 minutes from push to live deployment.

---

## 7. Deliverables Folder Organization

```
AI_Capstone/
├── config/
│   └── config.yaml                  # Central pipeline configuration
├── data/
│   ├── TCGA-COAD.star_counts.tsv    # Input RNA-seq data (gitignored)
│   └── ensembl_to_symbol_cache.tsv  # Gene mapping cache
├── docs/                            # Project documentation and reports
├── outputs/
│   ├── omics_evidence.csv           # Step 1: DE analysis (36,519 genes)
│   ├── ml_importance.csv            # Step 2: SHAP importance (36,519 genes)
│   ├── lit_evidence.csv             # Step 3: PubMed abstracts (~1,500)
│   ├── pathway_evidence.csv         # Step 4: Pathway annotations
│   ├── llm_scores.csv              # Step 5: LLM relevance scores (459 genes)
│   ├── llm_scores_cache.json        # Step 5: Content-hash LLM cache
│   ├── ranked_candidates.csv        # Step 6: Final rankings (36,519 genes)
│   ├── candidates.csv               # Top 500 candidates
│   ├── validation_report.csv        # Step 7: Per-biomarker report
│   └── validation_summary.json      # Step 7: Aggregate metrics
├── src/
│   ├── omics.py                     # Differential expression analysis
│   ├── ml_importance.py             # Random Forest + SHAP pipeline
│   ├── pubmed.py                    # PubMed literature retrieval
│   ├── pathway.py                   # g:Profiler pathway enrichment
│   ├── llm_scoring.py              # Claude-based literature scoring
│   ├── scoring.py                   # Weighted multi-evidence scoring
│   ├── validation.py                # Known biomarker validation
│   └── gene_mapping.py             # Ensembl ID to symbol mapping
├── web-app/
│   ├── main.py                      # FastAPI application
│   ├── api/
│   │   ├── biomarkers.py            # REST API endpoints
│   │   └── chat.py                  # LLM chat endpoint
│   ├── static/
│   │   ├── index.html               # Single-page web UI
│   │   └── styles.css               # Stylesheet
│   └── nginx/
│       └── nginx.conf               # Reverse proxy config
├── terraform/                       # Azure infrastructure as code
│   ├── main.tf                      # Provider and backend config
│   ├── webapp.tf                    # App Service + Web App
│   ├── acr.tf                       # Container Registry
│   ├── networking.tf                # VNet, subnets, private endpoints
│   ├── storage.tf                   # Azure Files storage
│   └── resource_group.tf            # Resource group
├── .github/workflows/
│   └── deploy.yml                   # CI/CD pipeline
├── Dockerfile                       # Multi-stage container build
├── docker-compose.yml               # Local development orchestration
├── run_pipeline.py                  # Pipeline CLI entrypoint
├── requirements.txt                 # Python dependencies
└── config/config.yaml               # Pipeline configuration
```

*[Screenshot of the deliverables folder structure should be inserted here]*

---

## 8. Installation and Execution

### 8.1 Running the Pipeline

```bash
# Clone the repository
git clone https://github.com/coughlka/AI_Capstone.git
cd AI_Capstone

# Install dependencies
pip install -r requirements.txt

# Set API key for LLM scoring
export ANTHROPIC_API_KEY=<your-key>

# Run the full pipeline
python run_pipeline.py --config config/config.yaml
```

**Reference:** `run_pipeline.py`, `requirements.txt`

### 8.2 Running the Web Application (Local)

```bash
# Option 1: Direct
cd web-app
python -m uvicorn main:app --host 0.0.0.0 --port 8000
# Open http://localhost:8000

# Option 2: Docker Compose
docker-compose up
# Open http://localhost:8080
```

**Reference:** `web-app/main.py`, `docker-compose.yml`

### 8.3 Deploying to Azure

```bash
# Initialize Terraform
cd terraform
terraform init

# Plan and apply infrastructure
terraform plan
terraform apply

# CI/CD handles Docker build + push automatically on merge to main
```

**Reference:** `terraform/`, `.github/workflows/deploy.yml`

---

## 9. Key Scripts and Notebooks Reference

| File | Purpose | How to Run |
|------|---------|------------|
| `run_pipeline.py` | Full 7-step pipeline orchestration | `python run_pipeline.py --config config/config.yaml` |
| `src/omics.py` | Welch's t-test DE on TCGA-COAD | Called by pipeline Step 1 |
| `src/ml_importance.py` | Random Forest + SHAP importance | Called by pipeline Step 2 |
| `src/pubmed.py` | PubMed abstract retrieval | Called by pipeline Step 3 |
| `src/pathway.py` | g:Profiler pathway enrichment | Called by pipeline Step 4 |
| `src/llm_scoring.py` | Claude-based literature scoring | Called by pipeline Step 5 |
| `src/scoring.py` | Multi-evidence weighted scoring | Called by pipeline Step 6 |
| `src/validation.py` | Known biomarker validation | Called by pipeline Step 7 |
| `web-app/main.py` | FastAPI web application | `uvicorn main:app --port 8000` |

---

## 10. Appendices

### Appendix A: Scoring Weight Configuration

Current weights in `config/config.yaml`:

```yaml
scoring:
  weights:
    omics: 0.40
    ml_importance: 0.10
    literature: 0.30
    pathway: 0.20
```

These weights were selected based on domain expertise and validated against the 33 known CRC biomarker panel. The omics stream receives the highest weight as differential expression provides the most direct molecular evidence. Literature receives the second-highest weight as published research contextualizes the biological relevance. Pathway enrichment captures systems-level evidence. ML importance provides a data-driven complement to the hypothesis-driven omics score.

### Appendix B: LLM Scoring Rubric

The Claude model scores each gene's literature evidence on a 0-100 scale:

| Range | Interpretation |
|-------|---------------|
| 80-100 | Direct, established CRC biomarker with clear clinical utility |
| 60-79 | Strong CRC evidence, emerging or not yet clinically validated |
| 40-59 | CRC mentioned but not primary focus; indirect evidence |
| 20-39 | Tangential CRC relevance only |
| 0-19 | No meaningful CRC biomarker evidence |

Criteria weights: Direct CRC evidence (40%), Biomarker potential (30%), Mechanistic insight (20%), Evidence quality (10%).

### Appendix C: Validation Biomarker Panel (80 markers, 3 tiers)

**Tier 1 - FDA/Guideline (14):** KRAS, NRAS, BRAF, EGFR, ERBB2, MLH1, MSH2, MSH6, PMS2, EPCAM, VEGFA, BMP3, NDRG4, DPYD

**Tier 2 - COSMIC CGC (25):** APC, TP53, PIK3CA, SMAD4, FBXW7, PTEN, CTNNB1, RNF43, POLE, TGFBR2, AXIN2, SOX9, TCF7L2, AMER1, ARID1A, SMAD2, ACVR2A, ATM, ZNRF3, BCL9L, RBM10, PCBP1, RPL22, PTPRT, FAM123B

**Tier 3 - Validated Literature (41):** VIM, CDH1, CEACAM5, MKI67, CDX2, DCC, MYC, MSH3, ERBB3, MET, VEGFB, FLT1, KDR, PDGFRA, IGF2, MAP2K1, TGFB1, SMAD3, NOTCH1, JAG1, DLL4, GSK3B, CSNK1A1, DVL2, MYB, ETV4, FOXQ1, CD274, PDCD1, CTLA4, LAG3, DNMT1, DNMT3B, TET2, KDM6A, MUTYH, CHEK2, BRCA1, BRCA2, PALB2, PIK3R1

**Functional Categories:** mutation/tumor-suppressor (27), signaling/pathway (31), MMR/MSI (6), expression/methylation (8), immune-checkpoint (4), epigenetic (4), pharmacogenomic (1)

---

*Note: Screenshots of the web application interface, Azure portal deployment, and GitHub Actions pipeline should be captured and inserted at the marked locations before final submission.*
