================================================================================
AI 894: Design & Implementation of AI Systems - Final Project
Spring 2026 - Team 3
================================================================================

PROJECT TITLE
-------------
Explainable AI for Colorectal Cancer Biomarker Evidence Synthesis and Ranking

TEAM
----
Ayan Choudhury
Keith Coughlin
Gabriel Saenz

ABSTRACT
--------
Biomarker prioritization for colorectal cancer (CRC) is a manual, fragmented
process in which researchers must weigh evidence from transcriptomic datasets,
PubMed literature, and biological pathway databases without a consistent,
auditable framework. This project delivers an explainable AI pipeline that
synthesizes four evidence streams - differential gene expression from
TCGA-COAD RNA-seq, PubMed abstract retrieval and LLM-based literature
scoring, multi-database pathway enrichment, and supervised ML feature
importance - into a single ranked list of CRC biomarker candidates. The
pipeline is validated against an 80-gene panel of known biomarkers drawn
from FDA guidelines, the COSMIC Cancer Gene Census, and clinical
literature; known biomarkers rank 1.46x better than chance (median rank
12,493 versus an expected 18,259 under a uniform null), and all 14
Tier-1 FDA / guideline biomarkers are recovered from the 36,519-gene
universe. A companion mismatch-repair (MMR) deficiency classifier
(LR / RF / XGBoost on TCGA-COAD IHC labels) is included as a second,
complementary analysis. All results are served through a FastAPI web
application with an interactive Claude chat endpoint, deployed to Azure
via Terraform-managed infrastructure.

================================================================================
INSTALLATION
================================================================================

REQUIREMENTS
------------
- Python 3.11 or later
- pip (or uv / conda)
- 8 GB RAM minimum (pipeline peak memory ~5 GB during DE analysis)
- ~3 GB free disk (expression matrix + caches + outputs)
- Optional: Docker + docker-compose for containerized web app
- Optional: Azure CLI + Terraform 1.6+ for cloud deployment

ENVIRONMENT VARIABLES
---------------------
Before running, create a .env file in the project root (or web-app/.env
for the web app) containing:

    ANTHROPIC_API_KEY=<your Claude API key>
    NCBI_API_KEY=<your NCBI E-Utilities key, optional but recommended>

A .env.example is provided as a template.

SETUP STEPS
-----------
1. Extract the Code.zip archive to a working directory.

2. (Optional but recommended) Create and activate a virtual environment:

       python -m venv .venv
       source .venv/bin/activate      # macOS / Linux
       .venv\Scripts\activate         # Windows

3. Install Python dependencies:

       pip install -r requirements.txt

4. Extract the dataset.zip archive into the data/ subdirectory. After
   extraction the following files must be present:

       data/TCGA-COAD.star_counts.tsv
       data/TCGA-COAD.GDC_phenotype.tsv
       data/ensembl_to_symbol_cache.tsv

================================================================================
EXECUTION
================================================================================

RUN THE FULL PIPELINE (STEPS 1-6)
---------------------------------
From the project root:

    python run_pipeline.py --config config/config.yaml

Execution order:
    Step 1: Omics differential expression   -> outputs/omics_evidence.csv
                                               outputs/candidates.csv
    Step 2: PubMed literature retrieval     -> outputs/lit_evidence.csv
    Step 3: Pathway enrichment (g:Profiler) -> outputs/pathway_evidence.csv
    Step 4: LLM scoring via Claude API      -> outputs/llm_scores.csv
    Step 5: Weighted scoring and ranking    -> outputs/ranked_candidates.csv
    Step 6: Validation vs known biomarkers  -> outputs/validation_report.csv
                                               outputs/validation_summary.json

Typical wall-clock: ~20 - 30 minutes on a modern laptop, dominated by
Step 2 (rate-limited NCBI calls) and Step 4 (LLM calls, ~400 genes).
LLM scores are content-hash cached; re-running with unchanged inputs
completes in under one minute.

RUN THE MMR DEFICIENCY CLASSIFIER
---------------------------------
A separate, standalone pipeline trains and evaluates LR / RF / XGBoost
classifiers on TCGA-COAD IHC-based MMR-deficiency labels:

    python run_mmr_classifier.py --config config/config.yaml

Outputs are written to outputs/coad_mmr_classifier/ and include CV
metrics, holdout results, SHAP feature importance, and an overlap
comparison against the main biomarker ranking.

RUN THE WEB APPLICATION (LOCAL)
-------------------------------
From the project root:

    cd web-app
    python -m uvicorn main:app --host 0.0.0.0 --port 8000

Then open http://localhost:8000 in a browser. The web app serves the
ranked candidate table, per-gene evidence details, validation dashboard,
and an LLM chat endpoint backed by the top 50 ranked genes.

Alternatively, with Docker:

    docker-compose up

DEPLOY TO AZURE
---------------
Infrastructure is defined in terraform/. From the project root:

    cd terraform
    terraform init
    terraform plan
    terraform apply

This provisions an Azure Container Registry, App Service (B1, Linux),
and associated networking. CI/CD is configured via
.github/workflows/deploy.yml; pushes to main trigger a Docker build,
push to ACR, and automated redeployment.

REPRODUCIBILITY NOTES
---------------------
- Random seeds are fixed in scoring.py and src/mmr/train_classifier.py
  (random_state = 42).
- LLM scoring is cached by content hash; identical gene inputs produce
  identical scores across runs.
- Pipeline outputs (outputs/*.csv, outputs/*.json) are produced
  deterministically from the inputs in data/ and the weights in
  config/config.yaml.

================================================================================
REPOSITORY LAYOUT (Code.zip)
================================================================================

    run_pipeline.py              Main pipeline entry point
    run_mmr_classifier.py        MMR classifier entry point
    config/config.yaml           All paths, weights, and hyperparameters
    src/                         Pipeline modules (omics, pubmed, pathway,
                                 llm_scoring, scoring, validation, mmr/)
    web-app/                     FastAPI application (api/, static/, main.py)
    terraform/                   Azure infrastructure as code
    .github/workflows/           CI/CD (deploy.yml)
    requirements.txt             Python dependencies
    Dockerfile, docker-compose.yml
    tests/                       Unit tests
    docs/                        Working drafts and supplementary notes

OUTPUTS (generated)
-------------------
    outputs/ranked_candidates.csv
    outputs/validation_summary.json
    outputs/coad_mmr_classifier/
    outputs/models/              (pickled scikit-learn pipelines)

================================================================================
CONTACT
================================================================================
Ayan Choudhury   - MMR deficiency classifier, PubMed module
Keith Coughlin   - Omics, LLM scoring, pipeline orchestration, ML classifier
Gabriel Saenz    - Pathway enrichment, web application, Azure deployment

Repository: https://github.com/coughlka/AI_Capstone
