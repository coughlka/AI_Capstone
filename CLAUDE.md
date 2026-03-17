# CLAUDE.md — AI Capstone: CRC Biomarker Discovery Pipeline

## Project Overview

An explainable AI pipeline for identifying and ranking colorectal cancer (CRC) biomarker candidates by synthesizing four evidence streams: omics (RNA-seq DE analysis), literature (PubMed), pathway enrichment (g:Profiler), and LLM scoring (Claude). Outputs a ranked list of 500 candidates with a web UI for exploration.

**GitHub:** https://github.com/coughlka/AI_Capstone
**Primary branch:** `main`
**Current working branch pattern:** `feature/<name>` or `gsaenz1/<name>`

---

## Team

| Member | GitHub | Modules Owned |
|--------|--------|---------------|
| Keith (you) | coughlka | `src/omics.py`, `src/scoring.py`, `src/llm_scoring.py`, `run_pipeline.py` |
| Ayan | — | `src/pubmed.py` |
| Gabriel | burrito227 (gsaenz1) | `src/pathway.py`, `web-app/`, Terraform/Azure deployment |

---

## Pipeline Architecture

### Execution Order (`run_pipeline.py`)

```
Step 1: Omics      → outputs/omics_evidence.csv + outputs/candidates.csv
Step 2: PubMed     → outputs/lit_evidence.csv
Step 3: Pathway    → outputs/pathway_evidence.csv
Step 4: LLM Score  → outputs/llm_scores.csv
Step 5: Scoring    → outputs/ranked_candidates.csv
Step 6: Validation → outputs/validation_report.csv + outputs/validation_summary.json
```

### Run the Pipeline

```bash
python run_pipeline.py --config config/config.yaml
```

### Launch the Web App (local)

```bash
cd web-app
python -m uvicorn main:app --host 0.0.0.0 --port 8000
# Open http://localhost:8000
```

---

## Key Source Files

| File | Purpose |
|------|---------|
| `run_pipeline.py` | CLI entrypoint, orchestrates all 6 steps |
| `config/config.yaml` | Central config — paths, weights, API keys, thresholds |
| `src/omics.py` | Welch's t-test DE analysis on TCGA-COAD, BH FDR correction |
| `src/pubmed.py` | NCBI E-Utilities esearch/efetch, rate-limited, up to 5 abstracts/gene |
| `src/pathway.py` | g:Profiler POST API, GO:BP + KEGG + Reactome + WikiPathways |
| `src/llm_scoring.py` | Claude API, scores each gene 0-100 on biomarker relevance |
| `src/scoring.py` | Min-max normalize each stream, weighted combination |
| `src/validation.py` | Enrichment ratio vs 33 known CRC biomarkers |
| `src/gene_mapping.py` | Ensembl ID → HGNC symbol via mygene.info (cached) |
| `web-app/main.py` | FastAPI app, CORS, static file serving |
| `web-app/api/biomarkers.py` | REST endpoints: /candidates, /genes/{id}, /stats, /validation |
| `web-app/api/chat.py` | Claude chat endpoint with tiered context (top 50 genes) |

---

## Scoring System

```
omics_signal    = |log2FC| × -log10(FDR + 1e-300)
omics_score     = min_max_normalize(omics_signal) → [0, 100]

literature_score = LLM score from Claude (0-100) if available, else abstract count normalized

pathway_score   = min_max_normalize(pathway_count) → [0, 100]

final_score     = 0.45 × omics_score + 0.35 × literature_score + 0.20 × pathway_score
```

Weights are configurable in `config/config.yaml` under `scoring.weights`.

---

## Data

- **Input:** `data/TCGA-COAD.star_counts.tsv` — 36,519 genes × 514 samples (473 tumor, 41 normal)
- **Download from:** UCSC Xena Browser (gitignored, not in repo)
- **Gene mapping cache:** `data/ensembl_to_symbol_cache.tsv`

---

## Current Results (as of March 2026)

- 25,800 significant genes (FDR < 0.05)
- 500 candidates selected by strongest DE signal
- ~395 genes with PubMed abstracts (~1,570 total)
- ~382 genes with pathway annotations
- LLM scores range 0-95, mean ~49
- **Validation enrichment: 2.07x** (known CRC biomarkers ranked 2.07x better than random)
- Best known biomarker: BMP3 at rank #23 (final score 60.52)

---

## Infrastructure

### Local Dev

```bash
docker-compose up   # FastAPI + outputs volume
```

### Azure Deployment (Gabriel owns this)

- Container registry: ACR (`psuai894acr`)
- App Service: `psuai894webapp` (B1 tier, Linux)
- Terraform state: `PSU-AI894-State-RG` / `psuai894storage`
- CI/CD: `.github/workflows/deploy.yml` builds and pushes Docker image to ACR on push to `main`
- **Known issue:** ACR webhook not provisioned — auto-redeploy won't fire without it

### Secrets Required

| Secret | Where Used |
|--------|-----------|
| `ANTHROPIC_API_KEY` | `web-app/.env` — LLM scoring + chat |
| `ACR_LOGIN_SERVER` | GitHub Actions secret |
| `ACR_USERNAME` | GitHub Actions secret |
| `ACR_PASSWORD` | GitHub Actions secret |
| PubMed API key | `config/config.yaml` (currently hardcoded — should move to env var) |

---

## Known Technical Debt

1. **PubMed API key hardcoded** in `config/config.yaml` — should be an env var
2. **ACR admin password** stored as plain-text app setting in `webapp.tf` — should use managed identity
3. **No ACR webhook** provisioned in Terraform — auto-deploy pipeline is broken
4. **No LLM scoring resume** — if interrupted, all 400 genes are re-scored from scratch (~15 min)
5. **Storage account** has public endpoint open despite private endpoint being configured
6. **Redundant data source** in `terraform/resource_group.tf` — never referenced

---

## Upcoming Work (Week of March 10)

| Task | Owner | Due | Notes |
|------|-------|-----|-------|
| Expand known biomarker labels (COSMIC/OncoKB/ClinVar, target 60-80+) | Ayan | 3/12 | |
| Train LR on scoring weights | Ayan | 3/14 | Compare to hand-tuned weights |
| Deploy web app to Azure | Gabriel | 3/12 | PR #11 open, needs webhook fix |
| Add pipeline staleness check | Gabriel | 3/14 | Content for Continual Learning writeup (Mar 29) |
| Add web app input sanitization | Gabriel | 3/15 | Content for Cybersecurity writeup (Apr 12) |
| Scoring weight ablation study | Keith | TBD | Complements Ayan's LR work |
| Precision/recall eval upgrade | Keith | TBD | Stronger than single enrichment ratio metric |

---

## Key Writeup Deadlines

| Writeup | Due |
|---------|-----|
| Continual Learning | Mar 29 |
| Cybersecurity | Apr 12 |

---

## PR Workflow

- Open PRs against `main`
- Keith is reviewer on Gabriel's PRs (requested reviewer)
- Use `gh pr` CLI for all PR operations
- Current open PR: #11 (Gabriel — Terraform deployment)
