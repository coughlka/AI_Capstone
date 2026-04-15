# Lesson 12: AI Cybersecurity, Ethics and Accountability

## Purpose

To detail the cybersecurity controls, ethical considerations, and accountability mechanisms applied to the CRC biomarker discovery pipeline and its deployed web application.

## Part 1: AI Cybersecurity (Week 11)

### Threat Model

The system has a narrow threat surface compared to a typical clinical or commercial application. It processes only public TCGA gene expression data, stores no patient-level information, and exposes a read-mostly API. The realistic threats are: leaked API credentials (Anthropic, NCBI, Azure), unauthorized access to the Azure deployment, prompt injection or abuse of the LLM chat endpoint, supply-chain compromise via third-party Python packages, and accidental commit of secrets to the public GitHub repository. The security controls described below address each of these threats directly.

### Credentials and Secrets Management

All API keys and access credentials are managed outside the source code repository. The Anthropic API key required for LLM scoring and the chat endpoint is loaded from a web-app/.env file that is listed in .gitignore and never committed. The .env file is generated locally on each developer machine and on the Azure App Service via app settings, but never enters version control. A second .gitignore entry covers any *.env variants to prevent accidental commits with alternative naming.

The PubMed (NCBI E-Utilities) API key was originally hardcoded in config/config.yaml during early development. This was identified as technical debt and migrated to an environment variable in a subsequent PR. The current pipeline reads NCBI_API_KEY from the environment, falling back to unauthenticated requests with the slower rate limit (1 request per second instead of 3) if the key is not provided.

CI/CD secrets used by the GitHub Actions workflow (ACR_LOGIN_SERVER, ACR_USERNAME, ACR_PASSWORD) are stored as encrypted GitHub Actions Secrets at the repository level. These secrets are injected into the workflow environment at build time and are never written to logs, build artifacts, or container layers. GitHub automatically masks secret values in workflow output, so even an accidental echo command would not leak them.

The Azure Container Registry admin password is currently stored as a plain-text app setting in webapp.tf, which was flagged as a known issue during the deployment review. Migration to Azure Managed Identity, which would eliminate this credential from Terraform state entirely, is identified as future work and is tracked in the technical debt log.

### Network Security

The Azure deployment applies defense-in-depth network controls through a Terraform-managed virtual network. The App Service is deployed with VNet integration, routing all outbound traffic through a dedicated web subnet (10.0.1.0/24). A separate private endpoint subnet (10.0.2.0/24) hosts the private endpoint for the Azure Files storage account, which holds pipeline output artifacts.

The Azure Files storage account is configured with a private DNS zone (privatelink.file.core.windows.net), preventing public internet access to the file share. Even an attacker who obtained the storage account access key could not connect from outside the VNet because the storage account's public network endpoint is disabled. This was a deliberate decision to ensure pipeline outputs never traverse the public internet.

The Azure App Service itself is publicly accessible at https://psuai894webapp.azurewebsites.net, which is required for the demo and reviewer access. Azure provides automatic TLS termination with managed certificates, so all traffic to the application is encrypted in transit. There is no plaintext HTTP fallback.

### Application-Level Security

The FastAPI web application applies several security controls at the application layer.

CORS is restricted via FastAPI's CORSMiddleware to localhost origins during development. Production deployment requires explicit allowlisting of the production domain before broader use. This prevents browsers from making cross-origin requests to the API from unauthorized sites, which is the primary defense against CSRF-style attacks on the read endpoints.

Input validation is enforced through Pydantic models at every API boundary. Query parameters are typed (int, str, Optional) with FastAPI's automatic validation rejecting malformed requests with 422 responses before they reach application logic. The chat endpoint validates the request body (ChatRequest model) for required fields and field types, preventing malformed JSON from reaching the Claude API call.

The API is read-only by design. Of the eight exposed endpoints, only two modify state: /api/reload (clears the in-memory cache) and /api/chat (records nothing persistently but does call the Anthropic API). All data retrieval endpoints use GET, so there is no risk of unintended writes via CSRF or accidental client behavior.

The Docker runtime image is built using a multi-stage Dockerfile that excludes build tools (gcc, build-essential, pip-compile) from the final image. The build stage compiles Python wheels and installs dependencies, then the runtime stage copies only the compiled artifacts and the application code. This reduces the attack surface by eliminating tools that could be exploited if a container were compromised. The runtime image is python:3.11-slim-bullseye with only curl added for health checks.

The local development reverse proxy uses nginx with HTTP basic authentication via .htpasswd, providing a minimal access control layer when running the stack on a developer workstation. The .htpasswd file is gitignored and generated per developer.

### LLM-Specific Security

The chat endpoint introduces a new class of risk because it sends user input to a third-party LLM API. Two specific concerns apply: prompt injection (a user crafts input designed to override the system prompt) and information disclosure (the model reveals information from its training data that the system is not authorized to share).

The current mitigations are partial. The system prompt provides strong context constraints, instructing Claude to discuss only the provided pipeline results, refuse questions about clinical decision-making, and decline to invent studies or biomarker claims. Conversation history is capped at the last 20 messages to limit context drift over long sessions. The Claude API call uses temperature=0 for the literature scorer and a low temperature for chat, reducing the variance of model outputs.

These mitigations do not constitute a formal prompt injection defense. An adversarial user could craft inputs designed to bypass the system prompt by inserting fake instructions or exploiting model behavior patterns. Server-side input sanitization (filtering for instruction-like patterns, length limits, and rate limiting per session) is identified as a known gap and remains in progress at the time of submission. The chat endpoint should be considered hardened against accidental misuse but not against a determined adversary.

### Supply Chain Security

The pipeline depends on a moderate number of third-party Python packages declared in requirements.txt. The most security-relevant dependencies are anthropic (LLM client), fastapi (web framework), shap (SHAP value computation), scikit-learn (Random Forest), and requests (HTTP client). All dependencies are pinned to specific versions in requirements.txt to prevent silent version drift, and the Docker build uses pip's standard wheel resolution.

GitHub Dependabot is enabled on the repository and provides automatic security advisories for known CVEs in the declared dependencies. As of the final report, no critical or high-severity advisories are outstanding. Routine dependency updates are handled via PR review.

The Docker base image (python:3.11-slim-bullseye) is pulled from the official Python image on Docker Hub. The slim-bullseye variant is a Debian-based minimal image that receives regular security updates. Image rebuilds via the CI/CD pipeline pick up the latest patched base image automatically.

### Known Cybersecurity Limitations

Three specific cybersecurity gaps are documented as known issues and tracked for future work:

The chat endpoint lacks formal prompt injection protection. The current system prompt and conversation history limit provide partial mitigation but do not prevent a determined adversary from crafting inputs that bypass the context constraints. A planned improvement is server-side input sanitization with pattern-based filtering, length limits, and per-session rate limiting.

The ACR admin password is stored as a plain-text app setting in webapp.tf rather than using Azure Managed Identity. Migration to managed identity would eliminate this credential from Terraform state and improve the secret rotation story.

The ACR webhook for automatic App Service redeployment was identified as missing during the Week 11 deployment work. PR #19 was opened to add the webhook to terraform/acr.tf but is pending review at the time of writing. Until the webhook is provisioned, deployments require a manual App Service restart in the Azure portal after each CI/CD push, which creates a window where the deployed code lags behind the latest commit on main.

## Part 2: AI Ethics and Accountability (Week 12)

### Responsible Use and Framing

The system is explicitly designed and labeled as a research decision support tool, not a clinical diagnostic. Every user-facing surface, including the web UI dashboard, API documentation, and the LLM chat system prompt, makes clear that the rankings represent computational prioritization of publicly available evidence and not clinical recommendations. No component of this system should be used to inform patient care without independent experimental validation and regulatory review.

The LLM chat system prompt explicitly instructs Claude to avoid clinical claims and to redirect users toward published literature when asked about treatment or diagnostic decisions. This guardrail is documented in web-app/api/chat.py and is applied to every conversation regardless of user input.

### Data Ethics and Privacy

All data sources used by the pipeline are publicly available under established data sharing frameworks. The TCGA-COAD RNA-seq data is distributed under the TCGA Data Use Agreement governing public-access tier data. PubMed abstracts are retrieved through the official NCBI E-Utilities API within its rate limits and terms of service (3 requests per second with an API key). Pathway data is retrieved from g:Profiler under its academic use terms. No proprietary clinical datasets, restricted-access genomic data, or patient identifiers are used at any point in the pipeline.

No patient-level or personally identifiable information is processed by this system. The TCGA sample barcodes used to label samples as tumor or normal are public identifiers that do not link to individual patients in any externally identifiable way. The system processes only aggregate statistical outputs (fold change, FDR, SHAP importance, pathway enrichment) derived from the source data. The web application's API exposes only gene-level rankings and evidence scores; no sample-level data is queryable through any endpoint. TCGA data access complies with the TCGA Data Use Agreement, which does not require IRB approval because the data has been pre-aggregated and de-identified by the TCGA consortium.

### Transparency and Explainability

Transparency is a core design principle of the pipeline. Every gene's final score is fully decomposable into its four component contributions (omics, ML importance, literature, pathway), each visible in the web UI gene detail view. A researcher can immediately see why a candidate was ranked highly, distinguishing genes elevated by a single strong signal from genes with convergent multi-stream evidence.

Three specific transparency mechanisms support auditability:

The LLM scoring module stores a 1-2 sentence rationale alongside each score in outputs/llm_scores.csv. When the web UI displays a gene's literature score of 75, it also displays the model's written justification, allowing users to evaluate whether the model's reasoning matches their domain expertise.

The ML evidence stream uses SHAP (SHapley Additive exPlanations) values, which are theoretically grounded in cooperative game theory and provide mathematically rigorous feature attribution. This was a deliberate replacement for an earlier supervised classifier that achieved AUC 1.0 by learning circular dependencies between features and labels. The SHAP redesign (PR #16) is documented in the repository commit history as a case study in catching and correcting an unsound modeling decision.

Scoring weights are stored in config/config.yaml and modifiable without code changes. The ablation study (src/ablation.py) provides quantitative justification for the current weight configuration by testing 246 alternatives. Users who disagree with the default weights can run their own configurations and immediately see the effect on the ranking.

### Bias Considerations

The pipeline has known systematic biases that are documented in the report and the codebase rather than obscured.

Literature bias: Genes that have been heavily studied have more PubMed abstracts and consequently higher LLM scores. Novel or understudied genes with genuine biomarker potential may be disadvantaged. This bias is partially mitigated by the 40 percent omics weight, which treats all genes equally regardless of prior study, and by the gene alias retry logic in src/pubmed.py, which queries PubMed using common alternative gene names to capture literature for genes with naming inconsistencies.

Expression bias: The pipeline cannot detect mutation-driven biomarkers because these genes do not necessarily show strong differential expression. The Tier 2 (COSMIC Cancer Gene Census) validation results make this limitation visible: zero Tier 2 markers appear in the top 500 candidates and the median rank for Tier 2 markers is 9,557. This is reported transparently rather than hidden.

Cohort bias: Results are based on a single TCGA-COAD cohort, which is predominantly European-ancestry and North American. Biomarker rankings may not generalize to other populations. Multi-cohort validation against TCGA-READ and independent GEO datasets is identified as the primary future work priority.

### Ethical Use of LLMs

The pipeline uses Claude in two distinct roles, each with appropriate constraints. As a literature scorer (src/llm_scoring.py), Claude is called with temperature=0 to ensure deterministic, reproducible scores. The system prompt restricts the model to evaluating provided abstracts, with explicit instructions not to invent studies or score genes based on training-data knowledge. Each score is bounded to a 0-100 integer with a written rationale, which is stored for auditability.

As an interactive chat assistant (web-app/api/chat.py), Claude is given a tiered context including the dataset summary, top 50 candidates, and optional gene-specific detail. The system prompt instructs the model to only discuss the provided pipeline results and to refuse questions about clinical decision-making. Conversation history is limited to the last 20 messages to prevent context drift.

### Version Control and Traceability

All pipeline code, configuration, and infrastructure definitions are version-controlled in the GitHub repository at https://github.com/coughlka/AI_Capstone. The team uses a feature branch workflow with pull request reviews. Notable PRs include PR #11 (Terraform deployment), PR #14 (initial supervised classifier), PR #16 (SHAP-based redesign after circularity was identified), PR #17 (PubMed alias retry), PR #18 (Docker outputs deployment fix), and PR #19 (ACR webhook).

The PR history provides a complete audit trail of every change, including the rationale for the ML approach redesign and the validation panel expansion from 33 to 80 known biomarkers. Pipeline outputs are committed to the repository at milestone points (commit 4e7d261) so reviewers can access the exact data used to generate the reported results without re-running the pipeline.

Docker images are tagged with both a rolling latest tag and an immutable git commit SHA tag (biomarker-cancer-web-app:{sha}) on every CI/CD build. This provides a complete chain of traceability from the deployed application version back to the source commit and pipeline configuration that produced the results. Rollback to any previous deployment is possible by redeploying the corresponding SHA tag.

### Performance Accountability via Validation

The validation framework (80 known CRC biomarkers organized in three confidence tiers) provides a quantitative performance baseline against which all pipeline changes are evaluated. The enrichment ratio, precision, recall, and per-tier metrics serve as the primary accountability mechanism: any modification to the pipeline that does not improve or maintain these metrics requires explicit justification before merging.

This was applied directly during the ML approach redesign. When the original supervised classifier reported AUC 1.0, the team did not accept the metric at face value because the validation enrichment did not improve in proportion. The investigation revealed circular feature construction, the approach was redesigned, and BMP3 (an FDA-approved CRC biomarker) rose from rank #163 to rank #14 after the SHAP-based replacement. The lesson is documented in the commit history and the final report: high model performance metrics must be interrogated, not just optimized.

The ablation study provides a similar accountability check on the scoring weights. The current configuration produces 1.462x enrichment, while the optimal configuration found by the ablation produces 1.537x. This 0.075x difference is documented and the decision to retain the current weights is justified explicitly, rather than tuning the weights silently to maximize a single metric.

### Logging and Monitoring

Each pipeline module uses Python's standard logging library to emit structured INFO and WARNING messages throughout execution. Log output includes step completion times, gene counts at each processing stage, API retry events, cache hit and miss rates, and SHAP fold-level metrics. This enables post-hoc execution tracing when diagnosing unexpected results and provides a written record of what the pipeline did during any given run.

The web application uses Uvicorn access logging and emits prefixed log messages from each API endpoint. The chat endpoint logs the user message length and the model's response token count (but not the full response content) for monitoring purposes. GitHub Actions logs every CI/CD build with full audit history, including the commit SHA, the triggering user, the build duration, and the success or failure status.

### Reproducibility as Accountability

The pipeline is designed to be deterministic given identical inputs. Three specific mechanisms enforce reproducibility:

The LLM scoring module uses temperature=0 and content-hash caching (MD5 of the abstract list per gene) so that re-running the pipeline against an unchanged literature corpus produces byte-identical scores in seconds rather than 15 minutes of API calls. This prevents accidental score drift between runs and ensures that reported metrics can be exactly reproduced.

The Random Forest classifier uses random_state=42 + fold_index for each cross-validation fold, ensuring SHAP importance scores are reproducible across runs given the same input expression matrix.

The full pipeline configuration is centralized in config/config.yaml. All scoring weights, FDR thresholds, top-N selections, and API parameters live in this file. To reproduce a previous run's results exactly, a researcher needs only the source code at a specific git commit and the corresponding config file.

### Team Accountability

Responsibility for each pipeline module is documented in the project README. Keith Coughlin owns the omics, ML importance, LLM scoring, scoring, and pipeline orchestration modules. Ayan Choudhury owns PubMed retrieval, the validation framework, and the MMR deficiency classifier. Gabriel Saenz owns pathway enrichment, the web application, Azure infrastructure, Terraform, and CI/CD.

Pull requests require review by at least one other team member before merging to main. Keith is the designated reviewer for Gabriel's infrastructure changes. This provides accountability for deployment decisions that affect the production Azure environment and ensures that no single team member can unilaterally change the production system.

## Summary

The cybersecurity, ethics, and accountability posture of this project reflects the realistic threat model of an academic research prototype handling public data. Credentials are managed outside source control, network access to sensitive resources is restricted by VNet integration and private endpoints, application inputs are validated by typed Pydantic models, and the Docker runtime image excludes build tools to minimize attack surface. The known security gaps (prompt injection, ACR password, ACR webhook) are documented rather than hidden and are tracked for future hardening.

On the ethics side, the system does not claim to solve CRC biomarker discovery; it claims to make biomarker prioritization more systematic, transparent, and auditable than the manual literature review process it supplements. Every limitation, including the expression-based bias, the literature bias toward well-studied genes, and the single-cohort validation, is documented in the report and visible in the validation results. Every modeling decision, including the ML approach redesign and the choice to retain current scoring weights despite the ablation findings, is recorded in the git history with explicit justification.

This approach trades some apparent performance for trust. The system reports a 1.46x enrichment ratio rather than tuning weights to claim a higher number; it acknowledges that mutation-driven biomarkers are invisible to expression-based analysis rather than hiding the Tier 2 results; it documents the circular classifier failure rather than quietly removing the original PR. For a research decision support tool intended for use by domain experts, transparency about limitations is more valuable than inflated performance claims.
