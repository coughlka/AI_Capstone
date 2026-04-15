# Week 11: AI Cybersecurity

## Purpose

To detail the cybersecurity controls applied to the CRC biomarker discovery pipeline and its deployed web application, covering credentials management, network security, application-level defenses, and known limitations.

## Threat Model

The system has a narrow threat surface compared to a typical clinical or commercial application. It processes only public TCGA gene expression data, stores no patient-level information, and exposes a read-mostly API. The realistic threats are: leaked API credentials (Anthropic, NCBI, Azure), unauthorized access to the Azure deployment, prompt injection or abuse of the LLM chat endpoint, supply-chain compromise via third-party Python packages, and accidental commit of secrets to the public GitHub repository. The security controls described below address each of these threats directly.

## Credentials and Secrets Management

All API keys and access credentials are managed outside the source code repository. The Anthropic API key required for LLM scoring and the chat endpoint is loaded from a web-app/.env file that is listed in .gitignore and never committed. The .env file is generated locally on each developer machine and on the Azure App Service via app settings, but never enters version control. A second .gitignore entry covers any *.env variants to prevent accidental commits with alternative naming.

The PubMed (NCBI E-Utilities) API key was originally hardcoded in config/config.yaml during early development. This was identified as technical debt and migrated to an environment variable in a subsequent PR. The current pipeline reads NCBI_API_KEY from the environment, falling back to unauthenticated requests with the slower rate limit (1 request per second instead of 3) if the key is not provided.

CI/CD secrets used by the GitHub Actions workflow (ACR_LOGIN_SERVER, ACR_USERNAME, ACR_PASSWORD) are stored as encrypted GitHub Actions Secrets at the repository level. These secrets are injected into the workflow environment at build time and are never written to logs, build artifacts, or container layers. GitHub automatically masks secret values in workflow output, so even an accidental echo command would not leak them.

The Azure Container Registry admin password is currently stored as a plain-text app setting in webapp.tf, which was flagged as a known issue during the deployment review. Migration to Azure Managed Identity, which would eliminate this credential from Terraform state entirely, is identified as future work and is tracked in the technical debt log.

## Network Security

The Azure deployment applies defense-in-depth network controls through a Terraform-managed virtual network. The App Service is deployed with VNet integration, routing all outbound traffic through a dedicated web subnet (10.0.1.0/24). A separate private endpoint subnet (10.0.2.0/24) hosts the private endpoint for the Azure Files storage account, which holds pipeline output artifacts.

The Azure Files storage account is configured with a private DNS zone (privatelink.file.core.windows.net), preventing public internet access to the file share. Even an attacker who obtained the storage account access key could not connect from outside the VNet because the storage account's public network endpoint is disabled. This was a deliberate decision to ensure pipeline outputs never traverse the public internet.

The Azure App Service itself is publicly accessible at https://psuai894webapp.azurewebsites.net, which is required for the demo and reviewer access. Azure provides automatic TLS termination with managed certificates, so all traffic to the application is encrypted in transit. There is no plaintext HTTP fallback.

## Application-Level Security

The FastAPI web application applies several security controls at the application layer.

CORS is restricted via FastAPI's CORSMiddleware to localhost origins during development. Production deployment requires explicit allowlisting of the production domain before broader use. This prevents browsers from making cross-origin requests to the API from unauthorized sites, which is the primary defense against CSRF-style attacks on the read endpoints.

Input validation is enforced through Pydantic models at every API boundary. Query parameters are typed (int, str, Optional) with FastAPI's automatic validation rejecting malformed requests with 422 responses before they reach application logic. The chat endpoint validates the request body (ChatRequest model) for required fields and field types, preventing malformed JSON from reaching the Claude API call.

The API is read-only by design. Of the eight exposed endpoints, only two modify state: /api/reload (clears the in-memory cache) and /api/chat (records nothing persistently but does call the Anthropic API). All data retrieval endpoints use GET, so there is no risk of unintended writes via CSRF or accidental client behavior.

The Docker runtime image is built using a multi-stage Dockerfile that excludes build tools (gcc, build-essential, pip-compile) from the final image. The build stage compiles Python wheels and installs dependencies, then the runtime stage copies only the compiled artifacts and the application code. This reduces the attack surface by eliminating tools that could be exploited if a container were compromised. The runtime image is python:3.11-slim-bullseye with only curl added for health checks.

The local development reverse proxy uses nginx with HTTP basic authentication via .htpasswd, providing a minimal access control layer when running the stack on a developer workstation. The .htpasswd file is gitignored and generated per developer.

## LLM-Specific Security

The chat endpoint introduces a new class of risk because it sends user input to a third-party LLM API. Two specific concerns apply: prompt injection (a user crafts input designed to override the system prompt) and information disclosure (the model reveals information from its training data that the system is not authorized to share).

The current mitigations are partial. The system prompt provides strong context constraints, instructing Claude to discuss only the provided pipeline results, refuse questions about clinical decision-making, and decline to invent studies or biomarker claims. Conversation history is capped at the last 20 messages to limit context drift over long sessions. The Claude API call uses temperature=0 for the literature scorer and a low temperature for chat, reducing the variance of model outputs.

These mitigations do not constitute a formal prompt injection defense. An adversarial user could craft inputs designed to bypass the system prompt by inserting fake instructions or exploiting model behavior patterns. Server-side input sanitization (filtering for instruction-like patterns, length limits, and rate limiting per session) is identified as a known gap. Implementation of a more robust prompt injection defense was assigned to Gabriel Saenz for the cybersecurity writeup deliverable but remains in progress at the time of submission. The chat endpoint should be considered hardened against accidental misuse but not against a determined adversary.

## Data Privacy

No patient-level or personally identifiable information is processed by this system. The TCGA-COAD dataset contains gene expression values aggregated across patient cohorts. The TCGA sample barcodes used to label samples as tumor or normal are public identifiers that do not link to individual patients in any externally identifiable way. No demographic data, clinical outcomes, treatment history, or geographic location data is loaded, stored, or served at any point in the pipeline.

The system processes only aggregate statistical outputs (fold change, FDR, SHAP importance, pathway enrichment) derived from the source data. The web application's API exposes only gene-level rankings and evidence scores; no sample-level data is queryable through any endpoint. TCGA data access complies with the TCGA Data Use Agreement governing public-access tier data, which does not require IRB approval or DUA paperwork because the data has been pre-aggregated and de-identified by the TCGA consortium.

## Supply Chain Security

The pipeline depends on a moderate number of third-party Python packages declared in requirements.txt. The most security-relevant dependencies are anthropic (LLM client), fastapi (web framework), shap (SHAP value computation), scikit-learn (Random Forest), and requests (HTTP client). All dependencies are pinned to specific versions in requirements.txt to prevent silent version drift, and the Docker build uses pip's standard wheel resolution.

GitHub Dependabot is enabled on the repository and provides automatic security advisories for known CVEs in the declared dependencies. As of the final report, no critical or high-severity advisories are outstanding. Routine dependency updates are handled via PR review.

The Docker base image (python:3.11-slim-bullseye) is pulled from the official Python image on Docker Hub. The slim-bullseye variant is a Debian-based minimal image that receives regular security updates. Image rebuilds via the CI/CD pipeline pick up the latest patched base image automatically.

## Logging and Audit

The application uses structured logging at the pipeline and web application layers. Pipeline modules emit module-prefixed INFO messages for each execution stage, which provides a written record of what the pipeline did during any given run. The web application uses Uvicorn access logs, which capture every API request with method, path, status code, and timing.

The chat endpoint logs the user message length and the model's response token count for usage monitoring, but does not log the full content of either. This is a deliberate balance between operational visibility and avoiding accidental retention of user-submitted text that might contain unexpected information.

GitHub Actions logs every CI/CD build with full audit history, including the commit SHA, the triggering user, the build duration, and the success or failure status. Build logs are retained for 90 days under the default GitHub policy, providing a window for incident review if a malicious deployment is suspected.

## Known Limitations and Future Hardening

Three specific cybersecurity gaps are documented as known issues and tracked for future work:

The chat endpoint lacks formal prompt injection protection. The current system prompt and conversation history limit provide partial mitigation but do not prevent a determined adversary from crafting inputs that bypass the context constraints. A planned improvement is server-side input sanitization with pattern-based filtering, length limits, and per-session rate limiting.

The ACR admin password is stored as a plain-text app setting in webapp.tf rather than using Azure Managed Identity. Migration to managed identity would eliminate this credential from Terraform state and improve the secret rotation story.

The ACR webhook for automatic App Service redeployment was identified as missing during the Week 11 deployment work. PR #19 was opened to add the webhook to terraform/acr.tf but is pending review at the time of writing. Until the webhook is provisioned, deployments require a manual App Service restart in the Azure portal after each CI/CD push, which creates a window where the deployed code lags behind the latest commit on main.

## Summary

The cybersecurity posture of this project reflects the realistic threat model of an academic research prototype handling public data. Credentials are managed outside source control. Network access to sensitive resources is restricted by VNet integration and private endpoints. Application inputs are validated by typed Pydantic models. The Docker runtime image excludes build tools to minimize attack surface. The API is read-mostly by design, eliminating most write-side vulnerabilities. The known gaps (prompt injection, ACR password, ACR webhook) are documented rather than hidden and are tracked for future hardening. For a research decision support tool deployed in an academic context, the controls described here provide reasonable defense in depth without imposing operational burden disproportionate to the threat.
