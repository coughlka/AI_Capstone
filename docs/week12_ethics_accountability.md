# Week 12: AI Ethics and Accountability

## Purpose

To detail the AI ethics and accountability considerations applied to the CRC biomarker discovery pipeline, building on the cybersecurity controls documented in Week 11.

## AI Ethics

### Responsible Use and Framing

The system is explicitly designed and labeled as a research decision support tool, not a clinical diagnostic. Every user-facing surface, including the web UI dashboard, API documentation, and the LLM chat system prompt, makes clear that the rankings represent computational prioritization of publicly available evidence and not clinical recommendations. No component of this system should be used to inform patient care without independent experimental validation and regulatory review.

The LLM chat system prompt explicitly instructs Claude to avoid clinical claims and to redirect users toward published literature when asked about treatment or diagnostic decisions. This guardrail is documented in web-app/api/chat.py and is applied to every conversation regardless of user input.

### Data Ethics

All data sources used by the pipeline are publicly available under established data sharing frameworks. The TCGA-COAD RNA-seq data is distributed under the TCGA Data Use Agreement governing public-access tier data. PubMed abstracts are retrieved through the official NCBI E-Utilities API within its rate limits and terms of service (3 requests per second with an API key). Pathway data is retrieved from g:Profiler under its academic use terms. No proprietary clinical datasets, restricted-access genomic data, or patient identifiers are used at any point in the pipeline.

The system processes only aggregate gene expression statistics. Individual patient barcodes from the TCGA dataset are used to label samples as tumor or normal during preprocessing, but no patient demographic, clinical outcome, or identifying information is stored, logged, or served through any API endpoint.

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

## Accountability

### Version Control and Traceability

All pipeline code, configuration, and infrastructure definitions are version-controlled in the GitHub repository at https://github.com/coughlka/AI_Capstone. The team uses a feature branch workflow with pull request reviews. Notable PRs include PR #11 (Terraform deployment), PR #14 (initial supervised classifier), PR #16 (SHAP-based redesign after circularity was identified), PR #17 (PubMed alias retry), PR #18 (Docker outputs deployment fix), and PR #19 (ACR webhook).

The PR history provides a complete audit trail of every change, including the rationale for the ML approach redesign and the validation panel expansion from 33 to 80 known biomarkers. Pipeline outputs are committed to the repository at milestone points (commit 4e7d261) so reviewers can access the exact data used to generate the reported results without re-running the pipeline.

Docker images are tagged with both a rolling latest tag and an immutable git commit SHA tag (biomarker-cancer-web-app:{sha}) on every CI/CD build. This provides a complete chain of traceability from the deployed application version back to the source commit and pipeline configuration that produced the results. Rollback to any previous deployment is possible by redeploying the corresponding SHA tag.

### Performance Accountability via Validation

The validation framework (80 known CRC biomarkers organized in three confidence tiers) provides a quantitative performance baseline against which all pipeline changes are evaluated. The enrichment ratio, precision, recall, and per-tier metrics serve as the primary accountability mechanism: any modification to the pipeline that does not improve or maintain these metrics requires explicit justification before merging.

This was applied directly during the ML approach redesign. When the original supervised classifier reported AUC 1.0, the team did not accept the metric at face value because the validation enrichment did not improve in proportion. The investigation revealed circular feature construction, the approach was redesigned, and BMP3 (an FDA-approved CRC biomarker) rose from rank #163 to rank #14 after the SHAP-based replacement. The lesson is documented in the commit history and the final report: high model performance metrics must be interrogated, not just optimized.

The ablation study (Section 4.4 of the final report) provides a similar accountability check on the scoring weights. The current configuration produces 1.462x enrichment, while the optimal configuration found by the ablation produces 1.537x. This 0.075x difference is documented and the decision to retain the current weights is justified explicitly, rather than tuning the weights silently to maximize a single metric.

### Logging and Monitoring

Each pipeline module uses Python's standard logging library to emit structured INFO and WARNING messages throughout execution. Log output includes step completion times, gene counts at each processing stage, API retry events, cache hit and miss rates, and SHAP fold-level metrics. This enables post-hoc execution tracing when diagnosing unexpected results and provides a written record of what the pipeline did during any given run.

The web application uses Uvicorn access logging and emits prefixed log messages from each API endpoint. The chat endpoint logs the user message and the model's response token count (but not the full response content) for monitoring purposes.

### Reproducibility as Accountability

The pipeline is designed to be deterministic given identical inputs. Three specific mechanisms enforce reproducibility:

The LLM scoring module uses temperature=0 and content-hash caching (MD5 of the abstract list per gene) so that re-running the pipeline against an unchanged literature corpus produces byte-identical scores in seconds rather than 15 minutes of API calls. This prevents accidental score drift between runs and ensures that reported metrics can be exactly reproduced.

The Random Forest classifier uses random_state=42 + fold_index for each cross-validation fold, ensuring SHAP importance scores are reproducible across runs given the same input expression matrix.

The full pipeline configuration is centralized in config/config.yaml. All scoring weights, FDR thresholds, top-N selections, and API parameters live in this file. To reproduce a previous run's results exactly, a researcher needs only the source code at a specific git commit and the corresponding config file.

### Team Accountability

Responsibility for each pipeline module is documented in the project README and CLAUDE.md. Keith Coughlin owns the omics, ML importance, LLM scoring, scoring, and pipeline orchestration modules. Ayan Choudhury owns PubMed retrieval, the validation framework, and the MMR deficiency classifier. Gabriel Saenz owns pathway enrichment, the web application, Azure infrastructure, Terraform, and CI/CD.

Pull requests require review by at least one other team member before merging to main. Keith is the designated reviewer for Gabriel's infrastructure changes. This provides accountability for deployment decisions that affect the production Azure environment and ensures that no single team member can unilaterally change the production system.

## Summary

The ethics and accountability posture of this project reflects the principle that responsible AI development requires documented decisions, decomposable outputs, and quantitative checks at every stage. The pipeline does not claim to solve CRC biomarker discovery; it claims to make biomarker prioritization more systematic, transparent, and auditable than the manual literature review process it supplements. Every limitation, including the expression-based bias, the literature bias toward well-studied genes, and the single-cohort validation, is documented in the report and visible in the validation results rather than obscured. Every modeling decision, including the ML approach redesign and the choice to retain current scoring weights despite the ablation findings, is recorded in the git history with explicit justification.

This approach trades some apparent performance for trust. The system reports a 1.46x enrichment ratio rather than tuning weights to claim a higher number; it acknowledges that mutation-driven biomarkers are invisible to expression-based analysis rather than hiding the Tier 2 results; it documents the circular classifier failure rather than quietly removing the original PR. For a research decision support tool intended for use by domain experts, transparency about limitations is more valuable than inflated performance claims.
