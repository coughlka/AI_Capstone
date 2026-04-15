# Discussion of Results

## Usefulness in Solving the Problem

The pipeline delivers a reproducible, multi-evidence ranking of colorectal cancer
(CRC) biomarker candidates that directly addresses the manual, fragmented
prioritization process described in our background section. Against an 80-gene
panel of known CRC biomarkers drawn from FDA/guideline sources, the COSMIC
Cancer Gene Census, and validated literature, 79 of 80 panel genes were
recovered from the TCGA-COAD expression matrix. Known biomarkers ranked
**1.46× better than chance** on the final score (median rank 12,493 vs. an
expected median of 18,259 under a uniform null), with **4 known biomarkers in
the top 500**. Tier-1 FDA/guideline markers performed best, with BMP3 ranking
at position 14 and all 14 Tier-1 markers successfully recovered. The
FOXQ1 gene led the literature tier at rank 53. These results demonstrate that
the weighted fusion of differential expression, literature signal, pathway
enrichment, LLM scoring, and supervised ML importance produces a ranking that
is both reproducible and interpretable — each candidate gene is backed by an
auditable evidence trail across all four streams.

The companion MMR-deficiency classifier answers a distinct biological question:
within the COAD cohort, can a supervised model distinguish mismatch repair
proficient from deficient tumors? Using 58 IHC-labeled samples (30 MSS / 28 MSI),
random forest achieved a tuned F1 of 0.67 (PR-AUC 0.72, threshold 0.1) on the
holdout set. SHAP feature attributions from the best model overlapped 99% with
our candidate pool at the gene-universe level but only 3% with the top-500
ranked biomarkers — a noteworthy finding that suggests the MMR-deficiency
signal is largely orthogonal to the generalized CRC biomarker ranking and
should be treated as a complementary, not redundant, evidence stream.

## Practical Implications

**Scientific.** The pipeline provides a template for evidence-triangulated
biomarker prioritization that can be redeployed to other cancer types by
swapping the TCGA source, the known-biomarker validation panel, and the
LLM prompt. The explicit, auditable scoring formula (omics 0.40 / ML importance
0.10 / literature 0.30 / pathway 0.20) allows collaborating labs to inspect
*why* any given gene was ranked where it was — a property that black-box
alternatives cannot offer.

**Commercial / Translational.** In a validation-limited environment where wet-lab
confirmation of a single biomarker candidate can cost tens of thousands of
dollars, a 1.46× concentration of known-good candidates toward the top of the
ranked list is economically meaningful: it increases the probability that
downstream validation efforts are spent on genes with real biological signal
rather than noise. The web application layer allows non-engineering stakeholders
(wet-lab collaborators, clinicians) to browse the ranked list and interact with
the LLM chat for context, lowering the barrier to adoption.

**Social / Clinical.** Faster, cheaper biomarker triage has the potential to
reduce the time between discovery-stage research and clinically validated
diagnostic panels. For CRC specifically — a cancer with well-documented
survival benefit from early detection — accelerating the biomarker-discovery
pipeline has direct patient-outcome implications if the downstream validation
and regulatory path can keep pace.

## Limitations

**Data coverage.** TCGA-COAD is a single cohort of 514 samples (473 tumor,
41 normal), heavily weighted toward U.S. patient populations. The ranking is
therefore biased toward variants and expression patterns well-represented in
this cohort and may underrepresent signal that is more prevalent in other
ancestral or geographic populations.

**MMR classifier sample size.** The IHC label column in the TCGA GDC phenotype
file yielded only 58 labeled samples after joining to expression. While class
balance is reasonable (30 / 28), the holdout set is only 12 samples, so point
estimates of F1, ROC-AUC, and PR-AUC carry wide confidence intervals. Attempting
to increase sample size by switching to the `microsatellite_instability` column
yielded 92 samples but with severe class imbalance (81 / 11), producing even
less trustworthy holdout metrics.

**Literature stream.** PubMed abstract retrieval is bounded by the NCBI rate
limit and the query design. A gene with strong but recent or non-English
literature may be under-scored. The LLM literature score, while auditable via
Claude API logs, is sensitive to prompt phrasing; the current rubric was
tuned by hand and is not backed by an inter-rater study against domain experts.

**Pathway stream.** g:Profiler enrichment is driven by gene-set annotations
that are themselves incomplete and biased toward well-studied pathways. Genes
with novel or under-characterized biological roles are systematically
disadvantaged in the pathway-score component.

**Single-snapshot validation.** Validation is against a fixed panel of 80 known
biomarkers. As new markers are added to FDA guidelines, COSMIC CGC, or clinical
literature, the panel — and therefore the reported enrichment ratio — may drift.
The pipeline re-run is cheap, but the validation dashboard numbers in this
report are point-in-time.

## Future Work

1. **Expand MMR sample size.** Replace the GDC IHC / MSI columns with the
   Liu et al. 2018 Pan-Cancer MSI calls (≈550 labeled COAD samples on
   Synapse), which should materially improve classifier reliability and allow
   external validation against an independently curated label set.

2. **Ablate LLM prompt variants.** Run an ablation study in which the LLM
   literature-scoring rubric is systematically varied (weighting of direct
   evidence vs. mechanistic vs. quality) and measure downstream impact on
   enrichment ratio. This quantifies the sensitivity of the final ranking to
   a human-design choice that is currently left implicit.

3. **Multi-cohort validation.** Port the pipeline to a second CRC cohort (e.g.,
   CPTAC or GTEx-normal comparison) and measure how stably the top-ranked
   candidates reproduce across cohorts. Stability across independent data is
   a stronger reliability signal than the current single-cohort enrichment.

4. **Clinical-outcome labels.** Where patient survival and treatment-response
   data are available in TCGA, extend the validation step from "is the gene a
   known biomarker" to "does elevated/depressed expression correlate with
   survival or response." This would connect the ranking to actionable
   clinical endpoints rather than a static annotation panel.

5. **Continuous learning integration.** Wire a scheduled job to refresh the
   PubMed and pathway streams monthly; when the ranking of a previously
   top-scored gene changes by more than a threshold, surface an alert in the
   web application so stakeholders can review what new evidence drove the change.

## Classifier Limits and Suggestions for Improvement

**Current limits.**

- **Random forest** is the best-performing MMR classifier on holdout F1, but its
  cross-validation ROC-AUC (0.52, std 0.17) is barely above chance. With only
  46 training samples after the stratified split, the model has limited capacity
  to learn stable gene-expression signatures and the variance across folds is
  high enough that the ranking of models (RF > LR > XGB) is not statistically
  reliable.
- **Logistic regression** with L1 regularization is the most interpretable
  option and produces a parsimonious coefficient vector, but its CV F1 (0.45,
  std 0.31) suggests that the true class boundary may be too non-linear for a
  linear model at this sample size.
- **XGBoost** performs worst on CV (ROC-AUC 0.47), consistent with tree-boosting
  methods generally requiring larger sample counts to avoid overfitting.

**Suggested improvements.**

1. **Bigger labeled sample.** The single most impactful change is increasing
   the labeled MMR-status set beyond 58 samples — see future work item 1.
2. **Gene-set-level features.** Collapse the ~60K-gene feature space into
   known MMR / DNA-damage-response pathway module scores (MSigDB Hallmark,
   Reactome MMR, KEGG mismatch repair). This would reduce the feature count
   from 60,660 to ~50 and sharply improve the sample-to-feature ratio.
3. **Hyperparameter search and calibration.** The current models use default
   scikit-learn / XGBoost parameters. A nested cross-validation hyperparameter
   search, plus probability calibration (isotonic or Platt), would produce
   more trustworthy class-probability outputs for downstream use.
4. **Uncertainty quantification.** Replace point estimates with bootstrap
   confidence intervals on all holdout metrics, so future readers understand
   the magnitude of the statistical noise in the current numbers rather than
   over-interpreting a single F1 of 0.67.
