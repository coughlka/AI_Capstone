"""Gold standard test cases for LLM-as-a-judge evaluation.

Each test case provides a gene symbol, realistic abstracts, expected score range,
and key terms the rationale should mention. Ground truth comes from published
CRC biomarker literature, not from our own pipeline outputs.
"""

# ---------------------------------------------------------------------------
# True Positives — well-established CRC biomarkers (expect 70-100)
# ---------------------------------------------------------------------------

TRUE_POSITIVES = [
    {
        "gene_symbol": "KRAS",
        "category": "true_positive",
        "expected_score_range": (70, 100),
        "rationale_should_mention": ["colorectal", "mutation", "EGFR", "resistance"],
        "abstracts": [
            {
                "title": "KRAS mutations predict resistance to anti-EGFR therapy in metastatic colorectal cancer",
                "snippet": "KRAS mutations are found in approximately 40% of colorectal cancers and are a well-established predictive biomarker for resistance to cetuximab and panitumumab. Patients with KRAS-mutant tumors show significantly worse progression-free survival when treated with anti-EGFR monoclonal antibodies.",
            },
            {
                "title": "Prognostic significance of KRAS mutation status in colorectal cancer: a meta-analysis",
                "snippet": "This meta-analysis of 27 studies comprising 4,853 patients demonstrates that KRAS mutations are associated with poorer overall survival in CRC (HR=1.45, 95% CI: 1.29-1.62). The association was strongest for codon 12 mutations in stage III-IV disease.",
            },
            {
                "title": "KRAS-driven MAPK signaling in colorectal tumorigenesis",
                "snippet": "Oncogenic KRAS activates the RAF-MEK-ERK signaling cascade, promoting cell proliferation, survival, and migration in colorectal epithelial cells. Constitutive MAPK activation is a hallmark of KRAS-mutant CRC and drives therapeutic resistance through feedback reactivation loops.",
            },
        ],
    },
    {
        "gene_symbol": "APC",
        "category": "true_positive",
        "expected_score_range": (70, 100),
        "rationale_should_mention": ["Wnt", "truncat", "colorectal", "tumor suppressor"],
        "abstracts": [
            {
                "title": "APC mutations initiate colorectal tumorigenesis through Wnt pathway activation",
                "snippet": "Truncating mutations in the APC tumor suppressor gene are the earliest and most frequent genetic event in colorectal cancer, occurring in over 80% of sporadic CRC cases. Loss of APC function leads to constitutive activation of the Wnt/beta-catenin signaling pathway, driving hyperproliferation of colonic epithelial cells.",
            },
            {
                "title": "Familial adenomatous polyposis and APC: from bench to bedside",
                "snippet": "Germline APC mutations cause familial adenomatous polyposis (FAP), an autosomal dominant syndrome characterized by hundreds to thousands of colorectal adenomas with near-100% progression to CRC by age 40. APC genotype-phenotype correlations guide clinical management including timing of prophylactic colectomy.",
            },
            {
                "title": "APC as a diagnostic biomarker: methylation-based detection in stool DNA",
                "snippet": "Hypermethylation of the APC promoter is detectable in stool DNA samples from CRC patients with 72% sensitivity and 94% specificity. Combined with other methylation markers, APC-based stool DNA testing offers a non-invasive screening alternative to colonoscopy for early CRC detection.",
            },
        ],
    },
    {
        "gene_symbol": "TP53",
        "category": "true_positive",
        "expected_score_range": (70, 100),
        "rationale_should_mention": ["tumor suppressor", "colorectal", "p53", "prognostic"],
        "abstracts": [
            {
                "title": "TP53 mutations in colorectal cancer: prognostic implications and therapeutic targets",
                "snippet": "TP53 mutations occur in approximately 60% of colorectal cancers and are associated with the adenoma-to-carcinoma progression sequence. Loss of p53 function impairs DNA damage response, apoptosis, and cell cycle arrest, accelerating genomic instability and tumor evolution in CRC.",
            },
            {
                "title": "Prognostic value of TP53 mutation status in stage II-III colorectal cancer",
                "snippet": "In a cohort of 1,248 stage II-III CRC patients, TP53 mutations were independently associated with reduced disease-free survival (HR=1.38, p=0.003). Gain-of-function TP53 mutations showed particularly poor prognosis compared to loss-of-function mutations.",
            },
        ],
    },
]

# ---------------------------------------------------------------------------
# Moderate Cases — relevant but not CRC-specific (expect 30-70)
# ---------------------------------------------------------------------------

MODERATE_CASES = [
    {
        "gene_symbol": "VEGFA",
        "category": "moderate",
        "expected_score_range": (30, 70),
        "rationale_should_mention": ["angiogenesis", "bevacizumab"],
        "abstracts": [
            {
                "title": "VEGFA and tumor angiogenesis: implications for anti-angiogenic therapy",
                "snippet": "Vascular endothelial growth factor A (VEGFA) is the primary mediator of tumor angiogenesis across multiple cancer types. Bevacizumab, a monoclonal antibody targeting VEGFA, has shown clinical benefit in metastatic CRC when combined with chemotherapy, though VEGFA expression alone is not a reliable predictive biomarker.",
            },
            {
                "title": "Circulating VEGF levels in cancer patients: a pan-cancer analysis",
                "snippet": "Elevated serum VEGF levels are observed across diverse malignancies including breast, lung, renal, and colorectal cancers. While VEGF correlates with tumor vascularity and stage, its lack of cancer-type specificity limits its utility as a standalone diagnostic biomarker for any individual cancer type.",
            },
        ],
    },
    {
        "gene_symbol": "MKI67",
        "category": "moderate",
        "expected_score_range": (30, 70),
        "rationale_should_mention": ["proliferation", "Ki-67"],
        "abstracts": [
            {
                "title": "Ki-67 proliferation index in colorectal cancer: systematic review",
                "snippet": "Ki-67 (MKI67) is a widely used proliferation marker across cancer types. In CRC, high Ki-67 expression has been variably associated with prognosis, with some studies showing correlation with advanced stage and others finding no independent prognostic value after multivariate adjustment.",
            },
            {
                "title": "Immunohistochemical markers in colorectal cancer: diagnostic and prognostic utility",
                "snippet": "MKI67 labeling index is routinely assessed in colorectal tumor pathology. While high proliferative activity indexed by Ki-67 indicates aggressive tumor biology, its prognostic utility in CRC remains controversial compared to well-established markers like microsatellite instability status.",
            },
        ],
    },
]

# ---------------------------------------------------------------------------
# True Negatives — should score low (expect 0-20)
# ---------------------------------------------------------------------------

TRUE_NEGATIVES = [
    {
        "gene_symbol": "GAPDH",
        "category": "true_negative",
        "expected_score_range": (0, 20),
        "rationale_should_mention": ["housekeeping", "reference", "normalization"],
        "abstracts": [
            {
                "title": "Selection of optimal reference genes for qRT-PCR normalization",
                "snippet": "GAPDH (glyceraldehyde-3-phosphate dehydrogenase) is among the most commonly used reference genes for quantitative RT-PCR normalization. While widely employed as a housekeeping control, GAPDH expression can vary under hypoxic conditions and in certain metabolic states, necessitating validation for each experimental context.",
            },
            {
                "title": "Stability assessment of reference genes across human tissues",
                "snippet": "Comprehensive analysis of reference gene stability across 72 human tissue types confirmed GAPDH as a reliable endogenous control for gene expression studies. No significant differential expression of GAPDH was observed between matched tumor and normal samples in the majority of tissue types examined.",
            },
        ],
    },
    {
        "gene_symbol": "BRCA1",
        "category": "true_negative",
        "expected_score_range": (0, 25),
        "rationale_should_mention": ["breast", "ovarian"],
        "abstracts": [
            {
                "title": "BRCA1 mutations and hereditary breast and ovarian cancer syndrome",
                "snippet": "Germline BRCA1 mutations confer a lifetime breast cancer risk of 55-72% and ovarian cancer risk of 39-44%. BRCA1 functions in homologous recombination DNA repair, and its loss leads to genomic instability primarily in breast and ovarian epithelial tissues.",
            },
            {
                "title": "PARP inhibitors in BRCA1-mutant cancers: from ovarian to pan-cancer applications",
                "snippet": "PARP inhibitor therapy has transformed the treatment of BRCA1-mutant ovarian and breast cancers. While BRCA1 mutations occur rarely in colorectal cancer (<1%), recent studies have explored PARP inhibitor sensitivity in microsatellite-stable CRC with homologous recombination deficiency from other causes.",
            },
        ],
    },
    {
        "gene_symbol": "FICTIONAL_GENE_X",
        "category": "true_negative",
        "expected_score_range": (0, 15),
        "rationale_should_mention": [],
        "abstracts": [
            {
                "title": "Characterization of zinc finger protein interactions in Drosophila wing development",
                "snippet": "We report the identification of novel zinc finger protein complexes regulating wing disc patterning in Drosophila melanogaster. These transcription factors coordinate hedgehog and decapentaplegic signaling during larval imaginal disc morphogenesis, with no orthologs detected in mammalian genomes.",
            },
        ],
    },
]

# ---------------------------------------------------------------------------
# Edge Cases — test boundary conditions
# ---------------------------------------------------------------------------

EDGE_CASES = [
    {
        "gene_symbol": "EMPTY_GENE",
        "category": "edge_no_abstracts",
        "expected_score_range": (0, 10),
        "rationale_should_mention": [],
        "abstracts": [],
    },
    {
        "gene_symbol": "CONFLICTED_GENE",
        "category": "edge_conflicting",
        "expected_score_range": (20, 60),
        "rationale_should_mention": ["conflict", "inconsistent"],
        "abstracts": [
            {
                "title": "CONFLICTED_GENE as a novel CRC biomarker: a pilot study",
                "snippet": "We report that CONFLICTED_GENE is significantly overexpressed in colorectal tumor tissue compared to adjacent normal mucosa (p<0.001). High expression was associated with advanced stage and poor prognosis, suggesting a role as a diagnostic and prognostic CRC biomarker.",
            },
            {
                "title": "Re-evaluation of CONFLICTED_GENE in colorectal cancer: no independent prognostic value",
                "snippet": "In our large multicenter validation study (n=2,450), we could not replicate the previously reported association between CONFLICTED_GENE expression and CRC outcomes. After adjustment for established prognostic factors, CONFLICTED_GENE showed no independent predictive value (HR=1.02, p=0.89).",
            },
        ],
    },
]


def get_all_test_cases():
    """Return all test cases as a flat list."""
    return TRUE_POSITIVES + MODERATE_CASES + TRUE_NEGATIVES + EDGE_CASES


def get_test_cases_by_category():
    """Return test cases organized by category."""
    return {
        "true_positive": TRUE_POSITIVES,
        "moderate": MODERATE_CASES,
        "true_negative": TRUE_NEGATIVES,
        "edge_cases": EDGE_CASES,
    }
