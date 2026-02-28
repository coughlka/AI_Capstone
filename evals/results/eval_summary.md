# LLM Evaluation Report
**Generated:** 2026-02-28T17:28:12.315661+00:00
**Model:** claude-sonnet-4-20250514
**Test Cases:** 10

---

## Score Accuracy
- **In expected range:** 9/10 (90%)

### By Category
- **true_positive**: 3 cases, mean score=87.3, in-range rate=100%, calibration=9.7/10
- **moderate**: 2 cases, mean score=65, in-range rate=100%, calibration=8.5/10
- **true_negative**: 3 cases, mean score=23.3, in-range rate=67%, calibration=7.3/10
- **edge_no_abstracts**: 1 cases, mean score=0, in-range rate=100%, calibration=10/10
- **edge_conflicting**: 1 cases, mean score=45, in-range rate=100%, calibration=9/10

## Judge Dimensions
- **Score Calibration:** 8.7/10 (range: 2-10)
- **Rubric Adherence:** 8.4/10 (range: 0-10)
- **Groundedness:** 9.8/10 (range: 9-10)
- **Format Compliance:** 100%

## Consistency
- **Mean score std dev:** 0.0
- **Max score std dev:** 0.0
- **All stable (std ≤ 15):** Yes

## Robustness
- **Pass rate:** 12/13 (92%)

## Detailed Results

| Gene | Category | Score | Expected | In Range | Calibration | Adherence | Grounded | Format |
|------|----------|-------|----------|----------|-------------|-----------|----------|--------|
| LY6G6D | true_positive | 92 | 70-100 | Yes | 10/10 | 10/10 | 10/10 | Yes |
| SCARA5 | true_positive | 85 | 70-100 | Yes | 9/10 | 10/10 | 10/10 | Yes |
| DPYD | true_positive | 85 | 70-100 | Yes | 10/10 | 10/10 | 10/10 | Yes |
| GSTM5 | moderate | 65 | 30-80 | Yes | 8/10 | 9/10 | 10/10 | Yes |
| HDAC9 | moderate | 65 | 30-80 | Yes | 9/10 | 9/10 | 9/10 | Yes |
| CNTN2 | true_negative | 25 | 0-25 | Yes | 10/10 | 8/10 | 10/10 | Yes |
| UCN3 | true_negative | 45 | 0-25 | **No** | 2/10 | 8/10 | 9/10 | Yes |
| FICTIONAL_GENE_X | true_negative | 0 | 0-15 | Yes | 10/10 | 10/10 | 10/10 | Yes |
| EMPTY_GENE | edge_no_abstracts | 0 | 0-10 | Yes | 10/10 | 0/10 | 10/10 | Yes |
| CONFLICTED_GENE | edge_conflicting | 45 | 20-60 | Yes | 9/10 | 10/10 | 10/10 | Yes |

## Flagged Issues

- **UCN3**: score 45 outside expected range 0-25