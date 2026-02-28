"""Consistency testing — same input should produce similar scores across runs.

For each test case, we run the scoring prompt N times with identical input
and measure score variance. High variance means the pipeline's rankings
are non-deterministic and unreliable.
"""

import statistics

from evals.judge import score_gene


def run_consistency_test(client, test_case: dict, n_runs: int = 3) -> dict:
    """Run the same scoring prompt multiple times and measure variance.

    Args:
        client: Anthropic client instance.
        test_case: Dict with gene_symbol, abstracts, etc.
        n_runs: Number of repeated runs.

    Returns:
        Dict with scores, std deviation, and stability assessment.
    """
    gene = test_case["gene_symbol"]
    abstracts = test_case["abstracts"]

    scores = []
    rationales = []

    for i in range(n_runs):
        result = score_gene(client, gene, abstracts)
        scores.append(result["parsed_score"])
        rationales.append(result["parsed_rationale"])

    std_dev = statistics.stdev(scores) if len(scores) > 1 else 0.0
    mean_score = statistics.mean(scores)
    score_range = max(scores) - min(scores)

    # Assess stability
    if std_dev <= 5:
        stability = "stable"
    elif std_dev <= 10:
        stability = "acceptable"
    elif std_dev <= 15:
        stability = "marginal"
    else:
        stability = "unstable"

    # Check rationale keyword overlap across runs
    keyword_sets = []
    for r in rationales:
        words = set(r.lower().split())
        # Filter to meaningful words (>4 chars) to avoid articles/prepositions
        keyword_sets.append({w for w in words if len(w) > 4})

    if len(keyword_sets) >= 2:
        # Jaccard similarity between all pairs
        overlaps = []
        for i in range(len(keyword_sets)):
            for j in range(i + 1, len(keyword_sets)):
                intersection = keyword_sets[i] & keyword_sets[j]
                union = keyword_sets[i] | keyword_sets[j]
                if union:
                    overlaps.append(len(intersection) / len(union))
        mean_overlap = statistics.mean(overlaps) if overlaps else 0.0
    else:
        mean_overlap = 1.0

    return {
        "gene_symbol": gene,
        "category": test_case["category"],
        "scores": scores,
        "mean_score": round(mean_score, 1),
        "std_dev": round(std_dev, 1),
        "score_range": score_range,
        "stability": stability,
        "rationale_overlap": round(mean_overlap, 3),
        "rationales": rationales,
    }


def run_all_consistency_tests(client, test_cases: list, n_runs: int = 3) -> dict:
    """Run consistency tests for all test cases.

    Returns:
        Summary dict with per-gene results and aggregate metrics.
    """
    results = []
    for tc in test_cases:
        print(f"  [consistency] Testing {tc['gene_symbol']} ({n_runs} runs)...")
        result = run_consistency_test(client, tc, n_runs=n_runs)
        results.append(result)
        print(f"    scores={result['scores']}, std={result['std_dev']}, stability={result['stability']}")

    # Aggregate
    all_stds = [r["std_dev"] for r in results]
    unstable = [r["gene_symbol"] for r in results if r["stability"] == "unstable"]

    return {
        "results": results,
        "mean_std": round(statistics.mean(all_stds), 1) if all_stds else 0.0,
        "max_std": round(max(all_stds), 1) if all_stds else 0.0,
        "unstable_genes": unstable,
        "all_stable": len(unstable) == 0,
    }
