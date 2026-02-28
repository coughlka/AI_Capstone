"""Robustness / perturbation testing.

Tests how the LLM scoring responds to systematic input perturbations.
A robust scorer should be resilient to irrelevant changes (abstract order,
duplication) and sensitive to meaningful changes (evidence removal).
"""

import copy
import random

from evals.judge import score_gene


def _run_perturbation(client, gene_symbol: str, abstracts: list[dict]) -> int:
    """Score a gene and return the parsed score."""
    result = score_gene(client, gene_symbol, abstracts)
    return result["parsed_score"]


def test_shuffle_order(client, test_case: dict, seed: int = 42) -> dict:
    """Shuffle abstract order. Score should be within ±10 of original.

    Tests for position bias — the LLM shouldn't weight the first
    abstract more heavily than the last.
    """
    gene = test_case["gene_symbol"]
    abstracts = test_case["abstracts"]

    if len(abstracts) <= 1:
        return {"skipped": True, "reason": "Need ≥2 abstracts to test shuffle"}

    # Original score
    original_score = _run_perturbation(client, gene, abstracts)

    # Shuffled score
    shuffled = copy.deepcopy(abstracts)
    rng = random.Random(seed)
    rng.shuffle(shuffled)
    shuffled_score = _run_perturbation(client, gene, shuffled)

    drift = abs(shuffled_score - original_score)

    return {
        "perturbation": "shuffle_order",
        "gene_symbol": gene,
        "original_score": original_score,
        "perturbed_score": shuffled_score,
        "drift": drift,
        "pass": drift <= 10,
        "threshold": 10,
    }


def test_remove_abstract(client, test_case: dict) -> dict:
    """Remove one abstract. Score should be within ±15 of original.

    Tests evidence sensitivity — removing evidence should cause a
    moderate score decrease, not a catastrophic one.
    """
    gene = test_case["gene_symbol"]
    abstracts = test_case["abstracts"]

    if len(abstracts) <= 1:
        return {"skipped": True, "reason": "Need ≥2 abstracts to test removal"}

    # Original score
    original_score = _run_perturbation(client, gene, abstracts)

    # Remove the first abstract
    reduced = abstracts[1:]
    reduced_score = _run_perturbation(client, gene, reduced)

    drift = abs(reduced_score - original_score)

    return {
        "perturbation": "remove_abstract",
        "gene_symbol": gene,
        "original_score": original_score,
        "perturbed_score": reduced_score,
        "drift": drift,
        "pass": drift <= 15,
        "threshold": 15,
    }


def test_no_abstracts(client, test_case: dict) -> dict:
    """Provide zero abstracts. Score should be ≤10.

    Tests graceful degradation — without evidence, the LLM should
    not rely on prior knowledge to score highly.
    """
    gene = test_case["gene_symbol"]

    score = _run_perturbation(client, gene, [])

    return {
        "perturbation": "no_abstracts",
        "gene_symbol": gene,
        "perturbed_score": score,
        "pass": score <= 10,
        "threshold": 10,
    }


def test_duplicate_abstracts(client, test_case: dict) -> dict:
    """Duplicate all abstracts. Score should be within ±10 of original.

    Tests inflation resistance — seeing the same evidence twice
    shouldn't convince the LLM that evidence is stronger.
    """
    gene = test_case["gene_symbol"]
    abstracts = test_case["abstracts"]

    if not abstracts:
        return {"skipped": True, "reason": "No abstracts to duplicate"}

    # Original score
    original_score = _run_perturbation(client, gene, abstracts)

    # Duplicated
    duplicated = abstracts + copy.deepcopy(abstracts)
    duplicated_score = _run_perturbation(client, gene, duplicated)

    drift = abs(duplicated_score - original_score)

    return {
        "perturbation": "duplicate_abstracts",
        "gene_symbol": gene,
        "original_score": original_score,
        "perturbed_score": duplicated_score,
        "drift": drift,
        "pass": drift <= 10,
        "threshold": 10,
    }


def test_gene_name_swap(client, source_test_case: dict, target_gene: str = "GAPDH") -> dict:
    """Give source gene's abstracts but label as target_gene.

    The most important robustness test: does the LLM score based on the
    abstract content (correct) or the gene name (incorrect)?

    If we give KRAS abstracts but call the gene "GAPDH", the score should
    still be high (70+) because the abstracts contain strong CRC evidence.
    """
    abstracts = source_test_case["abstracts"]

    if not abstracts:
        return {"skipped": True, "reason": "No abstracts for swap test"}

    # Score with original name
    original_name = source_test_case["gene_symbol"]
    original_score = _run_perturbation(client, original_name, abstracts)

    # Score with swapped name
    swapped_score = _run_perturbation(client, target_gene, abstracts)

    drift = abs(swapped_score - original_score)
    # The swapped score should still be high — evidence should dominate name
    followed_evidence = swapped_score >= 60

    return {
        "perturbation": "gene_name_swap",
        "original_gene": original_name,
        "swapped_to": target_gene,
        "original_score": original_score,
        "swapped_score": swapped_score,
        "drift": drift,
        "followed_evidence": followed_evidence,
        "pass": followed_evidence,
    }


def run_all_robustness_tests(client, representative_cases: list) -> dict:
    """Run all perturbation tests on representative genes.

    Args:
        client: Anthropic client.
        representative_cases: List of 3 test cases (high, medium, negative).

    Returns:
        Summary dict with all perturbation results.
    """
    all_results = []

    for tc in representative_cases:
        gene = tc["gene_symbol"]
        print(f"  [robustness] Testing {gene}...")

        # Run all perturbation types
        for test_fn, label in [
            (test_shuffle_order, "shuffle"),
            (test_remove_abstract, "remove"),
            (test_no_abstracts, "no_abstracts"),
            (test_duplicate_abstracts, "duplicate"),
        ]:
            print(f"    {label}...", end=" ")
            result = test_fn(client, tc)
            result["test_type"] = label
            all_results.append(result)
            if result.get("skipped"):
                print(f"skipped ({result['reason']})")
            else:
                status = "PASS" if result["pass"] else "FAIL"
                print(f"{status}")

    # Gene name swap: use the first case (should be high-confidence positive)
    if representative_cases:
        print(f"  [robustness] Gene name swap (KRAS abstracts → GAPDH label)...")
        swap_result = test_gene_name_swap(client, representative_cases[0], target_gene="GAPDH")
        swap_result["test_type"] = "gene_name_swap"
        all_results.append(swap_result)
        if not swap_result.get("skipped"):
            status = "PASS" if swap_result["pass"] else "FAIL"
            print(f"    {status} (original={swap_result['original_score']}, swapped={swap_result['swapped_score']})")

    # Aggregate
    non_skipped = [r for r in all_results if not r.get("skipped")]
    pass_count = sum(1 for r in non_skipped if r.get("pass"))
    total = len(non_skipped)

    return {
        "results": all_results,
        "pass_count": pass_count,
        "total_tests": total,
        "pass_rate": round(pass_count / total, 2) if total > 0 else 0.0,
    }
