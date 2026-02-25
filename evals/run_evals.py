"""CLI orchestrator for the full LLM evaluation suite.

Usage:
    python -m evals.run_evals                          # full eval
    python -m evals.run_evals --skip-consistency       # skip 3x repeat tests
    python -m evals.run_evals --skip-robustness        # skip perturbation tests
    python -m evals.run_evals --skip-consistency --skip-robustness  # judge only

Requires ANTHROPIC_API_KEY environment variable.
Estimated API calls: ~100-150 for full run, ~30 for judge-only.
"""

import argparse
import json
import os
import statistics
import sys
from datetime import datetime, timezone
from pathlib import Path

import anthropic

from evals.test_cases import get_all_test_cases, get_test_cases_by_category, TRUE_POSITIVES, MODERATE_CASES, TRUE_NEGATIVES
from evals.judge import evaluate_test_case, _get_client
from evals.consistency import run_all_consistency_tests
from evals.robustness import run_all_robustness_tests


RESULTS_DIR = Path(__file__).parent / "results"


def run_judge_evaluation(client, test_cases: list) -> list:
    """Run judge evaluation on all test cases."""
    results = []
    total = len(test_cases)

    for i, tc in enumerate(test_cases):
        gene = tc["gene_symbol"]
        print(f"  [{i+1}/{total}] Scoring + judging {gene} ({tc['category']})...")

        try:
            result = evaluate_test_case(client, tc)
            results.append(result)

            score = result["actual_score"]
            expected = result["expected_range"]
            in_range = result["in_expected_range"]
            cal = result["judge"]["score_calibration"]
            grnd = result["judge"]["groundedness"]
            status = "OK" if in_range else "OUT OF RANGE"

            print(f"    score={score} (expected {expected[0]}-{expected[1]}) [{status}]")
            print(f"    judge: calibration={cal}/10, groundedness={grnd}/10")

        except Exception as e:
            print(f"    ERROR: {e}")
            results.append({
                "gene_symbol": gene,
                "category": tc["category"],
                "error": str(e),
            })

    return results


def compute_summary(judge_results: list, consistency_results: dict = None, robustness_results: dict = None) -> dict:
    """Compute aggregate metrics from all evaluation results."""
    # Filter out errors
    valid = [r for r in judge_results if "error" not in r]

    # Accuracy: how many scores fell in expected range?
    in_range_count = sum(1 for r in valid if r.get("in_expected_range"))

    # Judge dimension means
    calibrations = [r["judge"]["score_calibration"] for r in valid if "judge" in r]
    adherences = [r["judge"]["rubric_adherence"] for r in valid if "judge" in r]
    groundedness_scores = [r["judge"]["groundedness"] for r in valid if "judge" in r]
    format_ok = [r["judge"]["format_ok"] for r in valid if "judge" in r]

    # By category breakdown
    by_category = {}
    for r in valid:
        cat = r["category"]
        if cat not in by_category:
            by_category[cat] = {"scores": [], "in_range": [], "calibrations": []}
        by_category[cat]["scores"].append(r["actual_score"])
        by_category[cat]["in_range"].append(r.get("in_expected_range", False))
        by_category[cat]["calibrations"].append(r["judge"]["score_calibration"])

    category_summary = {}
    for cat, data in by_category.items():
        category_summary[cat] = {
            "count": len(data["scores"]),
            "mean_score": round(statistics.mean(data["scores"]), 1),
            "in_range_rate": round(sum(data["in_range"]) / len(data["in_range"]), 2),
            "mean_calibration": round(statistics.mean(data["calibrations"]), 1),
        }

    summary = {
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "model": "claude-sonnet-4-20250514",
        "total_test_cases": len(judge_results),
        "errors": len(judge_results) - len(valid),
        "accuracy": {
            "in_range_count": in_range_count,
            "in_range_rate": round(in_range_count / len(valid), 2) if valid else 0,
            "by_category": category_summary,
        },
        "judge_scores": {
            "score_calibration": {
                "mean": round(statistics.mean(calibrations), 1) if calibrations else 0,
                "min": min(calibrations) if calibrations else 0,
                "max": max(calibrations) if calibrations else 0,
            },
            "rubric_adherence": {
                "mean": round(statistics.mean(adherences), 1) if adherences else 0,
                "min": min(adherences) if adherences else 0,
                "max": max(adherences) if adherences else 0,
            },
            "groundedness": {
                "mean": round(statistics.mean(groundedness_scores), 1) if groundedness_scores else 0,
                "min": min(groundedness_scores) if groundedness_scores else 0,
                "max": max(groundedness_scores) if groundedness_scores else 0,
            },
            "format_compliance": {
                "pass_rate": round(sum(format_ok) / len(format_ok), 2) if format_ok else 0,
            },
        },
    }

    if consistency_results:
        summary["consistency"] = {
            "mean_std": consistency_results["mean_std"],
            "max_std": consistency_results["max_std"],
            "unstable_genes": consistency_results["unstable_genes"],
            "all_stable": consistency_results["all_stable"],
        }

    if robustness_results:
        summary["robustness"] = {
            "pass_count": robustness_results["pass_count"],
            "total_tests": robustness_results["total_tests"],
            "pass_rate": robustness_results["pass_rate"],
        }

    return summary


def generate_markdown_report(summary: dict, judge_results: list, consistency_results: dict = None, robustness_results: dict = None) -> str:
    """Generate human-readable markdown report."""
    lines = [
        "# LLM Evaluation Report",
        f"**Generated:** {summary['timestamp']}",
        f"**Model:** {summary['model']}",
        f"**Test Cases:** {summary['total_test_cases']}",
        "",
        "---",
        "",
        "## Score Accuracy",
        f"- **In expected range:** {summary['accuracy']['in_range_count']}/{summary['total_test_cases'] - summary['errors']} ({summary['accuracy']['in_range_rate']:.0%})",
        "",
        "### By Category",
    ]

    for cat, data in summary["accuracy"]["by_category"].items():
        lines.append(f"- **{cat}**: {data['count']} cases, mean score={data['mean_score']}, "
                     f"in-range rate={data['in_range_rate']:.0%}, calibration={data['mean_calibration']}/10")

    lines.extend([
        "",
        "## Judge Dimensions",
        f"- **Score Calibration:** {summary['judge_scores']['score_calibration']['mean']}/10 "
        f"(range: {summary['judge_scores']['score_calibration']['min']}-{summary['judge_scores']['score_calibration']['max']})",
        f"- **Rubric Adherence:** {summary['judge_scores']['rubric_adherence']['mean']}/10 "
        f"(range: {summary['judge_scores']['rubric_adherence']['min']}-{summary['judge_scores']['rubric_adherence']['max']})",
        f"- **Groundedness:** {summary['judge_scores']['groundedness']['mean']}/10 "
        f"(range: {summary['judge_scores']['groundedness']['min']}-{summary['judge_scores']['groundedness']['max']})",
        f"- **Format Compliance:** {summary['judge_scores']['format_compliance']['pass_rate']:.0%}",
    ])

    if "consistency" in summary:
        c = summary["consistency"]
        lines.extend([
            "",
            "## Consistency",
            f"- **Mean score std dev:** {c['mean_std']}",
            f"- **Max score std dev:** {c['max_std']}",
            f"- **All stable (std ≤ 15):** {'Yes' if c['all_stable'] else 'No'}",
        ])
        if c["unstable_genes"]:
            lines.append(f"- **Unstable genes:** {', '.join(c['unstable_genes'])}")

    if "robustness" in summary:
        r = summary["robustness"]
        lines.extend([
            "",
            "## Robustness",
            f"- **Pass rate:** {r['pass_count']}/{r['total_tests']} ({r['pass_rate']:.0%})",
        ])

    # Detailed results table
    valid = [r for r in judge_results if "error" not in r]
    if valid:
        lines.extend([
            "",
            "## Detailed Results",
            "",
            "| Gene | Category | Score | Expected | In Range | Calibration | Adherence | Grounded | Format |",
            "|------|----------|-------|----------|----------|-------------|-----------|----------|--------|",
        ])
        for r in valid:
            j = r["judge"]
            in_range = "Yes" if r["in_expected_range"] else "**No**"
            fmt = "Yes" if j["format_ok"] else "**No**"
            lines.append(
                f"| {r['gene_symbol']} | {r['category']} | {r['actual_score']} | "
                f"{r['expected_range'][0]}-{r['expected_range'][1]} | {in_range} | "
                f"{j['score_calibration']}/10 | {j['rubric_adherence']}/10 | "
                f"{j['groundedness']}/10 | {fmt} |"
            )

    # Flagged issues
    issues = []
    for r in valid:
        if not r["in_expected_range"]:
            issues.append(f"- **{r['gene_symbol']}**: score {r['actual_score']} outside expected range "
                         f"{r['expected_range'][0]}-{r['expected_range'][1]}")
        if r["judge"]["groundedness"] < 5:
            issues.append(f"- **{r['gene_symbol']}**: low groundedness ({r['judge']['groundedness']}/10) — "
                         "rationale may contain hallucinated claims")
        if not r["judge"]["format_ok"]:
            issues.append(f"- **{r['gene_symbol']}**: format compliance failed")

    if issues:
        lines.extend(["", "## Flagged Issues", ""] + issues)
    else:
        lines.extend(["", "## Flagged Issues", "", "No issues flagged."])

    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description="Run LLM evaluation suite")
    parser.add_argument("--skip-consistency", action="store_true", help="Skip 3x repeat consistency tests")
    parser.add_argument("--skip-robustness", action="store_true", help="Skip perturbation robustness tests")
    args = parser.parse_args()

    # Verify API key
    if not os.getenv("ANTHROPIC_API_KEY"):
        print("ERROR: ANTHROPIC_API_KEY environment variable not set.")
        print("Set it with: export ANTHROPIC_API_KEY=your-key-here")
        sys.exit(1)

    client = _get_client()

    # Create results directory
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    all_test_cases = get_all_test_cases()

    print("=" * 60)
    print("LLM EVALUATION SUITE")
    print("=" * 60)
    print(f"Test cases: {len(all_test_cases)}")
    print(f"Consistency: {'skip' if args.skip_consistency else f'{len(all_test_cases)} genes x 3 runs'}")
    print(f"Robustness: {'skip' if args.skip_robustness else '3 genes x 5 perturbations'}")
    print()

    # Phase 1: Judge evaluation
    print("[Phase 1] Judge Evaluation")
    print("-" * 40)
    judge_results = run_judge_evaluation(client, all_test_cases)
    print()

    # Phase 2: Consistency
    consistency_results = None
    if not args.skip_consistency:
        print("[Phase 2] Consistency Testing (3 runs per gene)")
        print("-" * 40)
        consistency_results = run_all_consistency_tests(client, all_test_cases, n_runs=3)
        print()

    # Phase 3: Robustness
    robustness_results = None
    if not args.skip_robustness:
        print("[Phase 3] Robustness / Perturbation Testing")
        print("-" * 40)
        # Use 3 representative cases: high (KRAS), medium (VEGFA), negative (GAPDH)
        representative = []
        if TRUE_POSITIVES:
            representative.append(TRUE_POSITIVES[0])  # KRAS
        if MODERATE_CASES:
            representative.append(MODERATE_CASES[0])   # VEGFA
        if TRUE_NEGATIVES:
            representative.append(TRUE_NEGATIVES[0])   # GAPDH
        robustness_results = run_all_robustness_tests(client, representative)
        print()

    # Compute summary
    summary = compute_summary(judge_results, consistency_results, robustness_results)

    # Write JSON report
    report_path = RESULTS_DIR / "eval_report.json"
    full_report = {
        "summary": summary,
        "judge_results": judge_results,
    }
    if consistency_results:
        # Remove rationales from JSON to keep file size reasonable
        cr_clean = dict(consistency_results)
        for r in cr_clean.get("results", []):
            r.pop("rationales", None)
        full_report["consistency_results"] = cr_clean
    if robustness_results:
        full_report["robustness_results"] = robustness_results

    with open(report_path, "w") as f:
        json.dump(full_report, f, indent=2, default=str)
    print(f"JSON report: {report_path}")

    # Write markdown report
    md_report = generate_markdown_report(summary, judge_results, consistency_results, robustness_results)
    md_path = RESULTS_DIR / "eval_summary.md"
    with open(md_path, "w") as f:
        f.write(md_report)
    print(f"Markdown report: {md_path}")

    # Console summary
    print()
    print("=" * 60)
    print("SUMMARY")
    print("=" * 60)
    js = summary["judge_scores"]
    acc = summary["accuracy"]
    print(f"  Score accuracy:      {acc['in_range_rate']:.0%} in expected range")
    print(f"  Score calibration:   {js['score_calibration']['mean']}/10")
    print(f"  Rubric adherence:    {js['rubric_adherence']['mean']}/10")
    print(f"  Groundedness:        {js['groundedness']['mean']}/10")
    print(f"  Format compliance:   {js['format_compliance']['pass_rate']:.0%}")

    if consistency_results:
        print(f"  Consistency (std):   {summary['consistency']['mean_std']} mean, {summary['consistency']['max_std']} max")
    if robustness_results:
        print(f"  Robustness:          {summary['robustness']['pass_rate']:.0%} pass rate")

    print()
    print("Done.")


if __name__ == "__main__":
    main()
