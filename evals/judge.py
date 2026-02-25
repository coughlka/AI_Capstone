"""LLM-as-a-judge evaluation module.

Uses a separate Claude call to evaluate the quality of the pipeline's
LLM scoring responses across four dimensions:
- Score calibration: Is the numeric score appropriate for the evidence?
- Rubric adherence: Does the rationale address all 4 scoring criteria?
- Groundedness: Are claims supported by the provided abstracts only?
- Format compliance: Valid JSON with score and rationale?
"""

import json
import os

import anthropic

from src.llm_scoring import SYSTEM_PROMPT, _build_gene_prompt, _parse_llm_response


JUDGE_SYSTEM_PROMPT = """\
You are a strict evaluator assessing an AI biomarker scoring system. You will receive:

1. **SCORING INSTRUCTIONS** — the system prompt given to the AI scorer
2. **INPUT** — the gene name and abstracts the AI received
3. **RESPONSE** — the AI's scoring output (score + rationale)
4. **EXPECTED RANGE** — the expected score range based on external ground truth

Evaluate the response on these dimensions:

### Score Calibration (0-10)
Is the numeric score appropriate given ONLY the abstract content provided?
- 10: Score is well-calibrated within the expected range
- 7-9: Score is reasonable but slightly off from expected range
- 4-6: Score is noticeably miscalibrated (e.g., well-known CRC gene scored <50, or irrelevant gene scored >50)
- 0-3: Score is severely miscalibrated (e.g., strong CRC evidence scored <20, or no CRC evidence scored >80)

### Rubric Adherence (0-10)
The scoring instructions define 4 weighted criteria: Direct CRC evidence (40%), Biomarker potential (30%), \
Mechanistic insight (20%), Evidence quality (10%).
- 10: Rationale explicitly addresses all 4 criteria
- 7-9: Rationale addresses 3 of 4 criteria
- 4-6: Rationale addresses 2 of 4 criteria
- 0-3: Rationale addresses 0-1 criteria or is generic

### Groundedness (0-10)
Every factual claim in the rationale must be traceable to the provided abstracts.
- 10: All claims are directly supported by abstract content
- 7-9: Minor extrapolations but mostly grounded
- 4-6: Some claims not found in abstracts (possible hallucination from training data)
- 0-3: Major claims are unsupported by the provided abstracts

### Format Compliance (true/false)
Did the response contain valid JSON with "score" (integer 0-100) and "rationale" (string)?

Return ONLY a JSON object:
{"score_calibration": <int 0-10>, "rubric_adherence": <int 0-10>, "groundedness": <int 0-10>, \
"format_ok": <bool>, "reasoning": "<brief explanation of your evaluation>"}
"""


def _get_client():
    """Get Anthropic client."""
    api_key = os.getenv("ANTHROPIC_API_KEY")
    if not api_key:
        raise RuntimeError("ANTHROPIC_API_KEY environment variable not set")
    return anthropic.Anthropic(api_key=api_key)


def score_gene(client, gene_symbol: str, abstracts: list[dict]) -> dict:
    """Run the pipeline's LLM scoring on a gene and return raw response.

    Args:
        client: Anthropic client instance.
        gene_symbol: Gene symbol to score.
        abstracts: List of {"title": str, "snippet": str} dicts.

    Returns:
        Dict with "raw_text", "parsed_score", "parsed_rationale".
    """
    prompt = _build_gene_prompt(gene_symbol, abstracts)

    response = client.messages.create(
        model="claude-sonnet-4-20250514",
        max_tokens=256,
        system=SYSTEM_PROMPT,
        messages=[{"role": "user", "content": prompt}],
    )

    raw_text = response.content[0].text
    parsed = _parse_llm_response(raw_text)

    return {
        "raw_text": raw_text,
        "parsed_score": parsed["score"],
        "parsed_rationale": parsed["rationale"],
    }


def judge_response(
    client,
    gene_symbol: str,
    abstracts: list[dict],
    scorer_response: dict,
    expected_range: tuple[int, int],
) -> dict:
    """Have the judge evaluate a scorer's response.

    Args:
        client: Anthropic client instance.
        gene_symbol: Gene that was scored.
        abstracts: Abstracts that were provided to the scorer.
        scorer_response: Dict with "raw_text", "parsed_score", "parsed_rationale".
        expected_range: (min_score, max_score) expected for this gene.

    Returns:
        Dict with judge scores and reasoning.
    """
    # Build the input the scorer received
    scorer_input = _build_gene_prompt(gene_symbol, abstracts)

    user_prompt = f"""## SCORING INSTRUCTIONS (given to the AI scorer)
{SYSTEM_PROMPT}

## INPUT (gene and abstracts provided to the scorer)
{scorer_input}

## RESPONSE (the scorer's output)
{scorer_response['raw_text']}

## EXPECTED RANGE
Based on external ground truth, the expected score for {gene_symbol} is {expected_range[0]}-{expected_range[1]}.
The scorer gave: {scorer_response['parsed_score']}.

Evaluate this response strictly."""

    response = client.messages.create(
        model="claude-sonnet-4-20250514",
        max_tokens=512,
        system=JUDGE_SYSTEM_PROMPT,
        messages=[{"role": "user", "content": user_prompt}],
    )

    raw_judge = response.content[0].text

    # Parse judge response
    try:
        start = raw_judge.find("{")
        end = raw_judge.rfind("}") + 1
        if start >= 0 and end > start:
            result = json.loads(raw_judge[start:end])
            return {
                "score_calibration": int(result.get("score_calibration", 0)),
                "rubric_adherence": int(result.get("rubric_adherence", 0)),
                "groundedness": int(result.get("groundedness", 0)),
                "format_ok": bool(result.get("format_ok", False)),
                "reasoning": result.get("reasoning", ""),
                "raw_judge_response": raw_judge,
            }
    except (json.JSONDecodeError, ValueError):
        pass

    return {
        "score_calibration": 0,
        "rubric_adherence": 0,
        "groundedness": 0,
        "format_ok": False,
        "reasoning": f"Failed to parse judge response: {raw_judge[:200]}",
        "raw_judge_response": raw_judge,
    }


def evaluate_test_case(client, test_case: dict) -> dict:
    """Run full evaluation on a single test case: score + judge.

    Args:
        client: Anthropic client instance.
        test_case: Dict from test_cases.py with gene_symbol, abstracts, expected_score_range, etc.

    Returns:
        Complete evaluation result dict.
    """
    gene = test_case["gene_symbol"]
    abstracts = test_case["abstracts"]
    expected_range = test_case["expected_score_range"]

    # Step 1: Run the scorer
    scorer_result = score_gene(client, gene, abstracts)

    # Step 2: Run the judge
    judge_result = judge_response(client, gene, abstracts, scorer_result, expected_range)

    # Step 3: Check if score falls in expected range
    score = scorer_result["parsed_score"]
    in_expected_range = expected_range[0] <= score <= expected_range[1]

    return {
        "gene_symbol": gene,
        "category": test_case["category"],
        "expected_range": expected_range,
        "actual_score": score,
        "in_expected_range": in_expected_range,
        "rationale": scorer_result["parsed_rationale"],
        "judge": judge_result,
    }
