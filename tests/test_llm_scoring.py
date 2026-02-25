"""LLM response parsing and prompt construction tests.

Parsing failures silently assign score=0 to genes that might deserve high scores,
corrupting the entire ranking. Prompt construction errors mean the LLM scores
on incomplete or wrong evidence.
"""

import pytest

from src.llm_scoring import _parse_llm_response, _build_gene_prompt, SYSTEM_PROMPT


# ---------------------------------------------------------------------------
# Response Parsing
# ---------------------------------------------------------------------------

class TestParseLLMResponse:
    """Every parsing failure silently assigns score=0, potentially burying
    a real biomarker at the bottom of the ranking."""

    def test_valid_json(self):
        """Baseline: well-formed response should parse correctly."""
        result = _parse_llm_response('{"score": 75, "rationale": "Strong CRC evidence."}')
        assert result["score"] == 75
        assert result["rationale"] == "Strong CRC evidence."

    def test_json_with_surrounding_text(self):
        """Claude sometimes adds preamble before JSON. If we don't extract
        the JSON from within the text, we get score=0 for a valid response."""
        text = 'Based on my analysis:\n{"score": 85, "rationale": "Direct CRC role."}\nHope this helps!'
        result = _parse_llm_response(text)
        assert result["score"] == 85

    def test_json_with_markdown_code_block(self):
        """Claude frequently wraps JSON in ```json blocks. The current parser
        uses find('{') which should handle this, but verify explicitly."""
        text = '```json\n{"score": 90, "rationale": "Key biomarker."}\n```'
        result = _parse_llm_response(text)
        assert result["score"] == 90

    def test_score_clamped_above_100(self):
        """Score > 100 would break min-max normalization assumptions."""
        result = _parse_llm_response('{"score": 150, "rationale": "Very strong."}')
        assert result["score"] == 100

    def test_score_clamped_below_0(self):
        """Score < 0 is meaningless and would produce negative normalized values."""
        result = _parse_llm_response('{"score": -10, "rationale": "No evidence."}')
        assert result["score"] == 0

    def test_missing_score_key(self):
        """If LLM omits the score key, we should get 0, not crash."""
        result = _parse_llm_response('{"rationale": "Some text"}')
        assert result["score"] == 0

    def test_completely_invalid_json(self):
        """LLM returns prose. Should gracefully return score=0, not crash
        and kill the pipeline for the remaining 400+ genes."""
        result = _parse_llm_response("This gene appears to be relevant to CRC.")
        assert result["score"] == 0
        assert "Failed to parse" in result["rationale"]

    def test_empty_string(self):
        """Empty response (API timeout/error). Must not crash."""
        result = _parse_llm_response("")
        assert result["score"] == 0

    def test_float_score_converted_to_int(self):
        """LLM might return float like 72.5. Should be converted to int."""
        result = _parse_llm_response('{"score": 72.5, "rationale": "Moderate."}')
        assert isinstance(result["score"], int)
        assert result["score"] == 72


# ---------------------------------------------------------------------------
# Prompt Construction
# ---------------------------------------------------------------------------

class TestBuildGenePrompt:
    """Prompt construction errors mean the LLM scores on incomplete or
    incorrect evidence."""

    def test_gene_symbol_in_prompt(self):
        """If the gene symbol is missing, Claude can't identify what it's scoring."""
        prompt = _build_gene_prompt("KRAS", [{"title": "T1", "snippet": "S1"}])
        assert "KRAS" in prompt

    def test_all_abstracts_included(self):
        """If abstracts are silently dropped, the LLM scores on partial evidence."""
        abstracts = [
            {"title": f"Title {i}", "snippet": f"Snippet {i}"}
            for i in range(5)
        ]
        prompt = _build_gene_prompt("APC", abstracts)
        for i in range(5):
            assert f"Title {i}" in prompt
            assert f"Snippet {i}" in prompt

    def test_abstracts_are_numbered(self):
        """Numbered abstracts help the LLM reference specific evidence in rationale."""
        abstracts = [
            {"title": "First", "snippet": "S1"},
            {"title": "Second", "snippet": "S2"},
            {"title": "Third", "snippet": "S3"},
        ]
        prompt = _build_gene_prompt("TP53", abstracts)
        assert "[1]" in prompt
        assert "[2]" in prompt
        assert "[3]" in prompt

    def test_empty_abstracts_list(self):
        """No abstracts available for this gene. Prompt should still be valid
        and include the gene symbol — the LLM should return a low score."""
        prompt = _build_gene_prompt("RARE_GENE", [])
        assert "RARE_GENE" in prompt
        assert "Abstracts:" in prompt


# ---------------------------------------------------------------------------
# System Prompt Validation
# ---------------------------------------------------------------------------

class TestSystemPrompt:
    """The system prompt defines the scoring contract. Any drift from the
    intended criteria invalidates all LLM scores."""

    def test_all_four_criteria_present(self):
        """If a scoring criterion is missing from the prompt, the LLM won't
        consider it, biasing all scores toward the remaining criteria."""
        assert "Direct CRC evidence" in SYSTEM_PROMPT
        assert "Biomarker potential" in SYSTEM_PROMPT
        assert "Mechanistic insight" in SYSTEM_PROMPT
        assert "Evidence quality" in SYSTEM_PROMPT

    def test_criteria_weights_present(self):
        """Weights guide the LLM's scoring. Missing weights mean uncontrolled scoring."""
        assert "40%" in SYSTEM_PROMPT
        assert "30%" in SYSTEM_PROMPT
        assert "20%" in SYSTEM_PROMPT
        assert "10%" in SYSTEM_PROMPT

    def test_output_format_specified(self):
        """The LLM must know to return JSON. Without this, responses
        will be free-text and parsing will fail for every gene."""
        assert "JSON" in SYSTEM_PROMPT
        assert '"score"' in SYSTEM_PROMPT
        assert '"rationale"' in SYSTEM_PROMPT

    def test_evidence_only_instruction_present(self):
        """Without this instruction, the LLM uses prior knowledge about gene names
        instead of scoring based on the provided abstracts. This was the root cause
        of the gene name swap test failing (KRAS abstracts labeled GAPDH → scored 15)."""
        assert "ONLY" in SYSTEM_PROMPT
        assert "prior knowledge" in SYSTEM_PROMPT.lower()

    def test_calibration_anchors_present(self):
        """Without calibration anchors, the LLM inflated moderate-evidence genes
        (VEGFA, MKI67) to 75 when they should have been 30-70. The anchors define
        what each score range means concretely."""
        assert "80-100" in SYSTEM_PROMPT
        assert "60-79" in SYSTEM_PROMPT
        assert "40-59" in SYSTEM_PROMPT
        assert "20-39" in SYSTEM_PROMPT
        assert "0-19" in SYSTEM_PROMPT

    def test_no_abstracts_guidance_present(self):
        """Without this, KRAS with zero abstracts scored 95 purely from
        prior knowledge. The prompt must explicitly handle the empty case."""
        # Check that the prompt mentions what to do when no abstracts are provided
        prompt_lower = SYSTEM_PROMPT.lower()
        assert "no abstracts" in prompt_lower or "no crc-relevant content" in prompt_lower
