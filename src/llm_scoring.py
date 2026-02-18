"""LLM-based literature relevance scoring.

Uses Claude to analyze PubMed abstracts for each candidate gene and
score their relevance as colorectal cancer biomarkers.
"""

import json
import os
import time
from typing import Optional

import anthropic
import pandas as pd
from dotenv import load_dotenv

from src.utils import load_config, ensure_dirs, read_csv, write_csv

# Load .env from web-app directory (where the API key lives)
_PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
load_dotenv(os.path.join(_PROJECT_ROOT, "web-app", ".env"))

SYSTEM_PROMPT = """\
You are an expert biomedical researcher evaluating genes as potential \
colorectal cancer (CRC) biomarkers.

For each gene you will receive its name and PubMed abstract snippets. \
Score the gene's relevance as a CRC biomarker on a 0-100 scale based on:

- **Direct CRC evidence** (40%): Do the abstracts describe the gene's role \
specifically in colorectal/colon/rectal cancer?
- **Biomarker potential** (30%): Is there evidence the gene could serve as a \
diagnostic, prognostic, or therapeutic biomarker (e.g., differential expression, \
survival association, drug target)?
- **Mechanistic insight** (20%): Do the abstracts describe a biological mechanism \
linking the gene to CRC (e.g., Wnt signaling, EMT, immune evasion)?
- **Evidence quality** (10%): Are the abstracts from recent, peer-reviewed studies \
with clear methodology?

Return ONLY a JSON object with this exact format:
{"score": <integer 0-100>, "rationale": "<1-2 sentence explanation>"}
"""


def _build_gene_prompt(gene_symbol: str, abstracts: list[dict]) -> str:
    """Build the user prompt for a single gene."""
    parts = [f"Gene: {gene_symbol}\n\nAbstracts:"]
    for i, ab in enumerate(abstracts, 1):
        parts.append(f"\n[{i}] {ab['title']}\n{ab['snippet']}")
    return "\n".join(parts)


def _parse_llm_response(text: str) -> dict:
    """Extract score and rationale from LLM response."""
    text = text.strip()
    # Try to find JSON in the response
    start = text.find("{")
    end = text.rfind("}") + 1
    if start >= 0 and end > start:
        try:
            result = json.loads(text[start:end])
            score = int(result.get("score", 0))
            score = max(0, min(100, score))
            return {"score": score, "rationale": result.get("rationale", "")}
        except (json.JSONDecodeError, ValueError):
            pass
    return {"score": 0, "rationale": "Failed to parse LLM response"}


def run_llm_scoring(config_path: str) -> str:
    """Run LLM-based literature scoring for all candidate genes.

    Reads lit_evidence.csv, groups abstracts by gene, sends each gene's
    abstracts to Claude for relevance scoring, and writes llm_scores.csv.

    Args:
        config_path: Path to the configuration YAML file.

    Returns:
        Path to the output CSV file.
    """
    print("[llm_scoring] Loading configuration...")
    config = load_config(config_path)
    ensure_dirs(config)

    outputs_dir = config["paths"]["outputs_dir"]
    lit_path = os.path.join(outputs_dir, "lit_evidence.csv")
    omics_path = os.path.join(outputs_dir, "omics_evidence.csv")
    output_path = os.path.join(outputs_dir, "llm_scores.csv")

    # Load literature evidence
    if not os.path.exists(lit_path):
        print(f"[llm_scoring] No literature evidence found at {lit_path}")
        empty = pd.DataFrame(columns=["gene", "gene_symbol", "llm_score", "rationale"])
        write_csv(empty, output_path)
        return output_path

    lit_df = read_csv(lit_path)
    if lit_df.empty:
        print("[llm_scoring] Literature evidence is empty.")
        empty = pd.DataFrame(columns=["gene", "gene_symbol", "llm_score", "rationale"])
        write_csv(empty, output_path)
        return output_path

    # Load gene symbol mapping from omics
    symbol_map = {}
    if os.path.exists(omics_path):
        omics_df = read_csv(omics_path)
        if "gene_symbol" in omics_df.columns:
            symbol_map = dict(zip(omics_df["gene"], omics_df["gene_symbol"]))

    # Group abstracts by gene
    grouped = {}
    for _, row in lit_df.iterrows():
        gene = row["gene"]
        if gene not in grouped:
            grouped[gene] = []
        grouped[gene].append({
            "title": str(row.get("title", "")),
            "snippet": str(row.get("snippet", "")),
        })

    # Initialize Anthropic client
    api_key = os.getenv("ANTHROPIC_API_KEY")
    if not api_key:
        raise RuntimeError(
            "ANTHROPIC_API_KEY not set. Add it to web-app/.env or set as environment variable."
        )
    client = anthropic.Anthropic(api_key=api_key)

    print(f"[llm_scoring] Scoring {len(grouped)} genes with Claude...")
    rows = []
    total = len(grouped)

    for i, (gene, abstracts) in enumerate(grouped.items()):
        symbol = symbol_map.get(gene, gene)
        prompt = _build_gene_prompt(symbol, abstracts)

        try:
            response = client.messages.create(
                model="claude-sonnet-4-20250514",
                max_tokens=256,
                system=SYSTEM_PROMPT,
                messages=[{"role": "user", "content": prompt}],
            )
            result = _parse_llm_response(response.content[0].text)
        except anthropic.RateLimitError:
            print(f"[llm_scoring] Rate limited at gene {i+1}. Waiting 60s...")
            time.sleep(60)
            try:
                response = client.messages.create(
                    model="claude-sonnet-4-20250514",
                    max_tokens=256,
                    system=SYSTEM_PROMPT,
                    messages=[{"role": "user", "content": prompt}],
                )
                result = _parse_llm_response(response.content[0].text)
            except Exception as e:
                print(f"[llm_scoring] [{i+1}/{total}] {symbol}: retry failed ({e})")
                result = {"score": 0, "rationale": f"API error: {e}"}
        except Exception as e:
            print(f"[llm_scoring] [{i+1}/{total}] {symbol}: failed ({e})")
            result = {"score": 0, "rationale": f"Error: {e}"}

        rows.append({
            "gene": gene,
            "gene_symbol": symbol,
            "llm_score": result["score"],
            "rationale": result["rationale"],
        })

        if (i + 1) % 25 == 0 or i == 0:
            print(f"[llm_scoring] [{i+1}/{total}] {symbol}: score={result['score']}")

    scores_df = pd.DataFrame(rows)
    print(f"\n[llm_scoring] Score distribution:")
    print(f"  Mean: {scores_df['llm_score'].mean():.1f}")
    print(f"  Median: {scores_df['llm_score'].median():.1f}")
    print(f"  Min: {scores_df['llm_score'].min()}, Max: {scores_df['llm_score'].max()}")

    write_csv(scores_df, output_path)
    print(f"[llm_scoring] Wrote {len(scores_df)} scores to {output_path}")
    return output_path
