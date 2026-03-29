"""LLM-based literature relevance scoring.

Uses Claude to analyze PubMed abstracts for each candidate gene and
score their relevance as colorectal cancer biomarkers.
"""

import hashlib
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

Follow these steps IN ORDER:

STEP 1 — CHECK FOR ABSTRACTS: Look at the "Abstracts:" section in the user \
message. If it is empty or contains no abstracts, you MUST immediately return \
{"score": 0, "rationale": "No abstracts provided for evaluation."} and STOP. \
Do NOT use your prior knowledge about the gene to infer a score when no \
abstracts are provided. This applies regardless of how well-known the gene is.

STEP 2 — SCORE BASED ON ABSTRACT CONTENT ONLY: Read the provided abstracts \
and score the CRC biomarker evidence they contain. Base your score entirely \
on what the abstracts say about CRC relevance. Do NOT use prior knowledge \
about the gene.

Score on a 0-100 scale using these weighted criteria:

- **Direct CRC evidence** (40%): Do the abstracts describe a role \
specifically in colorectal/colon/rectal cancer?
- **Biomarker potential** (30%): Is there evidence of diagnostic, \
prognostic, or therapeutic biomarker utility (e.g., differential expression, \
survival association, drug target)?
- **Mechanistic insight** (20%): Do the abstracts describe a biological mechanism \
linking to CRC (e.g., Wnt signaling, EMT, immune evasion)?
- **Evidence quality** (10%): Are the abstracts from recent, peer-reviewed studies \
with clear methodology?

Score calibration guide:
- 80-100: Abstracts present direct, specific evidence of an established \
or validated CRC biomarker with clear diagnostic/prognostic/therapeutic utility.
- 60-79: Strong CRC-specific evidence but biomarker role is emerging, not yet \
validated, or limited to specific contexts.
- 40-59: CRC is mentioned but is not the primary focus, evidence is \
indirect, or findings are pan-cancer rather than CRC-specific.
- 20-39: Tangential CRC relevance only (e.g., mentioned in passing, general \
oncology context without CRC-specific data).
- 0-19: No meaningful CRC biomarker evidence in the provided abstracts.

Return ONLY a JSON object with this exact format:
{"score": <integer 0-100>, "rationale": "<1-2 sentence explanation>"}
"""


def _build_gene_prompt(gene_symbol: str, abstracts: list[dict]) -> str:
    """Build the user prompt for a single gene."""
    parts = [f"Gene: {gene_symbol}\n\nAbstracts:"]
    for i, ab in enumerate(abstracts, 1):
        parts.append(f"\n[{i}] {ab['title']}\n{ab['abstract']}")
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
        # Use 'abstract' column if available, fall back to 'snippet'
        abstract_text = str(row.get("abstract", "") or "")
        if not abstract_text.strip():
            abstract_text = str(row.get("snippet", "") or "")
        grouped[gene].append({
            "title": str(row.get("title", "")),
            "abstract": abstract_text,
        })

    # Build content hashes for each gene to detect changed abstracts
    def _gene_hash(abstracts: list[dict]) -> str:
        content = json.dumps({"abstracts": abstracts}, sort_keys=True)
        return hashlib.md5(content.encode()).hexdigest()

    gene_hashes = {gene: _gene_hash(abs_list) for gene, abs_list in grouped.items()}

    # Load existing cache to avoid re-scoring unchanged genes
    cached_scores = {}
    cache_path = os.path.join(outputs_dir, "llm_scores_cache.json")
    if os.path.exists(cache_path):
        try:
            with open(cache_path) as f:
                cached_scores = json.load(f)
            print(f"[llm_scoring] Loaded cache with {len(cached_scores)} existing scores")
        except (json.JSONDecodeError, IOError):
            cached_scores = {}

    # Determine which genes need (re-)scoring
    needs_scoring = []
    cached_rows = []
    for gene, abstracts in grouped.items():
        cache_entry = cached_scores.get(gene)
        if cache_entry and cache_entry.get("hash") == gene_hashes[gene]:
            cached_rows.append({
                "gene": gene,
                "gene_symbol": symbol_map.get(gene, gene),
                "llm_score": cache_entry["score"],
                "rationale": cache_entry["rationale"],
            })
        else:
            needs_scoring.append((gene, abstracts))

    print(f"[llm_scoring] {len(cached_rows)} genes cached, {len(needs_scoring)} need scoring")

    if not needs_scoring:
        print("[llm_scoring] All genes already scored — skipping API calls")
        scores_df = pd.DataFrame(cached_rows)
    else:
        # Initialize Anthropic client
        api_key = os.getenv("ANTHROPIC_API_KEY")
        if not api_key:
            raise RuntimeError(
                "ANTHROPIC_API_KEY not set. Add it to web-app/.env or set as environment variable."
            )
        client = anthropic.Anthropic(api_key=api_key)

        print(f"[llm_scoring] Scoring {len(needs_scoring)} genes with Claude...")
        new_rows = []
        total = len(needs_scoring)

        for i, (gene, abstracts) in enumerate(needs_scoring):
            symbol = symbol_map.get(gene, gene)
            prompt = _build_gene_prompt(symbol, abstracts)

            try:
                response = client.messages.create(
                    model="claude-sonnet-4-20250514",
                    max_tokens=512,
                    temperature=0,
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
                        max_tokens=512,
                        temperature=0,
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

            new_rows.append({
                "gene": gene,
                "gene_symbol": symbol,
                "llm_score": result["score"],
                "rationale": result["rationale"],
            })

            # Update cache immediately so progress is saved if interrupted
            cached_scores[gene] = {
                "hash": gene_hashes[gene],
                "score": result["score"],
                "rationale": result["rationale"],
            }

            if (i + 1) % 25 == 0 or i == 0:
                print(f"[llm_scoring] [{i+1}/{total}] {symbol}: score={result['score']}")

        # Save updated cache
        try:
            with open(cache_path, "w") as f:
                json.dump(cached_scores, f)
            print(f"[llm_scoring] Saved cache with {len(cached_scores)} entries")
        except IOError as e:
            print(f"[llm_scoring] Warning: failed to save cache ({e})")

        scores_df = pd.DataFrame(cached_rows + new_rows)

    print(f"\n[llm_scoring] Score distribution:")
    print(f"  Mean: {scores_df['llm_score'].mean():.1f}")
    print(f"  Median: {scores_df['llm_score'].median():.1f}")
    print(f"  Min: {scores_df['llm_score'].min()}, Max: {scores_df['llm_score'].max()}")

    write_csv(scores_df, output_path)
    print(f"[llm_scoring] Wrote {len(scores_df)} scores to {output_path}")
    return output_path
