"""LLM-powered chat endpoint for asking questions about biomarker data."""

import logging
import os
import re
import time
from collections import defaultdict
from typing import List, Optional

import anthropic
import pandas as pd
from dotenv import load_dotenv
from fastapi import APIRouter, HTTPException, Request
from pydantic import BaseModel

from api.biomarkers import _load_data

# Load .env from the web-app directory
load_dotenv(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), ".env"))

logger = logging.getLogger(__name__)

router = APIRouter(prefix="/api", tags=["Chat"])

# ---------------------------------------------------------------------------
# Input-validation constants
# ---------------------------------------------------------------------------
MAX_MESSAGE_LENGTH = 2000  # characters
MAX_HISTORY_MESSAGES = 20  # conversation turns sent to the LLM
# Ensembl gene ID pattern (e.g. ENSG00000141510)
_GENE_ID_RE = re.compile(r"^ENSG\d{11}$")

# ---------------------------------------------------------------------------
# Simple in-memory sliding-window rate limiter
# ---------------------------------------------------------------------------
_RATE_LIMIT_WINDOW = 60  # seconds
_RATE_LIMIT_MAX = 20     # max requests per window per IP
_request_log: dict[str, list[float]] = defaultdict(list)


def _check_rate_limit(client_ip: str) -> None:
    """Enforce a per-IP sliding-window rate limit.

    Raises HTTP 429 if the client has exceeded the allowed number of
    requests within the current window.
    """
    now = time.monotonic()
    timestamps = _request_log[client_ip]

    # Prune timestamps outside the window
    _request_log[client_ip] = [
        ts for ts in timestamps if now - ts < _RATE_LIMIT_WINDOW
    ]

    if len(_request_log[client_ip]) >= _RATE_LIMIT_MAX:
        logger.warning("Rate limit exceeded for %s", client_ip)
        raise HTTPException(
            status_code=429,
            detail="Rate limit exceeded. Please wait before sending more messages.",
        )

    _request_log[client_ip].append(now)


# ---------------------------------------------------------------------------
# Gene-context validation
# ---------------------------------------------------------------------------
def _validate_gene_context(gene_context: Optional[str]) -> Optional[str]:
    """Validate and sanitize the gene_context parameter.

    Returns the validated gene ID or None.  Rejects any value that does
    not match the expected Ensembl ID format *and* exist in the loaded
    dataset, preventing arbitrary text from being interpolated into the
    system prompt.
    """
    if gene_context is None:
        return None

    gene_context = gene_context.strip()
    if not gene_context:
        return None

    # Structural check — must look like an Ensembl gene ID
    if not _GENE_ID_RE.match(gene_context):
        logger.warning("Rejected invalid gene_context format: %s", gene_context[:60])
        return None

    # Existence check — must be present in our dataset
    try:
        data = _load_data()
        known_genes = set(data["ranked"]["gene"].tolist())
        if not data["omics"].empty:
            known_genes.update(data["omics"]["gene"].tolist())

        if gene_context not in known_genes:
            logger.warning("Rejected unknown gene_context: %s", gene_context)
            return None
    except Exception:
        # If data can't be loaded, disallow context injection entirely
        return None

    return gene_context


# ---------------------------------------------------------------------------
# Anthropic client (initialized lazily)
# ---------------------------------------------------------------------------
_client = None


def _get_client():
    """Get or create the Anthropic client."""
    global _client
    if _client is None:
        api_key = os.getenv("ANTHROPIC_API_KEY")
        if not api_key:
            raise HTTPException(
                status_code=500,
                detail="ANTHROPIC_API_KEY not set. Add it to web-app/.env"
            )
        _client = anthropic.Anthropic(api_key=api_key)
    return _client


class ChatMessage(BaseModel):
    role: str  # "user" or "assistant"
    content: str


class ChatRequest(BaseModel):
    message: str
    history: List[ChatMessage] = []
    gene_context: Optional[str] = None


def _build_system_prompt(gene_context: Optional[str] = None) -> str:
    """Build the system prompt with tiered data context.

    Tier 1: Static dataset summary (always included)
    Tier 2: Top 50 ranked candidates (always included)
    Tier 3: Gene-specific detail (included when gene_context is provided)
    """
    # Tier 1: Static summary
    system = (
        "You are a research assistant for a colorectal cancer (CRC) biomarker discovery project. "
        "You help users understand and interpret the pipeline results.\n\n"
        "## Dataset Background\n"
        "This project analyzes RNA sequencing data from TCGA-COAD (The Cancer Genome Atlas - "
        "Colon Adenocarcinoma). The dataset contains 473 tumor samples and 41 normal tissue samples.\n\n"
        "## Methodology\n"
        "- Differential expression analysis using Welch's t-test on log2-transformed expression values\n"
        "- Multiple testing correction using Benjamini-Hochberg FDR\n"
        "- Genes are scored using: omics_signal = |log2FC| * -log10(FDR)\n"
        "- Final scores combine three evidence streams with weights: "
        "omics (45%), literature (35%), pathway (20%)\n"
        "- Scores are normalized to a 0-100 scale using min-max normalization\n\n"
        "## Key Terms\n"
        "- log2FC: log2 fold change between tumor and normal expression. "
        "Positive = higher in tumor (upregulated), negative = lower in tumor (downregulated).\n"
        "- FDR: False Discovery Rate. Lower values mean higher statistical confidence.\n"
        "- Direction: 'up' means the gene is more active in tumors, 'down' means less active.\n\n"
        "## Instructions\n"
        "- Answer questions clearly and concisely\n"
        "- When discussing specific genes, reference their scores and evidence\n"
        "- If asked about something outside the dataset, say so\n"
        "- Use plain language accessible to someone without a medical background\n"
    )

    # Tier 2: Top 50 candidates
    try:
        data = _load_data()
        ranked = data["ranked"]
        omics = data["omics"]

        if not ranked.empty:
            top50 = ranked.head(50).copy()

            # Merge in extra columns from omics
            if not omics.empty:
                merge_cols = ["gene"]
                for col in ["gene_symbol", "direction", "log2fc", "fdr"]:
                    if col in omics.columns and col not in top50.columns:
                        merge_cols.append(col)
                if len(merge_cols) > 1:
                    top50 = top50.merge(omics[merge_cols], on="gene", how="left")

            system += "\n## Top 50 Ranked Candidates\n"
            system += top50.to_string(index=False, max_colwidth=20)
            system += "\n"
    except Exception:
        system += "\n(Pipeline data not available.)\n"

    # Tier 3: Gene-specific detail (gene_context has already been validated)
    if gene_context:
        try:
            data = _load_data()
            omics = data["omics"]
            ranked = data["ranked"]

            gene_info = []

            # Ranked scores
            row = ranked[ranked["gene"] == gene_context]
            if not row.empty:
                r = row.iloc[0]
                gene_info.append(
                    f"Scores: final={r.get('final_score', 'N/A'):.2f}, "
                    f"omics={r.get('omics_score', 'N/A'):.2f}, "
                    f"literature={r.get('literature_score', 'N/A'):.2f}, "
                    f"pathway={r.get('pathway_score', 'N/A'):.2f}"
                )

            # Omics evidence
            if not omics.empty:
                orow = omics[omics["gene"] == gene_context]
                if not orow.empty:
                    o = orow.iloc[0]
                    symbol = o.get("gene_symbol", "")
                    gene_info.append(f"Symbol: {symbol}")
                    gene_info.append(f"log2FC: {o.get('log2fc', 'N/A')}")
                    gene_info.append(f"FDR: {o.get('fdr', 'N/A')}")
                    gene_info.append(f"Direction: {o.get('direction', 'N/A')}")
                    gene_info.append(f"Tumor mean: {o.get('tumor_mean', 'N/A')}")
                    gene_info.append(f"Normal mean: {o.get('normal_mean', 'N/A')}")

            if gene_info:
                system += f"\n## Current Gene Context: {gene_context}\n"
                system += "The user is asking about this specific gene:\n"
                system += "\n".join(f"- {info}" for info in gene_info)
                system += "\n"
        except Exception:
            pass

    return system


@router.post("/chat")
async def chat(request: ChatRequest, raw_request: Request):
    """Send a message to the LLM with biomarker data context."""

    # --- Input validation & rate limiting ---
    client_ip = raw_request.client.host if raw_request.client else "unknown"
    _check_rate_limit(client_ip)

    # Enforce message length limit
    if len(request.message) > MAX_MESSAGE_LENGTH:
        raise HTTPException(
            status_code=400,
            detail=f"Message too long. Maximum length is {MAX_MESSAGE_LENGTH} characters.",
        )

    if not request.message.strip():
        raise HTTPException(status_code=400, detail="Message cannot be empty.")

    # Validate and sanitize gene_context before it reaches the system prompt
    validated_gene_context = _validate_gene_context(request.gene_context)

    client = _get_client()
    system_prompt = _build_system_prompt(validated_gene_context)

    # Build message history (last N messages + current message)
    messages = [
        {"role": msg.role, "content": msg.content}
        for msg in request.history[-MAX_HISTORY_MESSAGES:]
    ]
    messages.append({"role": "user", "content": request.message})

    try:
        response = client.messages.create(
            model="claude-sonnet-4-20250514",
            max_tokens=1024,
            system=system_prompt,
            messages=messages,
        )
        return {"response": response.content[0].text}
    except anthropic.APIError as e:
        raise HTTPException(status_code=502, detail=f"LLM API error: {str(e)}")
