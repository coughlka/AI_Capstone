"""Biomarker candidate API endpoints."""

import os
from typing import Optional

from fastapi import APIRouter, Query, HTTPException
import pandas as pd

router = APIRouter(prefix="/api", tags=["Biomarkers"])

# Data paths - relative to project root
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUTPUTS_DIR = os.path.join(PROJECT_ROOT, "outputs")

# Cache loaded dataframes
_data_cache = {}


def _load_data():
    """Load pipeline output files into memory."""
    if "ranked" not in _data_cache:
        ranked_path = os.path.join(OUTPUTS_DIR, "ranked_candidates.csv")
        omics_path = os.path.join(OUTPUTS_DIR, "omics_evidence.csv")
        pathway_path = os.path.join(OUTPUTS_DIR, "pathway_evidence.csv")

        if not os.path.exists(ranked_path):
            raise HTTPException(
                status_code=503,
                detail="Pipeline outputs not available. Run the pipeline first."
            )

        _data_cache["ranked"] = pd.read_csv(ranked_path)

        if os.path.exists(omics_path):
            _data_cache["omics"] = pd.read_csv(omics_path)
        else:
            _data_cache["omics"] = pd.DataFrame()

        if os.path.exists(pathway_path):
            _data_cache["pathway"] = pd.read_csv(pathway_path)
        else:
            _data_cache["pathway"] = pd.DataFrame()

        # Merge gene symbols from omics into ranked if available
        omics = _data_cache["omics"]
        if not omics.empty and "gene_symbol" in omics.columns:
            symbol_map = omics.set_index("gene")["gene_symbol"].to_dict()
            _data_cache["ranked"]["gene_symbol"] = _data_cache["ranked"]["gene"].map(symbol_map)

    return _data_cache


def _clear_cache():
    """Clear data cache (useful for reloading after pipeline runs)."""
    _data_cache.clear()


@router.get("/candidates")
async def get_candidates(
    page: int = Query(1, ge=1, description="Page number"),
    per_page: int = Query(50, ge=1, le=500, description="Results per page"),
    min_score: Optional[float] = Query(None, ge=0, le=100, description="Minimum final score"),
    direction: Optional[str] = Query(None, description="Filter by direction: 'up' or 'down'"),
    search: Optional[str] = Query(None, description="Search by gene ID or symbol"),
    sort_by: str = Query("final_score", description="Column to sort by"),
    sort_order: str = Query("desc", description="Sort order: 'asc' or 'desc'")
):
    """Get paginated list of ranked biomarker candidates."""
    data = _load_data()
    df = data["ranked"].copy()

    # Join with omics data for direction filtering and additional fields
    omics = data["omics"]
    if not omics.empty:
        merge_cols = ["gene"]
        for col in ["gene_symbol", "direction", "log2fc", "fdr"]:
            if col in omics.columns and col not in df.columns:
                merge_cols.append(col)
        if len(merge_cols) > 1:
            df = df.merge(omics[merge_cols], on="gene", how="left")

    # Apply filters
    if min_score is not None:
        df = df[df["final_score"] >= min_score]

    if direction and "direction" in df.columns:
        df = df[df["direction"] == direction]

    if search:
        search_lower = search.lower()
        mask = df["gene"].str.lower().str.contains(search_lower, na=False)
        if "gene_symbol" in df.columns:
            mask |= df["gene_symbol"].fillna("").str.lower().str.contains(search_lower, na=False)
        df = df[mask]

    # Sort
    ascending = sort_order.lower() == "asc"
    if sort_by in df.columns:
        df = df.sort_values(sort_by, ascending=ascending, na_position="last")

    # Paginate
    total = len(df)
    start = (page - 1) * per_page
    end = start + per_page
    page_data = df.iloc[start:end]

    # Convert to records, handling NaN values
    records = page_data.where(pd.notna(page_data), None).to_dict(orient="records")

    return {
        "total": total,
        "page": page,
        "per_page": per_page,
        "total_pages": (total + per_page - 1) // per_page,
        "candidates": records
    }


@router.get("/genes/{gene_id}")
async def get_gene_detail(gene_id: str):
    """Get detailed evidence for a specific gene."""
    data = _load_data()

    # Find in ranked
    ranked = data["ranked"]
    gene_row = ranked[ranked["gene"] == gene_id]
    if gene_row.empty:
        raise HTTPException(status_code=404, detail=f"Gene {gene_id} not found")

    gene_data = gene_row.iloc[0].to_dict()

    # Get omics evidence
    omics = data["omics"]
    omics_evidence = None
    if not omics.empty:
        omics_row = omics[omics["gene"] == gene_id]
        if not omics_row.empty:
            omics_evidence = omics_row.iloc[0].where(pd.notna(omics_row.iloc[0]), None).to_dict()

    # Get pathway evidence
    pathway = data["pathway"]
    pathway_evidence = None
    if not pathway.empty:
        pathway_row = pathway[pathway["gene"] == gene_id]
        if not pathway_row.empty:
            pathway_evidence = pathway_row.iloc[0].where(pd.notna(pathway_row.iloc[0]), None).to_dict()

    return {
        "gene": gene_id,
        "gene_symbol": omics_evidence.get("gene_symbol") if omics_evidence else None,
        "scores": {
            "final": gene_data.get("final_score"),
            "omics": gene_data.get("omics_score"),
            "literature": gene_data.get("literature_score"),
            "pathway": gene_data.get("pathway_score")
        },
        "omics_evidence": omics_evidence,
        "pathway_evidence": pathway_evidence
    }


@router.get("/stats")
async def get_stats():
    """Get pipeline summary statistics."""
    data = _load_data()
    ranked = data["ranked"]
    omics = data["omics"]

    stats = {
        "total_genes": len(ranked),
        "scoring_weights": {
            "omics": 0.45,
            "literature": 0.35,
            "pathway": 0.20
        },
        "score_range": {
            "min": float(ranked["final_score"].min()) if not ranked.empty else 0,
            "max": float(ranked["final_score"].max()) if not ranked.empty else 0,
            "mean": float(ranked["final_score"].mean()) if not ranked.empty else 0
        }
    }

    if not omics.empty:
        if "fdr" in omics.columns:
            sig = omics[omics["fdr"] < 0.05]
            stats["significant_genes"] = len(sig)

        if "direction" in omics.columns and "gene_symbol" in omics.columns:
            # Get top upregulated genes (lowest FDR among "up" direction)
            up = omics[omics["direction"] == "up"]
            if "fdr" in up.columns and not up.empty:
                up = up.nsmallest(5, "fdr")
                stats["top_upregulated"] = up["gene_symbol"].dropna().tolist()[:5]

            # Get top downregulated genes
            down = omics[omics["direction"] == "down"]
            if "fdr" in down.columns and not down.empty:
                down = down.nsmallest(5, "fdr")
                stats["top_downregulated"] = down["gene_symbol"].dropna().tolist()[:5]

    return stats


@router.post("/reload")
async def reload_data():
    """Reload pipeline data from disk (call after running pipeline)."""
    _clear_cache()
    _load_data()
    return {"status": "reloaded"}
