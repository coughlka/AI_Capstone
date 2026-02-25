"""Data contract tests for pipeline boundary validation.

Pipeline steps are loosely coupled through CSV files. If one step's output
schema changes, downstream steps silently break or produce wrong results.
These tests validate schemas and constraints at each boundary.
"""

import numpy as np
import pandas as pd
import pytest


# ---------------------------------------------------------------------------
# Expected schemas — the contract between pipeline stages
# ---------------------------------------------------------------------------

OMICS_REQUIRED_COLUMNS = {"gene", "gene_symbol", "log2fc", "p_value", "fdr", "direction", "tumor_mean", "normal_mean", "dataset"}
CANDIDATES_REQUIRED_COLUMNS = {"gene", "gene_symbol", "log2fc", "fdr", "direction", "rank"}
LIT_EVIDENCE_REQUIRED_COLUMNS = {"gene", "pmid", "year", "title", "abstract"}
LLM_SCORES_REQUIRED_COLUMNS = {"gene", "gene_symbol", "llm_score", "rationale"}
PATHWAY_REQUIRED_COLUMNS = {"gene", "pathway_count", "top_pathways"}
RANKED_REQUIRED_COLUMNS = {"gene", "final_score", "omics_score", "literature_score", "pathway_score"}


# ---------------------------------------------------------------------------
# Omics Output Contracts
# ---------------------------------------------------------------------------

class TestOmicsOutputContract:
    """Missing columns here cause silent fallbacks or crashes downstream."""

    def test_schema_has_required_columns(self, sample_omics_df):
        """Missing `fdr` would cause scoring.py to fall back to legacy
        mean_expr mode silently. Missing `gene_symbol` would break PubMed."""
        assert OMICS_REQUIRED_COLUMNS.issubset(set(sample_omics_df.columns))

    def test_fdr_in_valid_range(self, sample_omics_df):
        """FDR outside [0,1] is scientifically meaningless and would break
        the DE signal formula (-log10 of negative = NaN)."""
        assert (sample_omics_df["fdr"] >= 0).all(), "FDR values below 0"
        assert (sample_omics_df["fdr"] <= 1).all(), "FDR values above 1"

    def test_pvalue_in_valid_range(self, sample_omics_df):
        """Same as FDR — p-values outside [0,1] are invalid."""
        assert (sample_omics_df["p_value"] >= 0).all()
        assert (sample_omics_df["p_value"] <= 1).all()

    def test_direction_valid_values(self, sample_omics_df):
        """Direction must be 'up' or 'down'. Any other value would break
        the web app's direction filtering."""
        valid = {"up", "down"}
        assert set(sample_omics_df["direction"].unique()).issubset(valid)

    def test_log2fc_is_finite(self, sample_omics_df):
        """Infinite log2fc would produce infinite DE signal, causing one gene
        to dominate all rankings regardless of other evidence."""
        assert np.isfinite(sample_omics_df["log2fc"]).all()

    def test_direction_consistent_with_log2fc(self, sample_omics_df):
        """Direction should match the sign of log2fc. Inconsistency here means
        the web app shows 'upregulated' for genes that are actually downregulated."""
        for _, row in sample_omics_df.iterrows():
            if row["log2fc"] > 0:
                assert row["direction"] == "up", f"Gene {row['gene']}: log2fc={row['log2fc']} but direction={row['direction']}"
            elif row["log2fc"] < 0:
                assert row["direction"] == "down", f"Gene {row['gene']}: log2fc={row['log2fc']} but direction={row['direction']}"

    def test_no_duplicate_genes(self, sample_omics_df):
        """Duplicate gene IDs would cause ambiguous joins downstream,
        potentially doubling a gene's evidence scores."""
        assert not sample_omics_df["gene"].duplicated().any()


# ---------------------------------------------------------------------------
# LLM Scores Output Contract
# ---------------------------------------------------------------------------

class TestLLMScoresContract:

    def test_schema_has_required_columns(self, sample_llm_scores_df):
        assert LLM_SCORES_REQUIRED_COLUMNS.issubset(set(sample_llm_scores_df.columns))

    def test_scores_in_valid_range(self, sample_llm_scores_df):
        """Scores outside [0,100] break the normalization assumptions."""
        assert (sample_llm_scores_df["llm_score"] >= 0).all()
        assert (sample_llm_scores_df["llm_score"] <= 100).all()

    def test_no_nan_scores(self, sample_llm_scores_df):
        """NaN scores propagate through the weighted formula, silently
        dropping genes from the final ranking."""
        assert not sample_llm_scores_df["llm_score"].isna().any()

    def test_rationale_not_empty(self, sample_llm_scores_df):
        """Empty rationales mean the LLM failed to explain its scoring,
        making the results uninterpretable and unauditable."""
        for _, row in sample_llm_scores_df.iterrows():
            assert row["rationale"].strip(), f"Empty rationale for {row['gene']}"


# ---------------------------------------------------------------------------
# Pathway Output Contract
# ---------------------------------------------------------------------------

class TestPathwayContract:

    def test_schema_has_required_columns(self, sample_pathway_df):
        assert PATHWAY_REQUIRED_COLUMNS.issubset(set(sample_pathway_df.columns))

    def test_pathway_count_non_negative(self, sample_pathway_df):
        """Negative pathway counts are meaningless and would produce
        negative scores after normalization."""
        assert (sample_pathway_df["pathway_count"] >= 0).all()


# ---------------------------------------------------------------------------
# Ranked Output Contract
# ---------------------------------------------------------------------------

class TestRankedOutputContract:

    def test_schema_has_required_columns(self, sample_ranked_df):
        assert RANKED_REQUIRED_COLUMNS.issubset(set(sample_ranked_df.columns))

    def test_all_scores_in_valid_range(self, sample_ranked_df):
        """All individual and final scores should be in [0, 100]."""
        for col in ["final_score", "omics_score", "literature_score", "pathway_score"]:
            assert (sample_ranked_df[col] >= 0).all(), f"{col} has values below 0"
            assert (sample_ranked_df[col] <= 100).all(), f"{col} has values above 100"

    def test_sorted_by_final_score_descending(self, sample_ranked_df):
        """If not sorted, rank-based validation is meaningless — rank 1
        wouldn't be the top gene."""
        scores = sample_ranked_df["final_score"].tolist()
        assert scores == sorted(scores, reverse=True), "Not sorted by final_score descending"

    def test_no_duplicate_genes(self, sample_ranked_df):
        """Duplicate genes would inflate one gene's rank at the expense of others."""
        assert not sample_ranked_df["gene"].duplicated().any()

    def test_no_nan_scores(self, sample_ranked_df):
        """NaN in any score column means a gene can't be properly ranked."""
        for col in ["final_score", "omics_score", "literature_score", "pathway_score"]:
            assert not sample_ranked_df[col].isna().any(), f"NaN found in {col}"


# ---------------------------------------------------------------------------
# Cross-Stage Gene ID Consistency
# ---------------------------------------------------------------------------

class TestGeneIDConsistency:
    """Gene IDs must be consistent across pipeline stages. Mismatched formats
    (e.g., versioned vs unversioned Ensembl IDs) silently drop genes."""

    def test_omics_genes_present_in_ranked(self, sample_omics_df, sample_ranked_df):
        """Every gene in ranked output should trace back to omics evidence.
        If ranked contains IDs not in omics, something is fabricating genes."""
        ranked_genes = set(sample_ranked_df["gene"])
        omics_genes = set(sample_omics_df["gene"])
        assert ranked_genes.issubset(omics_genes), \
            f"Ranked output contains genes not in omics: {ranked_genes - omics_genes}"

    def test_llm_scored_genes_match_lit_evidence(self, sample_llm_scores_df, sample_lit_df):
        """Every gene with LLM scores should have literature evidence.
        LLM scores without source abstracts means the LLM was scoring blind."""
        scored_genes = set(sample_llm_scores_df["gene"])
        lit_genes = set(sample_lit_df["gene"])
        assert scored_genes.issubset(lit_genes), \
            f"LLM scored genes without literature: {scored_genes - lit_genes}"
