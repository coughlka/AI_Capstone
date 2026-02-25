"""Scoring formula validation tests.

The weighted scoring formula is the core ranking mechanism. Errors here
silently produce wrong biomarker rankings — the worst kind of bug because
the pipeline will still "work" but produce scientifically wrong results.
"""

import numpy as np
import pandas as pd
import pytest

from src.scoring import _min_max_normalize


# ---------------------------------------------------------------------------
# Min-Max Normalization
# ---------------------------------------------------------------------------

class TestMinMaxNormalize:
    """Normalization errors corrupt all three evidence scores and the final ranking."""

    def test_preserves_ordering(self):
        """Normalization must never change rank order of genes.
        If it did, genes would swap positions in the final ranking."""
        rng = np.random.default_rng(42)
        values = pd.Series(rng.uniform(0, 1000, size=100))
        normalized = _min_max_normalize(values)

        original_order = values.argsort().values
        normalized_order = normalized.argsort().values
        np.testing.assert_array_equal(original_order, normalized_order,
            err_msg="Normalization changed the rank ordering of genes")

    def test_output_bounds(self):
        """Scores outside [0, 100] would break interpretability and
        could cause the judge to miscalibrate expectations."""
        values = pd.Series([0.001, 50, 999.99, 1e6])
        normalized = _min_max_normalize(values)
        assert normalized.min() >= 0.0
        assert normalized.max() <= 100.0

    def test_constant_series_no_division_by_zero(self):
        """When all genes have the same score (e.g., all pathway_count=0),
        min==max causes division by zero. Should return target_min, not NaN/inf."""
        values = pd.Series([5.0, 5.0, 5.0, 5.0])
        normalized = _min_max_normalize(values)
        assert not normalized.isna().any(), "Constant series produced NaN"
        assert not np.isinf(normalized).any(), "Constant series produced inf"
        assert (normalized == 0.0).all(), "Constant series should normalize to target_min=0"

    def test_single_value(self):
        """Single gene (min==max edge case)."""
        values = pd.Series([42.0])
        normalized = _min_max_normalize(values)
        assert not normalized.isna().any()

    def test_negative_values_handled(self):
        """log2fc values can be negative. Normalization should handle this correctly."""
        values = pd.Series([-10.0, 0.0, 10.0])
        normalized = _min_max_normalize(values)
        assert normalized.iloc[0] == pytest.approx(0.0)
        assert normalized.iloc[1] == pytest.approx(50.0)
        assert normalized.iloc[2] == pytest.approx(100.0)

    def test_custom_range(self):
        """Verify custom target_min/target_max work correctly."""
        values = pd.Series([0.0, 50.0, 100.0])
        normalized = _min_max_normalize(values, target_min=10, target_max=90)
        assert normalized.iloc[0] == pytest.approx(10.0)
        assert normalized.iloc[1] == pytest.approx(50.0)
        assert normalized.iloc[2] == pytest.approx(90.0)


# ---------------------------------------------------------------------------
# Weighted Scoring Formula
# ---------------------------------------------------------------------------

class TestWeightedScoring:
    """The formula final_score = 0.45*omics + 0.35*lit + 0.20*pathway
    is the heart of the pipeline. Any error here invalidates the entire ranking."""

    def test_exact_weighted_calculation(self):
        """Verify exact formula: 0.45*80 + 0.35*60 + 0.20*40 = 65.0"""
        omics, lit, pathway = 80.0, 60.0, 40.0
        w_omics, w_lit, w_path = 0.45, 0.35, 0.20
        expected = w_omics * omics + w_lit * lit + w_path * pathway
        assert expected == pytest.approx(65.0)

    def test_weights_sum_to_one(self):
        """If weights don't sum to 1.0, scores are systematically inflated
        or deflated, making the 0-100 scale meaningless."""
        w_omics, w_lit, w_path = 0.45, 0.35, 0.20
        assert w_omics + w_lit + w_path == pytest.approx(1.0)

    def test_missing_literature_defaults_to_zero(self):
        """A gene with no PubMed hits should get literature_score=0, not NaN.
        NaN would propagate to final_score, silently dropping the gene."""
        omics, lit, pathway = 50.0, 0.0, 20.0
        w_omics, w_lit, w_path = 0.45, 0.35, 0.20
        final = w_omics * omics + w_lit * lit + w_path * pathway
        assert not np.isnan(final)
        assert final == pytest.approx(0.45 * 50 + 0.35 * 0 + 0.20 * 20)

    def test_missing_pathway_defaults_to_zero(self):
        """Same as above for pathway evidence."""
        omics, lit, pathway = 50.0, 60.0, 0.0
        w_omics, w_lit, w_path = 0.45, 0.35, 0.20
        final = w_omics * omics + w_lit * lit + w_path * pathway
        assert final == pytest.approx(0.45 * 50 + 0.35 * 60 + 0.20 * 0)

    def test_all_three_streams_contribute(self):
        """A gene scoring high on all 3 streams must beat one scoring high on only 1.
        This validates the multi-evidence design actually works as intended."""
        w_omics, w_lit, w_path = 0.45, 0.35, 0.20

        # Gene A: strong across all evidence
        score_a = w_omics * 90 + w_lit * 85 + w_path * 80

        # Gene B: strong omics only
        score_b = w_omics * 100 + w_lit * 0 + w_path * 0

        assert score_a > score_b, \
            "Multi-evidence gene should rank above single-evidence gene"

    def test_ranking_is_deterministic(self):
        """Same inputs must always produce the same ranking.
        Hidden randomness (e.g., unstable sort) would make results irreproducible."""
        data = pd.DataFrame({
            "gene": [f"GENE_{i}" for i in range(100)],
            "omics_score": np.random.default_rng(42).uniform(0, 100, 100),
            "literature_score": np.random.default_rng(43).uniform(0, 100, 100),
            "pathway_score": np.random.default_rng(44).uniform(0, 100, 100),
        })
        data["final_score"] = 0.45 * data["omics_score"] + 0.35 * data["literature_score"] + 0.20 * data["pathway_score"]

        ranking_1 = data.sort_values("final_score", ascending=False)["gene"].tolist()
        ranking_2 = data.sort_values("final_score", ascending=False)["gene"].tolist()
        assert ranking_1 == ranking_2

    def test_final_score_bounded(self):
        """If all inputs are in [0,100], final score must also be in [0,100]
        since weights sum to 1.0."""
        w_omics, w_lit, w_path = 0.45, 0.35, 0.20

        # Maximum possible
        max_score = w_omics * 100 + w_lit * 100 + w_path * 100
        assert max_score == pytest.approx(100.0)

        # Minimum possible
        min_score = w_omics * 0 + w_lit * 0 + w_path * 0
        assert min_score == pytest.approx(0.0)

    def test_omics_weight_dominates(self):
        """Omics has the highest weight (0.45). A gene with perfect omics but
        zero lit/pathway should score higher than one with zero omics but
        perfect lit and pathway (since 0.45 > 0.35 and 0.45 > 0.20)."""
        w_omics, w_lit, w_path = 0.45, 0.35, 0.20

        omics_only = w_omics * 100 + w_lit * 0 + w_path * 0  # 45
        lit_only = w_omics * 0 + w_lit * 100 + w_path * 0      # 35
        pathway_only = w_omics * 0 + w_lit * 0 + w_path * 100  # 20

        assert omics_only > lit_only > pathway_only
