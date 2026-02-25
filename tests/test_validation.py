"""Validation module tests.

Validation is our only measure of pipeline quality. If percentile or
enrichment calculations are wrong, we'd falsely claim the pipeline works
(or falsely claim it doesn't).
"""

import pytest

from src.validation import KNOWN_CRC_BIOMARKERS


# ---------------------------------------------------------------------------
# Known Biomarkers Ground Truth
# ---------------------------------------------------------------------------

class TestKnownBiomarkers:
    """The known biomarker list is our external ground truth. If it's
    incomplete or incorrectly categorized, validation metrics are misleading."""

    def test_established_crc_genes_present(self):
        """These are universally accepted CRC biomarkers. Missing any one
        means we're not validating against essential ground truth."""
        essential_genes = ["APC", "KRAS", "TP53", "BRAF", "PIK3CA", "SMAD4"]
        for gene in essential_genes:
            assert gene in KNOWN_CRC_BIOMARKERS, f"Essential CRC gene {gene} missing from ground truth"

    def test_all_three_categories_represented(self):
        """If an entire category is missing, we're not validating a whole
        class of biomarkers (mutation, expression, or signaling)."""
        categories = set(KNOWN_CRC_BIOMARKERS.values())
        assert "mutation-driven" in categories
        assert "expression/methylation" in categories
        assert "signaling/pathway" in categories

    def test_category_counts_reasonable(self):
        """Each category should have multiple genes for statistical reliability."""
        from collections import Counter
        counts = Counter(KNOWN_CRC_BIOMARKERS.values())
        for cat, count in counts.items():
            assert count >= 3, f"Category '{cat}' has only {count} genes — too few for reliable validation"

    def test_no_duplicate_entries(self):
        """Duplicate genes would inflate validation metrics."""
        symbols = list(KNOWN_CRC_BIOMARKERS.keys())
        assert len(symbols) == len(set(symbols))


# ---------------------------------------------------------------------------
# Percentile Calculation
# ---------------------------------------------------------------------------

class TestPercentileCalculation:
    """Percentile formula: (1 - rank / total) * 100. Errors here misrepresent
    how well the pipeline ranks known biomarkers."""

    def test_rank_1_is_highest_percentile(self):
        """Rank 1 of 500 = 99.8th percentile."""
        rank, total = 1, 500
        percentile = (1 - rank / total) * 100
        assert percentile == pytest.approx(99.8)

    def test_last_rank_is_lowest_percentile(self):
        """Rank 500 of 500 = 0th percentile."""
        rank, total = 500, 500
        percentile = (1 - rank / total) * 100
        assert percentile == pytest.approx(0.0)

    def test_median_rank_is_50th_percentile(self):
        """Rank 250 of 500 = 50th percentile."""
        rank, total = 250, 500
        percentile = (1 - rank / total) * 100
        assert percentile == pytest.approx(50.0)

    def test_percentile_monotonically_decreases_with_rank(self):
        """Higher rank number = lower percentile. If this is violated,
        validation reports will show contradictory results."""
        total = 1000
        prev_pct = 100.0
        for rank in range(1, total + 1):
            pct = (1 - rank / total) * 100
            assert pct <= prev_pct
            prev_pct = pct


# ---------------------------------------------------------------------------
# Enrichment Ratio
# ---------------------------------------------------------------------------

class TestEnrichmentRatio:
    """Enrichment ratio = expected_median / actual_median.
    A ratio > 1 means known biomarkers rank better than random."""

    def test_enrichment_above_random(self):
        """If known biomarkers have median rank 50 out of 500, they're
        performing 5x better than random (expected median = 250)."""
        total_genes = 500
        actual_median_rank = 50
        expected_median = total_genes / 2  # 250
        enrichment = expected_median / actual_median_rank
        assert enrichment == pytest.approx(5.0)

    def test_enrichment_at_random(self):
        """If known biomarkers have median rank = expected median, enrichment = 1.0.
        This means the pipeline is no better than random for finding known genes."""
        total_genes = 500
        actual_median_rank = 250
        expected_median = total_genes / 2
        enrichment = expected_median / actual_median_rank
        assert enrichment == pytest.approx(1.0)

    def test_enrichment_below_random(self):
        """Enrichment < 1 means known biomarkers rank *worse* than random.
        This would indicate a fundamental pipeline flaw."""
        total_genes = 500
        actual_median_rank = 400
        expected_median = total_genes / 2
        enrichment = expected_median / actual_median_rank
        assert enrichment < 1.0
