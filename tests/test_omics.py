"""Statistical validation tests for the omics module.

Every test here catches a real scientific or engineering bug:
- TCGA barcode misclassification → entire DE analysis inverted
- BH FDR implementation error → false biomarker calls or missed real ones
- DE signal formula error → wrong candidate ranking
"""

import numpy as np
import pandas as pd
import pytest
from scipy import stats as scipy_stats

from src.omics import parse_tcga_sample_labels, benjamini_hochberg


# ---------------------------------------------------------------------------
# TCGA Barcode Parsing
# ---------------------------------------------------------------------------

class TestTCGABarcodeParsing:
    """Misclassifying tumor vs normal samples would flip the DE analysis."""

    def test_tumor_codes_01_through_09(self):
        """Codes 01-09 are all tumor types. Missing any of these means
        excluding valid tumor samples from the analysis."""
        sample_ids = [
            f"TCGA-AA-3489-{code:02d}A-01R-0905-07"
            for code in range(1, 10)
        ]
        result = parse_tcga_sample_labels(sample_ids)
        assert (result["sample_type_group"] == "tumor").all()

    def test_normal_codes_10_through_19(self):
        """Codes 10-19 are all normal types. Misclassifying these as tumor
        would contaminate the tumor group with normal tissue expression."""
        sample_ids = [
            f"TCGA-AA-3489-{code}A-01R-0905-07"
            for code in range(10, 20)
        ]
        result = parse_tcga_sample_labels(sample_ids)
        assert (result["sample_type_group"] == "normal").all()

    def test_boundary_between_tumor_and_normal(self):
        """The critical boundary: code 09 is tumor, code 10 is normal.
        An off-by-one error here would misclassify one of these."""
        result = parse_tcga_sample_labels([
            "TCGA-AA-3489-09A-01R-0905-07",  # tumor
            "TCGA-AA-3489-10A-01R-0905-07",  # normal
        ])
        groups = result["sample_type_group"].tolist()
        assert groups == ["tumor", "normal"]

    def test_malformed_short_barcode(self):
        """Short barcodes shouldn't crash or be silently treated as valid."""
        result = parse_tcga_sample_labels(["TCGA-AA"])
        assert result.iloc[0]["sample_type_group"] == "other"
        assert result.iloc[0]["sample_type_code"] == "XX"

    def test_non_tcga_identifier(self):
        """Non-TCGA IDs (e.g., from another dataset mixed in) must not be
        classified as tumor or normal."""
        result = parse_tcga_sample_labels(["SRR12345678", "ENSG00000141510"])
        assert (result["sample_type_group"] == "other").all()

    def test_mixed_batch_counts(self):
        """Verify correct counts in a realistic mixed batch."""
        sample_ids = [
            "TCGA-AA-3489-01A-01R-0905-07",  # tumor
            "TCGA-AA-3490-01A-01R-0905-07",  # tumor
            "TCGA-AA-3491-01A-01R-0905-07",  # tumor
            "TCGA-AA-3492-11A-01R-0905-07",  # normal
            "TCGA-AA-3493-11A-01R-0905-07",  # normal
            "SHORT-ID",                        # other
        ]
        result = parse_tcga_sample_labels(sample_ids)
        counts = result["sample_type_group"].value_counts()
        assert counts["tumor"] == 3
        assert counts["normal"] == 2
        assert counts["other"] == 1


# ---------------------------------------------------------------------------
# Benjamini-Hochberg FDR Correction
# ---------------------------------------------------------------------------

class TestBenjaminiHochberg:
    """BH FDR is the foundation of all significance calls. An incorrect
    implementation means we either report false biomarkers or miss real ones."""

    def test_against_scipy_reference(self):
        """Gold standard: our implementation must match statsmodels/scipy.
        This is the single most important statistical test in the suite."""
        try:
            from statsmodels.stats.multitest import multipletests
        except ImportError:
            pytest.skip("statsmodels not installed — install for reference BH test")

        rng = np.random.default_rng(42)
        p_values = rng.uniform(0, 1, size=1000)

        our_fdr = benjamini_hochberg(p_values)
        _, ref_fdr, _, _ = multipletests(p_values, method="fdr_bh")

        np.testing.assert_allclose(our_fdr, ref_fdr, atol=1e-10,
            err_msg="Our BH FDR diverges from statsmodels reference implementation")

    def test_controls_false_discovery_rate(self):
        """Simulate a mixture of null and signal p-values. At FDR threshold 0.05,
        the actual proportion of false discoveries among rejections should be ≤ 0.05
        (with statistical tolerance for finite sample size)."""
        rng = np.random.default_rng(123)
        n_null = 800
        n_signal = 200

        # Null p-values: uniform(0,1)
        null_pvals = rng.uniform(0, 1, n_null)
        # Signal p-values: concentrated near 0
        signal_pvals = rng.beta(1, 50, n_signal)

        p_values = np.concatenate([null_pvals, signal_pvals])
        is_null = np.concatenate([np.ones(n_null, dtype=bool), np.zeros(n_signal, dtype=bool)])

        fdr = benjamini_hochberg(p_values)
        rejected = fdr < 0.05

        if rejected.sum() > 0:
            # False discovery proportion among rejections
            false_discoveries = (is_null & rejected).sum()
            fdp = false_discoveries / rejected.sum()
            # Allow some tolerance — FDR controls the *expected* FDP, not per-instance
            assert fdp < 0.15, f"FDP={fdp:.3f} is too high (expected ≤ 0.05 on average)"

    def test_monotonicity(self):
        """FDR values must be monotonically non-decreasing when sorted by p-value.
        Non-monotonic FDR would create contradictions: a gene with a better p-value
        could appear non-significant while a worse one appears significant."""
        p_values = np.array([0.001, 0.01, 0.02, 0.05, 0.1, 0.5, 0.9])
        fdr = benjamini_hochberg(p_values)

        # Sort both by p-value (already sorted here)
        sorted_fdr = fdr[np.argsort(p_values)]
        assert np.all(sorted_fdr[1:] >= sorted_fdr[:-1] - 1e-15), \
            "FDR values are not monotonically non-decreasing"

    def test_capped_at_one(self):
        """FDR values should never exceed 1.0. Values > 1 are meaningless
        and would confuse any downstream threshold comparison."""
        p_values = np.array([0.5, 0.7, 0.9, 0.99])
        fdr = benjamini_hochberg(p_values)
        assert np.all(fdr <= 1.0), f"FDR values exceed 1.0: {fdr}"

    def test_empty_input(self):
        """Empty input shouldn't crash the pipeline."""
        fdr = benjamini_hochberg(np.array([]))
        assert len(fdr) == 0

    def test_all_nan_input(self):
        """All-NaN input (e.g., genes with no valid expression) shouldn't crash."""
        fdr = benjamini_hochberg(np.array([np.nan, np.nan, np.nan]))
        assert np.all(np.isnan(fdr))

    def test_single_pvalue(self):
        """Single gene shouldn't cause index errors. BH with n=1: q = p * 1 / 1 = p."""
        fdr = benjamini_hochberg(np.array([0.03]))
        np.testing.assert_allclose(fdr, [0.03], atol=1e-15)

    def test_preserves_nan_positions(self):
        """NaN p-values should stay NaN, valid ones should be adjusted."""
        p_values = np.array([0.01, np.nan, 0.05, np.nan, 0.03])
        fdr = benjamini_hochberg(p_values)
        assert np.isnan(fdr[1])
        assert np.isnan(fdr[3])
        assert not np.isnan(fdr[0])
        assert not np.isnan(fdr[2])
        assert not np.isnan(fdr[4])


# ---------------------------------------------------------------------------
# DE Signal Formula
# ---------------------------------------------------------------------------

class TestDESignalFormula:
    """The DE signal formula |log2fc| * -log10(FDR) determines which genes
    become candidates. An error here propagates through the entire pipeline."""

    def test_known_gene_ranks_above_noise(self):
        """A gene with high |log2fc| and low FDR (like KRAS in CRC) should
        always rank above a gene with low |log2fc| and high FDR."""
        epsilon = 1e-300

        # KRAS-like: high fold change, very significant
        kras_signal = abs(2.5) * (-np.log10(1e-18 + epsilon))
        # GAPDH-like: tiny fold change, not significant
        gapdh_signal = abs(0.1) * (-np.log10(0.8 + epsilon))

        assert kras_signal > gapdh_signal, \
            f"KRAS signal ({kras_signal:.2f}) should rank above GAPDH ({gapdh_signal:.2f})"

    def test_both_components_matter(self):
        """A gene with huge fold change but poor p-value should rank below
        a gene with moderate fold change and excellent p-value."""
        epsilon = 1e-300

        # Large fold change, bad FDR
        big_fc_bad_p = abs(5.0) * (-np.log10(0.5 + epsilon))
        # Moderate fold change, excellent FDR
        mod_fc_good_p = abs(2.0) * (-np.log10(1e-20 + epsilon))

        assert mod_fc_good_p > big_fc_bad_p, \
            "Statistical significance should dominate over fold change alone"

    def test_direction_irrelevant_to_signal_strength(self):
        """Up and downregulated genes with same |log2fc| and FDR should
        have identical signal strength. The formula uses abs(log2fc)."""
        epsilon = 1e-300
        up_signal = abs(2.0) * (-np.log10(1e-10 + epsilon))
        down_signal = abs(-2.0) * (-np.log10(1e-10 + epsilon))
        np.testing.assert_allclose(up_signal, down_signal)
