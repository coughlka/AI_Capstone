#!/usr/bin/env python
"""CLI entrypoint for the within-COAD MMR deficiency classifier.

Original modules: Ayan Choudhury (train_coad_mmr_classifier.py,
validate_coad_mmr_labels.py, compare_shap_vs_biomarkers.py).
Wired into the project pipeline by Keith Coughlin.

Runs the three-step MMR analysis as a standalone pipeline alongside the
main biomarker scoring pipeline (run_pipeline.py):

  Step 1: Validate phenotype labels against the expression matrix
  Step 2: Train and evaluate LR / RF / XGBoost classifiers (CV + holdout, SHAP)
  Step 3: Compare SHAP top features against ranked biomarker candidates

This is intentionally separate from run_pipeline.py because it answers a
different biological question (within-cancer MMR subtype) and depends on
ranked_candidates.csv being already present for Step 3.
"""

import argparse
import sys

from src.mmr.validate_labels import run_mmr_label_validation
from src.mmr.train_classifier import run_mmr_classifier
from src.mmr.compare_shap import run_mmr_shap_vs_biomarkers


def main():
    parser = argparse.ArgumentParser(
        description="Run the COAD MMR deficiency classifier pipeline."
    )
    parser.add_argument(
        "--config",
        type=str,
        default="config/config.yaml",
        help="Path to the configuration YAML file (default: config/config.yaml)",
    )
    parser.add_argument(
        "--skip-validate",
        action="store_true",
        help="Skip the phenotype label validation step.",
    )
    parser.add_argument(
        "--skip-compare",
        action="store_true",
        help="Skip the SHAP-vs-biomarker overlap step.",
    )
    args = parser.parse_args()

    config_path = args.config
    print(f"Running MMR classifier pipeline with config: {config_path}")
    print("=" * 60)

    try:
        if not args.skip_validate:
            print("\n[Step 1/3] Validating MMR phenotype labels...")
            report_path = run_mmr_label_validation(config_path)
            print(f"  Report: {report_path}")
        else:
            print("\n[Step 1/3] Skipped (--skip-validate)")

        print("\n[Step 2/3] Training MMR classifier...")
        summary = run_mmr_classifier(config_path)
        print(f"  Best model: {summary['best_model']}")
        print(f"  Class counts: {summary['class_counts']}")

        if not args.skip_compare:
            print("\n[Step 3/3] Comparing SHAP top features vs ranked biomarkers...")
            summary_path = run_mmr_shap_vs_biomarkers(config_path)
            print(f"  Summary: {summary_path}")
        else:
            print("\n[Step 3/3] Skipped (--skip-compare)")

        print("\n" + "=" * 60)
        print("MMR classifier pipeline completed successfully.")

    except FileNotFoundError as e:
        print(f"\nError: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"\nUnexpected error: {e}", file=sys.stderr)
        raise


if __name__ == "__main__":
    main()
