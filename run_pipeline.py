#!/usr/bin/env python
"""Main CLI entrypoint for running the gene prioritization pipeline."""

import argparse
import sys

from src.omics import run_omics
from src.ml_importance import run_ml_importance
from src.pubmed import run_pubmed
from src.pathway import run_pathway
from src.llm_scoring import run_llm_scoring
from src.scoring import run_scoring
from src.validation import run_validation


def main():
    """Run the full pipeline: omics -> ml_importance -> pubmed -> pathway -> llm_scoring -> scoring -> validation."""
    parser = argparse.ArgumentParser(
        description="Run the gene prioritization pipeline."
    )
    parser.add_argument(
        '--config',
        type=str,
        default='config/config.yaml',
        help='Path to the configuration YAML file (default: config/config.yaml)'
    )
    args = parser.parse_args()

    config_path = args.config
    print(f"Running pipeline with config: {config_path}")
    print("=" * 60)

    try:
        # Step 1: Omics
        print("\n[Step 1/7] Running omics module...")
        omics_output = run_omics(config_path)
        print(f"  Output: {omics_output}")

        # Step 2: ML Importance (Random Forest + SHAP)
        print("\n[Step 2/7] Running ML importance (Random Forest + SHAP)...")
        ml_output = run_ml_importance(config_path)
        print(f"  Output: {ml_output}")

        # Step 3: PubMed
        print("\n[Step 3/7] Running pubmed module...")
        pubmed_output = run_pubmed(config_path)
        print(f"  Output: {pubmed_output}")

        # Step 4: Pathway
        print("\n[Step 4/7] Running pathway module...")
        pathway_output = run_pathway(config_path)
        print(f"  Output: {pathway_output}")

        # Step 5: LLM Scoring
        print("\n[Step 5/7] Running LLM literature scoring...")
        llm_output = run_llm_scoring(config_path)
        print(f"  Output: {llm_output}")

        # Step 6: Scoring
        print("\n[Step 6/7] Running scoring module...")
        scoring_output = run_scoring(config_path)
        print(f"  Output: {scoring_output}")

        # Step 7: Validation
        print("\n[Step 7/7] Running known biomarker validation...")
        validation_output = run_validation(config_path)
        print(f"  Output: {validation_output}")

        print("\n" + "=" * 60)
        print("Pipeline completed successfully!")
        print("\nOutput files:")
        print(f"  - {omics_output}")
        print(f"  - {ml_output}")
        print(f"  - {pubmed_output}")
        print(f"  - {pathway_output}")
        print(f"  - {llm_output}")
        print(f"  - {scoring_output}")
        print(f"  - {validation_output}")

    except FileNotFoundError as e:
        print(f"\nError: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"\nUnexpected error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
