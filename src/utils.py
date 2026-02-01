"""Utility functions for the pipeline."""

import os
from pathlib import Path

import pandas as pd
import yaml


def load_config(config_path: str) -> dict:
    """Load and parse a YAML configuration file.

    Args:
        config_path: Path to the YAML config file.

    Returns:
        Dictionary containing the parsed configuration.

    Raises:
        FileNotFoundError: If the config file does not exist.
        ValueError: If the YAML is invalid or cannot be parsed.
    """
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Config file not found: {config_path}")

    try:
        with open(config_path, 'r') as f:
            config = yaml.safe_load(f)
    except yaml.YAMLError as e:
        raise ValueError(f"Invalid YAML in config file '{config_path}': {e}")

    if config is None:
        raise ValueError(f"Config file '{config_path}' is empty or contains only comments.")

    return config


def ensure_dirs(config: dict) -> None:
    """Ensure the outputs directory exists.

    Args:
        config: Configuration dictionary with 'paths' section.
    """
    outputs_dir = config.get('paths', {}).get('outputs_dir', 'outputs')
    Path(outputs_dir).mkdir(parents=True, exist_ok=True)


def write_csv(df: pd.DataFrame, path: str) -> None:
    """Write a DataFrame to CSV, creating parent directories if needed.

    Args:
        df: DataFrame to write.
        path: Output file path.
    """
    parent = Path(path).parent
    if parent and not parent.exists():
        parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)


def read_csv(path: str) -> pd.DataFrame:
    """Read a CSV file into a DataFrame.

    Args:
        path: Path to the CSV file.

    Returns:
        DataFrame containing the CSV data.

    Raises:
        FileNotFoundError: If the file does not exist.
    """
    if not os.path.exists(path):
        raise FileNotFoundError(f"CSV file not found: {path}. Please ensure the file exists.")
    return pd.read_csv(path)
