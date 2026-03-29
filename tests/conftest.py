"""Shared fixtures for tests."""

import os
import tempfile

import pytest


@pytest.fixture
def tmp_dir():
    """Create a temporary directory for test outputs."""
    with tempfile.TemporaryDirectory() as d:
        yield d


@pytest.fixture
def sample_cache_3col(tmp_dir):
    """Create a 3-column (legacy) cache file."""
    path = os.path.join(tmp_dir, "cache_3col.tsv")
    with open(path, "w") as f:
        f.write("ensembl_id\tgene_symbol\ttype_of_gene\n")
        f.write("ENSG00000141510\tTP53\tprotein-coding\n")
        f.write("ENSG00000012048\tBRCA1\tprotein-coding\n")
        f.write("ENSG00000257151\tLINC02163\tncRNA\n")
    return path


@pytest.fixture
def sample_cache_4col(tmp_dir):
    """Create a 4-column (extended) cache file with aliases."""
    path = os.path.join(tmp_dir, "cache_4col.tsv")
    with open(path, "w") as f:
        f.write("ensembl_id\tgene_symbol\taliases\ttype_of_gene\n")
        f.write("ENSG00000141510\tTP53\tLFS1;TRP53;p53\tprotein-coding\n")
        f.write("ENSG00000012048\tBRCA1\tRNF53;BRCAI\tprotein-coding\n")
        f.write("ENSG00000257151\tLINC02163\tNIHCOLE\tncRNA\n")
        f.write("ENSG00000099999\tFAKEGENE\t\tprotein-coding\n")
    return path
