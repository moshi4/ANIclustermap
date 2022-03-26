from pathlib import Path

import pytest


@pytest.fixture(scope="session")
def data_dir() -> Path:
    """Data directory fixture"""
    return Path(__file__).parent / "data"


@pytest.fixture(scope="session")
def genome_fasta_dir(data_dir: Path) -> Path:
    """Genome fasta direcotry"""
    return data_dir / "genome_fasta"
