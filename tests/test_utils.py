"""Unit tests for pdb2reaction.utils helpers."""

from __future__ import annotations

import numpy as np


def test_pretty_block_with_numpy_scalars() -> None:
    """pretty_block should normalize NumPy scalars for YAML-safe dumping."""
    from pdb2reaction.utils import pretty_block

    payload = {
        "freeze_atoms": [np.int64(0), np.int64(3)],
        "ratio": np.float64(1.25),
    }
    text = pretty_block("freeze_atoms (effective)", payload)

    assert "freeze_atoms" in text
    assert "- 0" in text
    assert "- 3" in text
    assert "ratio: 1.25" in text
