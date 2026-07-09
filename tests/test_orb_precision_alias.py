"""Pins ORB's precision-string normalization.

orb_models validates the matmul-precision string strictly: the bare "float32"
is not a valid value and silently demotes to the slow "highest" mode, which also
turns on the donated-buffer aot_autograd path that breaks double-backward
Hessians. `_normalize_orb_precision` rewrites the legacy / unified tokens to the
canonical orb_models strings before they reach the loader, so a --config YAML
carrying calc.precision: float32 (or fp32/fp64) cannot slip a demoted precision
through. The module-level helper imports no torch / orb_models, so it is testable
without the ORB backend installed.
"""

from __future__ import annotations

import pytest

from pdb2reaction.backends.orb import _normalize_orb_precision


@pytest.mark.parametrize(
    "raw, expected",
    [
        ("float32", "float32-high"),   # invalid orb string -> canonical fast fp32
        ("fp32", "float32-high"),       # unified token
        ("FP32", "float32-high"),       # case-insensitive
        ("fp64", "float64"),            # unified token
        ("FP64", "float64"),
        ("float64", "float64"),         # already canonical -> untouched
        ("float32-high", "float32-high"),
        ("float32-highest", "float32-highest"),  # unknown-but-valid -> passthrough
    ],
)
def test_normalize_orb_precision(raw, expected) -> None:
    assert _normalize_orb_precision(raw) == expected
