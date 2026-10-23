"""Pins ORB's precision-string normalization.

orb_models validates the matmul-precision string strictly ("float32-high" =
TF32 matmul, "float32-highest" = IEEE FP32, "float64"). `_normalize_orb_precision`
rewrites the unified tokens and the bare "float32" to the canonical orb_models
strings before they reach the loader, so a --config YAML carrying
calc.precision: fp32 means IEEE FP32 as on every other backend. The module-level
helper imports no torch / orb_models, so it is testable without the ORB backend
installed.
"""

from __future__ import annotations

import inspect

import pytest

from pdb2reaction.backends.orb import (
    OrbASECalculator,
    OrbCalculator,
    _normalize_orb_precision,
)


@pytest.mark.parametrize(
    "raw, expected",
    [
        ("float32", "float32-highest"),  # bare orb string -> IEEE FP32
        ("fp32", "float32-highest"),     # unified token
        ("FP32", "float32-highest"),     # case-insensitive
        ("tf32", "float32-high"),        # unified token for TF32 matmul
        ("fp64", "float64"),            # unified token
        ("FP64", "float64"),
        ("float64", "float64"),         # already canonical -> untouched
        ("float32-high", "float32-high"),
        ("float32-highest", "float32-highest"),  # already canonical -> untouched
    ],
)
def test_normalize_orb_precision(raw, expected) -> None:
    assert _normalize_orb_precision(raw) == expected


def test_direct_orb_api_defaults_to_float64() -> None:
    """Direct API and ASE/DMF factory must match the CLI's ORB default.

    Loading ORB is intentionally unnecessary here: the public signatures are
    the fallback used by ``create_calculator`` / ``create_ase_calculator`` when
    callers do not pass a precision keyword.
    """
    calc_default = inspect.signature(OrbCalculator).parameters["precision"].default
    ase_default = inspect.signature(OrbASECalculator.__new__).parameters[
        "precision"
    ].default
    assert calc_default == "float64"
    assert ase_default == "float64"
