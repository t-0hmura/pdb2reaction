"""Regression tests for incompatible UMA Hessian/worker settings."""

import pytest

from pdb2reaction.backends.base import BackendError
from pdb2reaction.backends.uma import UMACalculator


def test_workers_gt_one_with_analytical_hessian_is_an_error():
    with pytest.raises(BackendError, match=r"workers>1"):
        UMACalculator(workers=2, hessian_calc_mode="Analytical")


def test_workers_gt_one_with_finite_difference_is_allowed():
    calc = UMACalculator(workers=2, hessian_calc_mode="FiniteDifference")
    assert calc._core_kw["workers"] == 2
