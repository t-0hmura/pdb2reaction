"""Regression tests for incompatible UMA Hessian/worker settings."""

import importlib.util
import sys
from types import ModuleType, SimpleNamespace

import pytest


# The CPU-only CI environment intentionally omits the optional FAIR-Chem runtime.
# UMACalculator validates this option pair before it loads a model, so minimal
# import stubs are sufficient to exercise that dependency-independent guard.
if importlib.util.find_spec("fairchem") is None:
    fairchem = ModuleType("fairchem")
    fairchem.__path__ = []
    core = ModuleType("fairchem.core")
    core.__path__ = []
    core.pretrained_mlip = SimpleNamespace()
    core.FAIRChemCalculator = type("FAIRChemCalculator", (), {})
    datasets = ModuleType("fairchem.core.datasets")
    datasets.__path__ = []
    datasets.data_list_collater = lambda *args, **kwargs: None
    atomic_data = ModuleType("fairchem.core.datasets.atomic_data")
    atomic_data.AtomicData = type("AtomicData", (), {})
    sys.modules.update(
        {
            "fairchem": fairchem,
            "fairchem.core": core,
            "fairchem.core.datasets": datasets,
            "fairchem.core.datasets.atomic_data": atomic_data,
        }
    )

from pdb2reaction.backends.base import BackendError
from pdb2reaction.backends.uma import UMACalculator


def test_workers_gt_one_with_analytical_hessian_is_an_error():
    with pytest.raises(BackendError, match=r"workers\s*>\s*1"):
        UMACalculator(workers=2, hessian_calc_mode="Analytical")


def test_workers_gt_one_with_finite_difference_is_allowed():
    calc = UMACalculator(workers=2, hessian_calc_mode="FiniteDifference")
    assert calc._core_kw["workers"] == 2
