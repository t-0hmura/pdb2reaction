"""Regression tests for incompatible UMA Hessian/worker settings."""

import importlib.util
import sys
from types import ModuleType, SimpleNamespace

import numpy as np
import pytest
from ase import Atoms


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
from pdb2reaction.backends import uma as uma_module
from pdb2reaction.backends.uma import UMAcore, UMACalculator, _positive_worker_count


def test_workers_gt_one_with_analytical_hessian_is_an_error():
    with pytest.raises(BackendError, match=r"workers\s*>\s*1"):
        UMACalculator(workers=2, hessian_calc_mode="Analytical")


def test_workers_gt_one_with_finite_difference_is_allowed():
    calc = UMACalculator(workers=2, hessian_calc_mode="FiniteDifference")
    assert calc._core_kw["workers"] == 2


def test_analytical_mode_requests_differentiable_inference_settings():
    analytical = UMACalculator(hessian_calc_mode="Analytical")
    finite_difference = UMACalculator(hessian_calc_mode="FiniteDifference")

    assert analytical._core_kw["analytical_hessian"] is True
    assert finite_difference._core_kw["analytical_hessian"] is False


@pytest.mark.parametrize(
    ("precision", "expected_dtype"),
    [("fp32", "float32"), ("fp64", "float64")],
)
def test_analytical_mode_uses_noncompiled_precision_matched_settings(
    monkeypatch, precision, expected_dtype
):
    captured = {}

    class FakeSettings:
        def __init__(self, **kwargs):
            self.kwargs = kwargs

    def fake_get_predict_unit(_model, **kwargs):
        captured.update(kwargs)
        return SimpleNamespace()

    monkeypatch.setattr(uma_module, "_UMAInferenceSettings", FakeSettings)
    monkeypatch.setattr(
        uma_module.pretrained_mlip,
        "get_predict_unit",
        fake_get_predict_unit,
        raising=False,
    )

    UMAcore(
        ["H"],
        device="cpu",
        precision=precision,
        analytical_hessian=True,
    )

    settings = captured["inference_settings"]
    assert settings.kwargs == {
        "compile": False,
        "base_precision_dtype": expected_dtype,
    }


def test_serial_uma_batch_stays_on_cpu_for_fairchem_lazy_initialization():
    class FakeData:
        dataset = None

    class FakeAtomicData:
        @staticmethod
        def from_ase(*_args, **_kwargs):
            return FakeData()

    class FakeBatch:
        def to(self, _device):
            pytest.fail("pdb2reaction must leave FAIR-Chem input device transfer to FAIR-Chem")

    core = object.__new__(UMAcore)
    core.has_torch_model = False
    core._AtomicData = FakeAtomicData
    core._collater = lambda *_args, **_kwargs: FakeBatch()
    core.parallel_predict = False
    core.elem = ["H"]
    core.charge = 0
    core.spin = 1
    core.task_name = "omol"
    core.precision = "fp32"
    core._max_neigh = None
    core._radius = None
    core._r_edges = False

    batch = core._ase_to_batch(Atoms("H", positions=[[0.0, 0.0, 0.0]]))

    assert isinstance(batch, FakeBatch)


@pytest.mark.parametrize("value", [0, -1, 1.5, True])
def test_invalid_worker_counts_are_rejected(value):
    with pytest.raises(BackendError, match="positive integer"):
        _positive_worker_count(value, "workers")


@pytest.mark.parametrize(("value", "expected"), [(None, 1), (1, 1), ("2", 2)])
def test_positive_worker_counts_are_preserved(value, expected):
    assert _positive_worker_count(value, "workers") == expected


def test_partial_fd_hessian_without_frozen_atoms_keeps_full_allocation(monkeypatch):
    torch = pytest.importorskip("torch")
    from pysisyphus import _array

    class FakeCore:
        device = torch.device("cpu")
        has_torch_model = False
        parallel_predict = False

        def compute(self, coordinates, *, forces, hessian):
            return {
                "energy": 0.0,
                "forces": -2.0 * np.asarray(coordinates, dtype=float),
            }

    fake = SimpleNamespace(
        _core=FakeCore(),
        _ensure_core=lambda elements: None,
        freeze_atoms=[],
        return_partial_hessian=True,
        hessian_double=True,
    )
    monkeypatch.setattr(
        _array,
        "active_square",
        lambda *args, **kwargs: pytest.fail("full active square copy is unnecessary"),
    )

    result = UMACalculator._build_fd_hessian_gpu(
        fake,
        ["H"],
        np.zeros((1, 3), dtype=float),
    )

    assert tuple(result["hessian"].shape) == (1, 3, 1, 3)


def test_uma_core_is_rebuilt_when_composition_changes(monkeypatch):
    from pdb2reaction.backends import uma

    built = []

    class FakeCore:
        def __init__(self, elements, **kwargs):
            self.elem = [element.capitalize() for element in elements]
            built.append(self.elem)

    monkeypatch.setattr(uma, "UMAcore", FakeCore)
    calculator = UMACalculator()

    calculator._ensure_core(["H", "H"])
    first = calculator._core
    calculator._ensure_core(["H", "H"])
    assert calculator._core is first

    calculator._ensure_core(["N", "N"])
    assert calculator._core is not first
    assert built == [["H", "H"], ["N", "N"]]
