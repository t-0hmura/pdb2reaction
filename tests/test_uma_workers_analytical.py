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

# UMA analytical Hessians need the EFS force graph, not backbone training.
import torch


def _scope_model(wrapped=True):
    inner = torch.nn.Module()
    inner.backbone = torch.nn.Sequential(torch.nn.Linear(3, 3), torch.nn.Dropout(0.2))
    inner.backbone.composition_dropout = 0.1
    head = torch.nn.Sequential(torch.nn.Linear(3, 1), torch.nn.Dropout(0.3))
    entry = torch.nn.Module()
    entry.head = head
    inner.output_heads = torch.nn.ModuleDict({
        "energyandforcehead": entry,
        "other": torch.nn.Linear(3, 1),
    })
    model = torch.nn.Module() if wrapped else inner
    if wrapped:
        model.module = inner
    model.eval()
    return model, inner, head


def _scope_state(model):
    return (
        tuple((name, module.training, getattr(module, "p", None))
              for name, module in model.named_modules()),
        tuple((name, parameter.requires_grad)
              for name, parameter in model.named_parameters()),
    )


@pytest.mark.parametrize("wrapped", [False, True])
@pytest.mark.parametrize("fail", [False, True])
def test_uma_head_scope_restores_mixed_state(wrapped, fail):
    model, inner, head = _scope_model(wrapped)
    # A saved training backbone is made eval during derivatives, then restored.
    model.train()
    inner.backbone[0].eval()
    head[0].eval()
    next(model.parameters()).requires_grad_(False)
    before = _scope_state(model)
    try:
        with _uma_analytical_head_scope(model):
            assert not model.training
            assert all(not module.training for module in inner.backbone.modules())
            assert not inner.output_heads["other"].training
            assert head.training and head[0].training
            assert not head[1].training and head[1].p == 0.0
            assert inner.backbone.composition_dropout == 0.1
            assert inner.backbone[1].p == 0.2
            assert not any(parameter.requires_grad for parameter in model.parameters())
            if fail:
                raise ValueError("injected derivative failure")
    except ValueError as exc:
        assert fail and str(exc) == "injected derivative failure"
    assert _scope_state(model) == before
    assert inner.backbone.composition_dropout == 0.1


@pytest.mark.parametrize("layout", ["missing", "not_module", "backbone_alias"])
def test_uma_head_scope_rejects_unknown_layout_without_mutation(layout):
    model, inner, head = _scope_model()
    if layout == "missing":
        del inner.output_heads["energyandforcehead"]
    elif layout == "not_module":
        del inner.output_heads["energyandforcehead"].head
        inner.output_heads["energyandforcehead"].head = object()
    else:
        inner.output_heads["energyandforcehead"].head = inner.backbone
    before = _scope_state(model)
    with pytest.raises(RuntimeError, match="UMA analytical Hessian requires"):
        with _uma_analytical_head_scope(model):
            pytest.fail("Unsupported UMA layout must not fall back to global training.")
    assert _scope_state(model) == before


def test_uma_head_scope_restores_after_partial_preparation(monkeypatch):
    model, _, head = _scope_model()
    before = _scope_state(model)
    original_train = head.train

    def interrupted_train(mode=True):
        original_train(mode)
        if mode:
            raise ValueError("injected head preparation failure")
        return head

    monkeypatch.setattr(head, "train", interrupted_train)
    with pytest.raises(ValueError, match="head preparation failure"):
        with _uma_analytical_head_scope(model):
            pytest.fail("Preparation did not fail.")
    assert _scope_state(model) == before


class _ScopePredictor:
    def __init__(self, model, inner, head, fail=False):
        self.model, self.inner, self.head = model, inner, head
        self.fail = fail
        self.derivative_calls = 0

    def predict(self, batch):
        if self.head.training:
            self.derivative_calls += 1
            assert not self.model.training
            assert all(not module.training for module in self.inner.backbone.modules())
            assert not self.inner.output_heads["other"].training
            assert self.inner.backbone.composition_dropout == 0.1
            assert not any(p.requires_grad for p in self.model.parameters())
            if self.fail:
                raise ValueError("injected predictor failure")
        # Deterministic stand-in for UMA's training-only functional routing.
        coefficient = 11.0 if self.inner.backbone.training else 1.0
        energy = coefficient * (batch.pos ** 2).sum()
        # Like UMA, the ordinary force evaluation consumes the energy graph;
        # EFS-head training must retain it for the subsequent Hessian.
        forces = -torch.autograd.grad(
            energy, batch.pos,
            create_graph=self.head.training, retain_graph=self.head.training,
        )[0]
        return {"energy": energy.reshape(1), "forces": forces}


@pytest.mark.parametrize("fail", [False, True])
def test_uma_analytical_owner_uses_head_scope_and_restores(fail):
    from types import SimpleNamespace

    model, inner, head = _scope_model()
    next(model.parameters()).requires_grad_(False)
    before = _scope_state(model)
    predictor = _ScopePredictor(model, inner, head, fail=fail)
    batch = SimpleNamespace(pos=torch.tensor([[0.2, 0.3, 0.4]], dtype=torch.float64,
                                            requires_grad=True))
    compute = _scope_owner_call(predictor, batch)
    if fail:
        with pytest.raises(ValueError, match="injected predictor failure"):
            compute()
    else:
        hessian = compute()
        torch.testing.assert_close(hessian.reshape(3, 3), 2.0 * torch.eye(3, dtype=torch.float64))
    assert predictor.derivative_calls == 1
    assert _scope_state(model) == before

from pdb2reaction.backends.uma import _uma_analytical_head_scope


def _scope_owner_call(predictor, batch):
    core = object.__new__(UMAcore)
    core.elem = ["H"]
    core.predict = predictor
    core.device = torch.device("cpu")
    core.parallel_predict = False
    core.has_torch_model = True
    core._ase_to_batch = lambda atoms: batch
    return lambda: core.compute(np.array([[0.2, 0.3, 0.4]]), hessian=True)["hessian"]


def test_uma_ordinary_prediction_does_not_require_efs_head():
    from types import SimpleNamespace

    model = torch.nn.Linear(3, 1).eval()
    predictor = SimpleNamespace(model=model, predict=lambda batch: {
        "energy": torch.tensor([2.0]), "forces": torch.zeros((1, 3)),
    })
    core = object.__new__(UMAcore)
    core.elem, core.predict = ["H"], predictor
    core.device = torch.device("cpu")
    core.parallel_predict, core.has_torch_model = False, True
    core._ase_to_batch = lambda atoms: SimpleNamespace(pos=torch.zeros((1, 3)))
    before = _scope_state(model)
    result = core.compute(np.zeros((1, 3)), forces=True)
    assert result["energy"] == 2.0 and result["hessian"] is None
    assert _scope_state(model) == before



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
