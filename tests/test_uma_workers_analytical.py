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


def _lazy_preparation_core(monkeypatch, *, device, precision, workers=1):
    """Model the PR298 preparation order with CPU tensors and logical devices."""
    dtype = torch.float64 if precision == "fp64" else torch.float32
    events = []

    class Batch:
        def __init__(self, data):
            self.pos = data.pos
            self.charge, self.spin = data.charge, data.spin
            self.device = self.charge_device = self.spin_device = "cpu"

        def to(self, target):
            # No CUDA allocation: these labels track the ordering contract.
            self.device = self.charge_device = self.spin_device = str(target)
            events.append(("batch_to", str(target)))
            return self

    class LazyPredictor:
        def __init__(self):
            if workers == 1:
                self.model = torch.nn.Module()
                self.model.backbone = torch.nn.Module()
            self.model_device = "cpu"
            self.initialized = False

        def predict(self, batch):
            if not self.initialized:
                assert (batch.device == batch.charge_device == batch.spin_device
                        == self.model_device), "lazy preparation device mismatch"
                events.append(("prepare", self.model_device))
                self.model_device = device
                self.initialized = True
            assert batch.device == "cpu", "the predictor must receive a host batch"
            assert (batch.charge, batch.spin) == (-1, 2)
            assert batch.pos.dtype == dtype and batch.pos.device.type == "cpu"
            batch.to(device)  # FAIR-Chem's transfer follows its lazy preparation.
            assert batch.device == self.model_device
            events.append(("predict", batch.device))
            return {"energy": (batch.pos ** 2).sum().reshape(1),
                    "forces": -2 * batch.pos}

    predictor = LazyPredictor()  # Deliberately exposes no move_to_device hook.

    class AtomicData:
        @staticmethod
        def from_ase(atoms, **kwargs):
            keys = kwargs["r_data_keys"]
            return SimpleNamespace(
                pos=torch.tensor(atoms.positions, dtype=kwargs["target_dtype"]),
                charge=atoms.info.get("charge", 0) if "charge" in keys else 0,
                spin=atoms.info.get("spin", 0) if "spin" in keys else 0,
                dataset=None,
            )

    def collate(data, **kwargs):
        assert len(data) == 1 and data[0].dataset == "omol"
        assert kwargs == {"otf_graph": True}
        return Batch(data[0])

    def serial_factory(model, **kwargs):
        assert model == "uma-s-1p2" and kwargs["device"] == device
        assert kwargs["workers"] == 1
        events.append(("construct", "serial"))
        return predictor

    def parallel_factory(**kwargs):
        assert kwargs["device"] == device and kwargs["num_workers"] == workers
        events.append(("construct", "parallel"))
        return predictor

    monkeypatch.setattr(uma_module, "AtomicData", AtomicData)
    monkeypatch.setattr(uma_module, "data_list_collater", collate)
    monkeypatch.setattr(uma_module, "_UMAInferenceSettings", SimpleNamespace)
    monkeypatch.setattr(uma_module.pretrained_mlip, "get_predict_unit", serial_factory, raising=False)
    monkeypatch.setattr(uma_module.pretrained_mlip, "pretrained_checkpoint_path_from_name",
                        lambda _model: "/unused-checkpoint", raising=False)
    monkeypatch.setattr(uma_module.pretrained_mlip, "get_reference_energies",
                        lambda *_args, **_kwargs: {}, raising=False)
    monkeypatch.setattr(uma_module, "ParallelMLIPPredictUnit", parallel_factory)
    monkeypatch.setattr(uma_module, "guess_inference_settings", lambda name: name)
    core = UMAcore(["O", "H", "H"], model="uma-s-1p2", charge=-1, spin=2,
                   device=device, precision=precision, workers=workers)
    assert events == [("construct", "parallel" if workers > 1 else "serial")]
    assert not predictor.initialized
    return core, predictor, events


@pytest.mark.parametrize("entry", ["compute", "forces_tensor"])
@pytest.mark.parametrize("device", ["cpu", "cuda"])
@pytest.mark.parametrize("precision", ["fp32", "fp64"])
def test_fresh_uma_entries_prepare_before_predictor_owned_transfer(
    monkeypatch, entry, device, precision,
):
    core, predictor, events = _lazy_preparation_core(
        monkeypatch, device=device, precision=precision,
    )
    positions = np.array([[0., 0., 0.], [.757, .586, 0.], [-.757, .586, 0.]])
    for displacement in (0., .001):
        current = positions.copy()
        current[1, 0] += displacement
        if entry == "compute":
            result = core.compute(current, forces=True)
            assert result["energy"] == pytest.approx((current ** 2).sum(), rel=1e-6)
            force = result["forces"]
        else:
            result = core.forces_tensor(current)
            assert not result.requires_grad
            assert result.dtype == (torch.float64 if precision == "fp64" else torch.float32)
            force = result.numpy()
        np.testing.assert_allclose(force, -2 * current, rtol=1e-6)
    assert predictor.initialized
    assert events[1:] == [("prepare", "cpu"), ("batch_to", device), ("predict", device),
                          ("batch_to", device), ("predict", device)]


def test_parallel_uma_constructor_keeps_host_handoff(monkeypatch):
    core, _, events = _lazy_preparation_core(
        monkeypatch, device="cuda", precision="fp32", workers=2,
    )
    assert core.parallel_predict and not core.has_torch_model
    result = core.compute(np.ones((3, 3)), forces=True)
    np.testing.assert_array_equal(result["forces"], -2 * np.ones((3, 3)))
    assert events == [("construct", "parallel"), ("prepare", "cpu"),
                      ("batch_to", "cuda"), ("predict", "cuda")]


def test_lazy_preparation_control_detects_premature_device_transfer(monkeypatch):
    core, predictor, _ = _lazy_preparation_core(
        monkeypatch, device="cuda", precision="fp32",
    )
    batch = core._ase_to_batch(Atoms("OHH", positions=np.zeros((3, 3))))
    batch.to("cuda")  # Recreate the reported old caller's ordering defect.
    with pytest.raises(AssertionError, match="lazy preparation device mismatch"):
        predictor.predict(batch)
    assert not predictor.initialized
