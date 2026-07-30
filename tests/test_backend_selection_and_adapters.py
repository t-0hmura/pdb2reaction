"""Backend selection, state restoration, and ASE adapter contracts."""

from __future__ import annotations

from pathlib import Path
import sys
from types import ModuleType

import numpy as np
import pytest


def test_auto_backend_checks_runtime_dependency(monkeypatch) -> None:
    import pdb2reaction.backends as backends

    imported = []

    def fake_import(name):
        imported.append(name)
        if name in {
            "fairchem.core",
            "orb_models.forcefield.pretrained",
        }:
            raise ImportError(name)
        return object()

    monkeypatch.setattr(backends.importlib, "import_module", fake_import)
    monkeypatch.setattr(backends, "_import_cls", lambda backend, key: object())

    assert backends.resolve_backend("auto") == "mace"
    assert "orb_models.forcefield.pretrained" in imported
    assert "mace.calculators" in imported


def test_auto_backend_provenance_uses_resolved_name(monkeypatch) -> None:
    import pdb2reaction.backends as backends
    from pdb2reaction.core.utils import calculator_provenance

    monkeypatch.setattr(backends, "resolve_backend", lambda backend: "mace")

    provenance = calculator_provenance({"backend": "auto"})

    assert provenance["mlip_backend"] == "mace"
    assert provenance["mlip_precision"] == "fp64"


def test_autograd_hessian_restores_each_module_state() -> None:
    torch = pytest.importorskip("torch")
    from pdb2reaction.backends.base import (
        _prepare_model_for_autograd_hessian,
        _restore_model_after_autograd_hessian,
    )

    model = torch.nn.Sequential(
        torch.nn.Linear(2, 2),
        torch.nn.Sequential(torch.nn.Dropout(p=0.25), torch.nn.Linear(2, 1)),
    )
    model.train(False)
    model[0].train(True)
    model[1].train(True)
    model[1][0].train(False)
    model[1][1].train(True)
    parameters = list(model.parameters())
    parameters[0].requires_grad_(False)
    expected_training = [module.training for module in model.modules()]
    expected_grad = [parameter.requires_grad for parameter in parameters]

    state = _prepare_model_for_autograd_hessian(model, torch)
    _restore_model_after_autograd_hessian(model, state)

    assert [module.training for module in model.modules()] == expected_training
    assert [parameter.requires_grad for parameter in parameters] == expected_grad
    assert model[1][0].p == pytest.approx(0.25)


def test_freeze_atoms_reject_negative_and_out_of_range_indices() -> None:
    from pdb2reaction.backends.base import MLIPCalculator

    with pytest.raises(ValueError, match="non-negative"):
        MLIPCalculator(freeze_atoms=[-1])

    calculator = MLIPCalculator(freeze_atoms=[2])
    with pytest.raises(ValueError, match="outside a 2-atom geometry"):
        calculator._active_and_frozen_dof_idx(2)


def test_aimnet_explicit_device_is_not_silently_dropped(monkeypatch) -> None:
    from pdb2reaction.backends import aimnet2
    from pdb2reaction.backends.base import BackendError

    package = ModuleType("aimnet")
    package.__path__ = []
    calculators = ModuleType("aimnet.calculators")

    class DeviceRejectingCalculator:
        def __init__(self, *args, **kwargs):
            if "device" in kwargs:
                raise TypeError("device is unsupported")

    calculators.AIMNet2Calculator = DeviceRejectingCalculator
    monkeypatch.setitem(sys.modules, "aimnet", package)
    monkeypatch.setitem(sys.modules, "aimnet.calculators", calculators)

    with pytest.raises(BackendError):
        aimnet2.AIMNet2Calculator(device="cuda")

    calculator = aimnet2.AIMNet2Calculator(device="auto")
    assert isinstance(calculator._calculator, DeviceRejectingCalculator)


def test_aimnet_ase_adapter_forwards_charge_spin_and_results(monkeypatch) -> None:
    from ase import Atoms
    from ase.calculators.calculator import Calculator
    from pdb2reaction.backends import aimnet2

    captured = {}

    class FakeBackend:
        def __init__(self, **kwargs):
            captured.update(kwargs)

        def _compute_energy_forces_ev(self, symbols, positions):
            captured["symbols"] = list(symbols)
            captured["positions"] = np.asarray(positions)
            return 1.25, np.full((len(symbols), 3), -0.5)

    monkeypatch.setattr(aimnet2, "AIMNet2Calculator", FakeBackend)
    calculator = aimnet2.AIMNet2ASECalculator(
        model="sentinel",
        device="cpu",
        charge=-1,
        spin=2,
    )
    atoms = Atoms("OH", positions=[[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    atoms.calc = calculator

    assert isinstance(calculator, Calculator)
    assert atoms.get_potential_energy() == pytest.approx(1.25)
    np.testing.assert_allclose(atoms.get_forces(), -0.5)
    assert captured["charge"] == -1
    assert captured["spin"] == 2


@pytest.mark.parametrize("dtype", ["float32", "float64"])
def test_mace_anicc_uses_supported_factory_signature(monkeypatch, dtype) -> None:
    from pdb2reaction.backends.mace import MACECalculator

    package = ModuleType("mace")
    package.__path__ = []
    calculators = ModuleType("mace.calculators")
    captured = {}

    def mace_anicc(**kwargs):
        captured["factory"] = kwargs
        return "raw-model" if kwargs.get("return_raw_model") else "calculator"

    class FakeMACECalculator:
        def __init__(self, **kwargs):
            captured["calculator"] = kwargs

    calculators.mace_anicc = mace_anicc
    calculators.mace_mp = lambda **kwargs: object()
    calculators.mace_off = lambda **kwargs: object()
    calculators.mace_omol = lambda **kwargs: object()
    calculators.MACECalculator = FakeMACECalculator
    mace_module = ModuleType("mace.calculators.mace")
    mace_module.MACECalculator = FakeMACECalculator
    foundations = ModuleType("mace.calculators.foundations_models")
    foundations.mace_mp_urls = {}
    monkeypatch.setitem(sys.modules, "mace", package)
    monkeypatch.setitem(sys.modules, "mace.calculators", calculators)
    monkeypatch.setitem(sys.modules, "mace.calculators.mace", mace_module)
    monkeypatch.setitem(sys.modules, "mace.calculators.foundations_models", foundations)

    calculator = MACECalculator.__new__(MACECalculator)
    calculator.device_str = "cpu"
    calculator.default_dtype = dtype

    calculator._build_calc("anicc")

    if dtype == "float64":
        assert captured == {"factory": {"device": "cpu"}}
    else:
        assert captured == {
            "factory": {"device": "cpu", "return_raw_model": True},
            "calculator": {
                "models": "raw-model",
                "device": "cpu",
                "default_dtype": "float32",
            },
        }


@pytest.mark.parametrize("size", ["small", "medium", "large"])
def test_mace_off23_display_alias_maps_to_upstream_size(
    monkeypatch, size
) -> None:
    from pdb2reaction.backends.mace import MACECalculator

    package = ModuleType("mace")
    package.__path__ = []
    calculators = ModuleType("mace.calculators")
    captured = {}
    calculators.mace_anicc = lambda **kwargs: object()
    calculators.mace_mp = lambda **kwargs: object()
    calculators.mace_omol = lambda **kwargs: object()
    calculators.mace_off = lambda **kwargs: captured.update(kwargs) or object()
    mace_module = ModuleType("mace.calculators.mace")
    mace_module.MACECalculator = object
    foundations = ModuleType("mace.calculators.foundations_models")
    foundations.mace_mp_urls = {}
    monkeypatch.setitem(sys.modules, "mace", package)
    monkeypatch.setitem(sys.modules, "mace.calculators", calculators)
    monkeypatch.setitem(sys.modules, "mace.calculators.mace", mace_module)
    monkeypatch.setitem(sys.modules, "mace.calculators.foundations_models", foundations)

    calculator = MACECalculator.__new__(MACECalculator)
    calculator.device_str = "cpu"
    calculator.default_dtype = "float64"
    calculator._build_calc(f"MACE-OFF23_{size}")

    assert captured == {
        "model": size,
        "device": "cpu",
        "default_dtype": "float64",
    }


def test_mace_download_cache_is_url_specific_and_atomic(
    tmp_path: Path, monkeypatch
) -> None:
    from pdb2reaction.backends import mace

    published = []

    def fake_retrieve(url, destination):
        Path(destination).write_text(url, encoding="utf-8")

    real_replace = mace.os.replace

    def record_replace(source, destination):
        assert Path(source).exists()
        real_replace(source, destination)
        published.append(Path(destination))

    monkeypatch.setattr(mace.tempfile, "gettempdir", lambda: str(tmp_path))
    monkeypatch.setattr(mace.urllib.request, "urlretrieve", fake_retrieve)
    monkeypatch.setattr(mace.os, "replace", record_replace)

    first = Path(mace.MACECalculator._download_to_tmp("https://one.test/model.pt"))
    second = Path(mace.MACECalculator._download_to_tmp("https://two.test/model.pt"))

    assert first != second
    assert first.read_text(encoding="utf-8") == "https://one.test/model.pt"
    assert second.read_text(encoding="utf-8") == "https://two.test/model.pt"
    assert published == [first, second]
