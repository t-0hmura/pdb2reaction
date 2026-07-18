"""Final representation and constraint contract for solvent corrections."""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest
import torch

from pdb2reaction.backends.base import BackendError
from pdb2reaction.backends.solvent import SolventCorrectedCalculator
import pdb2reaction.backends.solvent as solvent_module


class _Base:
    charge = 0
    mult = 1
    device_str = "cpu"
    hessian_calc_mode = "FiniteDifference"
    return_partial_hessian = False
    hessian_double = True
    out_hess_torch = True
    print_timing = False

    def __init__(self, result, *, freeze_atoms=(0,)):
        self.result = result
        self.freeze_atoms = list(freeze_atoms)

    def get_energy(self, elem, coords):
        return self.result

    def get_forces(self, elem, coords):
        return self.result

    def get_hessian(self, elem, coords):
        return self.result


@pytest.fixture(autouse=True)
def _stub_xtb_module(monkeypatch):
    stub = SimpleNamespace(
        solvent_correction_enabled=lambda value: str(value).lower() != "none",
        normalize_solvent_model=lambda value: str(value),
    )
    monkeypatch.setattr(solvent_module, "_xtb_mod", stub)


@pytest.mark.parametrize("as_torch", [False, True])
def test_solvent_full_hessian_preserves_representation_and_masks_last(
    as_torch: bool,
) -> None:
    base_hessian = np.eye(6, dtype=np.float32)
    if as_torch:
        base_hessian = torch.as_tensor(base_hessian)
    base_result = {
        "energy": 1.5,
        "forces": np.zeros(6, dtype=np.float32),
        "hessian": base_hessian,
    }
    base = _Base(base_result)
    base.out_hess_torch = as_torch
    calc = SolventCorrectedCalculator(base, solvent="water")
    calc._solvent_delta = lambda *args, **kwargs: (
        2.0,
        np.ones((2, 3), dtype=float),
        np.ones((6, 6), dtype=float),
    )

    force_result = calc.get_forces(["H", "H"], np.zeros(6))
    hessian_result = calc.get_hessian(["H", "H"], np.zeros(6))

    assert calc.out_hess_torch is as_torch
    assert isinstance(hessian_result["hessian"], torch.Tensor) is as_torch
    assert hessian_result["hessian"].dtype == base_hessian.dtype
    hessian_np = (
        hessian_result["hessian"].detach().cpu().numpy()
        if as_torch
        else hessian_result["hessian"]
    )
    assert np.array_equal(force_result["forces"][:3], np.zeros(3))
    assert np.array_equal(hessian_result["forces"][:3], np.zeros(3))
    assert np.array_equal(hessian_np[:3, :], np.zeros((3, 6)))
    assert np.array_equal(hessian_np[:, :3], np.zeros((6, 3)))
    assert np.any(hessian_np[3:, 3:] != np.eye(6, dtype=np.float32)[3:, 3:])
    assert base_result["energy"] == 1.5
    assert np.array_equal(base_result["forces"], np.zeros(6, dtype=np.float32))
    original_hessian = (
        base_hessian.detach().cpu().numpy() if as_torch else base_hessian
    )
    assert np.array_equal(original_hessian, np.eye(6, dtype=np.float32))


def test_solvent_partial_hessian_uses_existing_active_mapping() -> None:
    metadata = {
        "active_atoms": np.array([1]),
        "active_dofs": np.array([3, 4, 5]),
        "active_n_dof": 3,
        "full_n_dof": 6,
    }
    base_result = {
        "energy": 0.5,
        "forces": np.zeros(6),
        "hessian": torch.eye(3, dtype=torch.float64),
        "within_partial_hessian": metadata,
    }
    base = _Base(base_result)
    base.return_partial_hessian = True
    calc = SolventCorrectedCalculator(base, solvent="water")
    calc._solvent_delta = lambda *args, **kwargs: (
        0.0,
        np.ones((2, 3), dtype=float),
        np.arange(36, dtype=float).reshape(6, 6),
    )

    result = calc.get_hessian(["H", "H"], np.zeros(6))

    assert result["hessian"].shape == (3, 3)
    assert result["hessian"].dtype == torch.float64
    assert result["within_partial_hessian"] is metadata
    assert np.array_equal(result["forces"][:3], np.zeros(3))


def test_disabled_solvent_preserves_values_and_clones_buffers() -> None:
    metadata = {"active_dofs": np.arange(6)}
    base_result = {
        "energy": 0.25,
        "forces": np.arange(6, dtype=float),
        "hessian": torch.arange(36, dtype=torch.float64).reshape(6, 6),
        "within_partial_hessian": metadata,
    }
    calc = SolventCorrectedCalculator(
        _Base(base_result, freeze_atoms=()), solvent="none"
    )

    result = calc.get_hessian(["H", "H"], np.zeros(6))

    assert result["energy"] == base_result["energy"]
    assert np.array_equal(result["forces"], base_result["forces"])
    assert torch.equal(result["hessian"], base_result["hessian"])
    assert result["forces"] is not base_result["forces"]
    assert result["hessian"] is not base_result["hessian"]
    assert result["within_partial_hessian"] is metadata


def test_solvent_energy_is_added_once_per_public_evaluation() -> None:
    base_result = {
        "energy": 4.0,
        "forces": np.zeros(6),
        "hessian": np.eye(6),
    }
    calc = SolventCorrectedCalculator(
        _Base(base_result, freeze_atoms=()), solvent="water"
    )
    calls = []

    def correction(*args, **kwargs):
        calls.append((kwargs.get("need_forces"), kwargs.get("need_hessian")))
        return 2.0, np.zeros((2, 3)), np.zeros((6, 6))

    calc._solvent_delta = correction

    energy = calc.get_energy(["H", "H"], np.zeros(6))["energy"]
    force_energy = calc.get_forces(["H", "H"], np.zeros(6))["energy"]
    hessian_energy = calc.get_hessian(["H", "H"], np.zeros(6))["energy"]

    assert energy == force_energy == hessian_energy
    assert energy > base_result["energy"]
    assert base_result["energy"] == 4.0
    assert calls == [(False, None), (True, None), (True, True)]


def test_base_hessian_error_propagates_without_solvent_evaluation() -> None:
    class FailingBase(_Base):
        def get_hessian(self, elem, coords):
            raise BackendError("analytical unavailable")

    calc = SolventCorrectedCalculator(
        FailingBase({"energy": 0.0}, freeze_atoms=()), solvent="water"
    )
    called = False

    def correction(*args, **kwargs):
        nonlocal called
        called = True
        return 0.0, None, None

    calc._solvent_delta = correction

    with pytest.raises(BackendError, match="analytical unavailable"):
        calc.get_hessian(["H"], np.zeros(3))
    assert called is False
