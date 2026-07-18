"""Energy/force/Hessian consistency for harmonic distance restraints."""

from __future__ import annotations

import numpy as np
import pytest
import torch

from pysisyphus.constants import ANG2BOHR

from pdb2reaction.workflows.freq import _calc_full_hessian_torch
from pdb2reaction.workflows.opt import (
    _flatten_all_imag_modes_for_geom,
    _seed_rfo_initial_hessian,
)
from pdb2reaction.workflows.restraints import (
    HarmonicBiasCalculator,
    harmonic_pair_energy_forces_hessian,
)


def _pair_system():
    coords = np.array(
        [[0.0, 0.0, 0.0], [1.35 * ANG2BOHR, 0.21 * ANG2BOHR, 0.0]],
        dtype=float,
    )
    return coords, [(0, 1, 1.0)]


def test_harmonic_pair_hessian_matches_force_finite_difference() -> None:
    coords, pairs = _pair_system()
    k = 0.37
    _, _, analytical = harmonic_pair_energy_forces_hessian(
        coords, k, pairs, need_hessian=True
    )
    assert analytical is not None

    eps = 1.0e-5
    numerical = np.zeros_like(analytical)
    for column in range(coords.size):
        plus = coords.copy().reshape(-1)
        minus = coords.copy().reshape(-1)
        plus[column] += eps
        minus[column] -= eps
        _, force_plus, _ = harmonic_pair_energy_forces_hessian(
            plus, k, pairs, need_hessian=False
        )
        _, force_minus, _ = harmonic_pair_energy_forces_hessian(
            minus, k, pairs, need_hessian=False
        )
        numerical[:, column] = -(force_plus - force_minus) / (2.0 * eps)

    assert np.allclose(analytical, numerical, atol=2.0e-10, rtol=2.0e-9)
    assert np.allclose(analytical, analytical.T, atol=0.0, rtol=0.0)


def test_harmonic_pair_force_is_negative_energy_gradient() -> None:
    coords, pairs = _pair_system()
    k = 0.37
    _, force, _ = harmonic_pair_energy_forces_hessian(
        coords, k, pairs, need_hessian=False
    )

    eps = 1.0e-5
    gradient = np.zeros(coords.size, dtype=float)
    for column in range(coords.size):
        plus = coords.copy().reshape(-1)
        minus = coords.copy().reshape(-1)
        plus[column] += eps
        minus[column] -= eps
        e_plus, _, _ = harmonic_pair_energy_forces_hessian(
            plus, k, pairs, need_hessian=False
        )
        e_minus, _, _ = harmonic_pair_energy_forces_hessian(
            minus, k, pairs, need_hessian=False
        )
        gradient[column] = (e_plus - e_minus) / (2.0 * eps)

    assert np.allclose(force, -gradient, atol=2.0e-10, rtol=2.0e-9)


class _StaticBase:
    def __init__(self, result, *, freeze_atoms=()):
        self.result = result
        self.freeze_atoms = list(freeze_atoms)
        self.hessian_calls = 0

    def get_energy(self, elem, coords):
        return self.result

    def get_forces(self, elem, coords):
        return self.result

    def get_hessian(self, elem, coords):
        self.hessian_calls += 1
        return self.result


@pytest.mark.parametrize("use_torch", [False, True])
def test_restrained_hessian_preserves_full_representation_and_final_mask(
    use_torch: bool,
) -> None:
    coords, pairs = _pair_system()
    hessian = np.eye(6, dtype=np.float32)
    if use_torch:
        hessian = torch.as_tensor(hessian)
    base_result = {
        "energy": 1.0,
        "forces": np.ones(6, dtype=np.float32),
        "hessian": hessian,
    }
    base = _StaticBase(base_result, freeze_atoms=[0])
    wrapper = HarmonicBiasCalculator(base, k=2.0, pairs=pairs)

    result = wrapper.get_hessian(["H", "H"], coords)

    assert isinstance(result["hessian"], torch.Tensor) is use_torch
    assert result["hessian"].dtype == hessian.dtype
    hessian_np = (
        result["hessian"].detach().cpu().numpy()
        if use_torch
        else result["hessian"]
    )
    assert np.array_equal(hessian_np[:3, :], np.zeros((3, 6)))
    assert np.array_equal(hessian_np[:, :3], np.zeros((6, 3)))
    assert np.array_equal(result["forces"][:3], np.zeros(3))
    assert base_result["energy"] == 1.0
    assert np.array_equal(base_result["forces"], np.ones(6, dtype=np.float32))
    assert np.array_equal(
        hessian.detach().cpu().numpy() if use_torch else hessian,
        np.eye(6, dtype=np.float32),
    )


def test_restrained_hessian_maps_through_partial_metadata() -> None:
    coords, pairs = _pair_system()
    metadata = {
        "active_atoms": np.array([1]),
        "active_dofs": np.array([3, 4, 5]),
        "active_n_dof": 3,
        "full_n_dof": 6,
    }
    base_result = {
        "energy": 1.0,
        "forces": np.ones(6, dtype=np.float64),
        "hessian": torch.eye(3, dtype=torch.float64),
        "within_partial_hessian": metadata,
    }
    wrapper = HarmonicBiasCalculator(
        _StaticBase(base_result, freeze_atoms=[0]), k=2.0, pairs=pairs
    )

    result = wrapper.get_hessian(["H", "H"], coords)

    assert result["hessian"].shape == (3, 3)
    assert result["hessian"].dtype == torch.float64
    assert result["within_partial_hessian"] is metadata
    assert np.array_equal(result["forces"][:3], np.zeros(3))


def test_empty_restraints_match_base_values_and_clone_buffers() -> None:
    coords, _ = _pair_system()
    base_result = {
        "energy": 1.25,
        "forces": np.arange(6, dtype=float),
        "hessian": torch.arange(36, dtype=torch.float64).reshape(6, 6),
    }
    wrapper = HarmonicBiasCalculator(_StaticBase(base_result), k=2.0, pairs=[])

    energy_result = wrapper.get_energy(["H", "H"], coords)
    force_result = wrapper.get_forces(["H", "H"], coords)
    result = wrapper.get_hessian(["H", "H"], coords)

    assert energy_result["energy"] == base_result["energy"]
    assert np.array_equal(force_result["forces"], base_result["forces"])
    assert result["energy"] == base_result["energy"]
    assert np.array_equal(result["forces"], base_result["forces"])
    assert torch.equal(result["hessian"], base_result["hessian"])
    assert result["forces"] is not base_result["forces"]
    assert result["hessian"] is not base_result["hessian"]


class _CachedGeometry:
    atoms = ["H", "H"]
    cart_coords = np.zeros(6, dtype=float)

    def __init__(self):
        self._hessian = torch.full((6, 6), 99.0)
        self.results = None

    def set_results(self, results):
        self.results = results


def test_explicit_evaluator_bypasses_coordinate_only_hessian_cache() -> None:
    geometry = _CachedGeometry()
    expected = torch.eye(6, dtype=torch.float64) * 2.0
    evaluator = _StaticBase(
        {
            "energy": 0.0,
            "forces": np.zeros(6),
            "hessian": expected,
        }
    )

    result = _calc_full_hessian_torch(
        geometry,
        {},
        torch.device("cpu"),
        calculator=evaluator,
    )

    assert evaluator.hessian_calls == 1
    assert torch.equal(result, expected)
    assert geometry.results["hessian"] is expected


def test_restrained_rfo_seed_uses_exact_wrapper_and_never_reads_irc_cache(
    monkeypatch,
) -> None:
    class Geometry:
        cart_coords = np.zeros(6, dtype=float)
        cart_hessian = None

    geometry = Geometry()
    evaluator = object()
    expected = torch.eye(6, dtype=torch.float64) * 7.0
    seen = []

    def fake_hessian(geom, cfg, device, *, calculator=None, **kwargs):
        seen.append((geom, cfg, device, calculator))
        return expected

    def forbidden_cache(*args, **kwargs):
        raise AssertionError("restrained RFO must not read the IRC cache")

    monkeypatch.setattr(
        "pdb2reaction.workflows.opt._calc_full_hessian_torch", fake_hessian
    )
    monkeypatch.setattr(
        "pdb2reaction.io.hessian_cache.load", forbidden_cache
    )

    source = _seed_rfo_initial_hessian(
        geometry,
        {"device": "cpu", "backend": "sentinel"},
        evaluator,
        restraints_active=True,
    )

    assert source == "restrained"
    assert geometry.cart_hessian is expected
    assert seen[0][0] is geometry
    assert seen[0][1] == {"device": "cpu", "backend": "sentinel"}
    assert seen[0][2] == torch.device("cpu")
    assert seen[0][3] is evaluator


def test_flatten_energy_probes_use_exact_active_evaluator(monkeypatch) -> None:
    class Geometry:
        cart_coords = np.zeros(6, dtype=float)

        @property
        def coords(self):
            return self.cart_coords

        @coords.setter
        def coords(self, value):
            self.cart_coords = np.asarray(value, dtype=float).copy()

    geometry = Geometry()
    evaluator = object()
    seen = []

    def fake_energy(_geom, _kwargs, calc=None):
        seen.append(calc)
        return float(np.sum(_geom.cart_coords * _geom.cart_coords))

    monkeypatch.setattr(
        "pdb2reaction.workflows.opt._calc_energy", fake_energy
    )

    did_flatten = _flatten_all_imag_modes_for_geom(
        geometry,
        np.array([1.0, 1.0]),
        {"backend": "sentinel"},
        np.array([-100.0]),
        torch.tensor([[1.0, 0.0, 0.0, -1.0, 0.0, 0.0]]),
        5.0,
        0.1,
        calculator=evaluator,
    )

    assert did_flatten is True
    assert seen == [evaluator, evaluator, evaluator]
