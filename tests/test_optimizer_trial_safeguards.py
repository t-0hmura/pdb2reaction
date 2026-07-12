"""Regression tests for shared rollback and first-order-saddle safeguards."""

from __future__ import annotations

import numpy as np
import pytest
import torch

from pysisyphus.Geometry import Geometry
from pysisyphus.calculators.Calculator import Calculator
from pysisyphus.optimizers.LBFGS import LBFGS
from pysisyphus.tsoptimizers.RSIRFOptimizer import RSIRFOptimizer
from pysisyphus.tsoptimizers.RSPRFOptimizer import RSPRFOptimizer
from pysisyphus.tsoptimizers.TRIM import TRIM


class _QuadraticCalculator(Calculator):
    def __init__(self, out_dir):
        super().__init__(out_dir=out_dir, check_mem=False)

    @staticmethod
    def _results(coords):
        coords = np.asarray(coords, dtype=float)
        return {
            "energy": float(coords @ coords),
            "forces": -2.0 * coords,
            "hessian": 2.0 * np.eye(len(coords)),
        }

    def get_energy(self, atoms, coords, **prepare_kwargs):
        return self._results(coords)

    def get_forces(self, atoms, coords, **prepare_kwargs):
        return self._results(coords)

    def get_hessian(self, atoms, coords, **prepare_kwargs):
        return self._results(coords)


class _DoubleWellCalculator(Calculator):
    """One saddle coordinate plus two stiff minimizing coordinates."""

    def __init__(self, out_dir):
        super().__init__(out_dir=out_dir, check_mem=False)

    @staticmethod
    def _results(coords):
        x, y, z = np.asarray(coords, dtype=float)
        return {
            "energy": float((x * x - 1.0) ** 2 + 5.0 * y * y + 5.0 * z * z),
            "forces": np.array([-4.0 * x * (x * x - 1.0), -10.0 * y, -10.0 * z]),
            "hessian": np.diag([12.0 * x * x - 4.0, 10.0, 10.0]),
        }

    def get_energy(self, atoms, coords, **prepare_kwargs):
        return self._results(coords)

    def get_forces(self, atoms, coords, **prepare_kwargs):
        return self._results(coords)

    def get_hessian(self, atoms, coords, **prepare_kwargs):
        return self._results(coords)


class _TranslationalArtifactCalculator(Calculator):
    """Stationary diatomic with a spurious negative translation only."""

    def __init__(self, out_dir):
        super().__init__(out_dir=out_dir, check_mem=False)

    @staticmethod
    def _results(coords):
        coords = np.asarray(coords, dtype=float)
        translation_x = np.array([1.0, 0.0, 0.0, 1.0, 0.0, 0.0])
        translation_x /= np.linalg.norm(translation_x)
        hessian = 2.0 * np.eye(6) - 3.0 * np.outer(
            translation_x, translation_x
        )
        return {
            "energy": 0.0,
            "forces": np.zeros_like(coords),
            "hessian": hessian,
        }

    def get_energy(self, atoms, coords, **prepare_kwargs):
        return self._results(coords)

    def get_forces(self, atoms, coords, **prepare_kwargs):
        return self._results(coords)

    def get_hessian(self, atoms, coords, **prepare_kwargs):
        return self._results(coords)


class _SpuriousTranslationDoubleWellCalculator(Calculator):
    """Diatomic double well whose most-negative raw root is always TR noise."""

    well_position = 0.1

    def __init__(self, out_dir):
        super().__init__(out_dir=out_dir, check_mem=False)

    @classmethod
    def _results(cls, coords):
        coords = np.asarray(coords, dtype=float)
        stretch = np.array([-1.0, 0.0, 0.0, 1.0, 0.0, 0.0])
        stretch /= np.linalg.norm(stretch)
        translation = np.array([1.0, 0.0, 0.0, 1.0, 0.0, 0.0])
        translation /= np.linalg.norm(translation)
        q = float(stretch @ coords)
        a = cls.well_position
        d_energy_dq = 4.0 * q * (q * q - a * a)
        d2_energy_dq2 = 12.0 * q * q - 4.0 * a * a

        gradient = d_energy_dq * stretch
        # Keep the four transverse directions stiff. The x translation is a
        # deliberately inconsistent negative numerical artifact.
        transverse = np.diag([0.0, 2.0, 2.0, 0.0, 2.0, 2.0])
        hessian = (
            d2_energy_dq2 * np.outer(stretch, stretch)
            - np.outer(translation, translation)
            + transverse
        )
        return {
            "energy": float((q * q - a * a) ** 2),
            "forces": -gradient,
            "hessian": hessian,
        }

    def get_energy(self, atoms, coords, **prepare_kwargs):
        return self._results(coords)

    def get_forces(self, atoms, coords, **prepare_kwargs):
        return self._results(coords)

    def get_hessian(self, atoms, coords, **prepare_kwargs):
        return self._results(coords)


class _TorchSpuriousTranslationDoubleWellCalculator(
    _SpuriousTranslationDoubleWellCalculator
):
    def get_hessian(self, atoms, coords, **prepare_kwargs):
        results = self._results(coords)
        results["hessian"] = torch.as_tensor(
            results["hessian"], dtype=torch.float64
        )
        return results


def test_lbfgs_rejects_uphill_trial_through_shared_rollback(tmp_path) -> None:
    geom = Geometry(["H"], np.array([1.0, 0.0, 0.0]), coord_type="cart")
    geom.set_calculator(_QuadraticCalculator(tmp_path))
    opt = LBFGS(
        geom,
        max_step=0.1,
        line_search=False,
        out_dir=tmp_path,
    )
    assert opt.reject_uphill is True

    accepted = geom.coords.copy()
    accepted_cart = geom.cart_coords.copy()
    accepted_energy = geom.energy
    accepted_force = geom.forces.copy()
    rejected_step = np.array([1.0, 0.0, 0.0])
    geom.coords = accepted + rejected_step
    opt.coords = [accepted.copy(), geom.coords.copy()]
    opt.cart_coords = [accepted_cart.copy(), geom.cart_coords.copy()]
    opt.energies = [accepted_energy]
    opt.forces = [accepted_force]
    opt.steps = [rejected_step.copy()]
    opt.image_inds = [[0], [0]]
    opt.image_nums = [1, 1]
    opt.cur_cycle = 1

    retry_step = opt.optimize()

    np.testing.assert_allclose(geom.coords, accepted)
    np.testing.assert_allclose(geom.cart_coords, accepted_cart)
    assert len(opt.coords) == len(opt.cart_coords) == 1
    assert len(opt.energies) == len(opt.forces) == 1
    assert opt.steps == []
    assert opt.rejected_uphill_steps == 1
    assert opt._trial_max_step == 0.05
    assert np.abs(retry_step).max() <= 0.05


def _ts_optimizer(tmp_path, x, **kwargs):
    geom = Geometry(["H"], np.array([x, 0.0, 0.0]), coord_type="cart")
    geom.set_calculator(_DoubleWellCalculator(tmp_path))
    opt_kwargs = dict(
        hessian_init="calc",
        hessian_recalc=1,
        trust_radius=0.1,
        trust_min=1e-5,
        trust_max=0.1,
        thresh="baker",
        check_eigval_structure=True,
        out_dir=tmp_path,
    )
    opt_kwargs.update(kwargs)
    opt = RSPRFOptimizer(geom, **opt_kwargs)
    return geom, opt


def test_ts_mode_loss_rejects_trial_and_restores_hessian(tmp_path) -> None:
    geom, opt = _ts_optimizer(tmp_path, 0.4)
    opt.prepare_opt()
    accepted = geom.coords.copy()
    accepted_cart = geom.cart_coords.copy()
    accepted_energy = geom.energy
    accepted_force = geom.forces.copy()
    accepted_hessian = opt.H.copy()

    # x=0.4 has negative curvature; x=0.8 is already inside the convex basin.
    trial_step = np.array([0.4, 0.0, 0.0])
    geom.coords = accepted + trial_step
    opt.coords = [accepted.copy(), geom.coords.copy()]
    opt.cart_coords = [accepted_cart.copy(), geom.cart_coords.copy()]
    opt.energies = [accepted_energy]
    opt.forces = [accepted_force]
    opt.steps = [trial_step.copy()]
    opt.predicted_energy_changes = [-0.1]
    opt.image_inds = [[0], [0]]
    opt.image_nums = [1, 1]
    opt.cur_cycle = 1
    opt.hessian_recalc_in = 0

    _, _, _, eigvals, _, resetted = opt.housekeeping()

    np.testing.assert_allclose(geom.coords, accepted)
    np.testing.assert_allclose(opt.H, accepted_hessian)
    assert eigvals[0] < 0.0
    assert resetted is True
    assert opt.rejected_mode_loss_steps == 1
    assert opt.trust_radius == 0.025
    assert opt.steps == []
    assert opt.predicted_energy_changes == []


def test_energy_plateau_cannot_bypass_ts_curvature_requirement(tmp_path) -> None:
    geom, opt = _ts_optimizer(tmp_path, 1.0, energy_plateau_window=5)
    opt.cur_cycle = 0
    opt.last_cycle = 0
    opt.forces = [np.zeros(3)]
    opt.energies = [0.0] * 5
    opt.ts_mode_eigvals = np.array([1.0])

    converged, conv_info = opt.check_convergence(np.zeros(3))

    assert converged is False
    assert conv_info.desired_eigval_structure is False


def test_verified_phva_mode_throttles_repeated_exact_hessians(tmp_path) -> None:
    _, opt = _ts_optimizer(tmp_path, 0.0, energy_plateau=False)
    opt._physical_ts_mode = np.array([1.0, 0.0, 0.0])
    opt._last_exact_n_imaginary = 1
    opt.forces = [np.zeros(3)]
    opt.energies = [0.0, 0.0]
    opt.steps = [np.array([0.01, 0.0, 0.0])]

    assert opt._near_terminal_without_eigval_check() is False

    opt.steps[-1] = np.zeros(3)
    assert opt._near_terminal_without_eigval_check() is True


def test_exact_hessian_recovery_displaces_back_along_incoming_mode(tmp_path) -> None:
    geom, opt = _ts_optimizer(
        tmp_path,
        1.0,
        hessian_recalc=None,
        energy_plateau_window=5,
    )
    # Simulate a stale quasi-Newton model that still claims negative curvature
    # at the true local minimum.  The exact Hessian must override it.
    opt.H = np.diag([-1.0, 10.0, 10.0])
    opt.negative_mode_seen = True
    opt.coords = [geom.coords.copy()]
    opt.cart_coords = [geom.cart_coords.copy()]
    opt.energies = [0.0] * 4
    opt.forces = [geom.forces.copy()] * 4
    opt.steps = [np.array([0.1, 0.0, 0.0])]
    opt.image_inds = [[0]]
    opt.image_nums = [1]
    opt.cur_cycle = 4

    step = opt.optimize()

    assert opt.exact_saddle_checks == 1
    assert opt._saddle_recovery_active is True
    assert opt.saddle_recovery_steps == 1
    assert opt.ts_mode_eigvals[0] > 0.0
    # Reverse the incoming motion toward x=+1, back toward the nearby saddle x=0.
    assert step[0] == -opt.saddle_recovery_step
    np.testing.assert_allclose(step[1:], 0.0)


def test_initial_near_stationary_minimum_gets_finite_recovery_step(tmp_path) -> None:
    geom, opt = _ts_optimizer(tmp_path, 1.0, hessian_recalc=None)
    opt.prepare_opt()
    opt.coords = [geom.coords.copy()]
    opt.cart_coords = [geom.cart_coords.copy()]
    opt.image_inds = [[0]]
    opt.image_nums = [1]
    opt.cur_cycle = opt.last_cycle

    step = opt.optimize()

    assert opt.exact_saddle_checks == 1
    assert opt._saddle_recovery_active is True
    assert np.linalg.norm(step) == opt.saddle_recovery_step

    # A minimizing proposal on the next cycle would point back toward the
    # minimum. Recovery must retain its chosen uphill direction instead of
    # oscillating across the minimum.
    next_step = opt.apply_saddle_recovery_step(-step)
    np.testing.assert_allclose(next_step, step)


def test_exact_phva_ignores_negative_translation_and_uses_vibration(tmp_path) -> None:
    coords = np.array([-0.7, 0.0, 0.0, 0.7, 0.0, 0.0])
    geom = Geometry(["H", "H"], coords, coord_type="cart")
    geom.set_calculator(_TranslationalArtifactCalculator(tmp_path))
    opt = RSPRFOptimizer(
        geom,
        hessian_init="calc",
        hessian_recalc=None,
        trust_radius=0.1,
        trust_max=0.1,
        thresh="baker",
        check_eigval_structure=True,
        out_dir=tmp_path,
    )
    opt.prepare_opt()
    opt.coords = [geom.coords.copy()]
    opt.cart_coords = [geom.cart_coords.copy()]
    opt.image_inds = [[0]]
    opt.image_nums = [1]
    opt.cur_cycle = opt.last_cycle

    step = opt.optimize()

    # The unprojected Hessian still contains the artificial negative root.
    assert opt.ts_mode_eigvals[0] < 0.0
    # PHVA removes that translation and correctly identifies n_imag=0.
    assert opt._last_exact_n_imaginary == 0
    assert opt._saddle_recovery_active is True

    translation_x = np.array([1.0, 0.0, 0.0, 1.0, 0.0, 0.0])
    translation_x /= np.linalg.norm(translation_x)
    bond_stretch = np.array([-1.0, 0.0, 0.0, 1.0, 0.0, 0.0])
    bond_stretch /= np.linalg.norm(bond_stretch)
    assert abs(float(step @ translation_x)) < 1e-12
    assert abs(float(step @ bond_stretch)) == opt.saddle_recovery_step


@pytest.mark.parametrize(
    ("optimizer_cls", "calculator_cls"),
    [
        (RSPRFOptimizer, _SpuriousTranslationDoubleWellCalculator),
        (RSIRFOptimizer, _TorchSpuriousTranslationDoubleWellCalculator),
        (TRIM, _TorchSpuriousTranslationDoubleWellCalculator),
    ],
)
def test_recovery_tracks_physical_mode_past_more_negative_tr_artifact(
    tmp_path, optimizer_cls, calculator_cls
) -> None:
    a = _SpuriousTranslationDoubleWellCalculator.well_position
    coords = np.array([-a / np.sqrt(2.0), 0.0, 0.0, a / np.sqrt(2.0), 0.0, 0.0])
    geom = Geometry(["H", "H"], coords, coord_type="cart")
    geom.set_calculator(calculator_cls(tmp_path))
    opt = optimizer_cls(
        geom,
        hessian_init="calc",
        hessian_recalc=1,
        trust_radius=0.05,
        trust_update=False,
        trust_min=1e-5,
        trust_max=0.05,
        max_cycles=40,
        thresh="baker",
        check_eigval_structure=True,
        out_dir=tmp_path,
    )

    opt.run()

    stretch = np.array([-1.0, 0.0, 0.0, 1.0, 0.0, 0.0])
    stretch /= np.linalg.norm(stretch)
    final_q = float(stretch @ geom.coords)
    assert opt.is_converged is True
    assert abs(final_q) < 1e-3
    assert opt.saddle_recovery_steps > 0
    assert opt._last_exact_n_imaginary == 1
    # Root 0 remains the artificial translation (-1); the physical reaction
    # vibration is root 1 and must be followed after recovery.
    assert int(opt.roots[0]) == 1
