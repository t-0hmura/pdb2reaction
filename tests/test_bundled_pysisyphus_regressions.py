"""Focused regressions for bundled pysisyphus correctness fixes."""

from __future__ import annotations

import json
import inspect
from types import SimpleNamespace
import weakref

import numpy as np
import pytest
import torch

from pysisyphus.intcoords.BondedFragment import BondedFragment
from pysisyphus.intcoords.Cartesian import CartesianX, CartesianY
from pysisyphus.intcoords.LinearDisplacement import LinearDisplacement
from pysisyphus.intcoords.PrimTypes import PrimTypes
from pysisyphus.intcoords.RedundantCoords import RedundantCoords
from pysisyphus.intcoords.Torsion import Torsion
from pysisyphus.intcoords.exceptions import NeedNewInternalsException
from pysisyphus.intcoords.update import transform_int_step, update_internals
from pysisyphus.io.qcschema import geom_from_qcschema
from pysisyphus.Geometry import Geometry
from pysisyphus.calculators.Dimer import Dimer
from pysisyphus.irc.EulerPC import EulerPC
from pysisyphus.irc.DWI import DWI
from pysisyphus.irc.IRC import IRC
from pysisyphus.irc.Instanton import Instanton
from pysisyphus.helpers import geom_loader
from pysisyphus.helpers_pure import molecular_volume
from pysisyphus.linalg import quaternion_to_rot_mat
from pysisyphus.cos.ChainOfStates import ChainOfStates
from pysisyphus.cos.GrowingChainOfStates import GrowingChainOfStates
from pysisyphus.cos.GrowingString import GrowingString
from pysisyphus.optimizers.Optimizer import Optimizer
from pysisyphus.optimizers.HessianOptimizer import HessianOptimizer
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.optimizers.StringOptimizer import StringOptimizer
from pysisyphus.optimizers.closures import bfgs_multiply
import pysisyphus.optimizers.gdiis as gdiis_module
from pysisyphus.optimizers.hessian_updates import bofill_update
from pysisyphus.tsoptimizers.RSPRFOptimizer import RSPRFOptimizer
from pysisyphus.tsoptimizers.TSHessianOptimizer import TSHessianOptimizer
from thermoanalysis.constants import AMU2KG, C, KB, PLANCK, R
from thermoanalysis.thermo import (
    chai_head_gordon_weights,
    qrrho_vibrational_part_func,
    vibrational_heat_capacity,
    vibrational_part_funcs,
)


def test_dimer_bias_force_norm_axis_error_keeps_translation_step(
    monkeypatch,
) -> None:
    dimer = object.__new__(Dimer)
    dimer.freeze_atoms = np.empty(0, dtype=int)
    dimer.rigid_basis = None
    dimer.rigid_basis_getter = None
    dimer.rotation_remove_trans = False
    dimer._coords0 = None
    dimer._energy0 = None
    dimer._f0 = None
    dimer._f1 = None
    dimer.force_evals = 0
    dimer.N = np.array([1.0, 0.0, 0.0])
    dimer.rotation_disable = True
    dimer.rotation_disable_pos_curv = True
    dimer.write_orientations = False
    dimer.length = 0.1
    dimer.bias_translation = True
    dimer.curvature = lambda *_args: 1.0
    dimer.calculator = SimpleNamespace(
        get_energy=lambda _atoms, _coords: {"energy": 0.0},
        get_forces=lambda _atoms, _coords: {
            "energy": 0.0,
            "forces": np.array([1.0, 0.0, 0.0]),
        },
    )
    dimer.gaussians = []
    dimer.get_gaussian_energies = lambda _coords: 0.0
    dimer.get_gaussian_forces = lambda _coords, sum_: np.zeros(3)
    dimer.trans_force_f_perp = True
    dimer.calc_counter = 0
    messages = []
    dimer.log = lambda message="": messages.append(message)
    dimer.make_fn = lambda name: name
    monkeypatch.setattr(np, "savetxt", lambda *args, **kwargs: None)

    result = dimer.get_forces(["H"], np.zeros(3))

    np.testing.assert_allclose(result["forces"], [-1.0, 0.0, 0.0])
    assert "Skipping calculation of norm(bias_forces)" in messages


def test_cartesian_geometry_rejects_coordinate_kwargs_without_exit_zero() -> None:
    with pytest.raises(ValueError, match="coord_kwargs were given"):
        Geometry(
            atoms=["He"],
            coords=np.zeros(3),
            coord_type="cart",
            coord_kwargs={"define_prims": []},
        )


def _qcschema_payload() -> dict:
    return {
        "molecule": {
            "symbols": ["H", "H"],
            "geometry": [0.0, 0.0, 0.0, 0.0, 0.0, 1.4],
        },
        "comment": "x" * 300,
    }


def test_geometry_hessian_ownership_transfer_keeps_other_results() -> None:
    geom = Geometry(("H",), np.zeros(3), coord_type="cart")
    hessian = torch.eye(3, dtype=torch.float64)
    geom.cart_hessian = hessian
    geom.results = {"energy": -1.0, "hessian": hessian}

    taken = geom.take_hessian()

    assert taken is hessian
    assert geom._hessian is None
    assert geom.results == {"energy": -1.0}


def test_generic_hessian_recalculation_releases_all_old_dense_aliases() -> None:
    old_hessian = np.eye(3)
    old_ref = weakref.ref(old_hessian)

    class ReplacementGeometry:
        @property
        def hessian(self):
            assert old_ref() is None
            return 2.0 * np.eye(3)

    optimizer = SimpleNamespace(
        forces=[np.ones(3)],
        adapt_norm=None,
        hessian_recalc_adapt=None,
        hessian_recalc_in=0,
        hessian_recalc=5,
        H=old_hessian,
        cur_H=old_hessian,
        hessian_xtb=False,
        geometry=ReplacementGeometry(),
        using_active_dofs=False,
        cur_cycle=1,
        log=lambda _message: None,
    )
    del old_hessian

    HessianOptimizer.update_hessian(optimizer)

    np.testing.assert_allclose(optimizer.H, 2.0 * np.eye(3))
    assert optimizer.cur_H is None


def test_qcschema_accepts_json_text_and_file_paths(tmp_path) -> None:
    payload = _qcschema_payload()
    text = json.dumps(payload)
    from_text = geom_from_qcschema(text)
    path = tmp_path / "molecule.json"
    path.write_text(text)
    from_path = geom_from_qcschema(path)

    assert from_text.atoms == from_path.atoms == ("H", "H")
    np.testing.assert_allclose(from_text.coords, from_path.coords)


def test_bonded_fragment_gradient_is_local_to_bond_endpoints() -> None:
    coords = np.array([
        [4.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 2.0, 0.0],
    ])
    value, gradient = BondedFragment._calculate(
        coords, (0, 1), gradient=True, bond_indices=(1, 2),
    )

    assert value == pytest.approx(2.0)
    np.testing.assert_allclose(
        gradient.reshape(-1, 3),
        [[0.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 1.0, 0.0]],
    )


@pytest.mark.parametrize("z", [-0.5, 0.0, 0.5])
def test_torsion_jacobian_matches_finite_difference(z) -> None:
    coords = np.array([
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [1.0, 1.0, z],
    ])
    indices = [0, 1, 2, 3]
    analytic = Torsion._jacobian(coords, indices).reshape(coords.size, coords.size)
    numeric = np.zeros_like(analytic)
    step = 1.0e-5
    for column in range(coords.size):
        plus = coords.copy().reshape(-1)
        minus = coords.copy().reshape(-1)
        plus[column] += step
        minus[column] -= step
        plus_gradient = Torsion._calculate(
            plus.reshape(-1, 3), indices, gradient=True,
        )[1]
        minus_gradient = Torsion._calculate(
            minus.reshape(-1, 3), indices, gradient=True,
        )[1]
        numeric[:, column] = (plus_gradient - minus_gradient) / (2 * step)
    np.testing.assert_allclose(analytic, numeric, atol=2.0e-6, rtol=2.0e-6)


def test_complement_linear_displacement_derivatives_match_scalar() -> None:
    coords = np.array(
        [
            [-1.1, 0.2, 0.1],
            [0.0, 0.0, 0.0],
            [1.3, 0.4, -0.2],
        ],
        dtype=float,
    )
    indices = [0, 1, 2]
    cross_vec = LinearDisplacement._get_cross_vec(coords, indices)
    _, analytic_gradient = LinearDisplacement._calculate(
        coords,
        indices,
        gradient=True,
        complement=True,
        cross_vec=cross_vec.copy(),
    )
    numeric_gradient = np.zeros(9, dtype=float)
    scalar_step = 1.0e-5
    for dof in range(coords.size):
        plus = coords.copy().reshape(-1)
        minus = coords.copy().reshape(-1)
        plus[dof] += scalar_step
        minus[dof] -= scalar_step
        value_plus = LinearDisplacement._calculate(
            plus.reshape(-1, 3),
            indices,
            complement=True,
            cross_vec=cross_vec.copy(),
        )
        value_minus = LinearDisplacement._calculate(
            minus.reshape(-1, 3),
            indices,
            complement=True,
            cross_vec=cross_vec.copy(),
        )
        numeric_gradient[dof] = (value_plus - value_minus) / (2 * scalar_step)
    np.testing.assert_allclose(
        analytic_gradient, numeric_gradient, atol=2.0e-8, rtol=2.0e-8,
    )

    analytic = LinearDisplacement._jacobian(
        coords,
        indices,
        complement=True,
        cross_vec=cross_vec.copy(),
    ).reshape(9, 9)
    numeric = np.zeros((9, 9), dtype=float)
    step = 1.0e-5
    for column in range(coords.size):
        plus = coords.copy().reshape(-1)
        minus = coords.copy().reshape(-1)
        plus[column] += step
        minus[column] -= step
        plus_gradient = LinearDisplacement._calculate(
            plus.reshape(-1, 3),
            indices,
            gradient=True,
            complement=True,
            cross_vec=cross_vec.copy(),
        )[1]
        minus_gradient = LinearDisplacement._calculate(
            minus.reshape(-1, 3),
            indices,
            gradient=True,
            complement=True,
            cross_vec=cross_vec.copy(),
        )[1]
        numeric[:, column] = (plus_gradient - minus_gradient) / (2 * step)
    np.testing.assert_allclose(analytic, numeric, atol=3.0e-6, rtol=3.0e-6)


def test_invalid_dihedral_and_bend_indices_are_unioned(monkeypatch) -> None:
    class FakePrimitive:
        def __init__(self, indices):
            self.indices = list(indices)

        def calculate(self, coords3d, gradient=False):
            return 0.0, np.zeros(coords3d.size)

    primitives = [
        FakePrimitive((0, 1, 2)),
        FakePrimitive((0,)),
        FakePrimitive((0, 1, 2, 3)),
    ]
    monkeypatch.setattr(
        "pysisyphus.intcoords.update.dihedral_valid", lambda *_: False,
    )
    monkeypatch.setattr(
        "pysisyphus.intcoords.update.bend_valid", lambda *_: True,
    )

    with pytest.raises(NeedNewInternalsException) as caught:
        update_internals(
            np.zeros((4, 3)),
            np.zeros(3),
            primitives,
            dihedral_inds=[2],
            rotation_inds=[],
            bend_inds=[0],
            check_dihedrals=True,
            check_bends=True,
        )
    assert caught.value.invalid_inds == [2]


def test_all_constraints_must_be_satisfied_after_backtransform() -> None:
    primitives = [CartesianX([0]), CartesianY([0])]
    _, cart_step, failed = transform_int_step(
        int_step=np.array([0.0, 1.0]),
        old_cart_coords=np.zeros(3),
        cur_internals=np.zeros(2),
        Bt_inv_prim=np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),
        primitives=primitives,
        dihedral_inds=[],
        rotation_inds=[],
        bend_inds=[],
        constrained_inds=[0, 1],
        update_constraints=False,
    )
    assert not failed
    np.testing.assert_allclose(cart_step, 0.0, atol=1.0e-12)


def test_redundant_coords_string_uses_existing_index_properties() -> None:
    coords = RedundantCoords.__new__(RedundantCoords)
    coords._bond_inds = [0]
    coords._bend_inds = [1, 2]
    coords._dihedral_inds = [3]
    assert str(coords) == "RedundantCoords(1 bonds, 2 bends, 1 dihedrals)"


def test_instanton_hessian_and_analytical_flag_follow_images() -> None:
    instanton = Instanton.__new__(Instanton)
    instanton.images = [
        SimpleNamespace(cart_hessian=np.eye(3), is_analytical_2d=True),
        SimpleNamespace(cart_hessian=2.0 * np.eye(3), is_analytical_2d=True),
    ]
    np.testing.assert_allclose(
        instanton.cart_hessian,
        np.diag([1.0, 1.0, 1.0, 2.0, 2.0, 2.0]),
    )
    assert instanton.is_analytical_2d is True


def test_quaternion_rotation_matrix_is_orthogonal() -> None:
    quaternion = np.array([0.5, -0.5, 0.5, 0.5])
    rotation = quaternion_to_rot_mat(quaternion)
    np.testing.assert_allclose(rotation.T @ rotation, np.eye(3), atol=1.0e-12)
    assert np.linalg.det(rotation) == pytest.approx(1.0)


def test_qrrho_partition_function_uses_frequency_in_effective_inertia() -> None:
    temperature = 298.15
    frequencies = np.array([20.0, 50.0, 1600.0]) * 100.0 * C
    q_vib, q_vib_v0 = qrrho_vibrational_part_func(
        temperature, frequencies, I_mean=10.0, cutoff=100.0, alpha=4,
    )
    harmonic, harmonic_v0 = vibrational_part_funcs(temperature, frequencies)
    weights = chai_head_gordon_weights(frequencies, 100.0, 4)
    mu = PLANCK / (8.0 * np.pi**2 * frequencies)
    inertia = 10.0 * 1.0e-20 * AMU2KG
    effective_inertia = mu * inertia / (mu + inertia)
    free_rotor = np.sqrt(
        8.0 * np.pi**3 * effective_inertia * KB * temperature / PLANCK**2
    )
    expected = np.exp(
        np.sum(weights * np.log(harmonic) + (1.0 - weights) * np.log(free_rotor))
    )
    expected_v0 = np.exp(
        np.sum(
            weights * np.log(harmonic_v0)
            + (1.0 - weights) * np.log(free_rotor)
        )
    )
    assert q_vib == pytest.approx(expected, rel=1.0e-13)
    assert q_vib_v0 == pytest.approx(expected_v0, rel=1.0e-13)


def test_vibrational_heat_capacity_is_stable_in_both_limits() -> None:
    high_frequency = np.array([3000.0]) * 100.0 * C
    assert vibrational_heat_capacity(1.0, high_frequency) == pytest.approx(0.0)
    low_frequency = np.array([1.0e-9])
    assert vibrational_heat_capacity(298.15, low_frequency) == pytest.approx(R)


def test_regularized_lbfgs_scales_the_secant_correction() -> None:
    result = bfgs_multiply(
        [np.array([1.0, -1.0])],
        [np.array([-2.0, 2.0])],
        np.array([0.3, -0.4]),
        gamma_mult=False,
        mu_reg=0.1,
    )
    np.testing.assert_allclose(result, [3.45, -3.55])


def test_gediis_returns_none_when_the_inner_solve_fails(monkeypatch) -> None:
    monkeypatch.setattr(
        gdiis_module,
        "minimize",
        lambda *_args, **_kwargs: SimpleNamespace(success=False),
    )
    result = gdiis_module.gediis(
        np.array([[0.0, 0.0], [0.1, 0.0]]),
        np.array([0.0, 0.1]),
        np.array([[0.1, 0.0], [0.2, 0.0]]),
    )
    assert result is None


def test_geom_loader_dispatches_trajectory_suffix_before_xyz(tmp_path) -> None:
    trajectory = tmp_path / "two_trj.xyz"
    trajectory.write_text(
        "1\nfirst\nH 0 0 0\n1\nsecond\nH 0 0 1\n",
        encoding="utf-8",
    )
    geometries = geom_loader(trajectory, iterable=True)
    assert len(geometries) == 2
    assert geom_loader(f"{trajectory}[1]").comment == "second"


@pytest.mark.parametrize(
    ("kind", "eigvals", "gradient"),
    [
        ("max", np.array([-0.2]), np.array([0.1])),
        ("min", np.array([0.3]), np.array([0.12])),
    ],
)
def test_rsprfo_partition_derivative_matches_secular_finite_difference(
    kind, eigvals, gradient,
) -> None:
    opt = RSPRFOptimizer.__new__(RSPRFOptimizer)
    opt.rfo_dict = {"max": (None, "max"), "min": (None, "min")}
    opt.log = lambda *_: None
    alpha = 10.0
    eps = 1.0e-5

    minus = opt.solve_rfo_secular(
        eigvals, gradient, alpha - eps, kind=kind,
    )
    center = opt.solve_rfo_secular(eigvals, gradient, alpha, kind=kind)
    plus = opt.solve_rfo_secular(
        eigvals, gradient, alpha + eps, kind=kind,
    )
    assert minus is not None and center is not None and plus is not None
    numeric = (
        np.dot(plus[0], plus[0]) - np.dot(minus[0], minus[0])
    ) / (2.0 * eps)
    analytic = opt._partition_dstep2_dalpha(
        alpha, center[1], center[0], eigvals, gradient,
    )
    assert analytic == pytest.approx(numeric, rel=2.0e-7, abs=1.0e-10)


def test_rfoptimizer_rejects_oversized_accelerated_displacement() -> None:
    opt = RFOptimizer.__new__(RFOptimizer)
    opt.trust_radius = 0.1
    opt.log = lambda *_: None
    ref_step = np.array([0.08, 0.0])
    result = opt._accept_accelerated_step(
        np.array([0.08, 0.0]), np.array([0.08, 0.0]), ref_step,
    )
    np.testing.assert_allclose(result, ref_step)

    accepted = opt._accept_accelerated_step(
        np.array([0.04, 0.0]), np.array([0.03, 0.0]), ref_step,
    )
    np.testing.assert_allclose(accepted, [0.07, 0.0])


def test_rfoptimizer_accelerated_step_keeps_reference_tensor_representation() -> None:
    opt = RFOptimizer.__new__(RFOptimizer)
    opt.trust_radius = 0.1
    opt.log = lambda *_: None
    ref_step = torch.tensor([0.08, 0.0], dtype=torch.float64)

    accepted = opt._accept_accelerated_step(
        np.array([0.04, 0.0]), torch.tensor([0.03, 0.0]), ref_step,
    )

    assert isinstance(accepted, torch.Tensor)
    assert accepted.dtype == ref_step.dtype
    assert accepted.device == ref_step.device
    torch.testing.assert_close(accepted, torch.tensor([0.07, 0.0], dtype=torch.float64))


def test_full_string_budget_requests_nonconverged_stop(monkeypatch) -> None:
    monkeypatch.setattr(
        Optimizer, "check_convergence", lambda *_args, **_kwargs: (False, "no"),
    )
    opt = StringOptimizer.__new__(StringOptimizer)
    opt.geometry = SimpleNamespace(fully_grown=True)
    opt.stop_in = 1
    opt.stop_in_when_full = 1
    opt.stop_requested = False
    opt.stop_reason = ""
    opt.log = lambda *_: None

    converged, _ = opt.check_convergence()
    assert converged is False
    assert opt.stop_requested is True
    assert opt.stop_reason == "full-string cycle budget exhausted"


def test_automatic_climbing_never_selects_fixed_endpoint() -> None:
    cos = SimpleNamespace(
        get_hei_index=lambda: 2,
        moving_indices=np.array([1]),
        climb="one",
        started_climbing=True,
        fixed_climb_indices=None,
        log=lambda *_: None,
    )
    assert ChainOfStates.get_climbing_indices(cos) == ()


def test_growing_string_defers_climbing_after_node_insertion() -> None:
    cos = ChainOfStates.__new__(ChainOfStates)
    cos.images = [SimpleNamespace(), SimpleNamespace(), SimpleNamespace()]
    cos.coords_length = 2
    cos.forces_list = [np.zeros(4)]

    assert cos.check_for_climbing_start(1.0) is False


def test_cos_climbing_rms_uses_moving_active_dofs_only() -> None:
    image = SimpleNamespace(
        coord_type="cart", active_dof_indices=np.array([0, 1, 2]),
    )
    cos = SimpleNamespace(
        images=[image, image, image],
        coords_length=6,
        moving_indices=[1],
        forces_list=[
            np.array(
                [10.0] * 6
                + [0.1, 0.1, 0.1, 10.0, 10.0, 10.0]
                + [10.0] * 6,
            ),
        ],
        rms=lambda arr: np.sqrt(np.mean(np.square(arr))),
        fully_grown=True,
    )

    assert ChainOfStates.check_for_climbing_start(cos, 0.2) is True


def test_molecular_volume_counts_overlapping_spheres_as_a_union(
    monkeypatch,
) -> None:
    coords = np.zeros((1, 3))
    radii = np.array([1.5])

    np.random.seed(7)
    one_sphere = molecular_volume(coords, radii, n_trial=2_000, offset=0.5)
    np.random.seed(7)
    coincident_spheres = molecular_volume(
        np.repeat(coords, 2, axis=0),
        np.repeat(radii, 2),
        n_trial=2_000,
        offset=0.5,
    )

    assert coincident_spheres == one_sphere

    monkeypatch.setattr(
        np.random,
        "rand",
        lambda n_trial, dimensions: np.full((n_trial, dimensions), 0.5),
    )
    box_volume = (2.0 * (radii[0] + 0.5)) ** 3
    center_sample = molecular_volume(coords, radii, n_trial=1, offset=0.5)
    assert center_sample[0] == pytest.approx(box_volume)


def test_growing_string_reparametrization_guards_zero_density() -> None:
    assert GrowingString._reparam_step_fraction(0.0, 0.0, 1.0e-3) is None
    with pytest.raises(ValueError, match="coincident parameter densities"):
        GrowingString._reparam_step_fraction(0.1, 0.0, 1.0e-3)

    class Image:
        coords = np.zeros(1)

        def copy(self, **_kwargs):
            return self

        def __sub__(self, _other):
            return np.zeros(1)

    string = SimpleNamespace(
        images=[Image(), Image()],
        lf_ind=0,
        sk=0.1,
        get_cur_param_density=lambda: np.zeros(2),
        reset_geometries=lambda _image: None,
    )
    with pytest.raises(ValueError, match="zero path density"):
        GrowingString.get_new_image(string, 0)


def test_max_line_search_projects_endpoint_gradients(monkeypatch) -> None:
    captured = {}

    def capture_fit(**kwargs):
        captured.update(kwargs)
        return None

    monkeypatch.setattr(
        "pysisyphus.tsoptimizers.TSHessianOptimizer.poly_fit.quartic_fit",
        capture_fit,
    )
    step = np.array([0.5, -0.25])
    g0 = np.array([2.0, 4.0])
    g1 = np.array([-3.0, 1.0])
    optimizer = SimpleNamespace(
        max_line_search=True,
        min_line_search=False,
        cur_cycle=1,
        energies=[0.0, 0.2],
        forces=[-g0, -g1],
        steps=[step],
        logger=None,
        do_line_search=TSHessianOptimizer.do_line_search,
    )

    TSHessianOptimizer.step_and_grad_from_line_search(
        optimizer,
        0.2,
        g1,
        np.eye(2),
        np.array([], dtype=int),
        np.array([0, 1]),
    )

    assert captured["g0"] == pytest.approx(step.dot(g0))
    assert captured["g1"] == pytest.approx(step.dot(g1))
    assert captured["maximize"] is True


def test_growing_string_rejects_fewer_than_two_nodes() -> None:
    with pytest.raises(ValueError, match="at least 2"):
        GrowingChainOfStates([], lambda: None, max_nodes=1)


def test_optimizer_convergence_vector_excludes_frozen_cartesian_dofs() -> None:
    geometry = SimpleNamespace(
        coord_type="cart",
        active_dof_indices=np.array([0, 1, 2]),
        cart_coords=np.zeros(12),
    )
    opt = RFOptimizer.__new__(RFOptimizer)
    opt.is_cos = False
    opt.geometry = geometry
    active = opt._active_convergence_vector(
        np.array([4.0e-4] * 3 + [0.0] * 9),
    )
    assert np.sqrt(np.mean(active**2)) == pytest.approx(4.0e-4)


def test_cos_convergence_vector_uses_moving_active_dofs_only() -> None:
    image = SimpleNamespace(
        coord_type="cart", active_dof_indices=np.array([0, 1, 2]),
    )
    opt = RFOptimizer.__new__(RFOptimizer)
    opt.is_cos = True
    opt.geometry = SimpleNamespace(
        coords_length=6,
        moving_indices=np.array([1]),
        images=[image, image, image],
    )
    vector = np.arange(18.0)
    np.testing.assert_allclose(
        opt._active_convergence_vector(vector), vector[[6, 7, 8]],
    )


def test_equal_energy_upwinding_tangent_falls_back_to_path_geometry() -> None:
    class Image:
        def __init__(self, coords):
            self.coords = np.asarray(coords, dtype=float)
            self.energy = 0.0

        def __sub__(self, other):
            return self.coords - other.coords

    cos = ChainOfStates.__new__(ChainOfStates)
    cos.images = [
        Image([0.0, 0.0]),
        Image([1.0, 0.5]),
        Image([2.0, 0.0]),
    ]
    cos.started_climbing_lanczos = False

    tangent = cos.get_tangent(1, kind="upwinding")

    np.testing.assert_allclose(tangent, [1.0, 0.0])


def test_cos_public_exports_are_bound() -> None:
    import pysisyphus.cos as cos

    assert cos.ChainOfStates is ChainOfStates
    assert cos.GrowingChainOfStates is GrowingChainOfStates
    assert cos.GrowingString.__name__ == "GrowingString"


def test_irc_rms_gradient_uses_integration_basis() -> None:
    irc = IRC.__new__(IRC)
    irc._act_dofs = np.array([0, 1, 2])
    full = np.array([1.2e-3] * 3 + [0.0] * 9)
    assert irc.active_rms_gradient(full) == pytest.approx(1.2e-3)


def test_irc_releases_dense_direction_state_before_next_branch(
    monkeypatch,
) -> None:
    irc = IRC.__new__(IRC)
    irc.mw_hessian = object()
    irc.dwi = SimpleNamespace(hessians=[object(), object()])
    emptied = []
    monkeypatch.setattr("torch.cuda.is_available", lambda: True)
    monkeypatch.setattr("torch.cuda.empty_cache", lambda: emptied.append(True))

    irc._release_direction_hessian_state()

    assert irc.mw_hessian is None
    assert irc.dwi is None
    assert emptied == [True]


def test_irc_mass_weights_hessian_without_mutating_seed() -> None:
    irc = IRC.__new__(IRC)
    irc.mm_inv2 = np.array([0.5, 1.5, 2.0])
    seed = torch.arange(9, dtype=torch.float64).reshape(3, 3)
    before = seed.clone()

    weighted = irc._mw_hessian_active(seed)

    expected = irc.mm_inv2[:, None] * before.numpy() * irc.mm_inv2[None, :]
    np.testing.assert_allclose(weighted.numpy(), expected)
    assert torch.equal(seed, before)
    assert weighted.data_ptr() != seed.data_ptr()


def test_irc_releases_initial_and_finished_hessian_owners(
    monkeypatch,
) -> None:
    monkeypatch.setattr("torch.cuda.is_available", lambda: False)
    irc = IRC.__new__(IRC)
    irc.init_hessian = torch.eye(6)
    irc._release_initial_hessian()
    assert irc.init_hessian is None
    assert irc.init_hessian_shape == (6, 6)

    irc.dwi = SimpleNamespace(hessians=[torch.eye(2), torch.eye(2)])
    irc.mw_hessian = torch.eye(2)
    irc.backward = False
    irc._release_finished_interpolation_state()
    assert irc.dwi is None
    assert irc.mw_hessian is None

    retained = torch.eye(2)
    irc.dwi = SimpleNamespace(hessians=[torch.eye(2), torch.eye(2)])
    irc.mw_hessian = retained
    irc.backward = True
    irc._release_finished_interpolation_state()
    assert irc.dwi is None
    assert irc.mw_hessian is retained


def test_euler_bofill_update_is_in_place_and_matches_dense_formula() -> None:
    irc = EulerPC.__new__(EulerPC)
    original = torch.tensor(
        [[2.0, 0.2, 0.1], [0.2, 1.5, -0.3], [0.1, -0.3, 1.1]],
        dtype=torch.float64,
    )
    dx = torch.tensor([0.4, -0.2, 0.3], dtype=torch.float64)
    dg = torch.tensor([0.1, 0.5, -0.4], dtype=torch.float64)
    dense_update, _ = bofill_update(original.clone(), dx, dg)
    expected = original + dense_update
    irc.mw_hessian = original.clone()
    irc.hessian_update_func = bofill_update
    storage = irc.mw_hessian.data_ptr()

    key = irc._apply_hessian_update(dx, dg)

    assert key == "Bofill"
    assert irc.mw_hessian.data_ptr() == storage
    torch.testing.assert_close(irc.mw_hessian, expected)


@pytest.mark.parametrize("roots", ([1], [0, 1]))
def test_ts_fixed_root_mode_preserves_configured_roots(roots) -> None:
    optimizer = RSPRFOptimizer.__new__(RSPRFOptimizer)
    optimizer.small_eigval_thresh = 1.0e-8
    optimizer.log_negative_eigenvalues = lambda *_args: None
    optimizer._physical_ts_mode = None
    optimizer.track_mode_by_overlap = False
    optimizer.roots = list(roots)

    eigvals = np.array([-2.0, -1.0, 0.5])
    eigvecs = np.eye(3)
    optimizer.update_ts_mode(eigvals, eigvecs)

    np.testing.assert_array_equal(optimizer.roots, roots)
    np.testing.assert_allclose(optimizer.ts_modes, eigvecs[:, roots].T)


@pytest.mark.parametrize("roots", ([0, 0], [3]))
def test_ts_fixed_root_mode_rejects_invalid_roots(roots) -> None:
    optimizer = RSPRFOptimizer.__new__(RSPRFOptimizer)
    optimizer.small_eigval_thresh = 1.0e-8
    optimizer.log_negative_eigenvalues = lambda *_args: None
    optimizer._physical_ts_mode = None
    optimizer.track_mode_by_overlap = False
    optimizer.roots = list(roots)

    with pytest.raises(ValueError):
        optimizer.update_ts_mode(np.array([-2.0, -1.0, 0.5]), np.eye(3))


def test_dwi_evicts_old_hessian_before_copying_replacement() -> None:
    dwi = DWI(maxlen=2)
    matrices = [torch.eye(3) * value for value in (1.0, 2.0, 3.0)]
    for value, matrix in enumerate(matrices):
        dwi.update(
            np.array([value]),
            float(value),
            np.array([value]),
            matrix,
            copy_hessian=True,
        )

    assert len(dwi.hessians) == 2
    torch.testing.assert_close(dwi.hessians[0], matrices[1])
    torch.testing.assert_close(dwi.hessians[1], matrices[2])
    assert dwi.hessians[1].data_ptr() != matrices[2].data_ptr()
    source = inspect.getsource(DWI.update)
    assert source.index("self.hessians.popleft()") < source.index(
        "hessian.detach().clone()"
    )


def test_euler_corrector_reaches_unweighted_target_for_heavy_mass() -> None:
    class ConstantDWI:
        @staticmethod
        def interpolate(coords, gradient=True):
            return 0.0, np.ones_like(coords)

    irc = EulerPC.__new__(EulerPC)
    irc._m_sqrt = np.array([4.0])  # oxygen-like mass sqrt
    irc._act_dofs = np.array([0])
    irc.log = lambda *_: None
    start = np.zeros(1)
    corrected = irc.corrector_step(start, 0.1, ConstantDWI())
    unweighted_length = np.linalg.norm((corrected - start) / irc._m_sqrt)
    assert unweighted_length == pytest.approx(0.1, abs=1.0e-4)


def test_euler_corrector_degrades_instead_of_aborting_on_oscillation(capsys) -> None:
    """An oscillating DWI must cost one corrector, not the whole IRC.

    The corrector descends the two-point DWI *interpolation*, not the real PES,
    so a reversal there is an interpolation artefact. Raising instead of
    returning the last non-oscillating point aborts a complete run from inside
    a healthy IRC, because this branch is also the integration loop's escape
    hatch.
    """

    class OscillatingDWI:
        """1-D well at the origin; fixed-length descent overshoots and flips."""

        @staticmethod
        def interpolate(coords, gradient=True):
            return 0.0, coords.copy()

    irc = EulerPC.__new__(EulerPC)
    irc._m_sqrt = np.array([1.0])
    irc._act_dofs = np.array([0])
    irc.log = lambda *_: None
    start = np.array([0.02])

    corrected = irc.corrector_step(start, 0.1, OscillatingDWI())

    assert np.all(np.isfinite(corrected))
    # Degraded, not aborted: short of the requested 0.1 but still advancing, so
    # the caller gets a usable geometry and the IRC keeps going.
    advance = float(np.linalg.norm(corrected - start))
    assert 0.0 < advance < 0.1
    assert "oscillated" in capsys.readouterr().out


def test_euler_corrector_rejects_incomplete_zero_gradient() -> None:
    class ZeroDWI:
        @staticmethod
        def interpolate(coords, gradient=True):
            return 0.0, np.zeros_like(coords)

    irc = EulerPC.__new__(EulerPC)
    irc._m_sqrt = np.array([1.0])
    irc._act_dofs = np.array([0])
    irc.log = lambda *_: None
    with pytest.raises(RuntimeError, match="zero or non-finite gradient"):
        irc.corrector_step(np.zeros(1), 0.1, ZeroDWI())


def test_directional_irc_trajectory_contains_terminal_frame(tmp_path) -> None:
    irc = IRC.__new__(IRC)
    irc.atoms = ("H",)
    irc._m_sqrt = np.ones(3)
    irc.get_path_for_fn = lambda filename: str(tmp_path / filename)
    # Forward data has already been reversed by IRC.irc() for stitched-path
    # order: endpoint -> TS-adjacent.
    irc.irc_coords = [
        np.array([1.0, 0.0, 0.0]),
        np.array([0.0, 0.0, 0.0]),
    ]
    irc.irc_gradients = [np.zeros(3), np.zeros(3)]
    irc.irc_mw_coords = [coords.copy() for coords in irc.irc_coords]
    irc.irc_mw_gradients = [np.zeros(3), np.zeros(3)]
    irc.irc_energies = [-1.1, -1.0]
    irc.all_coords = []
    irc.all_gradients = []
    irc.all_mw_coords = []
    irc.all_mw_gradients = []
    irc.all_energies = []
    irc.converged = True
    irc.integration_stop_reason = ""
    irc.energy_increased = False
    irc.energy_converged = True
    irc.never_stop = False
    irc.cur_cycle = 1

    irc.set_data("forward")

    atom_lines = [
        line
        for line in (tmp_path / "forward_irc_trj.xyz").read_text().splitlines()
        if line.strip().startswith("H ")
    ]
    assert len(atom_lines) == len(irc.forward_energies) == 2
    # Directional file is chronological TS -> endpoint and includes endpoint.
    assert float(atom_lines[0].split()[1]) == pytest.approx(0.0)
    assert float(atom_lines[-1].split()[1]) == pytest.approx(
        0.529177, rel=1e-6
    )


def test_full_hessian_normal_modes_honor_geometry_freezes() -> None:
    geometry = Geometry(
        ["C"] * 5,
        np.array([
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [0.0, 2.0, 0.0],
            [0.0, 0.0, 2.0],
            [2.0, 2.0, 2.0],
        ]).reshape(-1),
        coord_type="cart",
        freeze_atoms=[0, 1, 2],
    )
    _, _, _, cart_modes = geometry.get_normal_modes(np.eye(15), full=True)

    assert cart_modes.shape == (15, 6)
    np.testing.assert_allclose(cart_modes[:9], 0.0)
    assert geometry.within_partial_hessian is None


def test_final_hessian_thermochemistry_uses_keyword_contract(monkeypatch) -> None:
    import pysisyphus.helpers as helpers

    calls = []

    class FakeGeometry:
        cart_hessian = np.eye(3)

        @staticmethod
        def mass_weigh_hessian(hessian):
            return hessian

        @staticmethod
        def eckart_projection(hessian):
            return hessian

        @staticmethod
        def get_thermoanalysis(*, T, p):
            calls.append((T, p))
            return "thermo"

    monkeypatch.setattr(helpers, "report_isotopes", lambda *_args: None)
    result = helpers.do_final_hessian(FakeGeometry(), T=310.0, p=98_000.0)

    assert calls == [(310.0, 98_000.0)]
    assert result.thermo == "thermo"


@pytest.mark.parametrize(
    "update_name",
    ["ts_bfgs_update", "ts_bfgs_update_org", "ts_bfgs_update_revised"],
)
@pytest.mark.parametrize("device", ["cpu", "cuda"])
def test_ts_bfgs_updates_preserve_tensor_device_and_numpy_values(
    update_name: str, device: str
) -> None:
    if device == "cuda" and not torch.cuda.is_available():
        pytest.skip("CUDA is unavailable")
    from pysisyphus.optimizers import hessian_updates

    update = getattr(hessian_updates, update_name)
    hessian = np.array(
        [[-1.0, 0.1, 0.0], [0.1, 2.0, 0.2], [0.0, 0.2, 3.0]],
        dtype=np.float64,
    )
    step = np.array([0.2, -0.1, 0.3], dtype=np.float64)
    gradient_delta = hessian @ step + np.array([0.1, 0.05, -0.02])
    expected, _ = update(hessian, step, gradient_delta)

    actual, _ = update(
        torch.as_tensor(hessian, device=device),
        torch.as_tensor(step, device=device),
        torch.as_tensor(gradient_delta, device=device),
    )

    assert isinstance(actual, torch.Tensor)
    assert actual.device.type == device
    assert actual.dtype == torch.float64
    np.testing.assert_allclose(actual.detach().cpu().numpy(), expected)


@pytest.mark.parametrize("translation", [0.0, 25.0])
def test_legacy_trans_rot_vectors_keep_linear_rank_under_translation(
    translation: float,
) -> None:
    """A translated linear molecule keeps rigid rank five, not six."""
    from pysisyphus.Geometry import get_trans_rot_vectors

    coords3d = np.array(
        [[0.0, 0.0, -1.2], [0.0, 0.0, 0.0], [0.0, 0.0, 1.2]], dtype=np.float64
    )
    coords3d[:, 2] += translation
    masses = np.array([15.999, 12.011, 15.999], dtype=np.float64)

    tr_vecs = get_trans_rot_vectors(coords3d.reshape(-1), masses)
    assert tr_vecs.shape[0] == 5


def test_active_tr_basis_reports_frozen_rank_for_rank_zero_rigid_basis() -> None:
    """A rank-zero rigid basis still reports an initialized frozen rank."""
    from pysisyphus.tr_projection import active_tr_basis

    coords = torch.tensor(
        [[0.0, 0.0, 0.0], [0.0, 0.0, 1.4], [0.0, 1.1, 0.0]], dtype=torch.float64
    )
    masses = torch.tensor([12.0, 1.0, 1.0], dtype=torch.float64)

    basis, info = active_tr_basis(coords, masses, [1, 2], rtol=1.0)
    assert info.full_rigid_rank == 0
    assert info.frozen_constraint_rank == 0
    assert basis.shape[1] == info.effective_rank == 0


def test_get_imag_frequencies_rejects_small_positive_eigenvalues() -> None:
    """A small positive eigenvalue is not reported as an imaginary frequency."""
    geom = Geometry(("h", "h"), np.array([0.0, 0.0, 0.0, 0.0, 0.0, 1.4]))
    eigvals = np.array([-1.0e-5, 0.0, 5.0e-7, 1.0e-3])
    nus = np.array([-11.0, 0.0, 12.0, 500.0])
    geom.get_normal_modes = lambda hessian=None: (nus, eigvals, None, None)

    np.testing.assert_array_equal(geom.get_imag_frequencies(), np.array([-11.0]))


def test_normal_modes_remove_exactly_the_rigid_rank_and_keep_soft_roots() -> None:
    """Low positive and low negative roots survive the rigid-space removal."""
    from pysisyphus.normal_modes import _frequencies_cm_and_modes
    from pysisyphus.tr_projection import active_tr_basis

    atomic_numbers = [6, 1, 8, 7]
    coords_bohr = np.array(
        [[0.0, 0.0, 0.0], [2.1, 0.0, 0.0], [0.2, 2.0, 0.0], [0.3, 0.4, 2.2]],
        dtype=np.float64,
    )
    freeze_idx = [0]
    device = torch.device("cpu")

    masses = torch.as_tensor(
        np.array([12.011, 1.008, 15.999, 14.007]) * 1822.888486209,
        dtype=torch.float64,
    )
    _, info = active_tr_basis(
        torch.as_tensor(coords_bohr, dtype=torch.float64),
        masses,
        [1, 2, 3],
        mode="constrained",
    )

    gen = torch.Generator().manual_seed(7)
    raw = torch.randn((12, 12), generator=gen, dtype=torch.float64)
    hessian = (raw.T @ raw + 0.2 * torch.eye(12, dtype=torch.float64)) * 1.0e-9
    freqs_cm, modes = _frequencies_cm_and_modes(
        hessian.clone(),
        atomic_numbers,
        coords_bohr,
        device,
        freeze_idx=freeze_idx,
        frequency_zero_cutoff_cm=0.0,
    )

    # Exactly the constrained rigid rank is removed from the 9 active DOF.
    assert freqs_cm.size == 9 - info.effective_rank
    assert modes.shape == (freqs_cm.size, 12)
    # The scaled Hessian only produces sub-cm^-1 roots, none of which is dropped.
    assert np.max(np.abs(freqs_cm)) < 5.0


@pytest.mark.parametrize("update_name", ["damped_bfgs_update", "flowchart_update"])
@pytest.mark.parametrize("device", ["cpu", "cuda"])
def test_minimization_updates_preserve_tensor_device_and_numpy_values(
    update_name: str, device: str
) -> None:
    if device == "cuda" and not torch.cuda.is_available():
        pytest.skip("CUDA is unavailable")
    from pysisyphus.optimizers import hessian_updates

    update = getattr(hessian_updates, update_name)
    hessian = np.array(
        [[1.5, 0.1, 0.0], [0.1, 2.0, 0.2], [0.0, 0.2, 3.0]], dtype=np.float64
    )
    step = np.array([0.2, -0.1, 0.3], dtype=np.float64)
    gradient_delta = hessian @ step + np.array([0.1, 0.05, -0.02])
    expected, _ = update(hessian, step, gradient_delta)

    actual, _ = update(
        torch.as_tensor(hessian, device=device),
        torch.as_tensor(step, device=device),
        torch.as_tensor(gradient_delta, device=device),
    )

    assert isinstance(actual, torch.Tensor)
    assert actual.device.type == device
    assert actual.dtype == torch.float64
    np.testing.assert_allclose(actual.detach().cpu().numpy(), expected)


@pytest.mark.parametrize("device", ["cpu", "cuda"])
def test_disabled_hessian_update_preserves_the_tensor_backend(device: str) -> None:
    if device == "cuda" and not torch.cuda.is_available():
        pytest.skip("CUDA is unavailable")
    from pysisyphus.optimizers.HessianOptimizer import dummy_hessian_update

    hessian = torch.eye(3, dtype=torch.float64, device=device)
    update, key = dummy_hessian_update(hessian, None, None)

    assert key == "no"
    assert isinstance(update, torch.Tensor)
    assert update.device.type == device
    assert update.dtype == torch.float64
    assert isinstance((hessian + update), torch.Tensor)


# Two roots below the historical 1e-6 magnitude cutoff, one of each sign, plus
# four ordinary roots.  A rank-based rigid removal must return all six.
_COMPLEMENT_OMEGA2 = (-4.0e-7, 2.5e-7, 0.02, 0.05, 0.08, 0.11)

_NORMAL_MODE_COORDS_BOHR = np.array(
    [
        [0.0, 0.0, 0.0],
        [2.05, 0.0, 0.0],
        [-0.7, 1.95, 0.0],
        [-0.6, -0.7, 1.85],
    ],
    dtype=np.float64,
)
_NORMAL_MODE_NUMBERS = [6, 1, 8, 7]


def _hessian_with_complement_spectrum(atomic_numbers, coords_bohr, active_idx, omega2):
    """Cartesian Hessian whose rigid complement has exactly ``omega2``."""
    from pysisyphus.constants import AMU2AU
    from pysisyphus.normal_modes import _safe_masses_amu
    from pysisyphus.tr_projection import active_tr_basis

    masses_amu = _safe_masses_amu(atomic_numbers)
    coords = torch.as_tensor(coords_bohr.reshape(-1, 3), dtype=torch.float64)
    basis, info = active_tr_basis(
        coords,
        torch.as_tensor(masses_amu * AMU2AU, dtype=torch.float64),
        list(active_idx),
        mode="constrained",
    )
    complete, _ = torch.linalg.qr(basis, mode="complete")
    lift = complete[:, info.effective_rank :].T
    assert lift.shape[0] == len(omega2)
    hessian_mw = lift.T @ torch.diag(torch.as_tensor(omega2, dtype=torch.float64)) @ lift
    sqrt_mass = torch.sqrt(
        torch.repeat_interleave(
            torch.as_tensor(masses_amu[list(active_idx)], dtype=torch.float64), 3
        )
    )
    return hessian_mw * sqrt_mass.view(-1, 1) * sqrt_mass.view(1, -1), info


def _expected_freqs_cm(omega2):
    import ase.units as units

    from pysisyphus.constants import AU2EV, BOHR2ANG

    scale = (
        units._hbar
        * 1e10
        / np.sqrt(units._e * units._amu)
        * np.sqrt(AU2EV)
        / BOHR2ANG
    )
    values = np.asarray(omega2, dtype=np.float64)
    return np.sign(values) * scale * np.sqrt(np.abs(values)) / units.invcm


@pytest.mark.parametrize("hessian_kind", ["full", "active_block"])
def test_compact_normal_modes_keep_low_signed_roots_under_phva(hessian_kind: str) -> None:
    """PHVA keeps genuine low positive and low negative complement roots."""
    from pysisyphus.normal_modes import _frequencies_cm_and_modes

    freeze_idx = [0]
    active_idx = [1, 2, 3]
    active_block, info = _hessian_with_complement_spectrum(
        _NORMAL_MODE_NUMBERS,
        _NORMAL_MODE_COORDS_BOHR,
        active_idx,
        _COMPLEMENT_OMEGA2,
    )
    # One frozen anchor leaves the three rotations about it as the rigid rank.
    assert info.effective_rank == 3
    if hessian_kind == "full":
        hessian = torch.zeros((12, 12), dtype=torch.float64)
        # Frozen rows/columns must be discarded, so they carry values that would
        # visibly contaminate the spectrum if they leaked into the active block.
        hessian[:3, :3] = torch.eye(3, dtype=torch.float64) * 7.0
        active_dofs = torch.as_tensor(
            [3 * atom + axis for atom in active_idx for axis in range(3)],
            dtype=torch.long,
        )
        hessian[active_dofs[:, None], active_dofs[None, :]] = active_block
    else:
        hessian = active_block.clone()

    freqs_cm, modes = _frequencies_cm_and_modes(
        hessian,
        _NORMAL_MODE_NUMBERS,
        _NORMAL_MODE_COORDS_BOHR,
        torch.device("cpu"),
        freeze_idx=freeze_idx,
        frequency_zero_cutoff_cm=0.0,
    )

    # Exactly the rigid rank is removed from the nine active DOF.
    assert freqs_cm.size == 9 - info.effective_rank == len(_COMPLEMENT_OMEGA2)
    np.testing.assert_allclose(
        np.sort(freqs_cm), np.sort(_expected_freqs_cm(_COMPLEMENT_OMEGA2)), rtol=1e-8
    )
    # Both roots below the historical tolerance survive with their own sign.
    assert min(abs(value) for value in _COMPLEMENT_OMEGA2[:2]) < 1.0e-6
    assert (freqs_cm < 0.0).sum() == 1
    assert modes.shape == (len(_COMPLEMENT_OMEGA2), 12)
    # Frozen degrees of freedom stay exactly zero after embedding.
    assert torch.count_nonzero(modes[:, :3]) == 0


def test_compact_normal_modes_keep_low_signed_roots_without_frozen_atoms() -> None:
    """The unconstrained kernel retains both low-magnitude complement roots."""
    from pysisyphus.constants import AMU2AU
    from pysisyphus.normal_modes import _frequencies_cm_and_modes, _safe_masses_amu
    from pysisyphus.tr_projection import active_tr_basis

    active_idx = [0, 1, 2, 3]
    # A chemically ordinary imaginary mode beside the two low-magnitude roots, so
    # a saddle-order consumer must count two negative roots, not one.
    omega2 = (-0.02,) + _COMPLEMENT_OMEGA2[:2] + (0.05, 0.08, 0.15)
    hessian, info = _hessian_with_complement_spectrum(
        _NORMAL_MODE_NUMBERS, _NORMAL_MODE_COORDS_BOHR, active_idx, omega2
    )
    # No frozen anchor leaves the six full rigid motions as the rigid rank.
    assert info.effective_rank == 6

    freqs_cm, modes = _frequencies_cm_and_modes(
        hessian.clone(),
        _NORMAL_MODE_NUMBERS,
        _NORMAL_MODE_COORDS_BOHR,
        torch.device("cpu"),
        frequency_zero_cutoff_cm=0.0,
    )

    assert freqs_cm.size == 12 - info.effective_rank == len(omega2)
    np.testing.assert_allclose(
        np.sort(freqs_cm), np.sort(_expected_freqs_cm(omega2)), rtol=1e-8
    )
    assert min(abs(value) for value in omega2[1:3]) < 1.0e-6
    assert (freqs_cm < 0.0).sum() == 2
    # Retained modes stay orthogonal to the removed rigid space.
    basis, _ = active_tr_basis(
        torch.as_tensor(_NORMAL_MODE_COORDS_BOHR, dtype=torch.float64),
        torch.as_tensor(
            _safe_masses_amu(_NORMAL_MODE_NUMBERS) * AMU2AU, dtype=torch.float64
        ),
        active_idx,
        mode="constrained",
    )
    overlap = basis.T @ modes.T
    assert float(torch.max(torch.abs(overlap))) < 1.0e-10


def test_normalize_prim_input_expands_translation_and_rotation_shortcuts() -> None:
    """The advertised TRANSLATION/ROTATION shortcuts reach real primitives."""
    from pysisyphus.intcoords.PrimTypes import (
        PrimTypes,
        normalize_prim_input,
        prims_from_prim_inputs,
    )

    translation = normalize_prim_input(["TRANSLATION", 0, 1])
    rotation = normalize_prim_input(["ROTATION", 0, 1])

    assert [tp[0] for tp in translation] == [
        PrimTypes.TRANSLATION_X,
        PrimTypes.TRANSLATION_Y,
        PrimTypes.TRANSLATION_Z,
    ]
    assert [tp[0] for tp in rotation] == [
        PrimTypes.ROTATION_A,
        PrimTypes.ROTATION_B,
        PrimTypes.ROTATION_C,
    ]
    # The generic enum members have no PrimMap entry, so only the expanded
    # component types can be instantiated.
    assert len(prims_from_prim_inputs([["TRANSLATION", 0, 1]])) == 3


def test_normalize_prim_input_keeps_the_distance_function_coefficient() -> None:
    """The DIST_FUNC coefficient stays a float instead of truncating to zero."""
    from pysisyphus.intcoords.PrimTypes import PrimTypes, normalize_prim_input

    (typed_prim,) = normalize_prim_input(["DIST_FUNC", 0, 1, 2, 3, 0.5])

    assert typed_prim[0] == PrimTypes.DISTANCE_FUNCTION
    assert typed_prim[1:5] == (0, 1, 2, 3)
    assert typed_prim[5] == pytest.approx(0.5)


def test_normalize_prim_input_returns_a_hashable_tuple_for_enum_input() -> None:
    """An enum-headed list input is normalized to the documented tuple."""
    from pysisyphus.intcoords.PrimTypes import PrimTypes, normalize_prim_input

    (typed_prim,) = normalize_prim_input([PrimTypes.BOND, 0, 1])

    assert typed_prim == (PrimTypes.BOND, 0, 1)
    assert hash(typed_prim)


def test_setup_redundant_keeps_explicit_definitions_on_their_own_atoms() -> None:
    """Explicit definitions are not remapped twice when atoms are frozen."""
    from pysisyphus.intcoords.PrimTypes import PrimTypes
    from pysisyphus.intcoords.setup import setup_redundant

    atoms = ("H", "C", "C", "H")
    coords3d = np.array(
        [
            [-2.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [2.4, 0.0, 0.0],
            [4.4, 0.0, 0.0],
        ]
    )
    coord_info = setup_redundant(
        atoms,
        coords3d,
        freeze_atoms=[0],
        define_prims=[(PrimTypes.AUX_BOND, 1, 3)],
    )

    assert (PrimTypes.AUX_BOND, 1, 3) in coord_info.typed_prims
    # Frozen atom 0 is excluded from the setup, so naming it is rejected.
    from pysisyphus.intcoords.exceptions import PrimitiveNotDefinedException

    with pytest.raises(PrimitiveNotDefinedException):
        setup_redundant(
            atoms,
            coords3d,
            freeze_atoms=[0],
            define_prims=[(PrimTypes.AUX_BOND, 0, 3)],
        )


def test_find_bonds_keeps_the_bond_array_shape_without_bonds() -> None:
    """A valid bondless structure serializes instead of sorting a missing axis."""
    from pysisyphus.intcoords.setup_fast import find_bonds
    from pysisyphus.io.pdb import get_conect_lines

    coords = np.array([[0.0, 0.0, 0.0]])
    bonds = find_bonds(("Ne",), coords)

    assert bonds.shape == (0, 2)
    assert bonds.dtype.kind == "i"
    assert get_conect_lines(("Ne",), coords) == []


def test_find_bonds_for_geom_forwards_the_bond_factor() -> None:
    """The accepted bond_factor reaches the detector instead of being dropped."""
    from pysisyphus.intcoords.setup_fast import find_bonds_for_geom

    # Two carbons far enough apart that only a generous factor bonds them.
    geom = Geometry(("c", "c"), np.array([0.0, 0.0, 0.0, 0.0, 0.0, 4.5]))

    assert find_bonds_for_geom(geom, bond_factor=1.3).tolist() == []
    assert find_bonds_for_geom(geom, bond_factor=1.7).tolist() == [[0, 1]]


def test_redundant_coords_registers_constrained_primitive_indices() -> None:
    """Constrained primitives go through the typed-primitive setter."""
    from pysisyphus.intcoords.PrimTypes import PrimTypes

    atoms = ("O", "H", "H")
    coords3d = np.array(
        [[0.0, 0.0, 0.0], [1.8, 0.0, 0.0], [-0.5, 1.7, 0.0]]
    )
    constraint = (PrimTypes.BOND, 0, 1)
    red = RedundantCoords(atoms, coords3d, constrain_prims=[constraint])

    assert constraint in red.typed_prims
    index = red.typed_prims.index(constraint)
    assert red.constrained_indices == [index]
    # The per-category index lists cover every typed primitive, constraints
    # included, so the projector and the primitive list stay aligned.
    assert index in red.bond_indices
    assert len(red.primitives) == len(red.typed_prims)


def test_augment_bonds_preserves_geometry_state(monkeypatch) -> None:
    """Augmentation keeps frozen atoms, isotopes, and coordinate options."""
    from pysisyphus.intcoords import augment_bonds as augment_bonds_module

    atoms = ("O", "H", "H", "H")
    coords = np.array(
        [0.0, 0.0, 0.0, 1.8, 0.0, 0.0, -0.5, 1.7, 0.0, 6.0, 0.0, 0.0]
    )
    geom = Geometry(
        atoms,
        coords,
        coord_type="redund",
        freeze_atoms=[0],
        isotopes=((1, 2.0),),
    )
    geom.cart_hessian = np.eye(coords.size)
    monkeypatch.setattr(
        augment_bonds_module, "find_missing_strong_bonds", lambda *a, **k: [(1, 3)]
    )

    new_geom = augment_bonds_module.augment_bonds(geom)

    assert new_geom is not geom
    assert list(new_geom.freeze_atoms) == [0]
    assert new_geom.isotopes == geom.isotopes
    assert new_geom.coord_type == geom.coord_type
    assert (PrimTypes.AUX_BOND, 1, 3) in new_geom.internal.typed_prims


def test_lanczos_returns_the_current_ritz_pair_on_residual_breakdown() -> None:
    """An exact residual breakdown must not produce a NaN climbing mode."""
    from pysisyphus.modefollow.lanczos import lanczos

    # A rank-one Hessian exhausts its Krylov space after the first vector.
    hessian = np.zeros((3, 3))
    hessian[0, 0] = -0.4

    def grad_getter(coords):
        return hessian @ np.asarray(coords, dtype=float)

    eigval, eigvec = lanczos(
        np.zeros(3), grad_getter, guess=np.array([1.0, 0.0, 0.0]), max_cycles=5
    )

    assert np.isfinite(eigval)
    assert np.all(np.isfinite(eigvec))
    assert eigval == pytest.approx(-0.4, abs=1e-6)


@pytest.mark.parametrize("guess", [np.zeros(3), np.array([np.nan, 0.0, 0.0])])
def test_lanczos_rejects_a_degenerate_initial_guess(guess: np.ndarray) -> None:
    from pysisyphus.modefollow.lanczos import lanczos

    with pytest.raises(ValueError):
        lanczos(np.zeros(3), lambda coords: np.zeros(3), guess=guess)


def test_quartic_fit_reports_no_fit_for_a_degenerate_line() -> None:
    """Equal endpoint energies with zero gradients no longer divide by zero."""
    from pysisyphus.optimizers.poly_fit import quartic_fit

    assert quartic_fit(-1.0, -1.0, 0.0, 0.0) is None


def test_sim_gediis_weights_use_the_inverse_hessian_quadratic_form() -> None:
    """The GEDIIS energy model uses f^T H^-1 f, not a row-summed contraction."""
    forces = np.array([[0.3, -0.1], [0.05, 0.2]])
    hessian = np.array([[2.0, 0.7], [0.7, 3.0]])
    hessian_inv = np.linalg.pinv(hessian, rcond=1e-6)
    expected = np.array([row @ hessian_inv @ row for row in forces])

    actual = np.einsum("ki,ij,kj->k", forces, hessian_inv, forces)

    np.testing.assert_allclose(actual, expected)
    # The previous row-summed contraction differs for an off-diagonal Hessian.
    assert not np.allclose(
        np.einsum("ki,ji,ki->k", forces, hessian_inv, forces), expected
    )
    source = inspect.getsource(gdiis_module.gediis)
    assert 'einsum("ki,ij,kj->k"' in source
    assert 'einsum("ki,ji,ki->k"' not in source


def test_eulerpc_corrector_keeps_the_last_advancing_microstep() -> None:
    """An immediate DWI reversal returns the advanced point, not the start."""
    init = np.array([0.0, 0.0])
    step = np.array([1.0, 0.0])
    calls = {"n": 0}

    class _ReversingDWI:
        def interpolate(self, coords, gradient=False):
            calls["n"] += 1
            # First microstep advances along +x, the second reverses.
            sign = -1.0 if calls["n"] == 1 else 1.0
            return 0.0, sign * step

    irc = EulerPC.__new__(EulerPC)
    irc.log = lambda *_: None
    irc._m_sqrt = np.ones(2)
    irc._act_dofs = np.array([0, 1])
    irc.get_integration_length_func = lambda _init: (
        lambda coords: float(np.linalg.norm(coords - init))
    )

    corrected = EulerPC.corrector_step(irc, init.copy(), 10.0, _ReversingDWI())

    assert calls["n"] == 2
    assert np.linalg.norm(corrected - init) > 0.0


def test_irc_gives_integration_and_energy_rise_priority_over_convergence() -> None:
    """A converged gradient no longer masks an integration stop or energy rise."""
    source = inspect.getsource(IRC.irc)
    stop_at = source.index("if self.integration_stop_requested:")
    increase_at = source.index("elif energy_increase_msg:")
    converged_at = source.index("elif self._gradient_converged(rms_grad):")
    assert stop_at < increase_at < converged_at

    irc = IRC.__new__(IRC)
    irc.never_stop = False
    irc.energy_increased = True
    irc.energy_converged = True
    assert irc._energy_increase_stop_message() == "Energy increased!"
    irc.never_stop = True
    assert irc._energy_increase_stop_message() == ""
    assert irc._energy_stop_message() == ""


def test_irc_resolves_the_requested_hessian_init() -> None:
    """hessian_init is resolved through get_guess_hessian instead of ignored."""
    source = inspect.getsource(IRC.run)
    assert "get_guess_hessian(" in source, source
    assert "self.init_hessian = self.geometry.hessian" not in source
