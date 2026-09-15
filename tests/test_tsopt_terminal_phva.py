"""Terminal PHVA runs only after numerical convergence."""

from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
import torch

from pdb2reaction.workflows import tsopt


class _Geometry:
    atomic_numbers = np.array([1])
    atoms = ["H"]
    cart_coords = np.zeros(3)


def test_optimizer_terminal_phva_carries_the_exact_hessian_for_irc_cache():
    exact_hessian = torch.diag(torch.tensor([-1.0, 2.0, 3.0]))
    optimizer = SimpleNamespace(
        _last_exact_cart_coords=np.zeros(3),
        _last_exact_frequencies_cm=np.array([-100.0, 20.0, 30.0]),
        _last_exact_modes=torch.eye(3),
        _last_rigid_projection_info={
            "active_atoms": [0], "frozen_atoms": [],
            "treatment": "constrained", "frequency_zero_cutoff_cm": 5.0,
            "raw_mode_count": 3, "near_zero_frequencies_cm": [],
        },
        cur_H=exact_hessian,
    )

    reused = tsopt._optimizer_exact_frequency_data(optimizer, _Geometry())

    assert reused is not None
    assert len(reused) == 4
    cached_hessian = reused[3]
    assert torch.equal(cached_hessian, exact_hessian)
    assert cached_hessian.data_ptr() != exact_hessian.data_ptr()


def _two_atom_exact_cache():
    geometry = SimpleNamespace(
        cart_coords=np.zeros(9), freeze_atoms=[2], tr_projection="constrained",
    )
    optimizer = SimpleNamespace(
        _last_exact_cart_coords=np.zeros(9),
        _last_exact_frequencies_cm=np.array([-100., 20., 30.]),
        _last_exact_modes=torch.zeros((3, 9)),
        _last_rigid_projection_info={
            "active_atoms": [0, 1], "frozen_atoms": [2],
            "treatment": "constrained", "frequency_zero_cutoff_cm": 5.,
            "raw_mode_count": 3, "near_zero_frequencies_cm": [],
        },
        cur_H=torch.eye(6),
    )
    return optimizer, geometry


@pytest.mark.parametrize(
    "requested, reusable", [([0, 1], True), ([0, 2], False), ([0, 1, 2], False)]
)
def test_exact_frequency_cache_requires_matching_active_atoms(requested, reusable):
    optimizer, geometry = _two_atom_exact_cache()
    cached = tsopt._optimizer_exact_frequency_data(
        optimizer, geometry, requested_active_atoms=requested,
    )
    assert (cached is not None) is reusable


def test_exact_frequency_cache_detects_restored_full_geometry_mask():
    optimizer, geometry = _two_atom_exact_cache()
    geometry.freeze_atoms = []
    assert tsopt._optimizer_exact_frequency_data(optimizer, geometry) is None


@pytest.mark.parametrize("change", ["cutoff", "projection", "missing_metadata"])
def test_exact_frequency_cache_requires_matching_analysis_conditions(change):
    optimizer, geometry = _two_atom_exact_cache()
    kwargs = {}
    if change == "cutoff":
        kwargs["frequency_zero_cutoff_cm"] = 10.
    elif change == "projection":
        geometry.tr_projection = "different"
    else:
        optimizer._last_rigid_projection_info.clear()
    assert tsopt._optimizer_exact_frequency_data(optimizer, geometry, **kwargs) is None


def _runner(tmp_path, monkeypatch, *, stalled):
    runner = object.__new__(tsopt.HessianDimer)
    runner.geom = _Geometry()
    runner.dump = False
    runner.optim_all_path = tmp_path / "optimization_all_trj.xyz"
    runner.mode_path = tmp_path / ".dimer_mode.dat"
    runner.out_dir = tmp_path
    runner.vib_dir = tmp_path / "vib"
    runner.freeze_atoms = []
    runner.masses_au_t = torch.ones(1)
    runner.masses_amu = np.ones(1)
    runner.device = torch.device("cpu")
    runner.root = 0
    runner.thresh_loose = "baker"
    runner.thresh = "baker"
    runner.flatten_max_iter = 1
    runner.flatten_loop_bofill = False
    runner.neg_freq_thresh_cm = 5.0
    runner.tr_projection = "none"
    runner.rigid_projection_info = {}
    runner.max_total_cycles = 1
    runner._cycles_spent = 0
    runner.is_converged = False
    runner.is_stalled = False
    runner.stop_reason = ""
    runner.n_imaginary_modes = None
    runner.imaginary_frequencies_cm = []
    runner.saddle_order_verified = False
    runner.prepared_input = None
    runner.ref_pdb = None

    hessian_calls = []
    mode_exports = []

    def fake_hessian(*, allow_reuse):
        hessian_calls.append(allow_reuse)
        return torch.eye(3)

    def fake_loop(_threshold):
        runner._cycles_spent = runner.max_total_cycles
        if stalled:
            runner.is_stalled = True
            runner.stop_reason = "energy plateau"
        return 1, False, False

    runner._calc_full_hessian_cached = fake_hessian
    runner._dimer_loop = fake_loop
    monkeypatch.setattr(
        tsopt,
        "_mode_direction_by_root",
        lambda *args, **kwargs: (np.ones((1, 3)), -100.0),
    )
    monkeypatch.setattr(
        tsopt,
        "_frequencies_cm_and_modes",
        lambda *args, **kwargs: (
            np.array([-100.0, 20.0, 30.0]),
            torch.eye(3),
        ),
    )
    monkeypatch.setattr(
        tsopt,
        "_write_all_imaginary_modes",
        lambda *args, **kwargs: mode_exports.append(True),
    )
    monkeypatch.setattr(
        tsopt,
        "write",
        lambda path, _atoms: Path(path).write_text("final\n", encoding="utf-8"),
    )
    return runner, hessian_calls, mode_exports


def test_dimer_max_cycles_saves_final_structure_and_skips_phva(
    monkeypatch, tmp_path, capsys,
):
    runner, hessian_calls, mode_exports = _runner(
        tmp_path, monkeypatch, stalled=False
    )

    runner.run()

    assert len(hessian_calls) == 1
    assert mode_exports == []
    assert runner.n_imaginary_modes is None
    assert runner.hessian_status == "skipped"
    assert (tmp_path / "final_geometry.xyz").is_file()
    assert "ERROR: Not converged." not in capsys.readouterr().err


def test_dimer_plateau_saves_final_structure_and_skips_phva(
    monkeypatch, tmp_path, capsys,
):
    runner, hessian_calls, mode_exports = _runner(
        tmp_path, monkeypatch, stalled=True
    )

    runner.run()

    assert len(hessian_calls) == 1
    assert mode_exports == []
    assert runner.n_imaginary_modes is None
    assert runner.hessian_status == "skipped"
    assert (tmp_path / "final_geometry.xyz").is_file()
    assert "ERROR: Not converged." not in capsys.readouterr().err
    assert runner.is_stalled is True


@pytest.mark.parametrize("key,value", [
    ("raw_mode_count", 4), ("near_zero_frequencies_cm", [float("nan")]),
    ("near_zero_frequencies_cm", None),
])
def test_legacy_cache_without_complete_partition_is_not_reused(key, value):
    optimizer, geometry = _two_atom_exact_cache()
    optimizer._last_rigid_projection_info[key] = value
    assert tsopt._optimizer_exact_frequency_data(optimizer, geometry) is None


@pytest.mark.parametrize("near,strict,expected", [([-3.2], 2, False), ([3.2], 1, True)])
def test_final_dimer_verdict_uses_near_zero_partition(near, strict, expected):
    runner = SimpleNamespace(rigid_projection_info={
        "raw_mode_count": 3, "near_zero_frequencies_cm": near,
    })
    indices = tsopt._finalize_dimer_saddle_status(runner, np.array([-100., 12.]), 5.)
    assert runner.saddle_order_verified is expected
    assert runner.n_imaginary_modes == 1
    assert runner.n_negative_modes == strict
    assert indices.tolist() == [0]


@pytest.mark.parametrize("chosen,expected", [("primary", 2), ("alternate", 1)])
def test_selected_flatten_branch_owns_opposite_near_sign_metadata(chosen, expected):
    """Exercise the actual branch packet/restore expressions without running PES."""
    import ast
    from copy import deepcopy
    from pysisyphus.normal_modes import _strict_negative_count

    tree = ast.parse(Path(tsopt.__file__).read_text())
    branch = next(n for n in ast.walk(tree) if isinstance(n, ast.FunctionDef) and n.name == "_run_flatten_branch")
    packet_expr = next(n.value for n in branch.body if isinstance(n, ast.Return))
    live = {"raw_mode_count": 3, "near_zero_frequencies_cm": [-2.]}
    env = dict(deepcopy=deepcopy, rigid_projection_info=live, branch_ready=True,
               branch_label="primary", label="primary", branch_optimizer=None,
               geometry=SimpleNamespace(cart_coords=np.zeros(3)),
               branch_freqs=np.array([-100., 12.]), branch_modes=torch.zeros((2, 3)),
               branch_n_imag=1, branch_ims=[-100.], converged=True,
               safeguards={}, cycles=1, branch_micro_obj=None, branch_micro_cycles=0)
    expression = compile(ast.Expression(packet_expr), "<native branch packet>", "eval")
    primary = eval(expression, env)
    live["near_zero_frequencies_cm"][0] = 2.
    env.update(branch_label="alternate", label="alternate")
    alternate = eval(expression, env)
    selected = primary if chosen == "primary" else alternate
    update = next(
        n for n in ast.walk(tree) if isinstance(n, ast.Expr) and isinstance(n.value, ast.Call)
        and isinstance(n.value.func, ast.Attribute) and n.value.func.attr == "update"
        and ast.unparse(n.value.func.value) == "rigid_projection_info"
        and any(isinstance(child, ast.Subscript)
                and isinstance(child.value, ast.Name) and child.value.id == "selected_result"
                and isinstance(child.slice, ast.Constant) and child.slice.value == "projection"
                for child in ast.walk(n.value))
    )
    env["selected_result"] = selected
    live.clear()
    exec(compile(ast.Module(body=[update], type_ignores=[]), "<native branch restore>", "exec"), env)
    assert primary["projection"]["near_zero_frequencies_cm"] == [-2.]
    assert alternate["projection"]["near_zero_frequencies_cm"] == [2.]
    assert _strict_negative_count(selected["freqs"], live) == expected


@pytest.mark.parametrize("near", [[-2.], [2.]])
def test_legacy_packet_with_omitted_vectors_must_be_recomputed(near):
    optimizer, geometry = _two_atom_exact_cache()
    optimizer._last_rigid_projection_info.update({
        "raw_mode_count": 4, "near_zero_frequencies_cm": near,
    })
    assert tsopt._optimizer_exact_frequency_data(optimizer, geometry) is None


def test_complete_cache_reuses_soft_pairs_without_double_counting():
    from pysisyphus.normal_modes import frequency_partition_info
    optimizer, geometry = _two_atom_exact_cache()
    optimizer._last_exact_frequencies_cm = np.array([-100., -2., 2.])
    optimizer._last_rigid_projection_info.update(
        frequency_partition_info(optimizer._last_exact_frequencies_cm, 5.))
    reused = tsopt._optimizer_exact_frequency_data(optimizer, geometry)
    assert reused is not None
    np.testing.assert_array_equal(reused[0], [-100., -2., 2.])
    assert reused[1].shape == (3, 9)
    assert tsopt._strict_negative_count(reused[0], reused[2]) == 2
