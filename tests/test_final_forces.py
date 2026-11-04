"""result.json final forces come from the final geometry, the final force and
step fields agree with run.log, and convergence_thresholds are the thresholds
of the optimizer that ran last."""

from __future__ import annotations

import json
import re
from types import SimpleNamespace

import numpy as np
import torch
from click.testing import CliRunner

from pysisyphus.Geometry import Geometry
from pysisyphus.calculators.Calculator import Calculator
from pysisyphus.optimizers.LBFGS import LBFGS
from pysisyphus.optimizers.Optimizer import CONV_THRESHS
from pdb2reaction.cli import cli as root_cli
from pdb2reaction.core.utils import optimizer_final_forces
from pdb2reaction.workflows import tsopt


class _QuadraticCalculator(Calculator):
    def __init__(self, out_dir):
        super().__init__(out_dir=out_dir, check_mem=False)

    @staticmethod
    def _results(coords):
        coords = np.asarray(coords, dtype=float)
        return {"energy": float(coords @ coords), "forces": -2.0 * coords}

    def get_energy(self, atoms, coords, **prepare_kwargs):
        return self._results(coords)

    def get_forces(self, atoms, coords, **prepare_kwargs):
        return self._results(coords)


def _run_lbfgs(tmp_path, atoms, coords, freeze_atoms=None, **kwargs):
    geom = Geometry(
        atoms, np.asarray(coords, dtype=float), coord_type="cart",
        freeze_atoms=freeze_atoms or [],
    )
    geom.set_calculator(_QuadraticCalculator(tmp_path))
    opt = LBFGS(geom, line_search=False, out_dir=tmp_path, **kwargs)
    opt.run()
    return geom, opt


def test_max_cycles_stop_reports_final_geometry_forces(tmp_path, capsys) -> None:
    geom, opt = _run_lbfgs(
        tmp_path, ["H", "H"], [1.0, 0.0, 0.0, 0.0, 0.5, 0.0],
        max_cycles=1, max_step=0.1,
    )
    assert not opt.is_converged
    # A step followed the only evaluation, so the history is one geometry behind.
    assert not np.array_equal(geom.cart_coords, opt.cart_coords[-1])

    final_forces = -2.0 * geom.cart_coords
    expected = (
        float(np.abs(final_forces).max()),
        float(np.sqrt(np.mean(final_forces**2))),
    )
    assert optimizer_final_forces(opt) == expected
    assert expected[0] != opt.max_forces[-1]
    # The JSON values are the ones the final summary prints.
    out = capsys.readouterr().out
    assert f"max(forces,cartesian): {expected[0]:.6f} hartree/bohr" in out
    assert f"rms(forces,cartesian): {expected[1]:.6f} hartree/bohr" in out


def test_max_cycles_stop_keeps_active_dof_basis_with_frozen_atoms(tmp_path) -> None:
    geom, opt = _run_lbfgs(
        tmp_path, ["H", "H"], [1.0, 0.0, 0.0, 0.0, 0.5, 0.0],
        freeze_atoms=[1], max_cycles=1, max_step=0.1,
    )
    active_forces = (-2.0 * geom.cart_coords)[:3]
    assert optimizer_final_forces(opt) == (
        float(np.abs(active_forces).max()),
        float(np.sqrt(np.mean(active_forces**2))),
    )


def test_converged_run_keeps_history_values(tmp_path) -> None:
    geom, opt = _run_lbfgs(
        tmp_path, ["H"], [1.0e-4, 0.0, 0.0], thresh="gau_loose", max_cycles=5,
    )
    assert opt.is_converged
    assert opt.max_forces[-1] > 0.0
    assert optimizer_final_forces(opt) == (opt.max_forces[-1], opt.rms_forces[-1])


def test_final_forces_fall_back_to_history_without_final_summary() -> None:
    stopped = SimpleNamespace(
        final_max_force=None, final_rms_force=None,
        max_forces=[0.3, 0.2], rms_forces=[0.1, 0.05],
    )
    assert optimizer_final_forces(stopped) == (0.2, 0.05)
    assert optimizer_final_forces(SimpleNamespace(max_forces=[], rms_forces=[])) is None
    assert optimizer_final_forces(None) is None


# E = 1/2 sum_i k_i (x_i - x0_i)^2 (Hartree, Bohr) with one negative k: a saddle.
_X0 = np.array([
    [0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 2.0, 0.0],
    [0.0, 0.0, 2.0], [2.0, 2.0, 2.0], [3.0, 1.0, 2.0],
]).ravel()
_CURVATURES = np.arange(1.0, 19.0) * 0.05
_CURVATURES[9] *= -1.0


class _QuadraticSaddle:
    def get_energy(self, atoms, coords):
        d = np.asarray(coords, dtype=float).reshape(-1) - _X0
        return {"energy": float(0.5 * np.dot(_CURVATURES * d, d))}

    def get_forces(self, atoms, coords):
        d = np.asarray(coords, dtype=float).reshape(-1) - _X0
        return {**self.get_energy(atoms, coords), "forces": -_CURVATURES * d}


def _last_final_summary(log):
    """(max, rms) forces printed in the last Final summary."""
    return (
        re.findall(r"max\(forces,cartesian\): (\S+)", log)[-1],
        re.findall(r"rms\(forces,cartesian\): (\S+)", log)[-1],
    )


def _last_cycle_row(log):
    """(max force, rms force, max step, rms step) of the last cycle row."""
    rows = [
        line.split() for line in log.splitlines()
        if len(line.split()) == 7 and line.split()[0].isdigit()
    ]
    return tuple(value.rstrip("*") for value in rows[-1][2:6])


def test_dimer_result_json_final_fields_match_run_log(tmp_path, monkeypatch):
    mode = np.zeros((6, 3))
    mode[3, 0] = 1.0
    monkeypatch.setattr(tsopt, "create_calculator", lambda **_kwargs: _QuadraticSaddle())
    monkeypatch.setattr(
        tsopt, "_calc_full_hessian_torch",
        lambda *_args, **_kwargs: torch.diag(torch.as_tensor(_CURVATURES)),
    )
    monkeypatch.setattr(
        tsopt, "_mode_direction_by_root", lambda *_args, **_kwargs: (mode.copy(), -100.0)
    )
    monkeypatch.setattr(tsopt, "_calc_energy", lambda *_args, **_kwargs: -1.0)
    start = _X0.copy()
    start[[9, 10, 15, 16]] += [0.4, 0.5, -0.6, 0.7]
    source = tmp_path / "input.xyz"
    source.write_text(
        "6\nquadratic saddle\n" + "".join(
            f"H {x:.8f} {y:.8f} {z:.8f}\n"
            for x, y, z in start.reshape(-1, 3) * 0.529177210903
        ),
        encoding="utf-8",
    )
    out_dir = tmp_path / "out"

    result = CliRunner().invoke(root_cli, [
        "tsopt", "-i", str(source), "-q", "0", "-m", "1", "-o", str(out_dir),
        "--opt-mode", "dimer", "--max-cycles", "3", "--print-every", "1",
        "--no-dump", "--out-json",
    ])

    report = json.loads((out_dir / "result.json").read_text(encoding="utf-8"))
    log = (out_dir / "run.log").read_text(encoding="utf-8")
    assert result.exit_code == 1, result.output
    assert report["execution_status"] == "completed"
    assert report["optimization_status"] == "not_converged"
    assert "Number of cycles exceeded!" in log
    # Forces of the final geometry, as the last LBFGS segment's Final summary
    # prints them; the step that produced it, as the last cycle row prints it.
    assert (
        f"{report['final_max_force']:.6f}", f"{report['final_rms_force']:.6f}",
    ) == _last_final_summary(log)
    assert (
        f"{report['final_max_step']:.6f}", f"{report['final_rms_step']:.6f}",
    ) == _last_cycle_row(log)[2:]
    # The cycle budget ends in the loose loop, so its LBFGS thresholds were used.
    assert "Reached --max-cycles budget after loose loop" in log
    keys = ("max_force_thresh", "rms_force_thresh", "max_step_thresh", "rms_step_thresh")
    assert report["convergence_thresholds"] == {
        key: float(value) for key, value in zip(keys, CONV_THRESHS["gau_loose"])
    }
