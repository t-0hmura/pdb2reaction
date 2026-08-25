"""An energy-only plateau is an additive 'stalled' outcome,
never convergence, and it stops further retry work.

These bind to the production optimizer/renderer code paths (no reimplemented
gate logic): the bundled ``pysisyphus`` Optimizer/RFOptimizer/TSHessianOptimizer
state machine and the product-local status helpers in
``pdb2reaction.core.utils`` / ``pdb2reaction.workflows.tsopt``.
"""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from pysisyphus.Geometry import Geometry
from pysisyphus.calculators.Calculator import Calculator
from pysisyphus.optimizers.LBFGS import LBFGS
from pysisyphus.optimizers.RFOptimizer import RFOptimizer

from pdb2reaction.core.utils import (
    optimizer_terminal_status,
    emit_optimizer_terminal_status,
)
from pdb2reaction.workflows.tsopt import (
    HessianDimer,
    _hessian_postprocessing_is_ready,
    _hessian_result_status,
    _thresholded_reaction_mode_index,
    _tsopt_terminal_outcome_message,
    _tsopt_terminal_status,
)


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


class _ConstantPlateauCalculator(Calculator):
    """Constant energy (an energy plateau) with a constant large force so the
    configured force/step criteria can never be met."""

    def __init__(self, out_dir):
        super().__init__(out_dir=out_dir, check_mem=False)

    @staticmethod
    def _results(coords):
        coords = np.asarray(coords, dtype=float)
        return {
            "energy": 1.0,
            "forces": np.array([0.5, 0.0, 0.0]),
            "hessian": np.eye(len(coords)),
        }

    def get_energy(self, atoms, coords, **prepare_kwargs):
        return self._results(coords)

    def get_forces(self, atoms, coords, **prepare_kwargs):
        return self._results(coords)

    def get_hessian(self, atoms, coords, **prepare_kwargs):
        return self._results(coords)


def _lbfgs(tmp_path, **kwargs):
    geom = Geometry(["H"], np.zeros(3), coord_type="cart")
    geom.set_calculator(_QuadraticCalculator(tmp_path))
    opt_kwargs = dict(
        thresh="gau",
        energy_plateau=True,
        energy_plateau_window=3,
        out_dir=tmp_path,
    )
    opt_kwargs.update(kwargs)
    return geom, LBFGS(geom, **opt_kwargs)


def _seed(opt, *, force, step, energies):
    opt.cur_cycle = 5
    opt.last_cycle = 0
    opt.forces = [np.asarray(force, dtype=float)]
    opt.steps = [np.asarray(step, dtype=float)]
    opt.energies = list(energies)


_HIGH = np.full(3, 1.0)
_LOW = np.full(3, 1.0e-8)
_FLAT = [1.0, 1.0, 1.0]  # range 0 over window 3 -> a plateau


# ---- Base state machine: force/step truth table under a flat energy window ----

@pytest.mark.parametrize(
    "force, step, expect_converged",
    [
        (_HIGH, _HIGH, False),   # high force, high step   -> stalled
        (_HIGH, _LOW, False),    # high force, low step    -> stalled
        (_LOW, _HIGH, False),    # low force, high step    -> stalled
        (_LOW, _LOW, True),      # all below thresh        -> converged (wins)
    ],
)
def test_energy_plateau_truth_table(tmp_path, force, step, expect_converged):
    _, opt = _lbfgs(tmp_path)
    _seed(opt, force=force, step=step, energies=_FLAT)

    converged, conv_info = opt.check_convergence(step)

    assert bool(converged) is expect_converged
    if expect_converged:
        # A real convergence wins over the plateau; no stall. (check_convergence
        # returns the bool; the run loop is what assigns self.is_converged.)
        assert opt.is_stalled is False
    else:
        # A plateau without real convergence stalls, never converges.  The
        # stall side effect (request_stall) sets the committed terminal state.
        assert opt.is_stalled is True
        assert opt.is_converged is False
        assert opt.stop_requested is True
        assert opt.termination_status == "stalled"
        assert "energy plateau" in opt.stop_reason
        # A plateau does not alter the ConvInfo force/step fields.
        assert bool(conv_info.max_force_converged) == bool(np.max(np.abs(force)) < 1e-3)


def test_energy_range_just_above_threshold_neither_converges_nor_stalls(tmp_path):
    _, opt = _lbfgs(tmp_path)
    # range 2e-5 over the window > energy_plateau_thresh (1e-5): not a plateau.
    _seed(opt, force=_HIGH, step=_HIGH, energies=[1.0, 1.0, 1.0 + 2.0e-5])

    converged, _ = opt.check_convergence(_HIGH)

    assert converged is False
    assert opt.is_stalled is False
    assert opt.termination_status == "not_converged"


def test_energy_plateau_disabled_does_not_stall(tmp_path):
    _, opt = _lbfgs(tmp_path, energy_plateau=False)
    _seed(opt, force=_HIGH, step=_HIGH, energies=_FLAT)

    converged, _ = opt.check_convergence(_HIGH)

    assert converged is False
    assert opt.is_stalled is False


def test_too_few_energy_samples_does_not_stall(tmp_path):
    _, opt = _lbfgs(tmp_path, energy_plateau_window=5)
    _seed(opt, force=_HIGH, step=_HIGH, energies=[1.0, 1.0])  # < window

    converged, _ = opt.check_convergence(_HIGH)

    assert converged is False
    assert opt.is_stalled is False


def test_thresh_never_neither_converges_nor_stalls(tmp_path):
    _, opt = _lbfgs(tmp_path, thresh="never")
    _seed(opt, force=_LOW, step=_LOW, energies=_FLAT)

    converged, _ = opt.check_convergence(_LOW)

    assert converged is False
    assert opt.is_stalled is False


# ---- RFOptimizer provisional probe must not stall; the final check must ------

def test_rfo_provisional_probe_suppresses_stall_but_final_check_stalls(tmp_path):
    geom = Geometry(["H"], np.zeros(3), coord_type="cart")
    geom.set_calculator(_QuadraticCalculator(tmp_path))
    opt = RFOptimizer(
        geom,
        thresh="gau",
        energy_plateau=True,
        energy_plateau_window=3,
        out_dir=tmp_path,
    )
    _seed(opt, force=_HIGH, step=_HIGH, energies=_FLAT)

    # Provisional probe (RFOptimizer's ref_step probe uses this) must NOT stall.
    provisional, _ = opt.check_convergence(_HIGH, allow_stall=False)
    assert provisional is False
    assert opt.is_stalled is False

    # The subsequent run-loop final check (allow_stall default True) stalls.
    final, _ = opt.check_convergence(_HIGH)
    assert final is False
    assert opt.is_stalled is True
    assert opt.termination_status == "stalled"


# ---- Run-loop integration: stall, no further step, terminal says stalled -----

def test_run_loop_stalls_on_energy_plateau(tmp_path, capsys):
    geom = Geometry(["H"], np.zeros(3), coord_type="cart")
    geom.set_calculator(_ConstantPlateauCalculator(tmp_path))
    opt = LBFGS(
        geom,
        thresh="gau",
        energy_plateau=True,
        energy_plateau_window=2,
        max_cycles=20,
        out_dir=tmp_path,
    )
    opt.run()

    assert opt.is_converged is False
    assert opt.is_stalled is True
    assert opt.stop_requested is True
    assert opt.stopped is True
    assert opt.termination_status == "stalled"
    # Exited well before max_cycles (stalled after the window filled).
    assert opt.cur_cycle < 19
    out = capsys.readouterr().out
    assert "Stalled" in out
    assert "Converged!" not in out


# ---- Product-local status renderers ------------------------------------------

class _FakeOpt:
    def __init__(self, *, is_converged=False, is_stalled=False, stop_reason=""):
        self.is_converged = is_converged
        self.is_stalled = is_stalled
        self.stop_reason = stop_reason

    @property
    def termination_status(self):
        if self.is_stalled:
            return "stalled"
        if self.is_converged:
            return "converged"
        return "not_converged"


def test_optimizer_terminal_status_maps_stalled_and_converged():
    assert optimizer_terminal_status(_FakeOpt(is_stalled=True)) == "stalled"
    assert optimizer_terminal_status(_FakeOpt(is_converged=True)) == "converged"
    assert optimizer_terminal_status(_FakeOpt()) == "not_converged"


def test_stalled_optimizer_never_reports_converged():
    stalled = _FakeOpt(is_stalled=True, is_converged=False, stop_reason="energy plateau: ...")
    assert optimizer_terminal_status(stalled) == "stalled"
    assert optimizer_terminal_status(stalled) != "converged"
    # A stalled TS optimizer is never a converged saddle, even at n_imag==1.
    assert _tsopt_terminal_status(stalled, saddle_verified=True) == "stalled"


def test_tsopt_terminal_status_composition():
    converged = _FakeOpt(is_converged=True)
    assert _tsopt_terminal_status(converged, saddle_verified=True) == "converged"
    # Numerical convergence is independent of saddle order.
    assert _tsopt_terminal_status(converged, saddle_verified=False) == "converged"
    assert _tsopt_terminal_status(_FakeOpt(), saddle_verified=True) == "not_converged"


def test_tsopt_terminal_outcome_messages_separate_numerical_and_saddle_status():
    first_order_not_converged = _tsopt_terminal_outcome_message(
        numerically_converged=False,
        hessian_ready=True,
        n_imaginary_modes=1,
    )
    assert first_order_not_converged == "[tsopt] ERROR: Not converged."

    converged_higher_order = _tsopt_terminal_outcome_message(
        numerically_converged=True,
        hessian_ready=True,
        n_imaginary_modes=2,
    )
    assert converged_higher_order == (
        "[tsopt] WARNING: Higher-order stationary point (n_imag=2). "
        "Try --flatten or all --refine-path."
    )

    converged_minimum = _tsopt_terminal_outcome_message(
        numerically_converged=True,
        hessian_ready=True,
        n_imaginary_modes=0,
    )
    assert converged_minimum == (
        "[tsopt] No imaginary mode detected. Try all --refine-path."
    )

    validated = _tsopt_terminal_outcome_message(
        numerically_converged=True,
        hessian_ready=True,
        n_imaginary_modes=1,
    )
    assert validated == "[tsopt] Converged (n_imag=1)."

    unavailable = _tsopt_terminal_outcome_message(
        numerically_converged=True,
        hessian_ready=False,
        n_imaginary_modes=None,
    )
    assert unavailable == "[tsopt] ERROR: Failed to complete terminal PHVA."


def test_reaction_mode_selection_uses_the_configured_saddle_threshold():
    frequencies = np.array([-450.0, -3.2, 20.0])
    assert _thresholded_reaction_mode_index(frequencies, 5.0, 1) == 0
    assert _thresholded_reaction_mode_index(np.array([-3.2, 20.0]), 5.0) is None


def test_hessian_postprocessing_requires_numerical_convergence():
    assert _hessian_postprocessing_is_ready(None) is False
    assert _hessian_postprocessing_is_ready(_FakeOpt()) is False
    assert _hessian_postprocessing_is_ready(_FakeOpt(is_converged=True)) is True
    assert _hessian_postprocessing_is_ready(_FakeOpt(is_stalled=True)) is False

    already_failed = _FakeOpt()
    already_failed._last_exact_failure_reason = "RuntimeError: failed"
    assert _hessian_postprocessing_is_ready(already_failed) is False


def test_hessian_result_status_marks_nonconverged_phva_as_skipped():
    assert _hessian_result_status(
        n_imaginary=None,
        hessian_error=None,
        postprocessing_ready=False,
    ) == "skipped"
    assert _hessian_result_status(
        n_imaginary=None,
        hessian_error="RuntimeError: failed",
        postprocessing_ready=False,
    ) == "failed"
    assert _hessian_result_status(
        n_imaginary=1,
        hessian_error="RuntimeError: failed",
        postprocessing_ready=True,
    ) == "failed"
    assert _hessian_result_status(
        n_imaginary=1,
        hessian_error=None,
        postprocessing_ready=True,
    ) == "completed"


def test_emit_terminal_status_stalled_and_converged_are_distinct(capsys):
    emit_optimizer_terminal_status(
        "opt", converged=False, cycles=3, max_cycles=20,
        stalled=True, stop_reason="energy plateau: range=0.00e+00",
    )
    stalled_out = capsys.readouterr().out
    assert "Stalled" in stalled_out
    assert "Converged!" not in stalled_out

    emit_optimizer_terminal_status(
        "opt", converged=True, cycles=7, max_cycles=20,
        converged_message="Numerical optimization converged.",
    )
    conv_out = capsys.readouterr().out
    assert "Numerical optimization converged." in conv_out
    assert "Converged!" not in conv_out
    assert "Stalled" not in conv_out


# ---- HessianDimer wrapper: a stalled child stops all further work -------------

def test_hessian_dimer_stops_after_child_stall(tmp_path):
    """A stalled child LBFGS makes the runner stalled and stops the segment
    loop before any further segment or Hessian update."""
    runner = HessianDimer.__new__(HessianDimer)
    runner.max_total_cycles = 100
    runner._cycles_spent = 0
    runner.update_interval_hessian = 5
    runner.is_stalled = False
    runner.is_converged = False
    runner.stop_reason = ""

    class _Geom:
        cart_coords = np.zeros(6)

    runner.geom = _Geom()

    calls = {"segments": 0, "hessian": 0}

    def _fake_segment(threshold, n_steps):
        calls["segments"] += 1
        # The child LBFGS stalled on an energy plateau.
        runner.is_stalled = True
        runner.stop_reason = "energy plateau: range=0.00e+00 au over 2 steps"
        return 3, False  # (steps, converged=False)

    def _fake_hessian(allow_reuse):
        calls["hessian"] += 1
        raise AssertionError("no Hessian update should run after a stall")

    runner._dimer_segment = _fake_segment
    runner._calc_full_hessian_cached = _fake_hessian

    steps, zero_step_converged, loop_converged = runner._dimer_loop("gau")

    assert runner.is_stalled is True
    assert runner.termination_status == "stalled"
    assert loop_converged is False
    assert calls["segments"] == 1          # only the stalled segment ran
    assert calls["hessian"] == 0           # no further Hessian/segment work
    assert steps == 3
    # And the public status mapper reports stalled, never converged.
    assert _tsopt_terminal_status(runner, saddle_verified=True) == "stalled"


def test_terminal_saddle_certification_uses_magnitude_threshold():
    """Certification ignores negative roots softer than the configured gate."""
    from pdb2reaction.workflows.tsopt import (
        _certified_negative_frequencies,
        _certified_saddle_order,
        _finalize_dimer_saddle_status,
        _imaginary_mode_indices_and_values,
    )

    # One strong imaginary mode plus a numerically soft negative root.
    freqs_cm = np.array([-450.0, -3.2, 12.0, 640.0])
    reported_idx, reported_values = _imaginary_mode_indices_and_values(freqs_cm, 5.0)

    # Mode export, printed output, and certification use one threshold.
    assert len(reported_idx) == 1
    assert reported_values == [-450.0]
    assert _certified_saddle_order(freqs_cm, 5.0) == 1
    assert _certified_negative_frequencies(freqs_cm, 5.0) == [-450.0]

    runner = _FakeOpt()
    runner.is_converged = True
    export_idx = _finalize_dimer_saddle_status(runner, freqs_cm, 5.0)

    # The public result is self-consistent with the thresholded certified set.
    assert runner.n_imaginary_modes == 1
    assert runner.imaginary_frequencies_cm == [-450.0]
    assert runner.saddle_order_verified is True
    assert runner.is_converged is True
    assert _tsopt_terminal_status(runner, saddle_verified=True) == "converged"
    # Mode export still follows the threshold-filtered indices.
    assert export_idx.tolist() == reported_idx.tolist()

    # A genuine first-order saddle still certifies.
    single = _FakeOpt()
    single.is_converged = True
    _finalize_dimer_saddle_status(single, np.array([-450.0, 12.0, 640.0]), 5.0)
    assert single.n_imaginary_modes == 1
    assert single.imaginary_frequencies_cm == [-450.0]
    assert single.saddle_order_verified is True
    assert single.is_converged is True
    assert _tsopt_terminal_status(single, saddle_verified=True) == "converged"

    # A lone sub-threshold negative root does not certify a transition state.
    soft = _FakeOpt()
    soft.is_converged = True
    soft_export = _finalize_dimer_saddle_status(soft, np.array([-3.2, 12.0, 640.0]), 5.0)
    assert soft.n_imaginary_modes == 0
    assert soft.imaginary_frequencies_cm == []
    assert soft.saddle_order_verified is False
    assert soft.is_converged is True
    assert soft_export.tolist() == []


def test_exact_phva_validation_ignores_soft_negative_roots():
    """The bundled exact PHVA branch reports order 1 for [-450, -3.2, +12]."""
    from pysisyphus.tsoptimizers.RSIRFOptimizer import RSIRFOptimizer

    modes = np.eye(3)
    printed: list[str] = []
    frequencies = {"value": np.array([-450.0, -3.2, 12.0])}

    # TSHessianOptimizer is abstract; RSIRFO is the shipped concrete owner.
    optimizer = RSIRFOptimizer.__new__(RSIRFOptimizer)
    # The shipped magnitude threshold is also the certification threshold.
    optimizer.saddle_imaginary_threshold_cm = 5.0
    optimizer.small_eigval_thresh = 1e-8
    optimizer.roots = np.array([0])
    optimizer.reference_mode = None
    optimizer.cur_cycle = 7
    optimizer.higher_order_saddle_checks = 0
    optimizer.max_higher_order_checks = 99
    optimizer.forces = []
    optimizer.geometry = SimpleNamespace(cart_coords=np.zeros(3))
    optimizer.table = SimpleNamespace(print=printed.append)
    optimizer._mw_frequencies_and_modes = lambda: (frequencies["value"], modes)
    optimizer._recovery_mode_from_mw = lambda _modes, index: modes[:, int(index)]
    optimizer._record_exact_saddle_candidate = lambda: None
    optimizer.request_stop = lambda *_a, **_k: None

    has_saddle_modes, _mode, _has_mode = (
        optimizer._verify_exact_vibrational_structure(
            np.array([-0.1, -1.0e-7, 0.2]), np.eye(3)
        )
    )

    assert optimizer._last_exact_n_imaginary == 1
    assert optimizer._last_exact_saddle_verified is True
    assert has_saddle_modes is True
    assert any("n_imag=1" in message for message in printed)

    # A single negative root still certifies first order.
    printed.clear()
    optimizer.higher_order_saddle_checks = 0
    frequencies["value"] = np.array([-450.0, 12.0, 88.0])
    optimizer._verify_exact_vibrational_structure(
        np.array([-0.1, 0.1, 0.2]), np.eye(3)
    )
    assert optimizer._last_exact_n_imaginary == 1
    assert optimizer._last_exact_saddle_verified is True
    assert any("n_imag=1" in message for message in printed)


def test_dimer_final_message_separates_no_mode_from_write_failure():
    from pdb2reaction.workflows.tsopt import _no_exported_mode_message

    assert _no_exported_mode_message(0, 5.0, 12.0) == (
        "[tsopt] No imaginary mode detected. Try all --refine-path."
    )
    assert _no_exported_mode_message(1, 5.0, -100.0) == (
        "[tsopt] ERROR: Failed to write imaginary mode trajectory."
    )
