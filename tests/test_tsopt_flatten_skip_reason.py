"""tsopt result.json records why a requested --flatten loop stopped or never ran."""

import json
from pathlib import Path

import numpy as np
import pytest
import torch
from click.testing import CliRunner

from pdb2reaction.cli import cli as root_cli
from pdb2reaction.workflows import tsopt
from pysisyphus.helpers import geom_loader
from pysisyphus.normal_modes import DEFAULT_FREQUENCY_ZERO_CUTOFF_CM

TWO_IMAGINARY = (-100.0, -50.0, 30.0)
ONE_IMAGINARY = (-100.0, 20.0, 30.0)


class _Geometry:
    atomic_numbers = np.array([1])
    atoms = ["H"]
    cart_coords = np.zeros(3)


# ---------------------------------------------------------------- Dimer runner

def _dimer_runner(
    tmp_path, monkeypatch, *, max_total_cycles, flatten_max_iter=1,
    freqs=TWO_IMAGINARY, flatten_moves=True,
):
    """A HessianDimer whose loops plateau after one cycle each."""
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
    runner.flatten_max_iter = flatten_max_iter
    runner.flatten_loop_bofill = False
    runner.skip_final_freq = False
    runner.neg_freq_thresh_cm = 5.0
    runner.tr_projection = "none"
    runner.rigid_projection_info = {}
    runner.max_total_cycles = max_total_cycles
    runner._cycles_spent = 0
    runner.is_converged = False
    runner.is_stalled = False
    runner.stop_reason = ""
    runner.flatten_skip_reason = None
    runner.n_imaginary_modes = None
    runner.imaginary_frequencies_cm = []
    runner.saddle_order_verified = False
    runner.prepared_input = None
    runner.ref_pdb = None
    runner.initial_hessian = None
    runner._store_ts_hessian = lambda _H: None

    def stalled_loop(_threshold):
        runner._cycles_spent += 1
        runner.is_stalled = True
        runner.stop_reason = "energy plateau"
        return 1, False, False

    runner._dimer_loop = stalled_loop
    runner._calc_full_hessian_cached = lambda *, allow_reuse: torch.eye(3)
    runner._flatten_once_with_modes = lambda _freqs, _modes: flatten_moves
    direction = lambda *args, **kwargs: (np.ones((1, 3)), -100.0)
    monkeypatch.setattr(tsopt, "_mode_direction_by_root", direction)
    monkeypatch.setattr(tsopt, "_mode_direction_by_root_from_Hact", direction)
    monkeypatch.setattr(
        tsopt,
        "_frequencies_cm_and_modes",
        lambda *args, **kwargs: (np.array(freqs), torch.eye(3)),
    )
    monkeypatch.setattr(tsopt, "_write_all_imaginary_modes", lambda *a, **k: 1)
    monkeypatch.setattr(
        tsopt,
        "write",
        lambda path, _atoms: Path(path).write_text("final\n", encoding="utf-8"),
    )
    return runner


@pytest.mark.parametrize(
    "case, expected",
    [
        ("budget_before", "max-cycles budget exhausted before flattening"),
        ("budget_during", "max-cycles budget exhausted during flattening"),
        ("no_eligible", "no eligible extra imaginary modes"),
        ("single_imaginary", None),
        ("not_requested", None),
    ],
)
def test_dimer_records_flatten_skip_reason(monkeypatch, tmp_path, case, expected):
    options = {
        "budget_before": dict(max_total_cycles=1),
        "budget_during": dict(max_total_cycles=2),
        "no_eligible": dict(max_total_cycles=10, flatten_moves=False),
        "single_imaginary": dict(max_total_cycles=10, freqs=ONE_IMAGINARY),
        "not_requested": dict(max_total_cycles=1, flatten_max_iter=0),
    }[case]
    runner = _dimer_runner(tmp_path, monkeypatch, **options)

    runner.run()

    assert runner.flatten_skip_reason == expected


class _FakeDimer:
    """Stand-in runner: the CLI must copy its skip reason into result.json."""

    def __init__(self, *, fn, out_dir, neg_freq_thresh_cm, uma_kwargs, **_kwargs):
        self.geom = geom_loader(fn, coord_type="cart")
        self.out_dir = Path(out_dir)
        self.neg_freq_thresh_cm = float(neg_freq_thresh_cm)
        self.uma_kwargs = dict(uma_kwargs)
        self._cycles_spent = 1
        self.is_converged = False
        self.is_stalled = False
        self.stop_reason = ""
        self.flatten_skip_reason = None
        self.saddle_order_verified = False
        self.n_imaginary_modes = None
        self.n_negative_modes = None
        self.imaginary_frequencies_cm = []
        self.hessian_status = "skipped"
        self.hessian_error = None
        self.rigid_projection_info = {}

    def run(self):
        self.flatten_skip_reason = "max-cycles budget exhausted before flattening"
        (self.out_dir / "final_geometry.xyz").write_text(
            self.geom.as_xyz(), encoding="utf-8"
        )


def test_dimer_result_json_carries_the_runner_skip_reason(monkeypatch, tmp_path):
    monkeypatch.setattr(tsopt, "HessianDimer", _FakeDimer)
    monkeypatch.setattr(tsopt, "_calc_energy", lambda *args, **kwargs: -1.0)
    source = tmp_path / "input.xyz"
    source.write_text("1\nfake dimer input\nHe 0 0 0\n", encoding="utf-8")
    out_dir = tmp_path / "out"

    result = CliRunner().invoke(
        root_cli,
        [
            "tsopt", "-i", str(source), "-q", "0", "-m", "1", "-o", str(out_dir),
            "--opt-mode", "dimer", "--flatten", "--max-cycles", "1",
            "--no-dump", "--out-json",
        ],
    )

    payload = json.loads((out_dir / "result.json").read_text(encoding="utf-8"))
    assert payload["flatten_requested"] is True, result.output
    assert payload["flatten_enabled"] is True
    assert payload["flatten_skip_reason"] == (
        "max-cycles budget exhausted before flattening"
    )


# ------------------------------------------------ Hessian-family (RS-P-RFO)

def _fake_ts_optimizer(*, first_cycles, target_mode_is_negative=None):
    """Optimizer stand-in: converges after ``first_cycles``; retries use all
    of their remaining budget."""

    class FakeTSOptimizer:
        instances = []

        def __init__(self, geometry, **kwargs):
            self.geometry = geometry
            self.max_cycles = kwargs.get("max_cycles")
            self.final_fn = Path(kwargs["out_dir"]) / "final_geometry.xyz"
            self.saddle_imaginary_threshold_cm = DEFAULT_FREQUENCY_ZERO_CUTOFF_CM
            self.forces = []
            FakeTSOptimizer.instances.append(self)

        def run(self):
            first = len(FakeTSOptimizer.instances) == 1
            used = first_cycles if first else int(self.max_cycles)
            self.cur_cycle = used - 1
            self.is_converged = True
            self.is_stalled = False
            self.stop_reason = ""
            self._last_exact_target_mode_is_negative = target_mode_is_negative
            self.final_fn.write_text(self.geometry.as_xyz(), encoding="utf-8")

    return FakeTSOptimizer


def _run_hessian_family(
    monkeypatch, tmp_path, *, max_cycles, first_cycles, flatten=True,
    freqs=TWO_IMAGINARY, flatten_moves=True, ref_mode=False,
    target_mode_is_negative=None,
):
    import pdb2reaction.io.hessian_cache as hessian_cache

    optimizer_cls = _fake_ts_optimizer(
        first_cycles=first_cycles,
        target_mode_is_negative=target_mode_is_negative,
    )
    monkeypatch.setitem(tsopt.TSOPT_CLASS_MAP, "rsprfo", optimizer_cls)
    monkeypatch.setattr(tsopt, "create_calculator", lambda **_kwargs: object())
    monkeypatch.setattr(
        tsopt, "_calc_full_hessian_torch", lambda *args, **kwargs: torch.eye(3)
    )
    monkeypatch.setattr(
        tsopt,
        "_frequencies_cm_and_modes",
        lambda *args, **kwargs: (np.array(freqs), torch.eye(3)),
    )
    monkeypatch.setattr(
        tsopt,
        "_flatten_once_with_modes_for_geom",
        lambda *args, **kwargs: flatten_moves,
    )
    monkeypatch.setattr(tsopt, "_write_all_imaginary_modes", lambda *a, **k: 0)
    monkeypatch.setattr(tsopt, "_calc_energy", lambda *args, **kwargs: -1.0)
    monkeypatch.setattr(hessian_cache, "store", lambda *args, **kwargs: None)
    monkeypatch.setattr(
        hessian_cache, "identity_from_context", lambda *args, **kwargs: {}
    )

    source = tmp_path / "input.xyz"
    source.write_text("1\nfake RS-P-RFO input\nHe 0 0 0\n", encoding="utf-8")
    out_dir = tmp_path / "out"
    args = [
        "tsopt", "-i", str(source), "-q", "0", "-m", "1", "-o", str(out_dir),
        "--opt-mode", "rsprfo", "--max-cycles", str(max_cycles),
        "--flatten" if flatten else "--no-flatten", "--no-dump", "--out-json",
    ]
    if ref_mode:
        mode_path = tmp_path / "ref_mode.npy"
        np.save(mode_path, np.array([1.0, 0.0, 0.0]))
        args += ["--ref-mode", str(mode_path)]

    result = CliRunner().invoke(root_cli, args)
    result_json = out_dir / "result.json"
    assert result_json.is_file(), result.output
    return json.loads(result_json.read_text(encoding="utf-8")), result


@pytest.mark.parametrize(
    "options, expected",
    [
        (
            dict(max_cycles=3, first_cycles=3),
            "max-cycles budget exhausted before flattening",
        ),
        (
            dict(max_cycles=3, first_cycles=1),
            "max-cycles budget exhausted during flattening",
        ),
        (
            dict(max_cycles=10, first_cycles=1, flatten_moves=False),
            "no eligible extra imaginary modes",
        ),
        (
            dict(max_cycles=10, first_cycles=1, ref_mode=True),
            "target mode sign never determined",
        ),
        (
            dict(
                max_cycles=10, first_cycles=1, ref_mode=True,
                target_mode_is_negative=False,
            ),
            "target mode is not negative",
        ),
        (dict(max_cycles=10, first_cycles=1, freqs=ONE_IMAGINARY), None),
    ],
    ids=[
        "budget_before", "budget_during", "no_eligible",
        "sign_undetermined", "target_positive", "single_imaginary",
    ],
)
def test_hessian_family_result_json_records_flatten_skip_reason(
    monkeypatch, tmp_path, options, expected,
):
    payload, result = _run_hessian_family(monkeypatch, tmp_path, **options)

    assert payload["flatten_requested"] is True, result.output
    assert payload["flatten_enabled"] is True
    assert payload["flatten_skip_reason"] == expected, result.output
    if expected is not None and expected.startswith("target mode"):
        assert f"Skipping extra-mode flattening: {expected}." in result.output


def test_hessian_family_without_flatten_reports_no_skip_reason(
    monkeypatch, tmp_path,
):
    payload, result = _run_hessian_family(
        monkeypatch, tmp_path, max_cycles=3, first_cycles=3, flatten=False,
    )

    assert payload["flatten_requested"] is False, result.output
    assert payload["flatten_enabled"] is False
    assert payload["flatten_skip_reason"] is None
