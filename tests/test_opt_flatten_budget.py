"""opt --flatten shares the command-level --max-cycles budget with its retries."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest
import torch
from click.testing import CliRunner

from pdb2reaction.workflows import opt as opt_workflow


class _ZeroCalculator:
    def get_energy(self, atoms, coords):
        return {"energy": 0.0}

    def get_forces(self, atoms, coords):
        return {"energy": 0.0, "forces": np.zeros(np.asarray(coords).size)}


def _run_opt(tmp_path, monkeypatch, *, max_cycles, runs, imag_counts):
    """Run the opt CLI with scripted (cycles, stalled) runs and imaginary-mode counts."""
    runs, imag_counts = iter(runs), iter(imag_counts)
    budgets, hessians = [], []

    class ScriptedOptimizer:
        def __init__(self, geometry, **kwargs):
            self.max_cycles = kwargs["max_cycles"]
            out_dir = Path(kwargs["out_dir"])
            self.final_fn = out_dir / "final_geometry.xyz"
            self.get_path_for_fn = lambda name: out_dir / name
            self.stop_reason = ""

        def run(self):
            budgets.append(self.max_cycles)
            cycles, stalled = next(runs)
            assert cycles <= self.max_cycles
            self.cur_cycle = cycles - 1
            self.is_stalled = stalled
            self.is_converged = not stalled and cycles < self.max_cycles

    def fake_hessian(geometry, *_args, **_kwargs):
        hessians.append(True)
        return torch.zeros(geometry.cart_coords.size, geometry.cart_coords.size)

    def fake_modes(hessian, *_args, **_kwargs):
        n_imag = next(imag_counts)
        freqs = np.array([-300.0] * n_imag + [300.0] * (hessian.shape[0] - n_imag))
        return freqs, torch.zeros(freqs.size, hessian.shape[0])

    monkeypatch.setattr(opt_workflow, "create_calculator", lambda **_k: _ZeroCalculator())
    monkeypatch.setattr(opt_workflow, "LBFGS", ScriptedOptimizer)
    monkeypatch.setattr(opt_workflow, "_calc_full_hessian_torch", fake_hessian)
    monkeypatch.setattr(opt_workflow, "_frequencies_cm_and_modes", fake_modes)
    monkeypatch.setattr(opt_workflow, "_flatten_all_imag_modes_for_geom", lambda *_a, **_k: True)
    source = tmp_path / "input.xyz"
    source.write_text("2\n\nH 0.0 0.0 0.0\nH 0.0 0.0 0.74\n", encoding="utf-8")
    out_dir = tmp_path / "opt"
    result = CliRunner().invoke(opt_workflow.cli, [
        "-i", str(source), "-q", "0", "-m", "1", "--opt-mode", "grad",
        "--max-cycles", str(max_cycles), "--flatten", "--out-json",
        "--out-dir", str(out_dir),
    ])
    assert result.exit_code == 0, result.output
    report = json.loads((out_dir / "result.json").read_text())
    return result, report, budgets, hessians


def test_stalled_flatten_retry_continues_on_the_remaining_budget(tmp_path, monkeypatch) -> None:
    result, report, budgets, hessians = _run_opt(
        tmp_path, monkeypatch, max_cycles=10,
        runs=[(3, False), (2, True), (2, False)], imag_counts=[2, 1, 0],
    )

    assert budgets == [10, 7, 5]
    assert len(hessians) == 3
    assert report["n_opt_cycles"] == 7
    assert report["status"] == "converged"
    assert "Remaining imaginary modes" not in result.output


@pytest.mark.parametrize(
    ("runs", "imag_counts", "budgets", "message"),
    [
        ([(4, False)], [], [4], "skipping flatten loop"),
        ([(1, False), (3, False)], [2, 1], [4, 3], "stopping flatten loop"),
    ],
    ids=["main-run", "retry"],
)
def test_flatten_stops_when_the_budget_is_spent(
    tmp_path, monkeypatch, runs, imag_counts, budgets, message
) -> None:
    result, report, spent_budgets, hessians = _run_opt(
        tmp_path, monkeypatch, max_cycles=4, runs=runs, imag_counts=imag_counts,
    )

    assert spent_budgets == budgets
    assert len(hessians) == len(imag_counts)
    assert f"Reached --max-cycles budget; {message}." in result.output
    assert report["n_opt_cycles"] == 4
    assert report["status"] == "not_converged"
