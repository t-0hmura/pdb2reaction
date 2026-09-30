"""tsopt reports why --flatten did not run and which nested YAML sections it ignores."""

import json
from pathlib import Path

import numpy as np
import pytest
import torch
from click.testing import CliRunner

from pdb2reaction.cli import cli as root_cli
from pdb2reaction.core.utils import unused_nested_yaml_sections
from pdb2reaction.workflows import tsopt
from pysisyphus.normal_modes import DEFAULT_FREQUENCY_ZERO_CUTOFF_CM

TWO_IMAGINARY = (-100.0, -50.0, 30.0)
BUDGET_BEFORE = "max-cycles budget exhausted before flattening"
BUDGET_DURING = "max-cycles budget exhausted during flattening"
FINAL_FREQ_SKIPPED = "final Hessian skipped (--skip-final-freq)"


def _fake_ts_optimizer(runs):
    """Optimizer stand-in; ``runs[i] = (cycles, converged)`` for the i-th run,
    where ``cycles=None`` spends the whole ``max_cycles`` given to that run."""

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
            cycles, converged = runs[len(FakeTSOptimizer.instances) - 1]
            used = int(self.max_cycles) if cycles is None else cycles
            self.cur_cycle = used - 1
            self.is_converged = converged
            self.is_stalled = False
            self.stop_reason = ""
            self._last_exact_target_mode_is_negative = None
            self.final_fn.write_text(self.geometry.as_xyz(), encoding="utf-8")

    return FakeTSOptimizer


def _run_rsprfo(monkeypatch, tmp_path, *, runs, max_cycles, extra_args):
    import pdb2reaction.io.hessian_cache as hessian_cache

    phva_calls = []

    def frequencies(*args, **kwargs):
        phva_calls.append(1)
        return np.array(TWO_IMAGINARY), torch.eye(3)

    monkeypatch.setitem(tsopt.TSOPT_CLASS_MAP, "rsprfo", _fake_ts_optimizer(runs))
    monkeypatch.setattr(tsopt, "create_calculator", lambda **_kwargs: object())
    monkeypatch.setattr(
        tsopt, "_calc_full_hessian_torch", lambda *args, **kwargs: torch.eye(3)
    )
    monkeypatch.setattr(tsopt, "_frequencies_cm_and_modes", frequencies)
    monkeypatch.setattr(
        tsopt, "_flatten_once_with_modes_for_geom", lambda *args, **kwargs: True
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
    result = CliRunner().invoke(root_cli, [
        "tsopt", "-i", str(source), "-q", "0", "-m", "1", "-o", str(out_dir),
        "--opt-mode", "rsprfo", "--max-cycles", str(max_cycles),
        "--no-dump", "--out-json", *extra_args,
    ])
    result_json = out_dir / "result.json"
    assert result_json.is_file(), result.output
    return json.loads(result_json.read_text(encoding="utf-8")), result, phva_calls


def test_unconverged_max_cycles_stop_records_budget_reason(monkeypatch, tmp_path):
    payload, result, phva_calls = _run_rsprfo(
        monkeypatch, tmp_path, runs=[(None, False)], max_cycles=3,
        extra_args=["--flatten"],
    )

    assert payload["flatten_requested"] is True, result.output
    assert payload["flatten_skip_reason"] == BUDGET_BEFORE, result.output
    assert "Reached --max-cycles budget; skipping flatten loop." in result.output
    # S2: an unconverged max-cycles stop computes no Hessian.
    assert payload["hessian_status"] == "skipped"
    assert phva_calls == []


def test_skip_final_freq_records_skipped_hessian_reason(monkeypatch, tmp_path):
    payload, result, phva_calls = _run_rsprfo(
        monkeypatch, tmp_path, runs=[(1, True)], max_cycles=10,
        extra_args=["--flatten", "--skip-final-freq"],
    )

    assert payload["optimization_status"] == "converged", result.output
    assert payload["flatten_skip_reason"] == FINAL_FREQ_SKIPPED
    assert phva_calls == []


def test_unconverged_flatten_retry_records_budget_reason(monkeypatch, tmp_path):
    payload, result, phva_calls = _run_rsprfo(
        monkeypatch, tmp_path, runs=[(1, True), (None, False)], max_cycles=3,
        extra_args=["--flatten"],
    )

    assert payload["flatten_skip_reason"] == BUDGET_DURING, result.output
    assert len(phva_calls) == 1


@pytest.mark.parametrize(
    "runs, extra_args",
    [([(None, False)], []), ([(1, True)], ["--skip-final-freq"])],
    ids=["max_cycles", "skip_final_freq"],
)
def test_no_reason_without_flatten(monkeypatch, tmp_path, runs, extra_args):
    payload, result, _ = _run_rsprfo(
        monkeypatch, tmp_path, runs=runs, max_cycles=3,
        extra_args=["--no-flatten", *extra_args],
    )

    assert payload["flatten_requested"] is False, result.output
    assert payload["flatten_skip_reason"] is None


# ------------------------------------------------------------ YAML NOTE

def test_unused_nested_sections_are_listed():
    yaml_cfg = {
        "opt": {"max_cycles": 5, "lbfgs": {"max_step": 0.1}, "rfo": {"trust_radius": 0.2}},
        "freq": {"thermo": {"temperature": 310.0}},
    }
    assert unused_nested_yaml_sections(yaml_cfg) == ["opt.lbfgs", "opt.rfo", "freq.thermo"]


def test_read_or_empty_nested_sections_are_not_listed():
    yaml_cfg = {"opt": {"lbfgs": {"max_step": 0.1}, "rfo": None}, "freq": {"thermo": {}}}
    assert unused_nested_yaml_sections(yaml_cfg, read=(("opt", "lbfgs"),)) == []
    assert unused_nested_yaml_sections({"geom": {"coord_type": "cart"}}) == []


def _tsopt_dry_run(tmp_path, config):
    source = tmp_path / "atom.xyz"
    source.write_text("1\natom\nHe 0 0 0\n", encoding="utf-8")
    path = tmp_path / "config.yaml"
    path.write_text(json.dumps(config), encoding="utf-8")
    return CliRunner().invoke(root_cli, [
        "tsopt", "-i", str(source), "-q", "0", "-m", "1", "--opt-mode", "rsprfo",
        "--dry-run", "--config", str(path), "-o", str(tmp_path / "out"),
    ])


def test_tsopt_notes_ignored_nested_sections_once(tmp_path):
    result = _tsopt_dry_run(tmp_path, {
        "opt": {"max_cycles": 7, "lbfgs": {"max_step": 0.1}, "rfo": {"trust_radius": 0.2}},
        "freq": {"thermo": {"temperature": 310.0}},
    })

    assert result.exit_code == 0, result.output
    notes = [line for line in result.output.splitlines() if "NOTE: Ignoring YAML" in line]
    assert notes == [
        "[tsopt] NOTE: Ignoring YAML sections that tsopt does not use: "
        "opt.lbfgs, opt.rfo, freq.thermo."
    ]


def test_tsopt_prints_no_note_without_nested_sections(tmp_path):
    result = _tsopt_dry_run(tmp_path, {"opt": {"max_cycles": 7}})

    assert result.exit_code == 0, result.output
    assert "NOTE: Ignoring YAML" not in result.output
