"""--skip-final-freq (tsopt, all) and --hess-device (freq, irc)."""

import json
from pathlib import Path

import numpy as np
import pytest
import torch
import yaml
from click.testing import CliRunner

from pdb2reaction.cli import cli as root_cli
from pdb2reaction.workflows import freq as freq_workflow
from pdb2reaction.workflows import irc as irc_workflow
from pdb2reaction.workflows import tsopt
from pysisyphus.normal_modes import DEFAULT_FREQUENCY_ZERO_CUTOFF_CM

ONE_IMAGINARY = (-100.0, 20.0, 30.0)
SKIP_WARNING = "TS saddle-point order is not verified (--skip-final-freq)."


def _helium(tmp_path: Path, name: str = "input.xyz") -> Path:
    source = tmp_path / name
    source.write_text("1\nhelium\nHe 0 0 0\n", encoding="utf-8")
    return source


# ------------------------------------------------------------------ tsopt CLI

def test_tsopt_rejects_dump_hess_with_skip_final_freq(tmp_path):
    result = CliRunner().invoke(root_cli, [
        "tsopt", "-i", str(_helium(tmp_path)), "-q", "0", "-m", "1",
        "-o", str(tmp_path / "out"), "--skip-final-freq",
        "--dump-hess", str(tmp_path / "ts.npy"),
    ])

    assert result.exit_code == 2, result.output
    assert "--dump-hess needs the final Hessian; drop --skip-final-freq." in result.output


# --------------------------------------------------------------- Dimer runner

class _Geometry:
    atomic_numbers = np.array([1])
    atoms = ["H"]
    cart_coords = np.zeros(3)


def _dimer_runner(tmp_path, monkeypatch, *, stalled, skip_final_freq):
    """A HessianDimer whose single loop converges or stops on a plateau."""
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
    runner.flatten_max_iter = 0
    runner.flatten_loop_bofill = False
    runner.neg_freq_thresh_cm = 5.0
    runner.tr_projection = "none"
    runner.rigid_projection_info = {}
    runner.max_total_cycles = 10
    runner._cycles_spent = 0
    runner.is_converged = False
    runner.is_stalled = False
    runner.stop_reason = ""
    runner.flatten_skip_reason = None
    runner.saddle_order_verified = False
    runner.n_imaginary_modes = None
    runner.n_negative_modes = None
    runner.imaginary_frequencies_cm = []
    runner.hessian_status = "not_run"
    runner.hessian_error = None
    runner.prepared_input = None
    runner.ref_pdb = None
    runner.initial_hessian = None
    runner.skip_final_freq = skip_final_freq
    runner.stored_ts_hessians = []
    runner._store_ts_hessian = runner.stored_ts_hessians.append

    def loop(_threshold):
        runner._cycles_spent += 1
        if stalled:
            runner.is_stalled = True
            runner.stop_reason = "energy plateau"
            return 1, False, False
        return 1, False, True

    runner._dimer_loop = loop
    runner._calc_full_hessian_cached = lambda *, allow_reuse: torch.eye(3)
    direction = lambda *args, **kwargs: (np.ones((1, 3)), -100.0)
    monkeypatch.setattr(tsopt, "_mode_direction_by_root", direction)
    monkeypatch.setattr(tsopt, "_mode_direction_by_root_from_Hact", direction)
    monkeypatch.setattr(
        tsopt,
        "_frequencies_cm_and_modes",
        lambda *args, **kwargs: (np.array(ONE_IMAGINARY), torch.eye(3)),
    )
    monkeypatch.setattr(tsopt, "_write_all_imaginary_modes", lambda *a, **k: 1)
    monkeypatch.setattr(
        tsopt,
        "write",
        lambda path, _atoms: Path(path).write_text("final\n", encoding="utf-8"),
    )
    return runner


def test_dimer_skip_final_freq_skips_terminal_phva_after_convergence(
    monkeypatch, tmp_path, capsys,
):
    runner = _dimer_runner(tmp_path, monkeypatch, stalled=False, skip_final_freq=True)

    runner.run()

    assert runner.is_converged is True
    assert runner.hessian_status == "skipped"
    assert runner.n_imaginary_modes is None
    assert runner.imaginary_frequencies_cm == []
    assert runner.saddle_order_verified is False
    assert runner.stored_ts_hessians == []
    assert SKIP_WARNING in capsys.readouterr().err


@pytest.mark.parametrize("skip_final_freq", [True, False])
def test_dimer_plateau_stop_always_runs_terminal_phva(
    monkeypatch, tmp_path, skip_final_freq,
):
    runner = _dimer_runner(
        tmp_path, monkeypatch, stalled=True, skip_final_freq=skip_final_freq,
    )

    runner.run()

    assert runner.is_stalled is True
    assert runner.hessian_status == "completed"
    assert runner.n_imaginary_modes == 1
    assert len(runner.stored_ts_hessians) == 1


# ------------------------------------------- Hessian family (RS-P-RFO) via CLI

def _fake_ts_optimizer(*, stalled):
    class FakeTSOptimizer:
        def __init__(self, geometry, **kwargs):
            self.geometry = geometry
            self.final_fn = Path(kwargs["out_dir"]) / "final_geometry.xyz"
            self.saddle_imaginary_threshold_cm = DEFAULT_FREQUENCY_ZERO_CUTOFF_CM
            self.forces = []

        def run(self):
            self.cur_cycle = 0
            self.is_converged = not stalled
            self.is_stalled = stalled
            self.stop_reason = "energy plateau" if stalled else ""
            self.final_fn.write_text(self.geometry.as_xyz(), encoding="utf-8")

    return FakeTSOptimizer


def _run_rsprfo(monkeypatch, tmp_path, *, stalled, extra_args):
    import pdb2reaction.io.hessian_cache as hessian_cache

    phva_calls = []

    def frequencies(*args, **kwargs):
        phva_calls.append(1)
        return np.array(ONE_IMAGINARY), torch.eye(3)

    monkeypatch.setitem(tsopt.TSOPT_CLASS_MAP, "rsprfo", _fake_ts_optimizer(stalled=stalled))
    monkeypatch.setattr(tsopt, "create_calculator", lambda **_kwargs: object())
    monkeypatch.setattr(
        tsopt, "_calc_full_hessian_torch", lambda *args, **kwargs: torch.eye(3)
    )
    monkeypatch.setattr(tsopt, "_frequencies_cm_and_modes", frequencies)
    monkeypatch.setattr(tsopt, "_write_all_imaginary_modes", lambda *a, **k: 0)
    monkeypatch.setattr(tsopt, "_calc_energy", lambda *args, **kwargs: -1.0)
    monkeypatch.setattr(hessian_cache, "store", lambda *args, **kwargs: None)
    monkeypatch.setattr(
        hessian_cache, "identity_from_context", lambda *args, **kwargs: {}
    )

    out_dir = tmp_path / "out"
    result = CliRunner().invoke(root_cli, [
        "tsopt", "-i", str(_helium(tmp_path)), "-q", "0", "-m", "1",
        "-o", str(out_dir), "--opt-mode", "rsprfo", "--max-cycles", "10",
        "--no-dump", "--out-json", *extra_args,
    ])
    result_json = out_dir / "result.json"
    assert result_json.is_file(), result.output
    return json.loads(result_json.read_text(encoding="utf-8")), result, phva_calls


@pytest.mark.parametrize("spelling", [["--skip-final-freq"], ["--skip-final-freq", "true"]])
def test_hessian_family_skip_final_freq_records_skipped_hessian(
    monkeypatch, tmp_path, spelling,
):
    payload, result, phva_calls = _run_rsprfo(
        monkeypatch, tmp_path, stalled=False, extra_args=[*spelling, "--flatten"],
    )

    assert payload["optimization_status"] == "converged", result.output
    assert payload["hessian_status"] == "skipped"
    assert payload["n_imaginary_modes"] is None
    assert payload["imaginary_frequencies_cm"] == []
    assert payload["saddle_order_verified"] is False
    assert phva_calls == []
    assert SKIP_WARNING in result.output


def test_hessian_family_plateau_stop_keeps_phva_with_skip_final_freq(
    monkeypatch, tmp_path,
):
    payload, result, phva_calls = _run_rsprfo(
        monkeypatch, tmp_path, stalled=True, extra_args=["--skip-final-freq"],
    )

    assert payload["optimization_status"] == "stalled", result.output
    assert payload["hessian_status"] == "completed"
    assert payload["n_imaginary_modes"] == 1
    assert len(phva_calls) == 1
    assert SKIP_WARNING not in result.output


# ----------------------------------------------------------------------- all

def test_continuation_stops_before_irc_when_final_freq_was_skipped():
    from pdb2reaction.workflows.all import _tsopt_continuation_decision

    payload = {
        "optimization_status": "converged",
        "hessian_status": "skipped",
        "n_imaginary_modes": None,
    }
    skipped = _tsopt_continuation_decision(payload, skip_final_freq=True)
    assert skipped["continue_irc"] is False
    assert skipped["reason"] == "terminal_hessian_explicitly_skipped"
    assert skipped["skip_final_freq"] is True

    # A plateau stop is reported as such even with the flag set.
    stalled = _tsopt_continuation_decision(
        {**payload, "optimization_status": "stalled", "hessian_status": "completed",
         "n_imaginary_modes": 1},
        skip_final_freq=True,
    )
    assert stalled["reason"] == "ts_optimization_stalled"
    assert _tsopt_continuation_decision(payload)["skip_final_freq"] is False


def test_all_forwards_skip_final_freq_to_tsopt(tmp_path, monkeypatch):
    from pdb2reaction.workflows import all as all_workflow

    hei = tmp_path / "hei.xyz"
    hei.write_text("1\nHEI\nH 0.0 0.0 0.0\n", encoding="utf-8")
    captured_args = []

    def fake_run_cli_main(_name, _cli, args, **_kwargs):
        captured_args.extend(args)
        ts_dir = Path(args[args.index("--out-dir") + 1])
        ts_dir.mkdir(parents=True, exist_ok=True)
        (ts_dir / "final_geometry.xyz").write_text("1\nTS\nH 0.0 0.0 0.0\n", encoding="utf-8")
        (ts_dir / "result.json").write_text(json.dumps({
            "status": "converged",
            "optimization_status": "converged",
            "hessian_status": "skipped",
            "n_imaginary_modes": None,
            "imaginary_frequencies_cm": [],
        }), encoding="utf-8")

    monkeypatch.setattr(all_workflow, "_run_cli_main", fake_run_cli_main)
    monkeypatch.setattr(all_workflow, "_echo_detail", lambda *_a, **_k: None)

    _ts_path, ts_geom = all_workflow._run_tsopt_on_hei(
        hei,
        charge=0,
        spin=1,
        calc_cfg={"backend": "uma"},
        args_yaml=None,
        out_dir=tmp_path / "segment",
        freeze_links=False,
        opt_mode_default="hess",
        ref_pdb=None,
        convert_files=False,
        overrides={"skip_final_freq": True},
    )

    assert "--skip-final-freq" in captured_args
    assert ts_geom._tsopt_continuation["continue_irc"] is False
    assert ts_geom._tsopt_continuation["reason"] == "terminal_hessian_explicitly_skipped"
    assert ts_geom._tsopt_continuation["skip_final_freq"] is True


def test_all_records_explicit_skip_final_freq_in_tsopt_overrides(tmp_path):
    source = _helium(tmp_path, "ts.xyz")
    base = ["all", "-i", str(source), "-q", "0", "--tsopt", "--dry-run", "--show-config"]

    explicit = CliRunner().invoke(root_cli, [*base, "--skip-final-freq"])
    omitted = CliRunner().invoke(root_cli, base)

    assert explicit.exit_code == 0, explicit.output
    assert omitted.exit_code == 0, omitted.output
    assert "skip_final_freq: true" in explicit.output
    assert "skip_final_freq:" not in omitted.output


# ---------------------------------------------------------------------- freq

@pytest.mark.parametrize(
    "hess_device, expected", [("auto", "meta"), ("cpu", "cpu")],
)
def test_freq_hess_device_sets_hessian_placement_device(
    tmp_path, monkeypatch, hess_device, expected,
):
    # The calculator device is "meta" so auto and cpu differ on a CPU-only host.
    config = tmp_path / "config.yaml"
    config.write_text(yaml.safe_dump({"calc": {"device": "meta"}}), encoding="utf-8")
    seen = []

    def evaluate(geometry, calc_cfg, device, **kwargs):
        seen.append((str(calc_cfg.get("device")), device.type))
        raise RuntimeError("stop after placement")

    monkeypatch.setattr(freq_workflow, "_calc_full_hessian_torch", evaluate)
    import pdb2reaction.io.hessian_cache as hessian_cache
    monkeypatch.setattr(hessian_cache, "load_matching", lambda *args, **kwargs: None)

    result = CliRunner().invoke(root_cli, [
        "freq", "-i", str(_helium(tmp_path)), "-q", "0", "-m", "1",
        "--config", str(config), "-o", str(tmp_path / "out"),
        "--hess-device", hess_device,
    ])

    # The calculator keeps its own device; only the Hessian placement moves.
    assert seen == [("meta", expected)], result.output
    cpu_note = "Hessian placement and diagonalization will run on CPU after evaluation."
    assert (cpu_note in result.output) is (expected == "cpu")


# ----------------------------------------------------------------------- irc

def test_resolve_hessian_device_never_moves_explicit_cuda_to_cpu():
    resolve = irc_workflow._resolve_hessian_device
    assert resolve("cpu", True) == ("cpu", "explicit_cpu")
    assert resolve("cuda", True) == ("cuda", "explicit_cuda")
    assert resolve("auto", False) == ("cpu", "auto_no_cuda")
    with pytest.raises(ValueError, match="no CUDA device is available"):
        resolve("cuda", False)


@pytest.mark.parametrize("command", ["irc", "freq"])
def test_hess_device_cuda_without_cuda_is_rejected_before_dry_run(
    tmp_path, monkeypatch, command
):
    monkeypatch.setattr(torch.cuda, "is_available", lambda: False)

    result = CliRunner().invoke(root_cli, [
        command, "-i", str(_helium(tmp_path)), "-q", "0", "-m", "1",
        "-o", str(tmp_path / "out"), "--hess-device", "cuda", "--dry-run",
    ])

    assert result.exit_code == 1, result.output
    assert "--hess-device cuda was requested but no CUDA device is available" in result.output


class _FakeHessianCalculator:
    def __init__(self):
        self.hessian_calls = 0

    def get_hessian(self, atoms, coords):
        self.hessian_calls += 1
        return {"energy": -1.0, "hessian": np.eye(len(coords))}


@pytest.mark.parametrize("hess_device", ["auto", "cpu"])
def test_irc_hess_device_cpu_seeds_the_fresh_hessian_on_cpu(
    tmp_path, monkeypatch, hess_device,
):
    import pdb2reaction.io.hessian_cache as hessian_cache

    calculator = _FakeHessianCalculator()
    seeded = []

    class StopAtEulerPC:
        def __init__(self, geometry, **_kwargs):
            seeded.append(geometry._hessian)
            raise RuntimeError("stop before IRC")

    monkeypatch.setattr(irc_workflow, "create_calculator", lambda **_kw: calculator)
    monkeypatch.setattr(irc_workflow, "EulerPC", StopAtEulerPC)
    monkeypatch.setattr(hessian_cache, "load_matching", lambda *args, **kwargs: None)
    monkeypatch.setattr(hessian_cache, "identity_from_context", lambda *args, **kwargs: {})

    result = CliRunner().invoke(root_cli, [
        "irc", "-i", str(_helium(tmp_path)), "-q", "0", "-m", "1",
        "-o", str(tmp_path / "out"), "--hess-device", hess_device,
    ])

    assert len(seeded) == 1, result.output
    assert f"IRC Hessian device: requested={hess_device}" in result.output
    if hess_device == "cpu":
        assert isinstance(seeded[0], torch.Tensor)
        assert seeded[0].device.type == "cpu"
        assert calculator.hessian_calls == 1
        assert "Hessian operations will run on CPU." in result.output
    else:
        # auto leaves the fresh Hessian to EulerPC on the calculator's device.
        assert seeded[0] is None
        assert calculator.hessian_calls == 0


def test_irc_other_hessian_init_ignores_the_cached_ts_hessian(tmp_path, monkeypatch):
    import pdb2reaction.io.hessian_cache as hessian_cache

    calculator = _FakeHessianCalculator()
    seeded, cache_reads = [], []

    class StopAtEulerPC:
        def __init__(self, geometry, **kwargs):
            seeded.append((kwargs.get("hessian_init"), geometry._hessian))
            raise RuntimeError("stop before IRC")

    def load_matching(*args, **kwargs):
        cache_reads.append(args)
        return {"hessian": np.eye(3), "active_dofs": None}

    monkeypatch.setattr(irc_workflow, "create_calculator", lambda **_kw: calculator)
    monkeypatch.setattr(irc_workflow, "EulerPC", StopAtEulerPC)
    monkeypatch.setattr(hessian_cache, "load_matching", load_matching)
    monkeypatch.setattr(hessian_cache, "identity_from_context", lambda *args, **kwargs: {})
    config = tmp_path / "config.yaml"
    config.write_text(yaml.safe_dump({"irc": {"hessian_init": "unit"}}), encoding="utf-8")

    result = CliRunner().invoke(root_cli, [
        "irc", "-i", str(_helium(tmp_path)), "-q", "0", "-m", "1",
        "-o", str(tmp_path / "out"), "--hess-device", "cpu", "--config", str(config),
    ])

    assert seeded == [("unit", None)], result.output
    assert cache_reads == []
    assert calculator.hessian_calls == 0
    assert "IRC Hessian device" not in result.output
