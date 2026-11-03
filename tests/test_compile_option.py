"""``--compile/--no-compile``: UMA and ORB run uncompiled unless it is asked for.

The flag and ``calc.compile_model`` resolve as one key (both given and
different is an error), other backends reject ``--compile``, the analytical
Hessian is never compiled, and a compiled run without a compiler stops with
the install guidance before the model loads. No model is loaded here.
"""

from __future__ import annotations

import json
import sys
import types
from pathlib import Path

import click
import pytest
import yaml
from click.testing import CliRunner

from pdb2reaction.cli import cli as root_cli
from pdb2reaction.core import compile_policy
from pdb2reaction.core import prerequisites as pre
from pdb2reaction.core.compile_policy import (
    apply_compile_option,
    find_compiler,
    resolve_compile,
    runs_compiled,
)
from pdb2reaction.core.utils import calculator_provenance, calculator_run_label

HELP = (
    "torch.compile the UMA or ORB model: slower first call, slightly faster "
    "afterwards; worth it only for long runs."
)


@pytest.fixture(autouse=True)
def _fresh_state(monkeypatch):
    """Reset the once-per-process notes and run at the default verbosity with the console gate off."""
    from pdb2reaction.core import output, utils

    monkeypatch.setattr(compile_policy, "_analytical_noted", False)
    monkeypatch.setattr(compile_policy, "_uma_tf32_noted", False)
    for module, name, value in (
        (output, "_VERBOSE_LEVEL", 2),
        (utils, "_GATE_ACTIVE", False),
        (utils, "_PIPELINE_MODE", False),
        (utils, "_CHILD_MODE", False),
    ):
        monkeypatch.setattr(module, name, value)


def _which_from(available, probed=None):
    def which(name):
        if probed is not None:
            probed.append(name)
        return f"/usr/bin/{name}" if name in available else None

    return which


# ---------------------------------------------------------------------------
# find_compiler
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("device", "variable", "order"),
    [("cpu", "CXX", ["c++", "g++", "clang++"]), ("cuda", "CC", ["cc", "gcc", "clang"])],
)
def test_find_compiler_probes_the_variable_then_the_names(monkeypatch, device, variable, order):
    monkeypatch.delenv("CC", raising=False)
    monkeypatch.delenv("CXX", raising=False)
    probed = []
    monkeypatch.setattr(compile_policy.shutil, "which", _which_from(set(), probed))
    assert find_compiler(device) is None
    assert probed == order

    monkeypatch.setenv(variable, "my-compiler")
    probed.clear()
    monkeypatch.setattr(compile_policy.shutil, "which", _which_from({"my-compiler", order[-1]}, probed))
    assert find_compiler(device) == "/usr/bin/my-compiler"
    assert probed == ["my-compiler"]


@pytest.mark.parametrize(
    ("device", "found"), [("cpu", "clang++"), ("cuda:0", "clang")],
)
def test_find_compiler_takes_the_first_name_on_path(monkeypatch, device, found):
    monkeypatch.delenv("CC", raising=False)
    monkeypatch.delenv("CXX", raising=False)
    monkeypatch.setattr(compile_policy.shutil, "which", _which_from({found}))
    assert find_compiler(device) == f"/usr/bin/{found}"


# ---------------------------------------------------------------------------
# resolve_compile
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("requested", [None, False])
def test_unrequested_compile_is_off_for_any_backend(monkeypatch, requested):
    monkeypatch.setattr(compile_policy.shutil, "which", _which_from(set()))
    for backend in ("uma", "orb", "mace", "dft"):
        assert resolve_compile(requested, backend=backend, device="cpu", analytical_hessian=False) is False


@pytest.mark.parametrize("backend", ["mace", "aimnet2", "dft", "custom"])
def test_compile_for_another_backend_is_a_usage_error(backend):
    with pytest.raises(click.UsageError, match="^--compile applies to the uma and orb backends$"):
        resolve_compile(True, backend=backend, device="cpu", analytical_hessian=False)


@pytest.mark.parametrize("backend", ["uma", "orb"])
def test_analytical_hessian_is_not_compiled_and_says_so_once(monkeypatch, capsys, backend):
    monkeypatch.setattr(compile_policy.shutil, "which", _which_from(set()))
    for _ in range(2):
        assert resolve_compile(True, backend=backend, device="cpu", analytical_hessian=True) is False
    err = capsys.readouterr().err
    assert err.count("[compile] --compile is ignored for the analytical Hessian.") == 1


def test_uma_tf32_on_one_worker_is_not_compiled_and_says_so_once(monkeypatch, capsys):
    monkeypatch.setattr(compile_policy.shutil, "which", _which_from(set()))
    for _ in range(2):
        assert resolve_compile(
            True, backend="uma", device="cpu", analytical_hessian=False, uma_tf32_single_worker=True,
        ) is False
    err = capsys.readouterr().err
    assert err.count(
        "[compile] --compile is ignored with --precision tf32 on UMA (UMA Nitro does not use torch.compile)."
    ) == 1
    assert runs_compiled({"backend": "uma", "precision": "tf32", "compile_model": True}) is False
    assert runs_compiled({"backend": "uma", "precision": "tf32", "compile_model": True, "workers": 2}) is True


@pytest.mark.parametrize(("device", "name"), [("cpu", "g++"), ("cuda", "gcc")])
def test_compile_with_a_compiler_is_on(monkeypatch, device, name):
    monkeypatch.delenv("CC", raising=False)
    monkeypatch.delenv("CXX", raising=False)
    monkeypatch.setattr(compile_policy.shutil, "which", _which_from({name}))
    assert resolve_compile(True, backend="orb", device=device, analytical_hessian=False) is True


@pytest.mark.parametrize(
    ("device", "kind"), [("cpu", pre.KIND_CXX_COMPILER), ("cuda", pre.KIND_C_COMPILER)],
)
def test_compile_without_a_compiler_stops_with_the_guidance(monkeypatch, capsys, device, kind):
    monkeypatch.delenv("CC", raising=False)
    monkeypatch.delenv("CXX", raising=False)
    monkeypatch.setattr(compile_policy.shutil, "which", _which_from(set()))
    with pytest.raises(pre.MissingPrerequisiteError) as missing:
        resolve_compile(True, backend="uma", device=device, analytical_hessian=False)
    assert missing.value.kind == kind
    with pytest.raises(SystemExit) as stop:
        pre.exit_if_missing_prerequisite(missing.value)
    assert stop.value.code == 1
    assert capsys.readouterr().err.strip() == pre.guidance(kind)


def test_auto_device_follows_cuda_availability(monkeypatch):
    import torch

    monkeypatch.delenv("CC", raising=False)
    monkeypatch.delenv("CXX", raising=False)
    monkeypatch.setattr(compile_policy.shutil, "which", _which_from(set()))
    for cuda, kind in ((True, pre.KIND_C_COMPILER), (False, pre.KIND_CXX_COMPILER)):
        monkeypatch.setattr(torch.cuda, "is_available", lambda cuda=cuda: cuda)
        with pytest.raises(pre.MissingPrerequisiteError) as missing:
            resolve_compile(True, backend="uma", device="auto", analytical_hessian=False)
        assert missing.value.kind == kind


# ---------------------------------------------------------------------------
# apply_compile_option: CLI and config
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("cli", "config", "expected"),
    [
        (None, None, False),
        (True, None, True),
        (False, None, False),
        (None, True, True),
        (None, False, False),
        (True, True, True),
        (False, False, False),
    ],
)
def test_one_explicit_value_wins_and_unset_is_off(cli, config, expected):
    calc_cfg = {"backend": "uma"}
    if config is not None:
        calc_cfg["compile_model"] = config
    apply_compile_option(calc_cfg, cli)
    assert calc_cfg["compile_model"] is expected


@pytest.mark.parametrize(("cli", "config"), [(True, False), (False, True)])
def test_cli_and_config_that_differ_are_a_usage_error(cli, config):
    flag = "--compile" if cli else "--no-compile"
    with pytest.raises(click.UsageError, match=f"^{flag} and calc.compile_model: {str(config).lower()} conflict"):
        apply_compile_option({"backend": "orb", "compile_model": config}, cli)


@pytest.mark.parametrize("backend", ["mace", "aimnet2", "dft", "custom"])
def test_other_backends_reject_compile_and_accept_no_compile(backend):
    with pytest.raises(click.UsageError, match="^--compile applies to the uma and orb backends$"):
        apply_compile_option({"backend": backend}, True)
    with pytest.raises(click.UsageError, match="^--compile applies to the uma and orb backends$"):
        apply_compile_option({"backend": backend, "compile_model": True}, None)
    calc_cfg = {"backend": backend}
    apply_compile_option(calc_cfg, False)
    assert calc_cfg["compile_model"] is False


def test_config_value_must_be_a_bool():
    with pytest.raises(click.BadParameter, match="calc.compile_model must be true or false"):
        apply_compile_option({"backend": "uma", "compile_model": "yes"}, None)


# ---------------------------------------------------------------------------
# The flag on every command that has --precision
# ---------------------------------------------------------------------------


def _commands_with_precision():
    from pdb2reaction.cli.app import _LAZY_SUBCOMMANDS
    import importlib

    commands = {}
    for name, (module_name, attr, _help) in _LAZY_SUBCOMMANDS.items():
        command = getattr(importlib.import_module(module_name), attr)
        names = {param.name for param in command.params}
        if "precision" in names:
            commands[name] = command
    return commands


def test_every_command_with_precision_has_compile():
    commands = _commands_with_precision()
    assert {"opt", "tsopt", "freq", "irc", "sp", "scan", "scan2d", "scan3d",
            "path-opt", "path-search", "all", "trj2fig"} <= set(commands)
    for name, command in commands.items():
        option = next((p for p in command.params if p.name == "compile_model"), None)
        assert option is not None, name
        assert option.opts == ["--compile"] and option.secondary_opts == ["--no-compile"], name
        assert option.default is None and option.show_default == "off", name
        assert option.help == HELP, name


# ---------------------------------------------------------------------------
# What ran: the backend line and mlip_compiled
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("calc_cfg", "label", "compiled"),
    [
        ({"backend": "uma", "precision": "fp64", "compile_model": True},
         "UMA (UMA-S-1.2.1 (OMol), fp64, compiled)", True),
        ({"backend": "uma", "precision": "fp64", "compile_model": False},
         "UMA (UMA-S-1.2.1 (OMol), fp64)", False),
        ({"backend": "uma", "precision": "fp64"}, "UMA (UMA-S-1.2.1 (OMol), fp64)", False),
        ({"backend": "uma", "precision": "fp32", "compile_model": True,
          "hessian_calc_mode": "Analytical"}, "UMA (UMA-S-1.2.1 (OMol), fp32)", False),
        ({"backend": "orb", "precision": "float64", "compile_model": True},
         "ORB (OrbMol-v2, fp64, compiled)", True),
        ({"backend": "mace", "default_dtype": "float64", "compile_model": False},
         "MACE (MACE-OMOL-0, fp64)", False),
    ],
)
def test_backend_line_and_provenance_record_compilation(calc_cfg, label, compiled):
    assert calculator_run_label(calc_cfg) == label
    assert calculator_provenance(calc_cfg)["mlip_compiled"] is compiled
    assert runs_compiled(calc_cfg) is compiled


def test_custom_calculator_is_recorded_uncompiled(tmp_path):
    from pdb2reaction.backends import apply_calc_file_to_calc_cfg

    calc_file = tmp_path / "toy.py"
    calc_file.write_text("def get_calculator(**kwargs):\n    return None\n")
    calc_cfg = {"backend": "uma"}
    apply_calc_file_to_calc_cfg(calc_cfg, str(calc_file), None)
    assert calculator_provenance(calc_cfg)["mlip_compiled"] is False


def test_opt_dry_run_prints_the_compiled_backend_line(tmp_path):
    xyz = tmp_path / "h2.xyz"
    xyz.write_text("2\nH2\nH 0 0 0\nH 0 0 0.74\n", encoding="utf-8")
    base = ["opt", "-i", str(xyz), "-q", "0", "-m", "1", "--out-dir", str(tmp_path / "out"),
            "--precision", "fp64", "--dry-run"]

    compiled = CliRunner().invoke(root_cli, [*base, "--compile"])
    assert compiled.exit_code == 0, compiled.output
    assert "[backend] UMA (UMA-S-1.2.1 (OMol), fp64, compiled)" in compiled.output

    default = CliRunner().invoke(root_cli, base)
    assert default.exit_code == 0, default.output
    assert "[backend] UMA (UMA-S-1.2.1 (OMol), fp64)" in default.output
    assert "compiled" not in default.output


@pytest.mark.parametrize(
    ("args", "config", "message"),
    [
        (["--compile", "-b", "mace"], None, "--compile applies to the uma and orb backends"),
        (["--compile"], "calc:\n  compile_model: false\n", "--compile and calc.compile_model: false conflict."),
        (["--no-compile"], "calc:\n  compile_model: true\n", "--no-compile and calc.compile_model: true conflict."),
    ],
)
def test_opt_rejects_compile_misuse(tmp_path, args, config, message):
    xyz = tmp_path / "h2.xyz"
    xyz.write_text("2\nH2\nH 0 0 0\nH 0 0 0.74\n", encoding="utf-8")
    cmd = ["opt", "-i", str(xyz), "-q", "0", "-m", "1", "--out-dir", str(tmp_path / "out"),
           "--dry-run", *args]
    if config is not None:
        path = tmp_path / "config.yaml"
        path.write_text(config, encoding="utf-8")
        cmd += ["--config", str(path)]
    result = CliRunner().invoke(root_cli, cmd)
    assert result.exit_code == 2, result.output
    assert message in result.output


def test_no_compile_is_accepted_by_other_backends(tmp_path):
    xyz = tmp_path / "h2.xyz"
    xyz.write_text("2\nH2\nH 0 0 0\nH 0 0 0.74\n", encoding="utf-8")
    result = CliRunner().invoke(root_cli, [
        "opt", "-i", str(xyz), "-q", "0", "-m", "1", "--out-dir", str(tmp_path / "out"),
        "-b", "mace", "--no-compile", "--dry-run",
    ])
    assert result.exit_code == 0, result.output


def test_sp_with_compile_and_no_compiler_stops_before_the_model_loads(monkeypatch, tmp_path):
    # The real FAIR-Chem, not the import stub test_uma_workers_analytical.py leaves without it.
    pytest.importorskip("fairchem.core.units.mlip_unit.api.inference")
    import pdb2reaction.backends.uma as uma_module

    def no_load(*args, **kwargs):
        pytest.fail("the model must not load when the compiler is missing")

    monkeypatch.setattr(compile_policy.shutil, "which", _which_from(set()))
    monkeypatch.setattr(compile_policy, "_on_cuda", lambda device: False)
    monkeypatch.setattr(uma_module.pretrained_mlip, "get_predict_unit", no_load, raising=False)
    inp = tmp_path / "geom.xyz"
    inp.write_text("1\nframe\nC 0.0 0.0 0.0\n", encoding="utf-8")
    result = CliRunner().invoke(
        root_cli,
        ["sp", "-i", str(inp), "-q", "0", "-m", "1", "--compile", "-o", str(tmp_path / "out")],
    )
    assert result.exit_code == 1
    assert "Missing system prerequisite: C++ compiler (g++)." in result.output
    assert "Traceback" not in result.output and "Unhandled" not in result.output
    record = json.loads((tmp_path / "out" / "result.json").read_text())
    assert record["execution_status"] == "failed"


# ---------------------------------------------------------------------------
# Backends
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("compile_model", "mode", "expected"),
    [(False, "FiniteDifference", False), (True, "FiniteDifference", True), (True, "Analytical", False)],
)
def test_uma_calculator_passes_the_resolved_choice_to_its_core(monkeypatch, compile_model, mode, expected):
    # The real FAIR-Chem, not the import stub test_uma_workers_analytical.py leaves without it.
    pytest.importorskip("fairchem.core.units.mlip_unit.api.inference")
    from pdb2reaction.backends.uma import UMACalculator

    monkeypatch.setattr(compile_policy.shutil, "which", _which_from({"g++"}))
    monkeypatch.delenv("CXX", raising=False)
    calc = UMACalculator(device="cpu", compile_model=compile_model, hessian_calc_mode=mode)
    assert calc._core_kw["compile_model"] is expected


@pytest.mark.parametrize("compile_model", [False, True])
def test_uma_ase_calculator_builds_the_fast_path(monkeypatch, compile_model):
    # The real FAIR-Chem, not the import stub test_uma_workers_analytical.py leaves without it.
    pytest.importorskip("fairchem.core.units.mlip_unit.api.inference")
    import pdb2reaction.backends.uma as uma_module

    class Reached(Exception):
        pass

    captured = {}

    def fake_get_predict_unit(_model, **kwargs):
        captured.update(kwargs)
        raise Reached

    monkeypatch.setattr(compile_policy.shutil, "which", _which_from({"g++"}))
    monkeypatch.delenv("CXX", raising=False)
    monkeypatch.setattr(uma_module.pretrained_mlip, "get_predict_unit", fake_get_predict_unit, raising=False)
    with pytest.raises(Reached):
        uma_module.UMAASECalculator(device="cpu", precision="fp64", compile_model=compile_model)
    settings = captured["inference_settings"]
    assert settings.compile is compile_model
    assert settings.merge_mole is True and settings.activation_checkpointing is False


class _FakeOrbModels:
    """Just enough of orb_models for OrbCalculator.__init__ to load a model."""

    def __init__(self, monkeypatch):
        self.calls = []

        class ConservativeForcefieldRegressor:
            pass

        class ORBCalculator:
            def __init__(self, *args, **kwargs):
                pass

        def loader(**kwargs):
            self.calls.append(kwargs)
            return ConservativeForcefieldRegressor(), None

        pretrained = types.ModuleType("orb_models.forcefield.pretrained")
        pretrained.ORB_PRETRAINED_MODELS = {"orbmol-v2": loader}
        calculator = types.ModuleType("orb_models.forcefield.inference.calculator")
        calculator.ORBCalculator = ORBCalculator
        modules = {
            "orb_models": types.ModuleType("orb_models"),
            "orb_models.forcefield": types.ModuleType("orb_models.forcefield"),
            "orb_models.forcefield.pretrained": pretrained,
            "orb_models.forcefield.inference": types.ModuleType("orb_models.forcefield.inference"),
            "orb_models.forcefield.inference.calculator": calculator,
        }
        modules["orb_models"].forcefield = modules["orb_models.forcefield"]
        modules["orb_models.forcefield"].pretrained = pretrained
        modules["orb_models.forcefield"].inference = modules["orb_models.forcefield.inference"]
        modules["orb_models.forcefield.inference"].calculator = calculator
        for name, module in modules.items():
            monkeypatch.setitem(sys.modules, name, module)


@pytest.mark.parametrize(
    ("kwargs", "expected"),
    [
        ({}, False),
        ({"compile_model": False}, False),
        ({"compile_model": True}, True),
        ({"compile_model": True, "hessian_calc_mode": "Analytical"}, False),
    ],
)
def test_orb_receives_an_explicit_bool(monkeypatch, kwargs, expected):
    from pdb2reaction.backends.orb import OrbCalculator

    fake = _FakeOrbModels(monkeypatch)
    monkeypatch.setattr(compile_policy.shutil, "which", _which_from({"g++"}))
    monkeypatch.delenv("CXX", raising=False)
    calc = OrbCalculator(device="cpu", **kwargs)
    assert calc.compile_model is expected
    assert fake.calls[0]["compile"] is expected


# ---------------------------------------------------------------------------
# all
# ---------------------------------------------------------------------------


def test_all_child_config_carries_only_an_explicit_value():
    from pdb2reaction.workflows._run_session import RunSession
    from pdb2reaction.workflows.all import _write_args_yaml_with_freeze_atoms

    session = RunSession()
    assert _write_args_yaml_with_freeze_atoms(None, [], compile_model=None, session=session) is None
    for value in (True, False):
        path = _write_args_yaml_with_freeze_atoms(None, [], compile_model=value, session=session)
        payload = yaml.safe_load(Path(path).read_text(encoding="utf-8"))
        assert payload["calc"]["compile_model"] is value
    session.close()


@pytest.mark.parametrize("flag", [None, "--compile", "--no-compile"])
def test_all_forwards_an_explicit_value_to_its_children(tmp_path, monkeypatch, flag):
    from pdb2reaction.core.result_commit import RUN_ID_ENV
    from pdb2reaction.workflows import all as all_workflow
    from test_all_tsopt_bond_diagnostic import _pdb, _write

    class ReachedChild(BaseException):
        pass

    captured = []

    def first_child(name, _cli, args, **kwargs):
        args = list(args)
        config = {}
        if "--config" in args:
            config = yaml.safe_load(Path(args[args.index("--config") + 1]).read_text(encoding="utf-8"))
        captured.append((name, args, (config.get("calc") or {}).get("compile_model")))
        raise ReachedChild

    monkeypatch.setattr(all_workflow, "_run_cli_main", first_child)
    monkeypatch.delenv(RUN_ID_ENV, raising=False)
    inputs = [_write(tmp_path / f"input{i}.pdb", _pdb(i)) for i in (0, 2)]
    args = ["all", *[arg for path in inputs for arg in ("-i", str(path))]]
    args += ["-q", "0", "-m", "1", "--out-dir", str(tmp_path / "out"),
             "--no-preopt", "--no-freeze-links"]
    if flag is not None:
        args.append(flag)
    with pytest.raises(ReachedChild):
        CliRunner().invoke(root_cli, args, catch_exceptions=False)

    (name, child_args, forwarded), = captured
    assert name == "path-opt"
    assert "--compile" not in child_args and "--no-compile" not in child_args
    assert forwarded is {None: None, "--compile": True, "--no-compile": False}[flag]


def test_all_rejects_compile_with_another_backend_in_dry_run(tmp_path):
    xyz = tmp_path / "ts.xyz"
    xyz.write_text("2\nH2\nH 0 0 0\nH 0 0 0.74\n", encoding="utf-8")
    result = CliRunner().invoke(
        root_cli, ["all", "-i", str(xyz), "-q", "0", "--tsopt", "-b", "mace", "--compile", "--dry-run"]
    )
    assert result.exit_code == 2, result.output
    assert "--compile applies to the uma and orb backends" in result.output


def test_all_show_config_prints_mlip_compiled(tmp_path):
    xyz = tmp_path / "h2.xyz"
    xyz.write_text("2\nH2\nH 0 0 0\nH 0 0 0.74\n", encoding="utf-8")
    result = CliRunner().invoke(root_cli, [
        "all", "-i", str(xyz), "-i", str(xyz), "-q", "0", "--compile",
        "--show-config", "--dry-run", "--out-dir", str(tmp_path / "out"),
    ])
    assert result.exit_code == 0, result.output
    assert "mlip_compiled: true" in result.output


def test_all_summary_copies_mlip_compiled(tmp_path):
    from pdb2reaction.workflows import all as all_workflow

    for compiled in (True, False):
        summary = {"segments": []}
        provenance = calculator_provenance(
            {"backend": "uma", "precision": "fp32", "compile_model": compiled}
        )
        all_workflow._enrich_summary(
            summary, pipeline_mode="path-opt", out_dir=tmp_path,
            mlip_backend="uma", mlip_model="uma-s-1p2p1", mlip_precision="fp32",
            charge=0, spin=1, config={"tsopt": False, "thermo": False, "dft": False},
            calculator_provenance_data=provenance,
        )
        assert summary["mlip_compiled"] is compiled
