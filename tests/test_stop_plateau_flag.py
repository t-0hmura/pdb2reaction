"""Contracts for the opt-in energy-plateau stop (``--stop-plateau``).

A plateau stop reports ``stalled``; it never reports convergence and it never
replaces ``--max-cycles`` as the real bound.  Because a run that stops there
leaves the geometry unconverged, it is opt-in.  These tests pin:

1. the shipped default is off in every optimizer config block AND in the
   bundled ``Optimizer`` signature, so a construction site that passes explicit
   kwargs instead of a config block cannot silently enable it;
2. ``opt`` / ``tsopt`` / ``all`` all expose the toggle and its two value
   options, with the toggle defaulting to off;
3. every flag ``all`` forwards to the ``tsopt`` child exists on that child;
4. the ``all`` gate resolves to ``None`` unless the toggle is explicit, and an
   explicit toggle reaches the in-process endpoint RFO;
5. a shared ``opt`` setting reaches the RS-I-RFO/RS-P-RFO macro kwargs.
"""

from __future__ import annotations

import inspect
from pathlib import Path

import click
import pytest
from click.testing import CliRunner

import pdb2reaction.core.defaults as defaults_mod
import pdb2reaction.workflows.all as allmod
from pdb2reaction.cli import cli as root_cli
from pdb2reaction.workflows.tsopt import _build_rsirfo_kwargs


def test_shipped_default_is_plateau_off() -> None:
    """Sweep every config block, not just ``OPT_BASE_KW``: a derived block that
    re-declared the key would reintroduce the default-on behavior for exactly
    one optimizer."""
    seen = 0
    for name, value in vars(defaults_mod).items():
        if not isinstance(value, dict):
            continue
        for cfg_name, cfg in [(name, value)] + [
            (f"{name}[{k!r}]", v) for k, v in value.items() if isinstance(v, dict)
        ]:
            if "energy_plateau" not in cfg:
                continue
            seen += 1
            assert cfg["energy_plateau"] is False, cfg_name
    assert seen >= 2, "no optimizer config block declares energy_plateau"
    assert defaults_mod.OPT_BASE_KW["energy_plateau"] is False
    assert defaults_mod.OPT_BASE_KW["energy_plateau_thresh"] == 1e-4
    assert defaults_mod.OPT_BASE_KW["energy_plateau_window"] == 50


def test_bundled_optimizer_default_is_plateau_off() -> None:
    """``align_freeze`` builds LBFGS from explicit kwargs with no config block,
    so the library default is the effective default there."""
    from pysisyphus.optimizers.Optimizer import Optimizer

    default = inspect.signature(Optimizer.__init__).parameters["energy_plateau"].default
    assert default is False


@pytest.mark.parametrize("command", ["opt", "tsopt", "all"])
def test_toggle_is_exposed_with_default_off(command: str) -> None:
    ctx = click.Context(root_cli)
    cmd = root_cli.get_command(ctx, command)
    params = {p.name: p for p in cmd.params if isinstance(p, click.Option)}

    toggle = params.get("stop_plateau")
    assert toggle is not None, f"{command} is missing --stop-plateau"
    assert toggle.opts == ["--stop-plateau"]
    assert toggle.secondary_opts == ["--no-stop-plateau"]
    assert toggle.is_bool_flag is True
    assert toggle.default is False

    for name, type_name in (
        ("stop_plateau_thresh", "float"),
        ("stop_plateau_window", "integer"),
    ):
        option = params.get(name)
        assert option is not None, f"{command} is missing --{name.replace('_', '-')}"
        assert option.default is None
        assert option.type.name == type_name


def test_forwarded_flags_exist_on_the_tsopt_child() -> None:
    """``all`` forwards these by literal string; a typo would be silent."""
    ctx = click.Context(root_cli)
    cmd = root_cli.get_command(ctx, "tsopt")
    flags = {
        opt
        for param in cmd.params
        for opt in list(param.opts) + list(param.secondary_opts)
    }
    assert {
        "--stop-plateau",
        "--no-stop-plateau",
        "--stop-plateau-thresh",
        "--stop-plateau-window",
    } <= flags


@pytest.mark.parametrize(
    "extra, expected_eff",
    [([], None), (["--stop-plateau"], True), (["--no-stop-plateau"], False)],
    ids=["default-unchanged", "on-arm", "off-arm"],
)
def test_all_gate_resolves_stop_plateau(tmp_path: Path, extra, expected_eff) -> None:
    """Reproduce the resolution line ``_stop_plateau_eff = bool(stop_plateau)
    if cli_param_overridden(...) else None``: the default arm forwards no
    override at all, so every optimizer keeps its own default-off config."""
    from pdb2reaction.workflows.all import cli as all_cli
    from pdb2reaction.core.utils import cli_param_overridden

    pdb = tmp_path / "x.pdb"
    pdb.write_text("END\n", encoding="utf-8")
    ctx = all_cli.make_context("all", ["-i", str(pdb), "--tsopt", "True", *extra])
    gate = cli_param_overridden(ctx, "stop_plateau")
    eff = bool(ctx.params["stop_plateau"]) if gate else None
    assert eff is expected_eff


class _StopBeforeRun(Exception):
    """Raised by the stub optimizer to capture cfg without running an opt."""


class _GeomStub:
    calculator = None
    freeze_atoms: list[int] = []

    def set_calculator(self, _calc) -> None:  # noqa: D401 - test stub
        pass


@pytest.mark.parametrize(
    "stop_plateau, expected",
    [(None, False), (True, True), (False, False)],
    ids=["default-none", "explicit-on", "explicit-off"],
)
def test_endpoint_rfo_receives_stop_plateau(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, stop_plateau, expected
) -> None:
    captured: dict[str, object] = {}

    class _StubRFO:
        def __init__(self, geom, **cfg) -> None:
            captured["cfg"] = cfg
            raise _StopBeforeRun()

    monkeypatch.setattr(allmod, "RFOptimizer", _StubRFO)

    with pytest.raises(_StopBeforeRun):
        allmod._optimize_endpoint_geom(
            _GeomStub(),
            "hess",  # -> rfo branch
            tmp_path,
            "reactant",
            dump=False,
            thresh=None,
            calc_identity_cfg=None,  # no cached Hessian seed -> straight to opt
            stop_plateau=stop_plateau,
            stop_plateau_window=None if stop_plateau is None else 7,
        )

    cfg = captured["cfg"]
    assert cfg["energy_plateau"] is expected
    assert cfg["energy_plateau_window"] == (50 if stop_plateau is None else 7)


def test_shared_opt_setting_reaches_the_ts_macro_kwargs(tmp_path: Path) -> None:
    """The RS-I-RFO/RS-P-RFO kwargs are assembled from the shared ``opt`` block
    plus ``rsirfo.*``; an opt-in plateau must survive that merge."""
    from pdb2reaction.core.defaults import OPT_BASE_KW, RSIRFO_KW

    opt_cfg = {**OPT_BASE_KW, "energy_plateau": True, "energy_plateau_window": 7}
    kwargs = _build_rsirfo_kwargs(opt_cfg, dict(RSIRFO_KW), tmp_path)

    assert kwargs["energy_plateau"] is True
    assert kwargs["energy_plateau_window"] == 7


@pytest.mark.parametrize("command", ["opt", "tsopt"])
def test_cli_accepts_the_toggle(command: str) -> None:
    runner = CliRunner()
    with runner.isolated_filesystem():
        Path("s.pdb").write_text("END\n", encoding="utf-8")
        result = runner.invoke(
            root_cli,
            [command, "-i", "s.pdb", "-q", "0", "--dry-run",
             "--stop-plateau", "--stop-plateau-window", "7"],
        )
    assert result.exit_code == 0, result.output


def test_all_cli_accepts_the_toggle() -> None:
    runner = CliRunner()
    with runner.isolated_filesystem():
        Path("s.pdb").write_text("END\n", encoding="utf-8")
        result = runner.invoke(
            root_cli,
            ["all", "-i", "s.pdb", "-q", "0", "--tsopt", "True", "--dry-run", "True",
             "--no-stop-plateau"],
        )
    assert result.exit_code == 0, result.output
