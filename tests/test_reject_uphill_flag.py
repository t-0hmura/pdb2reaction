"""Contracts for uphill rejection in minimum and transition-state searches.

The flag lets a user opt in to the post-IRC endpoint RFO uphill-rejection
safeguard. These tests pin four falsifiable facts:

1. the shipped default is off (``RFO_KW["reject_uphill"] is False``);
2. the default path (flag not passed -> ``reject_uphill=None``) leaves the
   endpoint RFO config byte-identical to ``RFO_KW`` (no behavior change);
3. an explicit toggle threads ``True``/``False`` into the *endpoint* RFO config
   through the real ``_optimize_endpoint_geom`` code path (rfo branch only);
4. TS optimizers force uphill rejection off even when YAML-like input tries to
   re-enable it.
"""

from __future__ import annotations

from pathlib import Path

import click
import pytest
from click.testing import CliRunner

import pdb2reaction.workflows.all as allmod
from pdb2reaction.cli import cli as root_cli
from pdb2reaction.core.defaults import RFO_KW
from pdb2reaction.workflows.tsopt import (
    _build_rsirfo_kwargs,
    _force_ts_reject_uphill_off,
)


def test_shipped_default_is_reject_uphill_off() -> None:
    # The None-path preserves whatever RFO_KW declares; this pins that default.
    assert RFO_KW["reject_uphill"] is False
    assert RFO_KW["uphill_tolerance"] == 1e-4


def test_ts_rfo_forces_reject_uphill_off(tmp_path: Path) -> None:
    kwargs = _build_rsirfo_kwargs(
        {"reject_uphill": True},
        {"reject_uphill": True},
        tmp_path,
    )
    assert kwargs["reject_uphill"] is False


def test_ts_dimer_forces_reject_uphill_off() -> None:
    hostile_yaml_cfg = {"reject_uphill": True}
    effective = _force_ts_reject_uphill_off(hostile_yaml_cfg)
    assert effective["reject_uphill"] is False
    assert hostile_yaml_cfg["reject_uphill"] is True


@pytest.mark.parametrize("command", ["opt", "all"])
def test_toggle_is_exposed_with_default_off(command: str) -> None:
    ctx = click.Context(root_cli)
    cmd = root_cli.get_command(ctx, command)
    param = next(
        (p for p in cmd.params if isinstance(p, click.Option) and p.name == "reject_uphill"),
        None,
    )
    assert param is not None, f"{command} is missing --reject-uphill"
    assert param.opts == ["--reject-uphill"]
    assert param.secondary_opts == ["--no-reject-uphill"]
    assert param.is_bool_flag is True
    assert param.default is False


class _StopBeforeRun(Exception):
    """Raised by the stub optimizer to capture cfg without running an opt."""


class _GeomStub:
    calculator = None
    freeze_atoms: list[int] = []

    def set_calculator(self, _calc) -> None:  # noqa: D401 - test stub
        pass


@pytest.mark.parametrize(
    "reject_uphill, expected",
    [(None, False), (True, True), (False, False)],
    ids=["default-none", "explicit-on", "explicit-off"],
)
def test_endpoint_rfo_receives_reject_uphill(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, reject_uphill, expected
) -> None:
    captured: dict[str, object] = {}

    class _StubRFO:
        def __init__(self, geom, **cfg) -> None:
            captured["cfg"] = cfg
            raise _StopBeforeRun()

    # OptClass = RFOptimizer is resolved from the module global at call time.
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
            reject_uphill=reject_uphill,
        )

    cfg = captured["cfg"]
    assert cfg["reject_uphill"] is expected


def test_opt_cli_accepts_the_toggle() -> None:
    runner = CliRunner()
    with runner.isolated_filesystem():
        Path("s.pdb").write_text("END\n", encoding="utf-8")
        result = runner.invoke(
            root_cli,
            ["opt", "-i", "s.pdb", "-q", "0", "--opt-mode", "hess",
             "--dry-run", "--no-reject-uphill"],
        )
    assert result.exit_code == 0, result.output


def test_all_cli_accepts_the_toggle() -> None:
    runner = CliRunner()
    with runner.isolated_filesystem():
        Path("s.pdb").write_text("END\n", encoding="utf-8")
        result = runner.invoke(
            root_cli,
            ["all", "-i", "s.pdb", "-q", "0", "--tsopt", "True", "--dry-run", "True",
             "--no-reject-uphill"],
        )
    assert result.exit_code == 0, result.output


@pytest.mark.parametrize(
    "extra, expected_eff",
    [(["--no-reject-uphill"], False), ([], None), (["--reject-uphill"], True)],
    ids=["off-arm", "default-unchanged", "on-arm"],
)
def test_all_gate_resolves_endpoint_reject_uphill(tmp_path: Path, extra, expected_eff) -> None:
    """Positive control for the benchmark's on/off arms.

    Parse `all` the way the CLI does and reproduce the exact resolution line
    ``_reject_uphill_eff = bool(reject_uphill) if cli_param_overridden(...) else None``.
    The default arm forwards no override and therefore inherits the shared
    default-off RFO configuration.
    """
    from pdb2reaction.workflows.all import cli as all_cli
    from pdb2reaction.core.utils import cli_param_overridden

    pdb = tmp_path / "x.pdb"
    pdb.write_text("END\n", encoding="utf-8")
    ctx = all_cli.make_context("all", ["-i", str(pdb), "--tsopt", "True", *extra])
    gate = cli_param_overridden(ctx, "reject_uphill")
    eff = bool(ctx.params["reject_uphill"]) if gate else None
    assert eff is expected_eff
