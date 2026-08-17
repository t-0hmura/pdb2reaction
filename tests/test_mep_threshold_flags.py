"""MEP threshold surfaces: GSM preset, DMF tolerance, and parent forwarding."""

from __future__ import annotations

import ast
from pathlib import Path
from typing import Set

import click
import pytest
from click.testing import CliRunner

from pdb2reaction.cli import cli as root_cli
from pdb2reaction.workflows import all as all_workflow
from pdb2reaction.workflows.path_opt import resolve_dmf_solve_tol

MEP_FLAGS = ("--thresh-gsm", "--thresh-dmf")


def _declared_flags(cli: click.Command) -> Set[str]:
    flags: Set[str] = set()
    for param in cli.params:
        flags.update(param.opts)
    return flags


def _flags_forwarded_to(var_name: str) -> Set[str]:
    """Collect literal long options appended to ``var_name`` inside ``all.py``."""
    tree = ast.parse(Path(all_workflow.__file__).read_text(encoding="utf-8"))
    flags: Set[str] = set()

    def _literals(node: ast.AST) -> Set[str]:
        return {
            child.value.strip()
            for child in ast.walk(node)
            if isinstance(child, ast.Constant)
            and isinstance(child.value, str)
            and child.value.strip().startswith("--")
        }

    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        func = node.func
        if (
            isinstance(func, ast.Attribute)
            and func.attr in {"append", "extend"}
            and isinstance(func.value, ast.Name)
            and func.value.id == var_name
        ):
            for arg in node.args:
                flags.update(_literals(arg))
        if (
            isinstance(func, ast.Name)
            and func.id in {"_append_cli_arg", "_append_toggle_arg"}
            and len(node.args) >= 2
            and isinstance(node.args[0], ast.Name)
            and node.args[0].id == var_name
        ):
            flags.update(_literals(node.args[1]))
    return flags


def test_parent_and_mep_children_declare_the_same_flags() -> None:
    from pdb2reaction.workflows.path_opt import cli as path_opt_cli
    from pdb2reaction.workflows.path_search import cli as path_search_cli

    for flag in MEP_FLAGS:
        assert flag in _declared_flags(all_workflow.cli)
        assert flag in _declared_flags(path_opt_cli)
        assert flag in _declared_flags(path_search_cli)


@pytest.mark.parametrize("var_name", ["po_args", "ps_args"])
def test_all_forwards_the_mep_flags_to_its_path_children(var_name: str) -> None:
    forwarded = _flags_forwarded_to(var_name)
    for flag in MEP_FLAGS:
        assert flag in forwarded


def test_default_keeps_the_historical_tight_tolerance() -> None:
    assert resolve_dmf_solve_tol({}) == "tight"


@pytest.mark.parametrize("raw", ["tight", " Middle ", "LOOSE"])
def test_presets_are_normalized(raw: str) -> None:
    assert resolve_dmf_solve_tol({"tol": raw}) == raw.strip().lower()


def test_float_tolerance_is_passed_through() -> None:
    assert resolve_dmf_solve_tol({"tol": "0.08"}) == pytest.approx(0.08)
    assert resolve_dmf_solve_tol({"tol": 0.12}) == pytest.approx(0.12)


def test_pinned_dual_inf_tol_survives_when_no_preset_is_requested() -> None:
    cfg = {"ipopt_options": {"dual_inf_tol": 0.07}}
    # None leaves pydmf's solve() from overwriting the caller's IPOPT option.
    assert resolve_dmf_solve_tol(cfg) is None


def test_explicit_tolerance_wins_over_a_pinned_ipopt_option() -> None:
    cfg = {"tol": "loose", "ipopt_options": {"dual_inf_tol": 0.07}}
    assert resolve_dmf_solve_tol(cfg) == "loose"


def test_a_gaussian_preset_is_rejected_with_a_pointer_to_the_right_flag() -> None:
    with pytest.raises(click.ClickException) as excinfo:
        resolve_dmf_solve_tol({"tol": "gau_tight"})
    message = str(excinfo.value)
    assert "tight|middle|loose" in message
    assert "--thresh-gsm" in message


@pytest.mark.parametrize("raw", ["0", "-1", "nan", "inf", "-inf", True, False])
def test_a_nonfinite_or_non_positive_tolerance_is_rejected(raw: object) -> None:
    with pytest.raises(click.ClickException, match="positive float|finite positive"):
        resolve_dmf_solve_tol({"tol": raw})


@pytest.mark.parametrize("command", ["path-opt", "path-search"])
def test_gsm_ignores_a_dormant_dmf_tolerance(
    tmp_path: Path, command: str,
) -> None:
    config = tmp_path / "dormant-dmf.yaml"
    config.write_text("dmf:\n  tol: nan\n", encoding="utf-8")
    smoke = Path(__file__).resolve().parent / "smoke"
    result = CliRunner().invoke(
        root_cli,
        [
            command, "-i", str(smoke / "r.pdb"), str(smoke / "p.pdb"),
            "-q", "-1", "-m", "1", "--mep-mode", "gsm",
            "--max-cycles-gsm", "1", "--no-preopt", "--config", str(config),
            "--dry-run", "--out-dir", str(tmp_path / command),
        ],
    )
    assert result.exit_code == 0, result.output


@pytest.mark.parametrize("command", ["path-opt", "path-search"])
def test_gsm_rejects_an_explicit_invalid_dmf_tolerance(
    tmp_path: Path, command: str,
) -> None:
    smoke = Path(__file__).resolve().parent / "smoke"
    result = CliRunner().invoke(
        root_cli,
        [
            command, "-i", str(smoke / "r.pdb"), str(smoke / "p.pdb"),
            "-q", "-1", "-m", "1", "--mep-mode", "gsm",
            "--thresh-dmf", "nan", "--max-cycles-gsm", "1", "--no-preopt",
            "--dry-run", "--out-dir", str(tmp_path / command),
        ],
    )
    assert result.exit_code != 0
    assert "finite positive" in result.output


def test_all_dry_run_rejects_an_explicit_invalid_dmf_tolerance(
    tmp_path: Path,
) -> None:
    smoke = Path(__file__).resolve().parent / "smoke"
    result = CliRunner().invoke(
        root_cli,
        [
            "all", "-i", str(smoke / "r.pdb"), str(smoke / "p.pdb"),
            "-q", "-1", "-m", "1", "--mep-mode", "gsm",
            "--thresh-dmf", "nan", "--dry-run",
            "--out-dir", str(tmp_path / "all"),
        ],
    )
    assert result.exit_code != 0
    assert "finite positive" in result.output


def test_all_show_config_reports_each_threshold_owner(tmp_path: Path) -> None:
    smoke = Path(__file__).resolve().parent / "smoke"
    result = CliRunner().invoke(
        root_cli,
        [
            "all", "-i", str(smoke / "r.pdb"), str(smoke / "p.pdb"),
            "-q", "-1", "-m", "1", "--thresh", "gau",
            "--thresh-gsm", "gau_loose", "--thresh-dmf", "middle",
            "--show-config", "--dry-run", "--out-dir", str(tmp_path / "all"),
        ],
    )
    assert result.exit_code == 0, result.output
    assert "thresh: gau" in result.output
    assert "thresh_gsm: gau_loose" in result.output
    assert "thresh_dmf: middle" in result.output
