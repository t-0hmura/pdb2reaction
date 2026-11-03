"""MEP threshold surfaces: GSM preset, DMF tolerance, and parent forwarding."""

from __future__ import annotations

from pathlib import Path
from typing import Set

import click
import pytest
from click.testing import CliRunner

from pdb2reaction.cli import cli as root_cli
from pdb2reaction.workflows import all as all_workflow
from pdb2reaction.workflows.path_opt import resolve_dmf_solve_tol
from test_all_root_summary_golden import pdb_frame, write
from _repo_files import require_repo_files

MEP_FLAGS = ("--gsm-param", "--thresh-gsm", "--thresh-dmf")


def _declared_flags(cli: click.Command) -> Set[str]:
    flags: Set[str] = set()
    for param in cli.params:
        flags.update(param.opts)
    return flags


def test_parent_and_mep_children_declare_the_same_flags() -> None:
    from pdb2reaction.workflows.path_opt import cli as path_opt_cli
    from pdb2reaction.workflows.path_search import cli as path_search_cli

    for flag in MEP_FLAGS:
        assert flag in _declared_flags(all_workflow.cli)
        assert flag in _declared_flags(path_opt_cli)
        assert flag in _declared_flags(path_search_cli)


class ReachedChild(BaseException):
    """Escape the parent pipeline at the first child dispatch."""


@pytest.mark.parametrize(
    "child,extra", [("path-opt", []), ("path_search", ["--refine-path"])],
)
def test_all_forwards_the_mep_flags_to_its_path_children(
    tmp_path: Path, monkeypatch, child: str, extra: list,
) -> None:
    captured = []

    def first_child(name, _cli, args, **_kwargs):
        captured.append((name, list(args)))
        raise ReachedChild

    monkeypatch.setattr(all_workflow, "_run_cli_main", first_child)
    inputs = []
    for frame in (0, 2):
        inputs += ["-i", str(write(tmp_path / f"input{frame}.pdb", pdb_frame(frame)))]
    values = {"--gsm-param": "energy", "--thresh-gsm": "gau", "--thresh-dmf": "0.08"}
    args = ["all", *inputs, "-q", "0", "-m", "1", "--out-dir", str(tmp_path / "out"),
            "--no-preopt", "--no-freeze-links", *extra]
    for flag, value in values.items():
        args += [flag, value]
    with pytest.raises(ReachedChild):
        CliRunner().invoke(root_cli, args, catch_exceptions=False)

    ((name, child_args),) = captured
    assert name == child
    for flag, value in values.items():
        assert child_args[child_args.index(flag) + 1] == value


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
@pytest.mark.parametrize(
    ("mode", "nodes", "source", "cli_override", "expected_exit"),
    [
        ("gsm", -1, "cli", None, 2),
        ("gsm", 1, "cli", None, 2),
        ("gsm", 2, "cli", None, 0),
        ("dmf", 1, "cli", None, 0),
        ("gsm", 1, "config", None, 2),
        ("gsm", 1, "config", 2, 0),
        ("gsm", 2, "config", 1, 2),
    ],
)
def test_path_dry_run_validates_effective_gsm_nodes(
    tmp_path: Path, command: str, mode: str, nodes: int,
    source: str, cli_override: int | None, expected_exit: int,
) -> None:
    require_repo_files("tests/smoke")
    smoke = Path(__file__).resolve().parent / "smoke"
    args = [
        command, "-i", str(smoke / "r.pdb"), str(smoke / "p.pdb"),
        "-q", "-1", "-m", "1", "--mep-mode", mode,
        "--no-preopt", "--dry-run", "--out-dir", str(tmp_path / command),
    ]
    key = "gs.max_nodes" if command == "path-opt" else "search.max_nodes_segment"
    if source == "cli":
        args.extend(["--max-nodes", str(nodes)])
    else:
        section, field = key.split(".")
        config = tmp_path / "nodes.yaml"
        config.write_text(f"{section}:\n  {field}: {nodes}\n", encoding="utf-8")
        args.extend([f"--{source}", str(config)])
    if cli_override is not None:
        args.extend(["--max-nodes", str(cli_override)])
    result = CliRunner().invoke(root_cli, args)
    assert result.exit_code == expected_exit, result.output
    if expected_exit:
        assert "--max-nodes" in result.output
        assert key in result.output
        assert "at least 2" in result.output
        assert "Input command is valid" not in result.output
    else:
        assert "Input command is valid" in result.output


@pytest.mark.parametrize("command", ["path-opt", "path-search"])
def test_gsm_ignores_a_dormant_dmf_tolerance(
    tmp_path: Path, command: str,
) -> None:
    require_repo_files("tests/smoke")
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
    require_repo_files("tests/smoke")
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
    require_repo_files("tests/smoke")
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
    require_repo_files("tests/smoke")
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


@pytest.mark.parametrize(
    ("mode_args", "expected_path", "expected_post"),
    [
        ([], "grad", "hess"),
        (["--opt-mode", "grad"], "grad", "grad"),
        (["--opt-mode", "hess", "--opt-mode-post", "grad"], "hess", "grad"),
    ],
)
def test_all_show_config_reports_effective_optimizer_modes(
    tmp_path: Path,
    mode_args: list[str],
    expected_path: str,
    expected_post: str,
) -> None:
    require_repo_files("tests/smoke")
    smoke = Path(__file__).resolve().parent / "smoke"
    result = CliRunner().invoke(
        root_cli,
        [
            "all", "-i", str(smoke / "r.pdb"), str(smoke / "p.pdb"),
            "-q", "-1", "--show-config", "--dry-run",
            "--out-dir", str(tmp_path / "all"), *mode_args,
        ],
    )

    assert result.exit_code == 0, result.output
    assert f"path_opt_mode: {expected_path}" in result.output
    assert f"post_opt_mode: {expected_post}" in result.output
    assert f"ts_opt_mode: {expected_post}" in result.output
    assert f"endpoint_opt_mode: {expected_post}" in result.output
    assert f"opt_mode_post: {expected_post}" in result.output


@pytest.mark.parametrize("command", ["path-opt", "path-search"])
def test_explicit_gsm_param_overrides_yaml(tmp_path: Path, command: str) -> None:
    require_repo_files("tests/smoke")
    smoke = Path(__file__).resolve().parent / "smoke"
    config = tmp_path / "gsm.yaml"
    config.write_text("gs:\n  param: equi\n", encoding="utf-8")
    result = CliRunner().invoke(
        root_cli,
        [
            command,
            "-i",
            str(smoke / "r.pdb"),
            str(smoke / "p.pdb"),
            "-q",
            "-1",
            "--config",
            str(config),
            "--gsm-param",
            "energy",
            "--show-config",
            "--dry-run",
            "-v",
            "3",
            "--out-dir",
            str(tmp_path / command),
        ],
    )
    assert result.exit_code == 0, result.output
    assert "gsm_param: energy" in result.output
