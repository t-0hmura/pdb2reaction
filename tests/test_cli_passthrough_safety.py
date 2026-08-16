"""Regression tests for permissive legacy CLI grammars."""

from __future__ import annotations

from pathlib import Path

import pytest
from click import UsageError
from click.testing import CliRunner

from pdb2reaction.cli import cli as root_cli
from pdb2reaction.core.utils import collect_option_values, reject_option_like_extra_args
from pdb2reaction.workflows import all as all_workflow


def test_root_invocation_resets_charge_multiplicity_override() -> None:
    from pdb2reaction.core.utils import (
        set_allow_charge_mult_mismatch,
        validate_charge_spin,
    )

    set_allow_charge_mult_mismatch(True)
    result = CliRunner().invoke(root_cli, ["sp", "--help"])

    assert result.exit_code == 0, result.output
    with pytest.raises(ValueError, match="electron count inconsistent"):
        validate_charge_spin(["H"], 0, 1)


@pytest.mark.parametrize(
    ("command", "option"),
    [
        ("dft", "--convert-files"),
        ("sp", "--convert-files"),
        ("sp", "--print-every"),
    ],
)
def test_noop_options_are_removed(command: str, option: str) -> None:
    result = CliRunner().invoke(root_cli, [command, option])

    assert result.exit_code == 2
    assert f"No such option: {option}" in result.output


def _xyz(path: Path) -> Path:
    path.write_text("1\nH\nH 0.0 0.0 0.0\n", encoding="utf-8")
    return path


def test_all_calculator_cli_values_override_yaml() -> None:
    cfg = all_workflow._build_calc_cfg(
        charge=-1,
        spin=2,
        workers=8,
        workers_per_node=4,
        yaml_cfg={
            "calc": {
                "charge": 0,
                "spin": 1,
                "workers": 2,
                "workers_per_node": 2,
            }
        },
    )

    assert cfg["charge"] == -1
    assert cfg["spin"] == 2
    assert cfg["workers"] == 8
    assert cfg["workers_per_node"] == 4


def test_child_runtime_forwarding_is_explicit_only() -> None:
    argv: list[str] = []
    all_workflow._append_explicit_child_runtime_argv(argv)
    assert argv == []

    all_workflow._append_explicit_child_runtime_argv(
        argv,
        workers=1,
        workers_per_node=2,
        dmf_backend="CPU",
        backend="UMA",
        solvent="NONE",
        solvent_model="ALPB",
        max_nodes=12,
        preopt=False,
    )
    assert argv == [
        "--dmf-backend",
        "cpu",
        "--workers",
        "1",
        "--workers-per-node",
        "2",
        "--backend",
        "uma",
        "--solvent",
        "none",
        "--solvent-model",
        "alpb",
        "--max-nodes",
        "12",
        "--no-preopt",
    ]


def test_only_claimed_bare_legacy_values_remain_allowed() -> None:
    reject_option_like_extra_args(
        ["second.pdb", "third.pdb"],
        allowed_values=["second.pdb", "third.pdb"],
    )
    with pytest.raises(UsageError, match="Unexpected extra argument: orphan"):
        reject_option_like_extra_args(["orphan"])

    reject_option_like_extra_args(
        ["same.pdb"],
        allowed_values=["same.pdb", "same.pdb"],
        consumed_values=["same.pdb"],
    )
    with pytest.raises(UsageError, match="Unexpected extra argument: same.pdb"):
        reject_option_like_extra_args(
            ["same.pdb"],
            allowed_values=["same.pdb"],
            consumed_values=["same.pdb"],
        )


def test_variadic_collector_preserves_equals_attached_and_grouped_order() -> None:
    argv = [
        "--input=first.pdb", "-i", "second.pdb", "third.pdb",
        "--output=first.out", "-osecond.out", "--ref-pdb=template.pdb",
        "--input=fourth.pdb",
    ]
    assert collect_option_values(argv, ("-i", "--input")) == [
        "first.pdb", "second.pdb", "third.pdb", "fourth.pdb",
    ]
    assert collect_option_values(argv, ("-o", "--output")) == [
        "first.out", "second.out",
    ]
    assert collect_option_values(argv, ("--ref-pdb",)) == ["template.pdb"]


@pytest.mark.parametrize("repeated", [False, True])
def test_all_accepts_grouped_and_repeated_inputs(
    tmp_path: Path,
    repeated: bool,
) -> None:
    first = _xyz(tmp_path / "first.xyz")
    second = _xyz(tmp_path / "second.xyz")
    inputs = (
        ["-i", str(first), "-i", str(second)]
        if repeated
        else ["-i", str(first), str(second)]
    )

    result = CliRunner().invoke(
        root_cli,
        ["all", *inputs, "-q", "0", "-m", "2", "--dry-run"],
    )

    assert result.exit_code == 0, result.output
    assert "Planned stages: path_opt." in result.output


def test_all_rejects_trailing_orphan_equal_to_consumed_input(tmp_path: Path) -> None:
    first = _xyz(tmp_path / "first.xyz")
    result = CliRunner().invoke(
        root_cli,
        ["all", "-i", str(first), "-q", "0", "--tsopt", "--dry-run", str(first)],
    )

    assert result.exit_code == 2
    assert f"Unexpected extra argument: {first}" in result.output


@pytest.mark.parametrize("unknown", ["--typo-option", "-Z"])
@pytest.mark.parametrize("command", ["all", "extract", "path-search"])
def test_passthrough_commands_reject_unknown_options(
    tmp_path: Path,
    command: str,
    unknown: str,
) -> None:
    first = _xyz(tmp_path / "first.xyz")
    second = _xyz(tmp_path / "second.xyz")
    if command == "all":
        argv = ["all", "-i", str(first), "-q", "0", "--dry-run", unknown]
    elif command == "extract":
        argv = ["extract", "-i", str(first), "-c", "LIG", unknown]
    else:
        argv = [
            "path-search",
            "-i",
            str(first),
            "-i",
            str(second),
            "-q",
            "0",
            "--dry-run",
            unknown,
        ]

    result = CliRunner().invoke(root_cli, argv)

    assert result.exit_code == 2
    assert f"No such option: {unknown}" in result.output


@pytest.mark.parametrize("command", ["all", "extract", "path-search"])
def test_passthrough_commands_reject_unclaimed_bare_arguments(
    tmp_path: Path,
    command: str,
) -> None:
    first = _xyz(tmp_path / "first.xyz")
    second = _xyz(tmp_path / "second.xyz")
    if command == "all":
        argv = ["all", "-i", str(first), "-q", "0", "--dry-run", "orphan"]
    elif command == "extract":
        argv = ["extract", "-i", str(first), "-c", "LIG", "--out-json", "orphan"]
    else:
        argv = [
            "path-search", "-i", str(first), "-i", str(second),
            "-q", "0", "--dry-run", "orphan",
        ]

    result = CliRunner().invoke(root_cli, argv)

    assert result.exit_code == 2
    assert "Unexpected extra argument: orphan" in result.output
