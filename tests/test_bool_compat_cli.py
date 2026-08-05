"""CLI-level bool compatibility regression tests."""

from __future__ import annotations

import click
from click.testing import CliRunner

from pdb2reaction.cli import cli as root_cli
from pdb2reaction.cli.bool_compat import normalize_argv_option_names


def _is_unavailable_command(cmd: click.Command | None) -> bool:
    if cmd is None:
        return True
    help_text = (cmd.help or "").strip()
    return help_text.startswith("[Unavailable]")


def test_option_name_normalization_stops_at_end_of_options() -> None:
    argv = ["bond-summary", "--JSON", "--", "--MixedCase.xyz"]
    assert normalize_argv_option_names(argv) == [
        "bond-summary",
        "--json",
        "--",
        "--MixedCase.xyz",
    ]


def test_all_bool_options_accept_toggle_and_value_styles() -> None:
    runner = CliRunner()
    ctx = click.Context(root_cli)

    tested = 0
    for command_name in root_cli.list_commands(ctx):
        command = root_cli.get_command(ctx, command_name)
        if _is_unavailable_command(command):
            continue

        value_opts, toggle_opts, _aliases, single_opts = root_cli._resolve_bool_options(
            ctx, command_name
        )
        option_names = sorted(
            {
                *value_opts,
                *toggle_opts,
                *single_opts,
            }
        )
        for opt in option_names:
            if not opt.startswith("--"):
                continue

            no_opt = f"--no-{opt[2:]}"

            result_flag = runner.invoke(root_cli, [command_name, opt, "--help"])
            assert result_flag.exit_code == 0, (
                f"{command_name} should accept {opt}. Output:\n{result_flag.output}"
            )

            result_no_flag = runner.invoke(root_cli, [command_name, no_opt, "--help"])
            assert result_no_flag.exit_code == 0, (
                f"{command_name} should accept {no_opt}. Output:\n{result_no_flag.output}"
            )

            result_value_false = runner.invoke(
                root_cli, [command_name, opt, "False", "--help"]
            )
            assert result_value_false.exit_code == 0, (
                f"{command_name} should accept '{opt} False'. Output:\n"
                f"{result_value_false.output}"
            )

            # Positive value form: `--flag True` should also be accepted
            # (an earlier version only tested False but bool_compat handles both).
            result_value_true = runner.invoke(
                root_cli, [command_name, opt, "True", "--help"]
            )
            assert result_value_true.exit_code == 0, (
                f"{command_name} should accept '{opt} True'. Output:\n"
                f"{result_value_true.output}"
            )

            tested += 1

    assert tested > 0
