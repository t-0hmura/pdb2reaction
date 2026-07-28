"""M38: one product-local advanced-help implementation, shared by `all` and
lazily-loaded subcommands.

Falsifiers:
1. Default ``all --help`` shows the primary set only; ``--help-advanced`` shows
   every live (advanced) option.
2. Repeating basic/advanced/basic in one process leaves basic output and the
   hidden state byte-identical (the ``try/finally`` restores every option).
3. The single callback drives both the ``all`` command and a generic
   (lazy-subcommand-style) configured command via a plain Click context.
4. Exactly one advanced-help callback implementation exists in the product.
"""

from __future__ import annotations

import ast
import pathlib

import click
import click.exceptions
import pytest
from click.testing import CliRunner

from pdb2reaction.cli.help_pages import (
    _hide_advanced_options,
    _show_advanced_subcommand_help,
)
from pdb2reaction.workflows.all import cli as all_cli
from pdb2reaction.workflows.all import _ALL_PRIMARY_HELP_OPTIONS


def _run_callback(command: click.Command) -> None:
    """Invoke the advanced-help callback on ``command`` with value=True."""
    ctx = click.Context(command, info_name=command.name or "cmd")
    try:
        _show_advanced_subcommand_help(ctx, None, True)
    except (SystemExit, click.exceptions.Exit):
        return
    pytest.fail("advanced-help callback did not exit")


def test_basic_help_hides_advanced_options_that_advanced_help_shows():
    runner = CliRunner()
    hidden = all_cli._advanced_hidden_options
    assert hidden, "expected `all` to hide at least one advanced option"

    basic = runner.invoke(all_cli, ["--help"])
    advanced = runner.invoke(all_cli, ["--help-advanced"])
    assert basic.exit_code == 0, basic.output
    assert advanced.exit_code == 0, advanced.output

    # A representative primary option stays visible in the basic help.
    assert "--tsopt" in basic.output

    # Compare on each option's rendered help-record column (e.g. "-r, --radius
    # FLOAT"); it carries the metavar so it cannot collide with prose that merely
    # mentions a flag. Every hidden option's column is present exactly once in
    # advanced help and absent from basic help.
    ctx = click.Context(all_cli, info_name="all")
    for opt in hidden:
        assert opt.hidden
        was_hidden = opt.hidden
        opt.hidden = False
        try:
            record = opt.get_help_record(ctx)
        finally:
            opt.hidden = was_hidden
        assert record is not None
        col = record[0]
        assert advanced.output.count(col) == 1, f"{col!r} not shown exactly once in advanced help"
        assert col not in basic.output, f"{col!r} leaked into basic help"


def test_repeated_help_calls_preserve_hidden_state_and_output():
    runner = CliRunner()
    assert all(opt.hidden for opt in all_cli._advanced_hidden_options)

    basic1 = runner.invoke(all_cli, ["--help"])
    advanced = runner.invoke(all_cli, ["--help-advanced"])
    basic2 = runner.invoke(all_cli, ["--help"])

    assert basic1.exit_code == basic2.exit_code == 0
    assert advanced.exit_code == 0
    # The try/finally re-hides everything: state is unchanged after --help-advanced.
    assert all(opt.hidden for opt in all_cli._advanced_hidden_options)
    # Basic output is identical before and after the advanced call.
    assert basic1.output == basic2.output


def test_single_callback_drives_all_command_and_generic_subcommand():
    # (a) the `all` command
    _run_callback(all_cli)
    assert all(opt.hidden for opt in all_cli._advanced_hidden_options)

    # (b) a generic command configured exactly like a lazily-loaded subcommand
    demo = click.Command(
        "demo",
        params=[
            click.Option(["-i", "--input"]),
            click.Option(["--advanced-opt"]),
        ],
    )
    _hide_advanced_options(demo, frozenset({"-i", "--input"}))
    adv_opt = next(p for p in demo.params if "--advanced-opt" in p.opts)
    assert adv_opt.hidden
    _run_callback(demo)
    # Restored to hidden after the callback's finally block.
    assert adv_opt.hidden


def test_hide_advanced_options_is_idempotent():
    demo = click.Command("demo", params=[click.Option(["--advanced-opt"])])
    _hide_advanced_options(demo, frozenset({"-i"}))
    first = demo._advanced_hidden_options
    # A second call must not re-hide / grow the recorded set.
    _hide_advanced_options(demo, frozenset({"-i"}))
    assert demo._advanced_hidden_options is first


def test_exactly_one_advanced_help_callback_in_product():
    root = pathlib.Path(__file__).resolve().parents[1] / "pdb2reaction"
    impls: list[str] = []
    for path in root.rglob("*.py"):
        tree = ast.parse(path.read_text(encoding="utf-8"))
        for node in ast.walk(tree):
            if isinstance(node, ast.FunctionDef) and node.name.startswith("_show_advanced"):
                impls.append(f"{path.relative_to(root)}::{node.name}")
    assert impls == ["cli/help_pages.py::_show_advanced_subcommand_help"], impls
