"""Keep literal smoke-test options synchronized with the Click CLI."""

from __future__ import annotations

import re
import shlex
from pathlib import Path

import click

from pdb2reaction.cli import cli


_COMMAND_RE = re.compile(
    r"(?<![-\w])pdb2reaction\s+([a-z][a-z0-9-]*)\b"
)


def _unknown_long_options(script: str) -> list[str]:
    """Return one diagnostic message per literal long option absent from its
    named subcommand."""

    root_context = click.Context(cli)
    commands = {
        name: cli.get_command(root_context, name)
        for name in cli.list_commands(root_context)
    }
    errors: list[str] = []
    for line_number, line in enumerate(
        script.replace("\\\n", " ").splitlines(), start=1
    ):
        if line.lstrip().startswith("#"):
            continue
        match = _COMMAND_RE.search(line)
        if match is None:
            continue
        name = match.group(1)
        # Root-level implicit `all` and the explicit `all` adapter forward a
        # staged option surface and are covered by the generated help/CLI
        # contracts. This guard targets strict named subcommands.
        if name not in commands or name == "all":
            continue
        command = commands[name]
        assert command is not None
        allowed = {
            option
            for parameter in command.params
            for option in (
                *getattr(parameter, "opts", ()),
                *getattr(parameter, "secondary_opts", ()),
            )
        }
        tokens = shlex.split(line[match.start() :], comments=True)
        literal_long_options = {
            token.split("=", 1)[0]
            for token in tokens
            if token.startswith("--")
        }
        for option in sorted(literal_long_options - allowed):
            errors.append(f"line {line_number}: {name} does not accept {option}")
    return errors


def test_smoke_literal_options_match_click_commands() -> None:
    script = (
        Path(__file__).parents[1] / "tests" / "smoke" / "run.sh"
    ).read_text(encoding="utf-8")
    assert _unknown_long_options(script) == []


def test_smoke_contract_rejects_unknown_path_search_option() -> None:
    script = (
        "pdb2reaction path-search -i r.pdb p.pdb "
        "--no-preopt --no-endopt --out-dir out\n"
    )
    assert _unknown_long_options(script) == [
        "line 1: path-search does not accept --no-endopt"
    ]
