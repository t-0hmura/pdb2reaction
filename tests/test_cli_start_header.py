"""Command-header labeling contracts."""

from __future__ import annotations

import click
import pytest


def test_lazy_subcommand_header_uses_context_info_name(monkeypatch) -> None:
    from pdb2reaction.cli import app
    from pdb2reaction.core import utils

    messages = []

    @click.command(name="cli")
    def command():
        pass

    context = click.Context(command, info_name="sp")
    monkeypatch.setattr(app.sys, "argv", ["pdb2reaction", "sp"])
    monkeypatch.setattr(utils, "is_child_mode", lambda: False)
    monkeypatch.setattr(utils, "verbose_level", lambda: 2)
    monkeypatch.setattr(
        utils,
        "emit",
        lambda message, **kwargs: messages.append(str(message)),
    )

    app._emit_start_header(context)

    assert "[mode] sp" in messages
    assert "[mode] cli" not in messages


@pytest.mark.parametrize(
    ("argv", "expected"),
    [
        (["pdb2reaction", "bond-summary", "--json"], True),
        (["pdb2reaction", "bond-summary", "--json=true"], True),
        (["pdb2reaction", "bond-summary", "--json", "yes"], True),
        (["pdb2reaction", "bond-summary", "--json=false"], False),
        (["pdb2reaction", "bond-summary", "--json", "--no-json"], False),
        (["pdb2reaction", "sp", "--json=true"], False),
    ],
)
def test_json_stdout_detection_covers_legacy_boolean_forms(
    argv: list[str],
    expected: bool,
) -> None:
    from pdb2reaction.cli.app import _requests_stdout_json

    assert _requests_stdout_json(argv) is expected


@pytest.mark.parametrize(
    ("argv", "expected"),
    [
        (["pdb2reaction", "opt", "--dry-run"], True),
        (["pdb2reaction", "opt", "--dry-run=true"], True),
        (["pdb2reaction", "opt", "--dry-run", "True"], True),
        (["pdb2reaction", "opt", "--dry-run=false"], False),
        (["pdb2reaction", "opt", "--dry-run", "--no-dry-run"], False),
        (["pdb2reaction", "opt", "--no-dry-run"], False),
        (["pdb2reaction", "opt", "--out-dir", "/tmp/dry-run-dir"], False),
        (["pdb2reaction", "opt"], False),
    ],
)
def test_dry_run_detection_covers_legacy_boolean_forms(
    argv: list[str],
    expected: bool,
) -> None:
    from pdb2reaction.cli.app import _requests_dry_run

    assert _requests_dry_run(argv) is expected


@pytest.mark.parametrize("dry_run", [True, False])
def test_a_dry_run_keeps_the_version_line_without_the_artwork(
    monkeypatch, dry_run: bool,
) -> None:
    """The artwork is a run banner; a dry run only reports a plan."""
    from pdb2reaction.cli import app
    from pdb2reaction.core import utils

    messages = []

    @click.command(name="cli")
    def command():
        pass

    argv = ["pdb2reaction", "opt"] + (["--dry-run"] if dry_run else [])
    monkeypatch.setattr(app.sys, "argv", argv)
    monkeypatch.setattr(utils, "is_child_mode", lambda: False)
    monkeypatch.setattr(utils, "verbose_level", lambda: 2)
    monkeypatch.setattr(
        utils,
        "emit",
        lambda message, **kwargs: messages.append(str(message)),
    )

    app._emit_start_header(click.Context(command, info_name="opt"))

    joined = "\n".join(messages)
    assert "pdb2reaction ver." in joined
    assert (app._PDB2REACTION_BANNER.strip() in joined) is not dry_run
