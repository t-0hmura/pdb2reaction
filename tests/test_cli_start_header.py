"""Command-header labeling contracts."""

from __future__ import annotations

import click


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
