"""Dependency-light product output adapter and verbosity state."""

from __future__ import annotations

from typing import Any

import click


_TAG_AWARE_MARKER = "__pdb2reaction_private_echo_tags__"
_VERBOSE_LEVEL = 0


def set_verbose_level(level: int) -> None:
    """Record the process-wide CLI verbosity level."""
    global _VERBOSE_LEVEL
    _VERBOSE_LEVEL = max(0, min(3, int(level)))


def verbose_level() -> int:
    """Return the process-wide CLI verbosity level."""
    return _VERBOSE_LEVEL


def is_verbose() -> bool:
    """Return whether detail output is enabled."""
    return _VERBOSE_LEVEL >= 2


def emit(
    message: str = "",
    *,
    narrative: bool | None = None,
    detail: bool = False,
    force: bool = False,
    raw_path: bool = False,
    **kwargs: Any,
) -> None:
    """Echo one line while safely consuming private tags before bootstrap."""

    if narrative is None:
        narrative = not (detail or force or raw_path)

    if not getattr(click.echo, _TAG_AWARE_MARKER, False):
        click.echo(message, **kwargs)
        return
    click.echo(
        message,
        narrative=narrative,
        detail=detail,
        force=force,
        raw_path=raw_path,
        **kwargs,
    )
