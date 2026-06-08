"""Shared logging configurator for pdb2reaction subcommands.

The root CLI group counts -v / -vv occurrences and calls `setup_logging` to
configure stdlib `logging` once per invocation. Subcommands then emit through
`logging.getLogger(__name__)` and gain `-v` visibility automatically:

    # root group callback (cli.py):
    from pdb2reaction.core.logging import setup_logging
    @click.option("-v", "--verbose", count=True, ...)
    def cli(verbose):
        setup_logging(verbose)

    # subcommand module:
    import logging
    logger = logging.getLogger(__name__)
    logger.info("...")  # visible at -v; suppressed by default

Level mapping:
    no -v   -> WARNING  (default; stdout volume unchanged from click.echo)
    -v      -> INFO
    -vv     -> DEBUG

Existing click.echo() calls are unaffected; the logger is opt-in per call site.
"""
from __future__ import annotations

import logging
import sys

_DEFAULT_FORMAT = "%(asctime)s [%(levelname)s] %(name)s: %(message)s"
_VERBOSE_TO_LEVEL = {
    0: logging.WARNING,
    1: logging.INFO,
    2: logging.DEBUG,
}


def setup_logging(verbose: int = 0) -> None:
    """Configure stdlib logging according to -v count. Idempotent."""
    level = _VERBOSE_TO_LEVEL.get(min(verbose, 2), logging.DEBUG)
    root = logging.getLogger()
    for h in list(root.handlers):
        root.removeHandler(h)
    handler = logging.StreamHandler(stream=sys.stderr)
    handler.setFormatter(logging.Formatter(_DEFAULT_FORMAT))
    handler.setLevel(level)
    root.addHandler(handler)
    root.setLevel(level)
