#!/usr/bin/env python3
"""Live-CLI / authored-command contract for docs-quality checks.

It defines the small ``AuthoredCommand(path, line, text, executable)`` shape and
pure functions shared by the docs command smoke, the example-script contract,
and the canonical-boolean style checker:

* explicit authored/public roots and the exact advertised shell examples;
* Markdown fenced-command and exact public-shell extraction (with file+line);
* live Click command/option discovery;
* authored option-name validation that RETAINS bracket/placeholder commands for
  STATIC validation and classifies execution eligibility separately;
* live-derived canonical ``--flag`` / ``--no-flag`` boolean option discovery.

Product constants (``TOOL_NAME``, the root CLI import, the example paths) stay
local to this file.
"""

from __future__ import annotations

import re
import shlex
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import click

REPO_ROOT = Path(__file__).resolve().parents[2]
TOOL_NAME = "pdb2reaction"

sys.path.insert(0, str(REPO_ROOT))

from pdb2reaction.cli import cli as root_cli  # noqa: E402

# Markdown fence languages that carry runnable shell commands.
CODE_FENCE_LANGS = frozenset({"", "bash", "sh", "shell", "console"})

# The exact list of advertised, runnable example scripts. Making this a
# checked contract means a renamed/removed script fails CI instead of silently
# dropping out of docs-quality automation.
PUBLIC_SHELL_EXAMPLES: tuple[Path, ...] = (REPO_ROOT / "examples" / "run.sh",)

# Angle-bracket placeholders (``<input.pdb>``, ``<charge>``) are synopsis
# notation and are never directly runnable.
_ANGLE_PLACEHOLDER_RE = re.compile(r"<[^>]*>")
_SHELL_BREAKS = frozenset({"|", "||", "&&", ";"})

# Value-style boolean literals rejected in authored guidance.
BOOL_VALUE_LITERALS = ("true", "false", "yes", "no", "on", "off")


@dataclass(frozen=True)
class AuthoredCommand:
    """One authored CLI command extracted from a public doc/shell source.

    ``executable`` is False when the command carries synopsis placeholders
    (``<input.pdb>``, an optional ``[--flag]`` group) that cannot be run as
    written. Such commands are STILL retained for static option validation; a
    quoted list/tuple value such as ``-s '[("A","B",1.6)]'`` is data, not
    synopsis, and stays executable.
    """

    path: Path
    line: int
    text: str
    executable: bool

    @property
    def rel(self) -> str:
        return rel_to_repo(self.path)

    @property
    def location(self) -> str:
        return f"{self.rel}:{self.line}"


def rel_to_repo(path: Path) -> str:
    """Repo-relative display path; falls back to the absolute path if outside."""
    try:
        return str(path.relative_to(REPO_ROOT))
    except ValueError:
        return str(path)


def is_executable(text: str) -> bool:
    """Classify whether an authored command can be run verbatim."""
    if _ANGLE_PLACEHOLDER_RE.search(text):
        return False
    try:
        tokens = shlex.split(text)
    except ValueError:
        return False
    for tok in tokens:
        # ``[--flag]`` / ``[options]`` optional-synopsis notation is not
        # runnable; a bracketed list/tuple value (contains ``(``) is real data.
        if tok.startswith("[") and "(" not in tok:
            return False
        if "<" in tok or ">" in tok:
            return False
    return True


def _make(path: Path, line: int, text: str) -> AuthoredCommand:
    return AuthoredCommand(path=path, line=line, text=text, executable=is_executable(text))


# --------------------------------------------------------------------------- #
# Explicit authored / public roots
# --------------------------------------------------------------------------- #
def docs_markdown_paths() -> list[Path]:
    """All docs Markdown pages (source of the command smoke surface)."""
    docs = REPO_ROOT / "docs"
    if not docs.is_dir():
        return []
    return sorted(docs.rglob("*.md"))


def authored_bool_style_paths() -> list[Path]:
    """Hand-authored surfaces scanned for canonical-boolean style.

    Generated command-reference pages and Sphinx build output are excluded:
    they are regenerated from the live CLI, not hand-authored.
    """
    paths: list[Path] = [REPO_ROOT / "README.md", REPO_ROOT / "CONTRIBUTING.md"]
    docs = REPO_ROOT / "docs"
    if docs.is_dir():
        paths.extend(
            p
            for p in docs.rglob("*.md")
            if "reference/commands" not in p.as_posix()
            and "docs/_build" not in p.as_posix()
        )
    skills = REPO_ROOT / "skills"
    if skills.is_dir():
        paths.extend(skills.rglob("*.md"))
    smoke = REPO_ROOT / "tests" / "smoke"
    if smoke.is_dir():
        paths.extend(smoke.glob("*.sh"))
    return sorted({p for p in paths if p.exists()})


# --------------------------------------------------------------------------- #
# Markdown fenced-command and exact public-shell extraction
# --------------------------------------------------------------------------- #
def _iter_command_blocks(path: Path) -> Iterable[list[tuple[int, str]]]:
    """Yield each shell-code fenced block as a list of (lineno, raw_line)."""
    in_fence = False
    fence_lang = ""
    block: list[tuple[int, str]] = []
    for lineno, line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
        stripped = line.strip()
        if stripped.startswith("```"):
            marker = stripped[3:].strip().lower()
            if not in_fence:
                in_fence = True
                fence_lang = marker
                block = []
            else:
                if fence_lang in CODE_FENCE_LANGS:
                    yield block
                in_fence = False
                fence_lang = ""
                block = []
            continue
        if in_fence:
            block.append((lineno, line))


def _commands_from_lines(
    path: Path, numbered_lines: Iterable[tuple[int, str]]
) -> list[AuthoredCommand]:
    out: list[AuthoredCommand] = []
    current = ""
    start_line: int | None = None
    for lineno, raw in numbered_lines:
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if " #" in line:
            line = line.split(" #", 1)[0].rstrip()
            if not line:
                continue
        if line.startswith("$"):
            line = line[1:].strip()
        if not current:
            start_line = lineno
        current = f"{current} {line}".strip() if current else line
        if current.endswith("\\"):
            current = current[:-1].rstrip()
            continue
        assert start_line is not None
        out.append(_make(path, start_line, current))
        current = ""
        start_line = None
    if current and start_line is not None:
        out.append(_make(path, start_line, current))
    return [cmd for cmd in out if cmd.text.startswith(TOOL_NAME)]


def extract_markdown_commands(paths: Iterable[Path]) -> list[AuthoredCommand]:
    """Extract every ``TOOL_NAME`` command from fenced blocks in ``paths``."""
    result: list[AuthoredCommand] = []
    for path in paths:
        for block in _iter_command_blocks(path):
            result.extend(_commands_from_lines(path, block))
    return result


def extract_shell_commands(path: Path) -> list[AuthoredCommand]:
    """Extract every ``TOOL_NAME`` invocation from a public shell script."""
    numbered = list(enumerate(path.read_text(encoding="utf-8").splitlines(), start=1))
    return _commands_from_lines(path, numbered)


# --------------------------------------------------------------------------- #
# Live Click command / option discovery
# --------------------------------------------------------------------------- #
def subcommand_from_tokens(tokens: list[str]) -> str:
    if len(tokens) < 2 or tokens[1].startswith("-"):
        return "all"
    if tokens[1].startswith("<") or tokens[1].startswith("["):
        # Generic usage templates still participate in token validation; use
        # `all` as the representative command for shared help flags.
        return "all"
    return tokens[1]


def _allowed_options(ctx: click.Context, subcmd: str) -> set[str] | None:
    command = root_cli.get_command(ctx, subcmd)
    if command is None:
        return None
    allowed = {
        opt
        for param in command.params
        if hasattr(param, "opts")
        for opt in (*param.opts, *getattr(param, "secondary_opts", ()))
    }
    allowed.update({"-h", "--help", "--help-advanced", "--version"})
    return allowed


def validate_option_names(commands: Iterable[AuthoredCommand]) -> list[str]:
    """Return ``file:line`` errors for unknown authored options (empty = OK).

    Executing a subcommand's ``--help`` only proves that the command exists;
    Click's eager help can return before parsing an invalid option. This static
    pass checks each option token against the live command object and
    deliberately retains examples containing ``<placeholder>`` / ``[...]``
    notation.
    """
    root_ctx = root_cli.make_context(TOOL_NAME, [], resilient_parsing=True)
    errors: list[str] = []
    try:
        for cmd in commands:
            try:
                tokens = shlex.split(cmd.text)
            except ValueError as exc:
                errors.append(f"{cmd.location}: cannot parse shell words: {cmd.text!r}: {exc}")
                continue
            if not tokens or tokens[0] != TOOL_NAME:
                continue

            subcmd = subcommand_from_tokens(tokens)
            allowed = _allowed_options(root_ctx, subcmd)
            if allowed is None:
                errors.append(f"{cmd.location}: unknown subcommand {subcmd!r}: {cmd.text}")
                continue

            start = 2 if len(tokens) > 1 and tokens[1] == subcmd else 1
            for token in tokens[start:]:
                token = token.lstrip("[").rstrip("] ,")
                if token in _SHELL_BREAKS:
                    break
                if token == "--":
                    break
                if not token.startswith("-") or token == "-":
                    continue
                # Negative numeric option values are not option names.
                try:
                    float(token)
                    continue
                except ValueError:
                    pass
                # Synopsis blocks often compress aliases/pairs into one shell
                # word (`-b/--backend`, `--dump/--no-dump`, `--a|--b`); validate
                # each advertised spelling independently.
                names = [
                    part.split("=", 1)[0].removesuffix("...")
                    for part in re.split(r"[/|]", token)
                    if part.startswith("-")
                ]
                for name in names:
                    if name not in allowed:
                        errors.append(
                            f"{cmd.location}: unknown option {name!r} for {subcmd}: {cmd.text}"
                        )
    finally:
        root_ctx.close()
    return errors


def live_bool_options() -> set[str]:
    """Union of canonical boolean option names across every live command.

    Derived from ``DefaultGroup._resolve_bool_options`` (value + toggle +
    single-flag + negative aliases), so a newly added toggle is covered without
    editing a manual table. This reads the live command graph only; it does not
    alter ``bool_compat`` runtime parsing.
    """
    ctx = click.Context(root_cli)
    known: set[str] = set()
    resolve = getattr(root_cli, "_resolve_bool_options", None)
    if resolve is None:  # pragma: no cover - defensive
        return known
    for name in root_cli.list_commands(ctx):
        value_opts, toggle_opts, neg_aliases, single_opts = resolve(ctx, name)
        known |= set(value_opts)
        known |= set(toggle_opts)
        known |= set(single_opts)
        known |= set(neg_aliases.values())
    return known
