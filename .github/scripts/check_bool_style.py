#!/usr/bin/env python3
"""Reject legacy value-style booleans in hand-authored user guidance."""

from __future__ import annotations

import re
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from pdb2reaction.cli.app import (
    _COMMAND_BOOL_TOGGLE_OPTIONS,
    _COMMAND_BOOL_VALUE_OPTIONS,
)


VALUE_RE = re.compile(
    r"(?P<option>--[a-z0-9][a-z0-9-]*)(?:[ \t]+|=)"
    r"(?P<value>true|false|yes|no|on|off)\b",
    re.IGNORECASE,
)


def _authored_markdown() -> list[Path]:
    files = [REPO_ROOT / "README.md", REPO_ROOT / "CONTRIBUTING.md"]
    files.extend(
        path for path in (REPO_ROOT / "docs").rglob("*.md")
        if "reference/commands" not in path.as_posix()
        and "docs/_build" not in path.as_posix()
    )
    files.extend((REPO_ROOT / "skills").rglob("*.md"))
    files.extend((REPO_ROOT / "tests" / "smoke").glob("*.sh"))
    return sorted(set(files))


def main() -> int:
    known = {"--flag"}
    for table in (_COMMAND_BOOL_VALUE_OPTIONS, _COMMAND_BOOL_TOGGLE_OPTIONS):
        for options in table.values():
            known.update(options)

    errors: list[str] = []
    for path in _authored_markdown():
        text = path.read_text(encoding="utf-8")
        for match in VALUE_RE.finditer(text):
            if match.group("option").lower() not in known:
                continue
            line = text.count("\n", 0, match.start()) + 1
            errors.append(
                f"{path.relative_to(REPO_ROOT)}:{line}: use --flag / --no-flag; "
                f"found {match.group(0)!r}"
            )

    if errors:
        print("\n".join(errors))
        return 1
    print("Canonical boolean style OK in authored docs, skills, and smoke scripts.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
