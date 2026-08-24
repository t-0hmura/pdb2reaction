#!/usr/bin/env python3
"""Reject legacy value-style booleans in hand-authored user guidance.

The set of boolean options is discovered from the live CLI command graph
(``docs_command_contract.live_bool_options``) rather than a hand-maintained
table, so a newly added toggle is covered automatically. Runtime ``bool_compat``
parsing is unchanged: legacy ``--flag True/False`` still works at the CLI; only
authored guidance is required to use the canonical ``--flag`` / ``--no-flag``.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))

import docs_command_contract as contract  # noqa: E402


VALUE_RE = re.compile(
    r"(?P<option>--[a-z0-9][a-z0-9-]*)(?:[ \t]+|=)"
    r"(?P<value>true|false|yes|no|on|off)\b",
    re.IGNORECASE,
)


def main() -> int:
    known = set(contract.live_bool_options())
    known.add("--flag")

    errors: list[str] = []
    for path in contract.authored_bool_style_paths():
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
