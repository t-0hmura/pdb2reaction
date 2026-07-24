"""Regression tests for the skill-command validator."""

from __future__ import annotations

import importlib.util
from pathlib import Path


def _load_module():
    path = Path(__file__).parents[1] / ".github" / "scripts" / "check_skill_commands.py"
    spec = importlib.util.spec_from_file_location("check_skill_commands", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_unknown_subcommand_is_rejected() -> None:
    module = _load_module()
    flags = module._collect_subcommand_flags()
    assert module._check_command("pdb2reaction tsotp -i ts.xyz", flags) == [
        "<unknown-subcommand:tsotp>"
    ]


def test_root_all_alias_is_validated() -> None:
    module = _load_module()
    flags = module._collect_subcommand_flags()
    assert module._check_command(
        "pdb2reaction -i R.pdb P.pdb --refine-path", flags
    ) == []


def test_unknown_flag_is_rejected() -> None:
    module = _load_module()
    flags = module._collect_subcommand_flags()
    assert module._check_command(
        "pdb2reaction opt -i input.xyz --no-such-flag", flags,
    ) == ["--no-such-flag"]


def test_non_shell_fence_does_not_hide_following_shell_command() -> None:
    module = _load_module()
    text = """\
```python
pdb2reaction opt -i ignored.xyz --python-only
```
```bash
pdb2reaction opt -i checked.xyz --no-such-flag
```
    """
    assert list(module._iter_command_lines(text)) == [
        (5, "pdb2reaction opt -i checked.xyz --no-such-flag"),
    ]
