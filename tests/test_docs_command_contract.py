"""Contract tests for the live-CLI docs command/option validator (M64, P04).

Loads the local `.github/scripts/docs_command_contract.py` module and checks
that bracket-bearing commands are retained for static validation, execution
eligibility is classified separately, unknown options fail with file:line, and
the canonical-boolean option set is derived from the live command graph.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path


def _load(name: str):
    path = Path(__file__).parents[1] / ".github" / "scripts" / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module  # let dataclasses resolve the module namespace
    spec.loader.exec_module(module)
    return module


contract = _load("docs_command_contract")
bool_style = _load("check_bool_style")


def _cmd(text: str, *, rel: str = "docs/example.md", line: int = 1):
    return contract._make(contract.REPO_ROOT / rel, line, text)


def test_quoted_list_value_stays_executable_and_validates_scan() -> None:
    # A quoted list/tuple is DATA, not synopsis: the command remains
    # execution-eligible and its -s option validates against the live CLI.
    cmd = _cmd('pdb2reaction all -i R.pdb -s "[(12,45,2.20)]"')
    assert cmd.executable is True
    assert contract.validate_option_names([cmd]) == []


def test_angle_placeholder_is_static_only_but_validated() -> None:
    cmd = _cmd("pdb2reaction opt -i <input.pdb> -q <charge> --dry-run")
    assert cmd.executable is False  # synopsis placeholder, not runnable
    assert contract.validate_option_names([cmd]) == []  # still statically valid


def test_optional_synopsis_bracket_is_not_executable() -> None:
    cmd = _cmd("pdb2reaction all -i R.pdb [--tsopt]")
    assert cmd.executable is False


def test_invented_option_fails_with_file_and_line() -> None:
    cmd = _cmd(
        'pdb2reaction all -i R.pdb -s "[(1,2,3.0)]" --invented-option',
        rel="docs/scan.md",
        line=42,
    )
    errors = contract.validate_option_names([cmd])
    assert len(errors) == 1
    assert "docs/scan.md:42" in errors[0]
    assert "--invented-option" in errors[0]


def test_unknown_subcommand_is_reported() -> None:
    cmd = _cmd("pdb2reaction not-a-command --help")
    errors = contract.validate_option_names([cmd])
    assert errors and "unknown subcommand" in errors[0]


def test_live_bool_options_cover_known_toggles() -> None:
    known = contract.live_bool_options()
    for opt in ("--tsopt", "--thermo", "--dry-run", "--climb"):
        assert opt in known


def test_real_docs_commands_validate() -> None:
    commands = contract.extract_markdown_commands(contract.docs_markdown_paths())
    assert commands, "no docs commands extracted"
    assert contract.validate_option_names(commands) == []


def test_value_style_bool_checker_rejects_known_bool(
    tmp_path: Path, monkeypatch, capsys,
) -> None:
    guide = tmp_path / "guide.md"
    guide.write_text(
        "pdb2reaction all -i R.pdb --deterministic True\n",
        encoding="utf-8",
    )
    monkeypatch.setattr(bool_style, "REPO_ROOT", tmp_path)
    monkeypatch.setattr(
        bool_style.contract,
        "authored_bool_style_paths",
        lambda: [guide],
    )

    assert bool_style.main() == 1
    assert "--deterministic True" in capsys.readouterr().out


def test_value_style_bool_checker_ignores_non_bool_near_miss(
    tmp_path: Path, monkeypatch, capsys,
) -> None:
    guide = tmp_path / "guide.md"
    guide.write_text(
        "pdb2reaction all -i R.pdb --label True\n",
        encoding="utf-8",
    )
    monkeypatch.setattr(bool_style, "REPO_ROOT", tmp_path)
    monkeypatch.setattr(
        bool_style.contract,
        "authored_bool_style_paths",
        lambda: [guide],
    )

    assert bool_style.main() == 0
    assert "Canonical boolean style OK" in capsys.readouterr().out
