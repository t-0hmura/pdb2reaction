"""Regression tests for the documentation command validator."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest


def _load_smoke_module():
    path = Path(__file__).parents[1] / ".github" / "scripts" / "smoke_docs_commands.py"
    spec = importlib.util.spec_from_file_location("smoke_docs_commands", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_placeholder_command_options_are_validated() -> None:
    smoke = _load_smoke_module()
    smoke._validate_option_names(
        ["pdb2reaction opt -i <input.pdb> -q <charge> --dry-run"]
    )


def test_invalid_standalone_option_is_rejected() -> None:
    smoke = _load_smoke_module()
    with pytest.raises(RuntimeError, match="--definitely-not-an-option"):
        smoke._validate_option_names(
            ["pdb2reaction opt -i <input.pdb> --definitely-not-an-option"]
        )


def test_custom_xyz_examples_are_fenced_and_supply_charge_spin() -> None:
    smoke = _load_smoke_module()
    commands = smoke._extract_docs_commands()
    custom = [cmd for cmd in commands if "--calc-file my_calc.py" in cmd]

    assert custom
    for command in custom:
        assert "-q 0" in command
        assert "-m 1" in command


def test_dynamic_smoke_excludes_usage_synopses() -> None:
    smoke = _load_smoke_module()

    selected = smoke._select_executable_all_commands(
        [
            "pdb2reaction all [OPTIONS]...",
            "pdb2reaction all -i R.pdb P.pdb -c LIG --dry-run",
        ]
    )

    assert selected == ["pdb2reaction all -i R.pdb P.pdb -c LIG --dry-run"]
