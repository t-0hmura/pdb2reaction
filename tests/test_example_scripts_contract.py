"""Contract tests for the advertised example-script checker (M65)."""

from __future__ import annotations

import importlib.util
import subprocess
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


checker = _load("check_example_scripts")
contract = checker.contract  # same module object the checker uses


def test_public_shell_examples_exist_with_two_invocations() -> None:
    scripts = contract.PUBLIC_SHELL_EXAMPLES
    assert scripts, "no advertised example scripts declared"
    total = 0
    for path in scripts:
        assert path.exists(), f"missing advertised script: {path}"
        total += len(contract.extract_shell_commands(path))
    assert total == 2  # examples/run.sh: MEP + scan invocations


def test_checker_passes_on_real_repo() -> None:
    assert checker.main() == 0


def test_checker_fails_on_shell_syntax_error(tmp_path, monkeypatch) -> None:
    bad = tmp_path / "run.sh"
    bad.write_text("pdb2reaction all -i R.pdb 'unterminated\n", encoding="utf-8")
    # sanity: bash itself rejects it
    proc = subprocess.run(["bash", "-n", str(bad)], text=True, capture_output=True)
    assert proc.returncode != 0
    monkeypatch.setattr(contract, "PUBLIC_SHELL_EXAMPLES", (bad,))
    assert checker.main() == 1


def test_checker_fails_on_invented_option(tmp_path, monkeypatch) -> None:
    script = tmp_path / "run.sh"
    script.write_text("pdb2reaction all -i R.pdb --invented-option\n", encoding="utf-8")
    monkeypatch.setattr(contract, "PUBLIC_SHELL_EXAMPLES", (script,))
    assert checker.main() == 1
