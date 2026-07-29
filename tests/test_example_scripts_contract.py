"""Contract tests for the advertised example-script checker."""

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
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


checker = _load("check_example_scripts")
contract = checker.contract  # same module object the checker uses


def _write_valid_script(path: Path, command: str) -> None:
    path.write_text(
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n"
        'SCRIPT_DIR=\"$(cd \"$(dirname \"${BASH_SOURCE[0]}\")\" && pwd)\"\n'
        f"{command}\n",
        encoding="utf-8",
    )
    path.chmod(0o755)


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
    _write_valid_script(
        bad, "pdb2reaction all -i \"$SCRIPT_DIR/R.pdb\" 'unterminated",
    )
    # sanity: bash itself rejects it
    proc = subprocess.run(["bash", "-n", str(bad)], text=True, capture_output=True)
    assert proc.returncode != 0
    monkeypatch.setattr(contract, "PUBLIC_SHELL_EXAMPLES", (bad,))
    assert checker.main() == 1


def test_checker_fails_on_invented_option(tmp_path, monkeypatch) -> None:
    script = tmp_path / "run.sh"
    _write_valid_script(
        script,
        'pdb2reaction all -i "$SCRIPT_DIR/R.pdb" --invented-option',
    )
    monkeypatch.setattr(contract, "PUBLIC_SHELL_EXAMPLES", (script,))
    assert checker.main() == 1
