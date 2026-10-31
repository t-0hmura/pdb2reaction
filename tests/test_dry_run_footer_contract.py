"""Shared completion contract for every dry-run-capable workflow."""

from __future__ import annotations

import re
from pathlib import Path

from click.testing import CliRunner

from pdb2reaction.cli import cli as root_cli
from pdb2reaction.core.utils import (
    DRY_RUN_COMPLETE_MESSAGE,
    emit_dry_run_complete,
)


EXPECTED = "[Dry run] --dry-run completed. Input command is valid."

_SMOKE = Path(__file__).resolve().parent / "smoke"
_R, _P, _TS = (str(_SMOKE / name) for name in ("r.pdb", "p.pdb", "ts.gjf"))
DRY_RUN_COMMANDS = {
    "all": ["all", "-i", _R, _P, "-q", "-1"],
    "dft": ["dft", "-i", _TS],
    "freq": ["freq", "-i", _TS],
    "irc": ["irc", "-i", _TS],
    "opt": ["opt", "-i", _R, "-q", "-1"],
    "path-opt": ["path-opt", "-i", _R, _P, "-q", "-1"],
    "path-search": ["path-search", "-i", _R, "-i", _P, "-q", "-1"],
    "scan": ["scan", "-i", _R, "-q", "-1", "--scan-lists", "[(1,5,1.8)]"],
    "scan2d": [
        "scan2d", "-i", _R, "-q", "-1", "--scan-lists", "[(1,5,1.8,2.2),(1,6,1.7,2.1)]",
    ],
    "scan3d": [
        "scan3d", "-i", _R, "-q", "-1",
        "--scan-lists", "[(1,5,1.8,2.2),(1,6,1.7,2.1),(2,3,1.0,1.2)]",
    ],
    "sp": ["sp", "-i", _R, "-q", "-1"],
    "tsopt": ["tsopt", "-i", _TS],
}
WORKFLOWS_WITH_PLAN_BLOCK = {"dft", "freq", "irc", "opt", "path-opt", "path-search", "tsopt"}
# The block is a titled key/value table; the title also appears in the echoed command line
# when this test is selected by name, so match the underlined title only.
PLAN_BLOCK = re.compile(r"^dry_run_plan\n-+$", re.MULTILINE)


def test_canonical_dry_run_footer_is_exact(capsys) -> None:
    assert DRY_RUN_COMPLETE_MESSAGE == EXPECTED
    emit_dry_run_complete()
    assert capsys.readouterr().out.rstrip().splitlines()[-1] == EXPECTED


def test_every_dry_run_workflow_uses_the_canonical_terminal_footer(tmp_path, monkeypatch) -> None:
    monkeypatch.chdir(tmp_path)
    for name, argv in DRY_RUN_COMMANDS.items():
        result = CliRunner().invoke(root_cli, [*argv, "--dry-run"])
        assert result.exit_code == 0, (name, result.output)
        assert result.output.rstrip().splitlines()[-1] == EXPECTED, name
        assert "[dry-run] Validation complete" not in result.output, name
        assert "Elapsed" not in result.output, name


def test_dry_run_plan_block_is_shown_only_at_verbose_three(tmp_path, monkeypatch) -> None:
    monkeypatch.chdir(tmp_path)
    shown = set()
    for name, argv in DRY_RUN_COMMANDS.items():
        default = CliRunner().invoke(root_cli, [*argv, "--dry-run"])
        full = CliRunner().invoke(root_cli, [*argv, "--dry-run", "-v", "3"])
        assert default.exit_code == 0 and full.exit_code == 0, (name, full.output)
        assert not PLAN_BLOCK.search(default.output), name
        if PLAN_BLOCK.search(full.output):
            shown.add(name)
    assert shown == WORKFLOWS_WITH_PLAN_BLOCK
