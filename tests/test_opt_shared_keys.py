"""Precedence between the shared ``opt`` block and the lbfgs/rfo sections."""

from __future__ import annotations

from pathlib import Path

import pytest
from click.testing import CliRunner

from pdb2reaction.workflows import opt as opt_workflow


def _dry_run(tmp_path: Path, monkeypatch, yaml_text: str, *extra: str):
    xyz = tmp_path / "h2.xyz"
    xyz.write_text("2\nH2\nH 0 0 0\nH 0 0 0.74\n")
    config = tmp_path / "config.yaml"
    config.write_text(yaml_text)
    blocks = {}

    def capture(title, content, **_kwargs):
        blocks[title] = dict(content)
        return ""

    monkeypatch.setattr(opt_workflow, "pretty_block", capture)
    result = CliRunner().invoke(
        opt_workflow.cli,
        [
            "-i", str(xyz), "-q", "0", "--config", str(config), "--dry-run",
            "-o", str(tmp_path / "out"), *extra,
        ],
    )
    return result, blocks


def test_optimizer_section_values_apply_when_opt_keeps_defaults(tmp_path, monkeypatch):
    result, blocks = _dry_run(
        tmp_path,
        monkeypatch,
        "lbfgs:\n  max_cycles: 50\n  thresh: gau_tight\n  print_every: 7\n",
        "--opt-mode", "grad",
    )

    assert result.exit_code == 0, result.output
    assert blocks["opt"]["max_cycles"] == 50
    assert blocks["opt"]["thresh"] == "gau_tight"
    assert blocks["opt"]["print_every"] == 7


def test_unused_optimizer_section_does_not_leak(tmp_path, monkeypatch):
    result, blocks = _dry_run(
        tmp_path, monkeypatch, "rfo:\n  max_cycles: 50\n", "--opt-mode", "grad",
    )

    assert result.exit_code == 0, result.output
    assert blocks["opt"]["max_cycles"] == 100000


@pytest.mark.parametrize(
    ("yaml_text", "extra"),
    [
        ("opt:\n  max_cycles: 10\nrfo:\n  max_cycles: 50\n", ()),
        ("rfo:\n  max_cycles: 50\n", ("--max-cycles", "10")),
    ],
)
def test_explicit_conflict_is_an_error(tmp_path, monkeypatch, yaml_text, extra):
    result, _ = _dry_run(
        tmp_path, monkeypatch, yaml_text, "--opt-mode", "hess", *extra,
    )

    assert result.exit_code == 2, result.output
    assert "opt.max_cycles and rfo.max_cycles conflict" in result.output
