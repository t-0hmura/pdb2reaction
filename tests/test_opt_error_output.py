"""Optimization error envelopes follow the effective output directory."""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from click.testing import CliRunner

from pdb2reaction.cli import cli as root_cli


def test_opt_invalid_cartesian_coord_kwargs_exit_nonzero(tmp_path: Path) -> None:
    source = tmp_path / "atom.xyz"
    source.write_text("1\natom\nHe 0 0 0\n", encoding="utf-8")
    config = tmp_path / "config.yaml"
    config.write_text(
        "geom:\n  coord_type: cart\n  coord_kwargs:\n    define_prims: []\n",
        encoding="utf-8",
    )

    result = CliRunner().invoke(
        root_cli,
        [
            "opt", "-i", str(source), "-q", "0", "--dry-run",
            "--config", str(config), "-o", str(tmp_path / "out"),
        ],
    )

    assert result.exit_code == 1
    assert "coord_kwargs were given" in result.output


def test_opt_rejects_zero_cycle_budget(tmp_path: Path) -> None:
    source = tmp_path / "atom.xyz"
    source.write_text("1\natom\nHe 0 0 0\n", encoding="utf-8")

    result = CliRunner().invoke(
        root_cli,
        [
            "opt", "-i", str(source), "-q", "0", "--max-cycles", "0",
            "--dry-run", "-o", str(tmp_path / "out"),
        ],
    )

    assert result.exit_code != 0


def test_opt_rejects_negative_cycle_budget(tmp_path: Path) -> None:
    source = tmp_path / "atom.xyz"
    source.write_text("1\natom\nHe 0 0 0\n", encoding="utf-8")

    result = CliRunner().invoke(
        root_cli,
        [
            "opt", "-i", str(source), "-q", "0", "--max-cycles", "-1",
            "--dry-run", "-o", str(tmp_path / "out"),
        ],
    )

    assert result.exit_code == 2
    assert "Invalid value for '--max-cycles'" in result.output


def test_yaml_output_directory_receives_runtime_error_envelope(tmp_path) -> None:
    source = tmp_path / "atom.xyz"
    source.write_text("1\natom\nHe 0 0 0\n", encoding="utf-8")
    out_dir = tmp_path / "yaml-output"
    config = tmp_path / "config.yaml"
    config.write_text(
        f"opt:\n  out_dir: {out_dir}\n",
        encoding="utf-8",
    )
    calculator = tmp_path / "broken_calculator.py"
    calculator.write_text(
        "def get_calculator(**kwargs):\n"
        "    raise RuntimeError('calculator construction failed')\n",
        encoding="utf-8",
    )

    result = CliRunner().invoke(
        root_cli,
        [
            "opt",
            "-i",
            str(source),
            "-q",
            "0",
            "--config",
            str(config),
            "--calc-file",
            str(calculator),
        ],
    )

    assert result.exit_code == 1
    payload = json.loads((out_dir / "result.json").read_text(encoding="utf-8"))
    assert payload["status"] == "error"
    assert "calculator construction failed" in payload["error"]
