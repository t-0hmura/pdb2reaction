"""TS-optimization configuration validation."""

from __future__ import annotations

from click.testing import CliRunner

from pdb2reaction.cli import cli as root_cli


def test_dimer_hessian_update_interval_must_advance(tmp_path) -> None:
    source = tmp_path / "atom.xyz"
    source.write_text("1\natom\nHe 0 0 0\n", encoding="utf-8")
    config = tmp_path / "config.yaml"
    config.write_text(
        "hessian_dimer:\n  update_interval_hessian: 0\n",
        encoding="utf-8",
    )

    result = CliRunner().invoke(
        root_cli,
        [
            "tsopt",
            "-i",
            str(source),
            "-q",
            "0",
            "--opt-mode",
            "dimer",
            "--dry-run",
            "--config",
            str(config),
        ],
    )

    assert result.exit_code != 0
    assert "update_interval_hessian must be a positive integer" in result.output
