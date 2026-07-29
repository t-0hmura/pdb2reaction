"""IRC option-combination validation."""

from __future__ import annotations

from click.testing import CliRunner

from pdb2reaction.cli import cli as root_cli


def test_irc_requires_at_least_one_direction(tmp_path) -> None:
    source = tmp_path / "atom.xyz"
    source.write_text("1\natom\nHe 0 0 0\n", encoding="utf-8")

    result = CliRunner().invoke(
        root_cli,
        [
            "irc",
            "-i",
            str(source),
            "-q",
            "0",
            "--no-forward",
            "--no-backward",
            "--dry-run",
            "--out-dir",
            str(tmp_path / "out"),
        ],
    )

    assert result.exit_code != 0
    assert "requires at least one enabled direction" in result.output


def test_irc_rejects_hidden_downhill_mode(tmp_path) -> None:
    source = tmp_path / "atom.xyz"
    source.write_text("1\natom\nHe 0 0 0\n", encoding="utf-8")
    config = tmp_path / "config.yaml"
    config.write_text("irc:\n  downhill: true\n", encoding="utf-8")

    result = CliRunner().invoke(
        root_cli,
        [
            "irc",
            "-i",
            str(source),
            "-q",
            "0",
            "--dry-run",
            "--config",
            str(config),
            "--out-dir",
            str(tmp_path / "out"),
        ],
    )

    assert result.exit_code != 0
    assert "downhill is not supported" in result.output


def test_irc_yaml_prefix_rejects_path_components(tmp_path) -> None:
    source = tmp_path / "atom.xyz"
    source.write_text("1\natom\nHe 0 0 0\n", encoding="utf-8")
    config = tmp_path / "config.yaml"
    config.write_text("irc:\n  prefix: ../escape\n", encoding="utf-8")

    result = CliRunner().invoke(
        root_cli,
        [
            "irc",
            "-i",
            str(source),
            "-q",
            "0",
            "--dry-run",
            "--config",
            str(config),
            "--out-dir",
            str(tmp_path / "out"),
        ],
    )

    assert result.exit_code != 0
    assert "without path components" in result.output
