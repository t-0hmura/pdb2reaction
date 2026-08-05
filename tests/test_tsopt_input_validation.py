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


def test_tsopt_dry_run_rejects_invalid_configured_hessian_mode(tmp_path) -> None:
    source = tmp_path / "atom.xyz"
    source.write_text("1\natom\nHe 0 0 0\n", encoding="utf-8")
    config = tmp_path / "config.yaml"
    config.write_text("calc:\n  hessian_calc_mode: typo\n", encoding="utf-8")

    result = CliRunner().invoke(
        root_cli,
        [
            "tsopt", "-i", str(source), "-q", "0", "--dry-run",
            "--config", str(config),
        ],
    )

    assert result.exit_code == 1
    assert "Unsupported hessian_calc_mode 'typo'" in result.output


def test_tsopt_cli_print_every_overrides_nested_dimer_config(
    monkeypatch, tmp_path,
) -> None:
    from pdb2reaction.workflows import tsopt

    stripped = []
    original_strip = tsopt.strip_inherited_keys

    def capture_strip(values, *args, **kwargs):
        stripped.append(dict(values))
        return original_strip(values, *args, **kwargs)

    monkeypatch.setattr(tsopt, "strip_inherited_keys", capture_strip)
    source = tmp_path / "atom.xyz"
    source.write_text("1\natom\nHe 0 0 0\n", encoding="utf-8")
    config = tmp_path / "config.yaml"
    config.write_text(
        "hessian_dimer:\n  lbfgs:\n    print_every: 99\n",
        encoding="utf-8",
    )

    result = CliRunner().invoke(
        root_cli,
        [
            "tsopt", "-i", str(source), "-q", "0", "--opt-mode", "dimer",
            "--print-every", "7", "--show-config", "--dry-run",
            "--config", str(config),
        ],
    )

    assert result.exit_code == 0, result.output
    assert any(values.get("print_every") == 7 for values in stripped)


def test_tsopt_rejects_input_aliasing_final_geometry(tmp_path) -> None:
    out = tmp_path / "out"
    out.mkdir()
    source = out / "final_geometry.xyz"
    original = "1\natom\nHe 0 0 0\n"
    source.write_text(original, encoding="utf-8")

    result = CliRunner().invoke(
        root_cli,
        ["tsopt", "-i", str(source), "-q", "0", "-o", str(out)],
    )

    assert result.exit_code == 2
    assert "collides with a reserved TSOPT output" in result.output
    assert source.read_text(encoding="utf-8") == original
