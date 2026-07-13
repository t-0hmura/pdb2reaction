"""Smoke regressions for previously low-coverage utility subcommands."""

from __future__ import annotations

from pathlib import Path
import re

from click.testing import CliRunner
import pytest

import pdb2reaction.io.energy_diagram as energy_diagram
from pdb2reaction.cli import cli as root_cli


def _write_text(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


def test_add_elem_info_smoke(tmp_path: Path) -> None:
    in_pdb = _write_text(
        tmp_path / "input_no_elem.pdb",
        (
            "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 10.00\n"
            "END\n"
        ),
    )
    out_pdb = tmp_path / "output_add_elem.pdb"

    runner = CliRunner()
    result = runner.invoke(
        root_cli,
        ["add-elem-info", "-i", str(in_pdb), "-o", str(out_pdb)],
        catch_exceptions=False,
    )

    assert result.exit_code == 0, result.output
    assert out_pdb.exists()

    atom_line = next(line for line in out_pdb.read_text(encoding="utf-8").splitlines() if line.startswith("ATOM"))
    assert len(atom_line) >= 78
    assert atom_line[76:78].strip() == "C"


def test_fix_altloc_smoke(tmp_path: Path) -> None:
    in_pdb = _write_text(
        tmp_path / "altloc_input.pdb",
        (
            "ATOM      1  CA AALA A   1       0.000   0.000   0.000  0.60 10.00           C\n"
            "ATOM      2  CA BALA A   1       0.100   0.000   0.000  0.40 10.00           C\n"
            "END\n"
        ),
    )
    out_pdb = tmp_path / "altloc_clean.pdb"

    runner = CliRunner()
    result = runner.invoke(
        root_cli,
        ["fix-altloc", "-i", str(in_pdb), "-o", str(out_pdb)],
        catch_exceptions=False,
    )

    assert result.exit_code == 0, result.output
    assert out_pdb.exists()

    atom_lines = [line for line in out_pdb.read_text(encoding="utf-8").splitlines() if line.startswith("ATOM")]
    assert len(atom_lines) == 1
    assert atom_lines[0][16] == " "


def test_trj2fig_csv_smoke(tmp_path: Path) -> None:
    trj = _write_text(
        tmp_path / "traj.xyz",
        (
            "1\n"
            "0.000000\n"
            "H 0.0 0.0 0.0\n"
            "1\n"
            "0.500000\n"
            "H 0.0 0.0 0.1\n"
        ),
    )
    out_csv = tmp_path / "energy.csv"

    runner = CliRunner()
    result = runner.invoke(
        root_cli,
        ["trj2fig", "-i", str(trj), "-o", str(out_csv), "--unit", "kcal", "-r", "init"],
        catch_exceptions=False,
    )

    assert result.exit_code == 0, result.output
    assert out_csv.exists()

    lines = out_csv.read_text(encoding="utf-8").splitlines()
    assert lines[0] == "frame,energy_hartree,delta_kcal"
    assert len(lines) == 3


def test_energy_diagram_smoke(tmp_path: Path, monkeypatch) -> None:
    out_png = tmp_path / "energy_diagram.png"

    class _DummyFigure:
        def write_image(self, out_path: str, scale: int = 2) -> None:
            _ = scale
            Path(out_path).write_text("dummy-image", encoding="utf-8")

    monkeypatch.setattr(
        energy_diagram,
        "build_energy_diagram",
        lambda **_kwargs: _DummyFigure(),
    )

    runner = CliRunner()
    result = runner.invoke(
        root_cli,
        [
            "energy-diagram",
            "-i",
            "0",
            "-i",
            "1.2",
            "-o",
            str(out_png),
            "--label-x",
            "R",
            "--label-x",
            "TS",
        ],
        catch_exceptions=False,
    )

    assert result.exit_code == 0, result.output
    assert out_png.exists()
    assert out_png.read_text(encoding="utf-8") == "dummy-image"


def test_verbose_is_a_per_subcommand_option() -> None:
    """`-v/--verbose LEVEL` is injected into every subcommand; it is no longer a
    root-group option, so a root-placed `-v` is rejected."""
    runner = CliRunner()
    for name in (
        "opt", "tsopt", "freq", "irc", "sp", "scan",
        "path-search", "dft", "all", "extract",
    ):
        res = runner.invoke(root_cli, [name, "--help"])
        assert res.exit_code == 0, res.output
        assert "-v, --verbose" in res.output, f"{name} --help is missing -v"

    # Root-placed `-v` no longer exists (it moved onto the subcommands).
    root = runner.invoke(root_cli, ["-v", "2", "opt", "--help"])
    assert root.exit_code != 0
    assert "No such option" in root.output

    # The level is an IntRange(0, 3); 0/1/2/3 are accepted (default 2) and
    # out-of-range values are rejected.
    for cmd in (["opt", "-v", "4"], ["extract", "-v", "4"]):
        bad = runner.invoke(root_cli, cmd)
        assert bad.exit_code != 0, cmd
        assert "is not in the range" in bad.output, cmd


def test_tsopt_ref_mode_is_documented_as_advanced_all_handoff() -> None:
    runner = CliRunner()
    primary_help = runner.invoke(root_cli, ["tsopt", "--help"])
    advanced_help = runner.invoke(root_cli, ["tsopt", "--help-advanced"])
    old_name = runner.invoke(
        root_cli,
        ["tsopt", "--reference-mode", "unused.npy"],
    )

    assert primary_help.exit_code == 0, primary_help.output
    assert "--ref-mode" not in primary_help.output

    assert advanced_help.exit_code == 0, advanced_help.output
    assert "--ref-mode" in advanced_help.output
    assert "--reference-mode" not in advanced_help.output
    assert "all workflow supplies" in advanced_help.output
    assert "standalone tsopt" in advanced_help.output

    assert old_name.exit_code != 0
    assert "No such option: --reference-mode" in old_name.output


@pytest.mark.parametrize(
    ("command", "help_flag"),
    [("path-opt", "--help"), ("path-search", "--help"), ("all", "--help-advanced")],
)
def test_path_workflow_max_nodes_defaults_to_twenty(
    command: str, help_flag: str
) -> None:
    result = CliRunner().invoke(root_cli, [command, help_flag])

    assert result.exit_code == 0, result.output
    option_start = result.output.index("--max-nodes")
    option_help = result.output[option_start : option_start + 500]
    assert re.search(r"\[default:\s*20\]", option_help)
