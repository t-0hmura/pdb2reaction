"""Smoke regressions for previously low-coverage utility subcommands."""

from __future__ import annotations

import json
from pathlib import Path
import re

from click.testing import CliRunner
import pytest

import pdb2reaction.io.energy_diagram as energy_diagram
import pdb2reaction.io.trj2fig as trj2fig
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


def test_fix_altloc_selects_one_coherent_residue_conformer(tmp_path: Path) -> None:
    """Do not merge atoms unique to A and B into a hybrid residue."""
    in_pdb = _write_text(
        tmp_path / "altloc_conformers.pdb",
        (
            "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 10.00           N\n"
            "ATOM      2  CA AALA A   1       1.000   0.000   0.000  0.60 10.00           C\n"
            "ATOM      3  CB AALA A   1       1.500   1.000   0.000  0.60 10.00           C\n"
            "ATOM      4  CG AALA A   1       2.500   1.000   0.000  0.60 10.00           C\n"
            "ATOM      5  CA BALA A   1       1.100   0.000   0.000  0.40 10.00           C\n"
            "ATOM      6  CB BALA A   1       1.600  -1.000   0.000  0.40 10.00           C\n"
            "ATOM      7  CD BALA A   1       2.600  -1.000   0.000  0.40 10.00           C\n"
            "END\n"
        ),
    )
    out_pdb = tmp_path / "altloc_conformers_clean.pdb"

    result = CliRunner().invoke(
        root_cli,
        ["fix-altloc", "-i", str(in_pdb), "-o", str(out_pdb)],
        catch_exceptions=False,
    )

    assert result.exit_code == 0, result.output
    atom_lines = [
        line
        for line in out_pdb.read_text(encoding="utf-8").splitlines()
        if line.startswith("ATOM")
    ]
    assert [line[12:16].strip() for line in atom_lines] == ["N", "CA", "CB", "CG"]
    assert all(line[16] == " " for line in atom_lines)


def test_fix_altloc_uses_residue_level_mean_occupancy(tmp_path: Path) -> None:
    """Per-atom occupancy conflicts must not mix conformer labels."""
    in_pdb = _write_text(
        tmp_path / "altloc_mean.pdb",
        (
            "ATOM      1  CA AALA A   1       1.000   0.000   0.000  0.70 10.00           C\n"
            "ATOM      2  CB AALA A   1       1.500   1.000   0.000  0.50 10.00           C\n"
            "ATOM      3  CA BALA A   1       1.100   0.000   0.000  0.30 10.00           C\n"
            "ATOM      4  CB BALA A   1       1.600  -1.000   0.000  0.60 10.00           C\n"
            "END\n"
        ),
    )
    out_pdb = tmp_path / "altloc_mean_clean.pdb"

    result = CliRunner().invoke(
        root_cli,
        ["fix-altloc", "-i", str(in_pdb), "-o", str(out_pdb)],
        catch_exceptions=False,
    )

    assert result.exit_code == 0, result.output
    serials = [
        int(line[6:11])
        for line in out_pdb.read_text(encoding="utf-8").splitlines()
        if line.startswith("ATOM")
    ]
    assert serials == [1, 2]


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
        [
            "trj2fig", "-i", str(trj), "-o", str(out_csv),
            "--unit", "kcal", "-r", "init", "--out-json",
        ],
        catch_exceptions=False,
    )

    assert result.exit_code == 0, result.output
    assert out_csv.exists()

    lines = out_csv.read_text(encoding="utf-8").splitlines()
    assert lines[0] == "frame,energy_hartree,delta_kcal"
    assert len(lines) == 3
    payload = json.loads((tmp_path / "result.json").read_text(encoding="utf-8"))
    assert payload["energy_source"] == "trajectory_comment"
    assert payload["backend"] is None
    assert payload["charge"] is None
    assert payload["solvent"] is None


def test_trj2fig_recompute_json_records_actual_provenance(
    tmp_path: Path, monkeypatch
) -> None:
    trj = _write_text(
        tmp_path / "traj.xyz",
        "1\nno-energy-needed\nH 0.0 0.0 0.0\n",
    )
    out_csv = tmp_path / "energy.csv"
    seen: dict = {}

    def _fake_recompute(path, charge, multiplicity, **kwargs):
        seen.update(
            path=path, charge=charge, multiplicity=multiplicity, **kwargs
        )
        return [-0.5]

    monkeypatch.setattr(trj2fig, "recompute_energies", _fake_recompute)
    result = CliRunner().invoke(
        root_cli,
        [
            "trj2fig", "-i", str(trj), "-o", str(out_csv),
            "-q", "-1", "-m", "2", "-b", "orb",
            "--solvent", "water", "--out-json",
        ],
        catch_exceptions=False,
    )
    assert result.exit_code == 0, result.output
    payload = json.loads((tmp_path / "result.json").read_text(encoding="utf-8"))
    assert payload["energy_source"] == "mlip_recomputed"
    assert payload["backend"] == "orb"
    assert payload["charge"] == -1
    assert payload["multiplicity"] == 2
    assert payload["solvent"] == "water"
    assert seen["backend"] == "orb"
    assert seen["solvent"] == "water"


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


@pytest.mark.parametrize(
    "args",
    [
        ["energy-diagram", "-i", "0", "12.5", "4.3"],
        [
            "energy-diagram", "-i", "0", "-i", "12.5", "-i", "4.3",
            "--label-x", "R", "TS", "P",
        ],
    ],
)
def test_energy_diagram_rejects_unflagged_variadic_tokens(args: list[str]) -> None:
    result = CliRunner().invoke(root_cli, args)
    assert result.exit_code != 0
    assert "unexpected extra arguments" in result.output.lower()


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
    # Click 8.2 and 8.3 punctuate this diagnostic differently.  The contract
    # is that the unpublished long name is rejected, not the exact wording.
    assert "No such option" in old_name.output
    assert "--reference-mode" in old_name.output


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
    assert "DMF" in option_help
    assert "max_nodes+2" in option_help


def test_opt_accepts_yaml_only_charge_and_spin(tmp_path: Path) -> None:
    config = tmp_path / "calc.yaml"
    config.write_text("calc:\n  charge: -1\n  spin: 1\n", encoding="utf-8")
    input_path = Path(__file__).parent / "smoke" / "r.pdb"

    result = CliRunner().invoke(
        root_cli,
        ["opt", "-i", str(input_path), "--config", str(config), "--dry-run"],
    )

    assert result.exit_code == 0, result.output
    assert "Validation complete" in result.output


@pytest.mark.parametrize(
    ("extra", "expected"),
    [([], "legacy-active"), (["--tr-projection", "constrained"], "constrained")],
)
def test_opt_tr_projection_cli_overrides_yaml(
    tmp_path: Path, extra: list[str], expected: str,
) -> None:
    config = tmp_path / "projection.yaml"
    config.write_text(
        "calc:\n  charge: -1\n  spin: 1\n"
        "geom:\n  tr_projection: legacy-active\n",
        encoding="utf-8",
    )
    input_path = Path(__file__).parent / "smoke" / "r.pdb"
    result = CliRunner().invoke(
        root_cli,
        [
            "opt", "-i", str(input_path), "--config", str(config),
            "--dry-run", "-v", "3", *extra,
        ],
    )

    assert result.exit_code == 0, result.output
    assert f"tr_projection: {expected}" in result.output


@pytest.mark.parametrize(
    ("toggle", "expected"),
    [([], None), (["--flatten"], True), (["--no-flatten"], False)],
)
def test_all_only_forwards_explicit_flatten(
    tmp_path: Path, toggle: list[str], expected: bool | None,
) -> None:
    config = tmp_path / "all.yaml"
    config.write_text(
        "calc:\n  charge: -1\n  spin: 1\n"
        "hessian_dimer:\n  flatten_max_iter: 7\n",
        encoding="utf-8",
    )
    input_path = Path(__file__).parent / "smoke" / "r.pdb"
    result = CliRunner().invoke(
        root_cli,
        [
            "all", "-i", str(input_path), "-i", str(input_path),
            "--config", str(config), "--show-config", "--dry-run", *toggle,
        ],
    )

    assert result.exit_code == 0, result.output
    marker = "overrides:\n  tsopt:"
    assert marker in result.output
    tsopt_block = result.output.split(marker, 1)[1].split("  freq:", 1)[0]
    if expected is None:
        assert "flatten:" not in tsopt_block
    else:
        assert f"flatten: {str(expected).lower()}" in tsopt_block


@pytest.mark.parametrize(
    ("extra", "expected_cycles"),
    [([], 123), (["--preopt-max-cycles", "456"], 456)],
)
def test_path_opt_explicit_preopt_cycles_override_yaml(
    tmp_path: Path, extra: list[str], expected_cycles: int,
) -> None:
    config = tmp_path / "path-opt.yaml"
    config.write_text(
        "calc:\n  charge: -1\n  spin: 1\n"
        "lbfgs:\n  max_cycles: 123\n",
        encoding="utf-8",
    )
    input_path = Path(__file__).parent / "smoke" / "r.pdb"
    result = CliRunner().invoke(
        root_cli,
        [
            "path-opt", "-i", str(input_path), str(input_path),
            "--config", str(config), "--show-config", "--dry-run",
            "-v", "3", *extra,
        ],
    )

    assert result.exit_code == 0, result.output
    assert f"preopt_max_cycles: {expected_cycles}" in result.output


def test_scan_passes_merged_yaml_mapping_to_runtime(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Guard against unpacking the loader's mapping instead of passing it on."""
    import sys
    import pdb2reaction.workflows.scan as scan_workflow

    config = tmp_path / "scan.yaml"
    config.write_text(
        "calc:\n  charge: -1\n  spin: 1\n"
        "bias:\n  k: 17.5\n",
        encoding="utf-8",
    )
    input_path = Path(__file__).parent / "smoke" / "r.pdb"
    seen: dict[str, object] = {}

    def _capture_config(yaml_cfg, **_kwargs):
        seen["yaml_cfg"] = yaml_cfg
        raise RuntimeError("stop-after-config-capture")

    monkeypatch.setattr(scan_workflow, "build_scan_configs", _capture_config)
    args = [
        "scan", "-i", str(input_path),
        "--config", str(config),
        "--scan-lists", "[(1, 2, 1.0)]",
        "--out-dir", str(tmp_path / "out"),
    ]
    # scan deliberately re-reads raw argv to preserve repeated -s groups.
    monkeypatch.setattr(sys, "argv", ["pdb2reaction", *args])
    result = CliRunner().invoke(root_cli, args)

    assert result.exit_code != 0
    assert isinstance(seen.get("yaml_cfg"), dict)
    assert seen["yaml_cfg"]["bias"]["k"] == 17.5


@pytest.mark.parametrize(
    ("command", "scan_spec"),
    [
        ("scan2d", "[(1,2,0.7,1.0),(2,3,0.7,1.0)]"),
        ("scan3d", "[(1,2,0.7,1.0),(2,3,0.7,1.0),(1,3,0.7,1.0)]"),
    ],
)
def test_grid_scans_accept_yaml_only_charge_and_spin(
    tmp_path: Path, command: str, scan_spec: str,
) -> None:
    config = tmp_path / f"{command}.yaml"
    config.write_text("calc:\n  charge: -1\n  spin: 1\n", encoding="utf-8")
    input_path = Path(__file__).parent / "smoke" / "r.pdb"
    result = CliRunner().invoke(
        root_cli,
        [
            command, "-i", str(input_path),
            "--scan-lists", scan_spec,
            "--config", str(config), "--dry-run",
            "--out-dir", str(tmp_path / "out"),
        ],
    )

    assert result.exit_code == 0, result.output
    assert "resolved charge : -1" in result.output
