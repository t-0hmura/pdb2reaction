"""CLI documentation regression tests for ``all --dry-run``."""

import importlib

from click.testing import CliRunner

from pdb2reaction.cli.app import cli as root_cli
from pdb2reaction.workflows.all import cli


def test_all_dry_run_help_describes_temporary_extract_precheck():
    result = CliRunner().invoke(cli, ["--help"])
    assert result.exit_code == 0
    help_text = " ".join(result.output.split())
    assert "runs extraction in a temporary" in help_text
    assert "derived charge and electron parity" in help_text
    assert "no computational stage or persistent output" in help_text


def test_all_dry_run_plan_omits_unrequested_extraction(tmp_path):
    xyz = tmp_path / "ts.xyz"
    xyz.write_text("2\nH2\nH 0 0 0\nH 0 0 0.74\n", encoding="utf-8")

    result = CliRunner().invoke(
        root_cli, ["all", "-i", str(xyz), "-q", "0", "--tsopt", "--dry-run"]
    )

    assert result.exit_code == 0, result.output
    assert "extraction was not requested" in result.output
    assert "Planned stages: tsopt -> irc." in result.output
    assert "Planned stages: extract" not in result.output


def test_all_scan_list_dry_run_reports_scan_and_path(tmp_path):
    xyz = tmp_path / "reactant.xyz"
    xyz.write_text("2\nH2\nH 0 0 0\nH 0 0 0.74\n", encoding="utf-8")

    result = CliRunner().invoke(
        root_cli,
        [
            "all", "-i", str(xyz), "-q", "0",
            "--scan-lists", "[(1,2,1.0)]",
            "--dry-run",
        ],
    )

    assert result.exit_code == 0, result.output
    assert "Planned stages: scan -> path_opt." in result.output


def test_all_dry_run_cleans_extract_tempdir_on_failure(tmp_path, monkeypatch):
    all_workflow = importlib.import_module("pdb2reaction.workflows.all")
    pdb = tmp_path / "input.pdb"
    pdb.write_text(
        "HETATM    1  C1  SAM A   1       0.000   0.000   0.000  1.00  0.00           C  \nEND\n",
        encoding="utf-8",
    )
    dry_dir = tmp_path / "dry_extract"

    def fake_mkdtemp(*args, **kwargs):
        dry_dir.mkdir()
        return str(dry_dir)

    def fail_extract(**kwargs):
        raise RuntimeError("synthetic extract failure")

    monkeypatch.setattr(all_workflow.tempfile, "mkdtemp", fake_mkdtemp)
    monkeypatch.setattr(all_workflow, "extract_api", fail_extract)

    result = CliRunner().invoke(
        root_cli,
        ["all", "-i", str(pdb), "-c", "SAM", "--tsopt", "--dry-run"],
    )

    assert result.exit_code != 0
    assert "synthetic extract failure" in result.output
    assert not dry_dir.exists()
