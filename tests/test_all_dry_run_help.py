"""Regression tests for ``all --dry-run``: help text, stage planning and
temp-dir cleanup."""

import importlib

from click.testing import CliRunner

from pdb2reaction.cli.app import cli as root_cli
from pdb2reaction.workflows.all import _charge_override_message, cli


def test_all_dry_run_help_describes_temporary_extract_precheck():
    result = CliRunner().invoke(cli, ["--help"])
    assert result.exit_code == 0
    help_text = " ".join(result.output.split())
    assert "runs extraction in a temporary" in help_text
    assert "derived charge and electron parity" in help_text
    assert "no computational stage or persistent output" in help_text


def test_all_help_describes_charge_precedence_with_and_without_extraction():
    result = CliRunner().invoke(cli, ["--help"])
    assert result.exit_code == 0
    help_text = " ".join(result.output.split())
    assert "Without extraction, -q/--charge explicitly sets the total" in help_text
    assert "with extraction, it asserts the derived total" in help_text


def test_all_no_extract_charge_override_reports_the_actual_relation():
    differs = _charge_override_message("-q/--charge", 0, 1)
    matches = _charge_override_message("-q/--charge", 1, 1)

    assert "WARNING: -q/--charge supplied" in differs
    assert "(overrides workflow-derived +1)" in differs
    assert "(matches workflow-derived +1)" not in differs
    assert "(matches workflow-derived +1)" in matches


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


def test_all_xyz_ref_pdb_derives_residue_charge_before_dry_run(tmp_path):
    xyz = tmp_path / "ts.xyz"
    ref = tmp_path / "ts_ref.pdb"
    xyz.write_text(
        "2\nTS\nC 0.000 0.000 0.000\nO 1.200 0.000 0.000\n",
        encoding="utf-8",
    )
    ref.write_text(
        "HETATM    1  C1  LIG A   1       0.000   0.000   0.000  1.00  0.00           C  \n"
        "HETATM    2  O1  LIG A   1       1.200   0.000   0.000  1.00  0.00           O  \n"
        "END\n",
        encoding="utf-8",
    )

    result = CliRunner().invoke(
        root_cli,
        [
            "all", "-i", str(xyz), "--ref-pdb", str(ref),
            "-l", "LIG:0", "--tsopt", "--dry-run",
        ],
    )

    assert result.exit_code == 0, result.output
    assert "charge=+0, spin(multiplicity)=1" in result.output

    bad_xyz = tmp_path / "bad.xyz"
    bad_xyz.write_text("1\nbad\nC 0 0 0\n", encoding="utf-8")
    mismatch = CliRunner().invoke(
        root_cli,
        [
            "all", "-i", str(bad_xyz), "--ref-pdb", str(ref),
            "-l", "LIG:0", "--tsopt", "--dry-run",
        ],
    )
    assert mismatch.exit_code != 0
    assert "atom count" in mismatch.output.lower()


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


def test_all_rejects_repeated_scan_stage_flags(tmp_path):
    xyz = tmp_path / "reactant.xyz"
    xyz.write_text("2\nH2\nH 0 0 0\nH 0 0 0.74\n", encoding="utf-8")

    result = CliRunner().invoke(
        root_cli,
        [
            "all", "-i", str(xyz), "-q", "0",
            "--scan-lists", "[(1,2,1.0)]",
            "--scan-lists", "[(1,2,1.2)]",
            "--dry-run",
        ],
    )

    assert result.exit_code != 0
    assert "repeated flags are not accepted" in result.output


def test_all_rejects_scan_lists_with_multiple_inputs(tmp_path):
    reactant = tmp_path / "reactant.xyz"
    product = tmp_path / "product.xyz"
    reactant.write_text("2\nR\nH 0 0 0\nH 0 0 0.74\n", encoding="utf-8")
    product.write_text("2\nP\nH 0 0 0\nH 0 0 0.80\n", encoding="utf-8")

    result = CliRunner().invoke(
        root_cli,
        [
            "all",
            "-i",
            str(reactant),
            str(product),
            "-q",
            "0",
            "--scan-lists",
            "[(1,2,1.0)]",
            "--dry-run",
        ],
    )

    assert result.exit_code != 0
    assert "--scan-lists requires exactly one input structure" in result.output


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


def test_all_dry_run_cleans_session_owned_extract_tempdir_on_success(
    tmp_path, monkeypatch,
):
    all_workflow = importlib.import_module("pdb2reaction.workflows.all")
    pdb = tmp_path / "input.pdb"
    pdb.write_text(
        "HETATM    1  H1  SAM A   1       0.000   0.000   0.000  1.00  0.00           H  \n"
        "HETATM    2  H2  SAM A   1       0.000   0.000   0.740  1.00  0.00           H  \n"
        "END\n",
        encoding="utf-8",
    )
    dry_dir = tmp_path / "dry_extract"

    def fake_mkdtemp(*args, **kwargs):
        dry_dir.mkdir()
        return str(dry_dir)

    def successful_extract(**kwargs):
        output = kwargs["output"][0]
        all_workflow.shutil.copy2(pdb, output)
        return {"charge_summary": {"total_charge": 0.0}}

    monkeypatch.setattr(all_workflow.tempfile, "mkdtemp", fake_mkdtemp)
    monkeypatch.setattr(all_workflow, "extract_api", successful_extract)

    result = CliRunner().invoke(
        root_cli,
        ["all", "-i", str(pdb), "-c", "SAM", "-q", "0", "--tsopt", "--dry-run"],
    )

    assert result.exit_code == 0, result.output
    assert not dry_dir.exists()


def test_all_dry_run_does_not_report_prior_summary(tmp_path):
    xyz = tmp_path / "ts.xyz"
    xyz.write_text("2\nH2\nH 0 0 0\nH 0 0 0.74\n", encoding="utf-8")
    out_dir = tmp_path / "existing"
    out_dir.mkdir()
    (out_dir / "summary.json").write_text(
        '{"status": "STALE_SENTINEL", "n_segments_reactive": 99}',
        encoding="utf-8",
    )

    result = CliRunner().invoke(
        root_cli,
        [
            "all",
            "-i",
            str(xyz),
            "-q",
            "0",
            "--tsopt",
            "--dry-run",
            "--out-dir",
            str(out_dir),
        ],
    )

    assert result.exit_code == 0, result.output
    assert "STALE_SENTINEL" not in result.output
    assert "Reactive segments: 99" not in result.output
