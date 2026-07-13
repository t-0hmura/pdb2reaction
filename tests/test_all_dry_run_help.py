"""CLI documentation regression tests for ``all --dry-run``."""

from click.testing import CliRunner

from pdb2reaction.workflows.all import cli


def test_all_dry_run_help_describes_temporary_extract_precheck():
    result = CliRunner().invoke(cli, ["--help"])
    assert result.exit_code == 0
    help_text = " ".join(result.output.split())
    assert "runs extraction in a temporary" in help_text
    assert "derived charge and electron parity" in help_text
    assert "no computational stage or persistent output" in help_text
