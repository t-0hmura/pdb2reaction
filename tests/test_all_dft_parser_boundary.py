"""Exercise all's actual DFT argv against the real child parser, without SCF."""
import subprocess
from types import SimpleNamespace

import pytest
from click.testing import CliRunner

from pdb2reaction.workflows import all as workflow, dft


@pytest.mark.parametrize("convert_files", [True, False])
@pytest.mark.parametrize("invalid_option", [False, True])
def test_parent_dft_argv_and_parser_error_logs(tmp_path, monkeypatch, convert_files, invalid_option):
    geometry = tmp_path / "h2.xyz"
    geometry.write_text("2\nH2\nH 0 0 0\nH 0 0 0.74\n")
    calls, echoes = [], []

    def run(cmd, **kwargs):
        calls.append(cmd)
        args = cmd[4:] + ["--dry-run"]
        if invalid_option:
            args += ["--not-a-dft-option"]
        result = CliRunner().invoke(dft.cli, args)
        return SimpleNamespace(returncode=result.exit_code,
                               stdout=result.output if not result.exit_code else "",
                               stderr=result.output if result.exit_code else "")

    monkeypatch.setattr(subprocess, "run", run)
    monkeypatch.setattr(workflow, "_echo", lambda text, **kw: echoes.append((text, kw)))
    result = workflow._run_dft_for_state(
        geometry, 0, 1, tmp_path / "dft", None,
        func_basis="hf/sto-3g", engine="cpu", convert_files=convert_files,
        overrides={"max_cycle":40,"grid_level":0,"conv_tol":1e-5},
    )
    assert "--convert-files" not in calls[0] and "--no-convert-files" not in calls[0]
    assert result["_dft_returncode"] == (2 if invalid_option else 0)
    log = (tmp_path / "dft/run.log").read_text()
    if invalid_option:
        assert "No such option" in log and "--not-a-dft-option" in log
        assert "No such option" in result["_dft_stderr"]
        assert any("No such option" in text and kw.get("err") and kw.get("narrative") for text, kw in echoes)
        assert result["_dft_failed"] is True
    else:
        assert "dry-run" in log.lower()
