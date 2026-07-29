"""Invocation-level restoration tests for ``pdb2reaction all``."""

from __future__ import annotations

import os
import signal
from pathlib import Path

from click.testing import CliRunner
import pytest

from pdb2reaction.workflows import all as all_workflow
from pdb2reaction.workflows._run_session import RunSession
from pdb2reaction.workflows._run_session import declare_public_output
from pdb2reaction.core.result_commit import RUN_ID_ENV


def _write_h2(path: Path) -> None:
    path.write_text("2\nH2\nH 0 0 0\nH 0 0 0.74\n", encoding="utf-8")


def _seed_process_state(monkeypatch):
    from pdb2reaction.core import utils

    sentinel_handler = lambda *_args: None
    monkeypatch.setattr(signal, "getsignal", lambda _sig: sentinel_handler)
    restored_handlers = []
    monkeypatch.setattr(
        signal,
        "signal",
        lambda sig, handler: restored_handlers.append((sig, handler)),
    )
    utils.set_pipeline_mode(True)
    utils.set_child_mode(True)
    utils.set_convert_file_enabled(False)
    all_workflow._echo_state._started = True
    all_workflow._FREEZE_ATOMS_GLOBAL = [9, 3]
    all_workflow._FREEZE_ATOMS_YAML = [7]
    return utils, sentinel_handler, restored_handlers


@pytest.mark.parametrize("failure", [None, RuntimeError("boom"), SystemExit(130)])
@pytest.mark.parametrize("prior_run_id", [None, "prior-mcp-run"])
def test_all_restores_exact_process_state_and_prepared_resources(
    tmp_path: Path, monkeypatch, failure, prior_run_id,
) -> None:
    inputs = [tmp_path / "left.xyz", tmp_path / "right.xyz"]
    for path in inputs:
        _write_h2(path)
    utils, sentinel_handler, restored_handlers = _seed_process_state(monkeypatch)
    if prior_run_id is None:
        monkeypatch.delenv(RUN_ID_ENV, raising=False)
    else:
        monkeypatch.setenv(RUN_ID_ENV, prior_run_id)
    cleanup_order: list[str] = []

    class Prepared:
        def __init__(self, path: Path) -> None:
            self.source_path = path
            self.geom_path = path
            self.structure_template = None
            self.gjf_template = None

        def cleanup(self) -> None:
            cleanup_order.append(self.source_path.name)

    monkeypatch.setattr(
        all_workflow,
        "prepare_input_structure",
        lambda path: Prepared(Path(path)),
    )
    if failure is not None:
        def fail_after_preparation():
            raise failure
        monkeypatch.setattr(all_workflow, "verbose_level", fail_after_preparation)

    result = CliRunner().invoke(
        all_workflow.cli,
        ["-i", str(inputs[0]), "-i", str(inputs[1]), "-q", "0", "--dry-run"],
    )

    if failure is None:
        assert result.exit_code == 0, result.output
    else:
        assert result.exit_code != 0
    assert cleanup_order == ["right.xyz", "left.xyz"]
    assert utils.pipeline_mode_enabled() is True
    assert utils.is_child_mode() is True
    assert utils.is_convert_file_enabled() is False
    assert all_workflow._echo_state._started is True
    assert all_workflow._FREEZE_ATOMS_GLOBAL == [9, 3]
    assert all_workflow._FREEZE_ATOMS_YAML == [7]
    if prior_run_id is None:
        assert RUN_ID_ENV not in os.environ
    else:
        assert os.environ[RUN_ID_ENV] == prior_run_id
    assert restored_handlers[-1] == (signal.SIGINT, sentinel_handler)
    utils.set_pipeline_mode(False)
    utils.set_child_mode(False)
    utils.set_convert_file_enabled(True)
    all_workflow._echo_state._started = False
    all_workflow._FREEZE_ATOMS_GLOBAL = None
    all_workflow._FREEZE_ATOMS_YAML = None


def test_generated_yaml_is_session_owned_and_caller_yaml_is_preserved(
    tmp_path: Path, monkeypatch,
) -> None:
    caller = tmp_path / "caller.yaml"
    caller.write_text("geom:\n  coord_type: cart\n", encoding="utf-8")
    owned = tmp_path / "owned"
    monkeypatch.setattr(
        all_workflow.tempfile,
        "mkdtemp",
        lambda **_kwargs: str(owned.mkdir() or owned),
    )
    session = RunSession()

    generated = all_workflow._write_args_yaml_with_freeze_atoms(
        caller,
        [0],
        session=session,
    )
    assert generated is not None and generated.exists()
    session.close()

    assert not owned.exists()
    assert caller.read_text(encoding="utf-8") == "geom:\n  coord_type: cart\n"


def test_generated_yaml_directory_is_owned_before_serialization(
    tmp_path: Path, monkeypatch,
) -> None:
    owned = tmp_path / "owned_failure"
    monkeypatch.setattr(
        all_workflow.tempfile,
        "mkdtemp",
        lambda **_kwargs: str(owned.mkdir() or owned),
    )
    monkeypatch.setattr(
        all_workflow.yaml,
        "safe_dump",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(RuntimeError("dump failed")),
    )
    session = RunSession()

    with pytest.raises(RuntimeError, match="dump failed"):
        all_workflow._write_args_yaml_with_freeze_atoms(
            None,
            [0],
            session=session,
        )
    session.close()

    assert not owned.exists()


def test_key_outputs_exclude_stale_segments_and_aggregate_diagrams(
    tmp_path: Path,
) -> None:
    root = tmp_path / "out"
    current_segment = root / "segments" / "seg_01" / "structures" / "ts.pdb"
    stale_segment = root / "segments" / "seg_02" / "structures" / "ts.pdb"
    stale_diagram = root / "irc_plot_all.png"
    for path in (current_segment, stale_segment, stale_diagram):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("old", encoding="utf-8")
    manifest = all_workflow.InvocationManifest()
    replacement = root / ".current"
    summary = root / "summary.json"
    replacement.write_text("current", encoding="utf-8")
    declare_public_output(manifest, root, current_segment)
    declare_public_output(manifest, root, summary)
    replacement.replace(current_segment)
    summary.write_text("{}", encoding="utf-8")

    key_outputs = all_workflow._current_key_output_files(
        manifest, root
    )

    assert "seg_01" in key_outputs
    assert "structures/ts.pdb" in key_outputs["seg_01"]["files"]
    assert "seg_02" not in key_outputs
    assert "irc_plot_all.png" not in key_outputs
    assert "summary.json" in key_outputs


class _FakeCommand:
    def __init__(self, error=None) -> None:
        self.error = error

    def main(self, **_kwargs) -> None:
        if self.error is not None:
            raise self.error


@pytest.mark.parametrize("error", [None, RuntimeError("child failure")])
def test_run_cli_main_restores_prior_child_mode(monkeypatch, error) -> None:
    from pdb2reaction.core.utils import is_child_mode, set_child_mode

    monkeypatch.setattr(all_workflow.gc, "collect", lambda: 0)
    monkeypatch.setattr(all_workflow.torch.cuda, "is_available", lambda: False)
    set_child_mode(True)
    try:
        if error is None:
            all_workflow._run_cli_main("fake", _FakeCommand(), [])
        else:
            with pytest.raises(Exception, match="child failure"):
                all_workflow._run_cli_main("fake", _FakeCommand(error), [])
        assert is_child_mode() is True
    finally:
        set_child_mode(False)
