"""Unit tests for the p2r MCP runner envelope."""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

from pdb2reaction.mcp import _runner
from pdb2reaction.mcp._runner import (
    SubcmdResult,
    SubcmdResultDict,
    MCP_SUBCMD_RESULT_SCHEMA_VERSION,
    MCP_SUBCMD_RESULT_STATUSES,
)
from pdb2reaction.core.result_commit import RUN_ID_ENV


def test_subcmd_result_to_dict_carries_schema_version() -> None:
    r = SubcmdResult(status="ok", exit_code=0, argv=["pdb2reaction", "opt"])
    d = r.to_dict()
    assert d["schema_version"] == MCP_SUBCMD_RESULT_SCHEMA_VERSION
    assert d["status"] == "ok"
    assert d["exit_code"] == 0
    assert d["argv"] == ["pdb2reaction", "opt"]
    assert "run_id" in d


def test_subcmd_result_status_enum() -> None:
    assert MCP_SUBCMD_RESULT_STATUSES == (
        "ok",
        "failed",
        "summary_missing",
        "summary_parse_error",
        "summary_run_mismatch",
    )


def test_subcmd_result_dict_keys_match_to_dict() -> None:
    typed_keys = set(SubcmdResultDict.__annotations__.keys())
    r = SubcmdResult(status="ok", exit_code=0)
    runtime_keys = set(r.to_dict().keys())
    assert runtime_keys <= typed_keys


def test_timeout_hint_names_public_tool_argument(monkeypatch) -> None:
    def raise_timeout(*_args, **_kwargs):
        raise subprocess.TimeoutExpired(cmd=["pdb2reaction", "opt"], timeout=12.0)

    monkeypatch.setattr(_runner.subprocess, "run", raise_timeout)
    result = _runner.run_subcmd(
        ["pdb2reaction", "opt"],
        timeout=12.0,
    )

    assert result.status == "failed"
    assert result.exit_code == 124
    assert "timeout_seconds" in (result.hint or "")


def test_stale_summary_is_never_returned(
    tmp_path: Path, monkeypatch
) -> None:
    (tmp_path / "summary.json").write_text(
        json.dumps({"run_id": "old", "secret_old_payload": 1}),
        encoding="utf-8",
    )

    def no_write(argv, **_kwargs):
        return subprocess.CompletedProcess(argv, 0, stdout="", stderr="")

    monkeypatch.setattr(_runner.subprocess, "run", no_write)
    result = _runner.run_subcmd(
        ["pdb2reaction", "opt", "-i", "relative.xyz"],
        out_dir=tmp_path,
    )

    assert result.status == "summary_run_mismatch"
    assert result.summary == {}
    assert result.run_id != "old"
    assert result.argv[:3] == [sys.executable, "-m", "pdb2reaction"]


def test_missing_and_malformed_summaries_return_no_payload(
    tmp_path: Path, monkeypatch
) -> None:
    def no_write(argv, **_kwargs):
        return subprocess.CompletedProcess(argv, 0, stdout="", stderr="")

    monkeypatch.setattr(_runner.subprocess, "run", no_write)
    missing = _runner.run_subcmd(
        ["pdb2reaction", "all", "-i", "relative.xyz"],
        out_dir=tmp_path,
    )
    assert missing.status == "summary_missing"
    assert missing.summary == {}

    (tmp_path / "summary.json").write_text("{broken", encoding="utf-8")
    malformed = _runner.run_subcmd(
        ["pdb2reaction", "all", "-i", "relative.xyz"],
        out_dir=tmp_path,
    )
    assert malformed.status == "summary_parse_error"
    assert malformed.summary == {}


def test_current_summary_uses_bound_interpreter_source_path_and_unchanged_cwd(
    tmp_path: Path, monkeypatch
) -> None:
    captured: dict = {}

    def write_current(argv, **kwargs):
        captured["argv"] = list(argv)
        captured["kwargs"] = kwargs
        current = kwargs["env"][RUN_ID_ENV]
        payload = json.dumps({"run_id": current, "status": "success"})
        (tmp_path / "summary.json").write_text(payload, encoding="utf-8")
        (tmp_path / "result.json").write_text(payload, encoding="utf-8")
        return subprocess.CompletedProcess(argv, 0, stdout="ok", stderr="")

    monkeypatch.setattr(_runner.subprocess, "run", write_current)
    result = _runner.run_subcmd(
        ["pdb2reaction", "opt", "-i", "relative/input.xyz"],
        out_dir=tmp_path,
        env_overrides={"PYTHONPATH": "/external/source"},
    )

    assert result.status == "ok"
    assert result.summary["run_id"] == result.run_id
    assert result.argv == captured["argv"]
    assert result.argv == [
        sys.executable,
        "-m",
        "pdb2reaction",
        "opt",
        "-i",
        "relative/input.xyz",
    ]
    assert "cwd" not in captured["kwargs"]
    source_root = str(Path(_runner.__file__).resolve().parents[2])
    pythonpath = captured["kwargs"]["env"]["PYTHONPATH"].split(os.pathsep)
    assert pythonpath[0] == source_root
    assert "/external/source" in pythonpath[1:]


def test_path_shadow_cannot_replace_imported_cli(
    tmp_path: Path,
) -> None:
    sentinel = tmp_path / "shadow-ran"
    shadow = tmp_path / "pdb2reaction"
    shadow.write_text(
        f"#!/bin/sh\ntouch {sentinel}\nexit 99\n",
        encoding="utf-8",
    )
    shadow.chmod(0o755)

    result = _runner.run_subcmd(
        ["pdb2reaction", "--version"],
        env_overrides={"PATH": f"{tmp_path}{os.pathsep}{os.environ.get('PATH', '')}"},
    )

    assert result.status == "ok"
    assert result.exit_code == 0
    assert result.argv[:3] == [sys.executable, "-m", "pdb2reaction"]
    assert not sentinel.exists()
