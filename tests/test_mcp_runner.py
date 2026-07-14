"""Unit tests for the p2r MCP runner envelope."""

from __future__ import annotations

import subprocess

from pdb2reaction.mcp import _runner
from pdb2reaction.mcp._runner import (
    SubcmdResult,
    SubcmdResultDict,
    MCP_SUBCMD_RESULT_SCHEMA_VERSION,
    MCP_SUBCMD_RESULT_STATUSES,
)


def test_subcmd_result_to_dict_carries_schema_version() -> None:
    r = SubcmdResult(status="ok", exit_code=0, argv=["pdb2reaction", "opt"])
    d = r.to_dict()
    assert d["schema_version"] == MCP_SUBCMD_RESULT_SCHEMA_VERSION
    assert d["status"] == "ok"
    assert d["exit_code"] == 0
    assert d["argv"] == ["pdb2reaction", "opt"]


def test_subcmd_result_status_enum() -> None:
    assert MCP_SUBCMD_RESULT_STATUSES == (
        "ok",
        "failed",
        "summary_missing",
        "summary_parse_error",
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
