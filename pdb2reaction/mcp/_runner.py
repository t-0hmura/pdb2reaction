"""Shared subprocess-runner + summary-parser used by every MCP tool body.

Design intent (see docs/mcp_server.md):
- Each MCP tool body assembles a CLI argv list (`["pdb2reaction", "opt", ...]`)
  and calls `run_subcmd(argv, out_dir)`.
- The runner spawns the CLI, captures stdout / stderr, parses the produced
  `summary.json` if any, and returns a structured `SubcmdResult`.
- Errors (non-zero exit code, missing summary.json, parse error) are surfaced
  as structured fields so the calling agent sees a stable schema, not a
  Python exception traceback.
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
import uuid
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional, Sequence, TypedDict

from pdb2reaction.core.result_commit import RUN_ID_ENV


# Schema version for the MCP tool return envelope. Bump when the field
# set / value types in `SubcmdResultDict` change. 1.0 matches the
# baseline (status / exit_code / out_dir / summary / stderr_tail /
# stdout_tail / hint / argv / schema_version).
MCP_SUBCMD_RESULT_SCHEMA_VERSION = "1.1"

# Allowed values for the `status` field. Documented in docs/mcp_server.md.
MCP_SUBCMD_RESULT_STATUSES = (
    "ok",
    "failed",
    "summary_missing",
    "summary_parse_error",
    "summary_run_mismatch",
)


class SubcmdResultDict(TypedDict, total=False):
    """TypedDict shape returned by every MCP tool body.

    Mirrors :class:`SubcmdResult.to_dict`. ``total=False`` so consumers
    can omit optional fields (out_dir / summary / hint) on early-failure
    paths.
    """

    schema_version: str
    status: str
    exit_code: int
    out_dir: Optional[str]
    summary: dict[str, Any]
    stderr_tail: str
    stdout_tail: str
    hint: Optional[str]
    argv: list[str]
    run_id: str


@dataclass
class SubcmdResult:
    """Structured result of a single pdb2reaction subcmd invocation.

    `status` is one of :data:`MCP_SUBCMD_RESULT_STATUSES`. The envelope
    carries `schema_version = "1.1"` so MCP clients can pin the contract
    and migrate when the structure changes.
    """

    status: str  # see MCP_SUBCMD_RESULT_STATUSES
    exit_code: int
    out_dir: Optional[str] = None
    summary: dict[str, Any] = field(default_factory=dict)
    stderr_tail: str = ""
    stdout_tail: str = ""
    hint: Optional[str] = None
    argv: list[str] = field(default_factory=list)
    run_id: str = ""

    def to_dict(self) -> SubcmdResultDict:
        """Serialise to a plain dict the MCP framework can ship over JSON-RPC."""
        return {
            "schema_version": MCP_SUBCMD_RESULT_SCHEMA_VERSION,
            "status": self.status,
            "exit_code": self.exit_code,
            "out_dir": self.out_dir,
            "summary": self.summary,
            "stderr_tail": self.stderr_tail,
            "stdout_tail": self.stdout_tail,
            "hint": self.hint,
            "argv": self.argv,
            "run_id": self.run_id,
        }


_PRODUCT_ALIASES = {"pdb2reaction", "p2r"}


def _bind_server_argv(argv: Sequence[str]) -> list[str]:
    """Bind product console aliases to this interpreter and imported package."""

    requested = list(argv)
    if requested and Path(requested[0]).name in _PRODUCT_ALIASES:
        return [sys.executable, "-m", "pdb2reaction", *requested[1:]]
    return requested


def _child_env(
    *,
    run_id: str,
    env_overrides: Optional[dict[str, str]] = None,
) -> dict[str, str]:
    """Build a child environment tied to this source tree and invocation."""

    env = os.environ.copy()
    if env_overrides:
        env.update(env_overrides)

    source_root = str(Path(__file__).resolve().parents[2])
    existing = [part for part in env.get("PYTHONPATH", "").split(os.pathsep) if part]
    env["PYTHONPATH"] = os.pathsep.join(
        [source_root, *(part for part in existing if part != source_root)]
    )
    # The runner owns the identity even if an inherited/override environment
    # carried an older invocation's value.
    env[RUN_ID_ENV] = run_id
    return env


def _read_current_summary(
    summary_path: Path,
    *,
    run_id: str,
    expected_primary: bool = False,
) -> tuple[str, dict[str, Any], Optional[str]]:
    """Read a summary only when it belongs to the current invocation."""

    if not summary_path.exists():
        return "summary_missing", {}, None
    try:
        loaded = json.loads(summary_path.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError) as exc:
        return (
            "summary_parse_error",
            {},
            f"{summary_path.name} present but not valid JSON: {exc}",
        )
    if not isinstance(loaded, dict):
        return (
            "summary_parse_error",
            {},
            f"{summary_path.name} must contain a JSON object",
        )
    if loaded.get("run_id") != run_id:
        return (
            "summary_run_mismatch",
            {},
            f"{summary_path.name} does not belong to current run_id {run_id}",
        )
    # Leaf commands publish summary.json first and authoritative result.json
    # last.  If a primary replace failed, the files can be missing/mixed or
    # carry different generations.  ``expected_primary`` is derived from the
    # invoked command, so an unrelated stale result.json cannot invalidate a
    # legitimate one-file aggregate summary.
    primary_path = summary_path.with_name("result.json")
    if expected_primary:
        if not primary_path.exists():
            return (
                "summary_run_mismatch",
                {},
                "current leaf summary is missing authoritative result.json",
            )
        try:
            primary = json.loads(primary_path.read_text(encoding="utf-8"))
        except (json.JSONDecodeError, OSError) as exc:
            return (
                "summary_run_mismatch",
                {},
                f"result.json is not a valid current-run primary: {exc}",
            )
        if (
            not isinstance(primary, dict)
            or primary.get("run_id") != run_id
            or primary != loaded
        ):
            return (
                "summary_run_mismatch",
                {},
                "summary.json and authoritative result.json are different generations",
            )
    return "ok", loaded, None


_SUMMARY_ONLY_COMMANDS = {"all", "path-search"}


def _expects_result_primary(argv: Sequence[str]) -> bool:
    """Return whether a managed summary command publishes a result.json pair."""

    requested = list(argv)
    if not requested or Path(requested[0]).name not in _PRODUCT_ALIASES:
        return False
    if len(requested) < 2:
        return False
    return requested[1] not in _SUMMARY_ONLY_COMMANDS


def _tail(text: str, max_lines: int = 60) -> str:
    """Return at most the last `max_lines` lines of `text`."""
    if not text:
        return ""
    lines = text.rstrip().splitlines()
    if len(lines) <= max_lines:
        return text.rstrip()
    return "...\n" + "\n".join(lines[-max_lines:])


def _extract_hint(stderr: str) -> Optional[str]:
    """Extract the most recent `; recover: <hint>` suffix from stderr, if any.

    Subcommands emit recovery hints in this form so MCP clients can surface
    them to the agent without re-parsing the full error.
    """
    hint = None
    for line in stderr.splitlines():
        marker = "; recover:"
        if marker in line:
            hint = line.split(marker, 1)[1].strip()
    return hint


def run_subcmd(
    argv: Sequence[str],
    *,
    out_dir: Optional[Path] = None,
    timeout: Optional[float] = None,
    env_overrides: Optional[dict[str, str]] = None,
    summary_filename: str = "summary.json",
) -> SubcmdResult:
    """Spawn a pdb2reaction subcmd and collect a structured result.

    Parameters
    ----------
    argv
        Argv list (e.g. ``["pdb2reaction", "opt", "-i", "r.pdb", "-q", "-1"]``).
        The leading executable must already be on PATH inside the MCP server's
        environment.
    out_dir
        Expected output directory (must match the `--out-dir` argv entry).
        Used to locate `summary.json`. If None, the runner does not parse a
        summary and only reports exit code + stderr/stdout.
    timeout
        Subprocess timeout in seconds. None = no timeout.
    env_overrides
        Optional environment variables to set for the subprocess.
    summary_filename
        Override for the summary file (default ``summary.json``).
    """
    run_id = str(uuid.uuid4())
    executed_argv = _bind_server_argv(argv)
    env = _child_env(run_id=run_id, env_overrides=env_overrides)

    try:
        proc = subprocess.run(
            executed_argv,
            capture_output=True,
            text=True,
            timeout=timeout,
            env=env,
        )
    except FileNotFoundError as exc:
        return SubcmdResult(
            status="failed",
            exit_code=127,
            stderr_tail=str(exc),
            hint=(
                "The pdb2reaction CLI is not on PATH. Install pdb2reaction "
                "into the environment that hosts the MCP server."
            ),
            argv=executed_argv,
            run_id=run_id,
        )
    except subprocess.TimeoutExpired as exc:
        return SubcmdResult(
            status="failed",
            exit_code=124,
            stderr_tail=f"TIMEOUT after {timeout}s: {exc}",
            hint=(
                "Increase the `timeout_seconds` tool argument or rerun with a smaller "
                "system / fewer cycles."
            ),
            argv=executed_argv,
            run_id=run_id,
        )

    exit_code = proc.returncode
    stderr_tail = _tail(proc.stderr)
    stdout_tail = _tail(proc.stdout)
    hint = _extract_hint(proc.stderr)

    summary: dict[str, Any] = {}
    summary_status: str = "summary_missing"
    if out_dir is not None:
        summary_path = Path(out_dir) / summary_filename
        summary_status, summary, summary_hint = _read_current_summary(
            summary_path,
            run_id=run_id,
            expected_primary=_expects_result_primary(argv),
        )
        hint = hint or summary_hint

    if exit_code != 0:
        status = "failed"
    elif out_dir is not None and summary_status != "ok":
        status = summary_status
    else:
        status = "ok"

    return SubcmdResult(
        status=status,
        exit_code=exit_code,
        out_dir=str(out_dir) if out_dir else None,
        summary=summary,
        stderr_tail=stderr_tail,
        stdout_tail=stdout_tail,
        hint=hint,
        argv=executed_argv,
        run_id=run_id,
    )
