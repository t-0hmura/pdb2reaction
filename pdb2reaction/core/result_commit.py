"""Exact-path atomic commits for machine-readable result artifacts.

The public ``result.json``/``summary.json`` pair cannot be replaced as one
filesystem transaction.  This module therefore stages every destination,
publishes compatibility mirrors first, and publishes the authoritative
primary path last.  A successful return means every destination contains the
same serialized bytes.  If primary publication fails after a mirror was
published, consumers must use the current-run identity to reject the mixed
generation.
"""

from __future__ import annotations

import json
import os
import tempfile
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any


RUN_ID_ENV = "PDB2REACTION_RUN_ID"


class ResultCommitError(OSError):
    """Raised when serialization, staging, or publication cannot complete."""

    def __init__(self, phase: str, path: Path, cause: BaseException) -> None:
        self.phase = str(phase)
        self.path = Path(path)
        self.cause = cause
        super().__init__(f"result commit failed during {self.phase} for {self.path}: {cause}")


class RunIdentityError(ValueError):
    """Raised when a payload conflicts with the active MCP run identity."""


def apply_current_run_id(
    payload: Mapping[str, Any],
    *,
    environ: Mapping[str, str] | None = None,
) -> dict[str, Any]:
    """Return a copy carrying the active run ID, if one is configured.

    A caller-provided matching ID is retained.  A conflicting ID is rejected
    before any public file is staged or replaced.
    """

    result = dict(payload)
    env = os.environ if environ is None else environ
    active_run_id = env.get(RUN_ID_ENV)
    if not active_run_id:
        return result

    supplied = result.get("run_id")
    if supplied is not None and supplied != active_run_id:
        raise RunIdentityError(
            f"payload run_id {supplied!r} conflicts with active run_id {active_run_id!r}"
        )
    result["run_id"] = active_run_id
    return result


def serialize_json_bytes(payload: Mapping[str, Any]) -> bytes:
    """Serialize one immutable JSON payload using the legacy public format."""

    return json.dumps(payload, indent=2, ensure_ascii=False).encode("utf-8")


def stage_exact(path: Path, payload: bytes) -> Path:
    """Write and fsync *payload* to a temporary sibling of *path*."""

    destination = Path(path)
    staged: Path | None = None
    try:
        destination.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile(
            mode="wb",
            prefix=f".{destination.name}.",
            suffix=".tmp",
            dir=destination.parent,
            delete=False,
        ) as handle:
            staged = Path(handle.name)
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        return staged
    except Exception as exc:
        if staged is not None:
            try:
                staged.unlink(missing_ok=True)
            except OSError:
                pass
        if isinstance(exc, ResultCommitError):
            raise
        raise ResultCommitError("stage", destination, exc) from exc


def _replace_exact(staged: Path, destination: Path) -> None:
    """Publish one staged sibling; kept separate as a fault-injection seam."""

    os.replace(staged, destination)


def commit_exact(
    primary: Path,
    payload: bytes,
    *,
    mirrors: Sequence[Path] = (),
) -> Path:
    """Atomically publish exact bytes to mirrors first and primary last."""

    primary = Path(primary)
    destinations: list[Path] = []
    for item in mirrors:
        destination = Path(item)
        if destination != primary and destination not in destinations:
            destinations.append(destination)
    destinations.append(primary)

    staged: dict[Path, Path] = {}
    try:
        for destination in destinations:
            try:
                staged[destination] = stage_exact(destination, payload)
            except Exception as exc:
                if isinstance(exc, ResultCommitError):
                    raise
                raise ResultCommitError("stage", destination, exc) from exc
        for destination in destinations:
            try:
                _replace_exact(staged[destination], destination)
            except Exception as exc:
                if isinstance(exc, ResultCommitError):
                    raise
                raise ResultCommitError("replace", destination, exc) from exc
        return primary
    finally:
        for staged_path in staged.values():
            try:
                staged_path.unlink(missing_ok=True)
            except OSError:
                pass


def commit_payloads(
    primary: Path,
    payloads: Mapping[Path, bytes],
) -> Path:
    """Stage distinct payloads, then publish companions before *primary*.

    This is the heterogeneous counterpart of :func:`commit_exact`.  Staging
    every destination and its prior generation first guarantees that
    serialization, staging, or backup failure cannot replace a public
    artifact.  Publishing the authoritative primary last keeps companion
    failures away from it; if a later publication fails, already-published
    companions are atomically restored before the error is returned.
    """

    primary = Path(primary)
    normalized = {Path(path): bytes(payload) for path, payload in payloads.items()}
    if primary not in normalized:
        raise ValueError(f"primary destination {primary} is missing from payloads")
    destinations = [path for path in normalized if path != primary] + [primary]

    staged: dict[Path, Path] = {}
    prior: dict[Path, Path | None] = {}
    published: list[Path] = []
    try:
        for destination in destinations:
            staged[destination] = stage_exact(destination, normalized[destination])
        for destination in destinations:
            if destination.exists():
                try:
                    prior[destination] = stage_exact(
                        destination, destination.read_bytes()
                    )
                except Exception as exc:
                    if isinstance(exc, ResultCommitError):
                        raise
                    raise ResultCommitError("backup", destination, exc) from exc
            else:
                prior[destination] = None
        for destination in destinations:
            try:
                _replace_exact(staged[destination], destination)
            except Exception as exc:
                for published_destination in reversed(published):
                    backup = prior[published_destination]
                    try:
                        if backup is None:
                            published_destination.unlink(missing_ok=True)
                        else:
                            _replace_exact(backup, published_destination)
                            prior[published_destination] = None
                    except Exception as rollback_exc:
                        raise ResultCommitError(
                            "rollback", published_destination, rollback_exc
                        ) from exc
                raise ResultCommitError("replace", destination, exc) from exc
            published.append(destination)
        return primary
    finally:
        for staged_path in (*staged.values(), *prior.values()):
            if staged_path is None:
                continue
            try:
                staged_path.unlink(missing_ok=True)
            except OSError:
                pass


def commit_json(
    primary: Path,
    payload: Mapping[str, Any],
    *,
    mirrors: Sequence[Path] = (),
) -> Path:
    """Serialize once, then publish the identical bytes to every destination."""

    try:
        serialized = serialize_json_bytes(payload)
    except Exception as exc:
        if isinstance(exc, ResultCommitError):
            raise
        raise ResultCommitError("serialize", Path(primary), exc) from exc
    return commit_exact(Path(primary), serialized, mirrors=mirrors)
