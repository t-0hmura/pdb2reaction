"""Fault-injection tests for exact-path result publication."""

from __future__ import annotations

import json
import os
from pathlib import Path

import pytest

from pdb2reaction.core import result_commit
from pdb2reaction.core.result_commit import ResultCommitError, RUN_ID_ENV
from pdb2reaction.core.utils import write_result_json
from pdb2reaction.mcp._runner import _read_current_summary


def _seed_pair(directory: Path, run_id: str = "old") -> tuple[bytes, bytes]:
    result = directory / "result.json"
    summary = directory / "summary.json"
    result.write_text(json.dumps({"run_id": run_id, "generation": "result-old"}))
    summary.write_text(json.dumps({"run_id": run_id, "generation": "summary-old"}))
    return result.read_bytes(), summary.read_bytes()


def _assert_no_staged_siblings(directory: Path) -> None:
    assert not list(directory.glob(".*.tmp"))


def test_successful_pair_is_byte_identical_and_does_not_mutate_input(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv(RUN_ID_ENV, "current")
    supplied = {
        "status": "success",
        "run_id": "current",
        "nested": {"value": 1},
    }

    primary = write_result_json(tmp_path, supplied, command="opt")

    assert primary == tmp_path / "result.json"
    assert primary.read_bytes() == (tmp_path / "summary.json").read_bytes()
    assert json.loads(primary.read_text())["run_id"] == "current"
    assert supplied == {
        "status": "success",
        "run_id": "current",
        "nested": {"value": 1},
    }
    _assert_no_staged_siblings(tmp_path)


def test_stage_exact_preserves_existing_destination_mode(tmp_path: Path) -> None:
    destination = tmp_path / "result.json"
    destination.write_text("old")
    destination.chmod(0o640)

    staged = result_commit.stage_exact(destination, b"new")
    try:
        assert os.stat(staged).st_mode & 0o777 == 0o640
    finally:
        staged.unlink()


def test_stage_exact_uses_ordinary_umask_mode_for_new_destination(tmp_path: Path) -> None:
    destination = tmp_path / "new-result.json"
    previous_umask = os.umask(0o027)
    try:
        staged = result_commit.stage_exact(destination, b"new")
    finally:
        os.umask(previous_umask)
    try:
        assert os.stat(staged).st_mode & 0o777 == 0o640
    finally:
        staged.unlink()


def test_serialization_failure_preserves_old_pair(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    old_result, old_summary = _seed_pair(tmp_path)

    def fail_serialize(_payload):
        raise TypeError("injected serialization failure")

    monkeypatch.setattr(result_commit, "serialize_json_bytes", fail_serialize)
    with pytest.raises(ResultCommitError, match="serialize"):
        write_result_json(tmp_path, {"status": "success"}, command="opt")

    assert (tmp_path / "result.json").read_bytes() == old_result
    assert (tmp_path / "summary.json").read_bytes() == old_summary
    _assert_no_staged_siblings(tmp_path)


@pytest.mark.parametrize("failed_stage", [1, 2])
def test_staging_failure_preserves_old_pair(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    failed_stage: int,
) -> None:
    old_result, old_summary = _seed_pair(tmp_path)
    real_stage = result_commit.stage_exact
    calls = 0

    def fail_selected_stage(path: Path, payload: bytes) -> Path:
        nonlocal calls
        calls += 1
        if calls == failed_stage:
            raise OSError("injected stage failure")
        return real_stage(path, payload)

    monkeypatch.setattr(result_commit, "stage_exact", fail_selected_stage)
    with pytest.raises(ResultCommitError, match="stage"):
        write_result_json(tmp_path, {"status": "success"}, command="opt")

    assert (tmp_path / "result.json").read_bytes() == old_result
    assert (tmp_path / "summary.json").read_bytes() == old_summary
    _assert_no_staged_siblings(tmp_path)


def test_mirror_replace_failure_preserves_old_pair(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    old_result, old_summary = _seed_pair(tmp_path)
    real_replace = result_commit._replace_exact

    def fail_mirror(staged: Path, destination: Path) -> None:
        if destination.name == "summary.json":
            raise OSError("injected mirror replace failure")
        real_replace(staged, destination)

    monkeypatch.setattr(result_commit, "_replace_exact", fail_mirror)
    with pytest.raises(ResultCommitError, match="replace"):
        write_result_json(tmp_path, {"status": "success"}, command="opt")

    assert (tmp_path / "result.json").read_bytes() == old_result
    assert (tmp_path / "summary.json").read_bytes() == old_summary
    _assert_no_staged_siblings(tmp_path)


def test_primary_replace_failure_leaves_valid_pair_but_runner_rejects_mismatch(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    old_result, _old_summary = _seed_pair(tmp_path)
    monkeypatch.setenv(RUN_ID_ENV, "current")
    real_replace = result_commit._replace_exact

    def fail_primary(staged: Path, destination: Path) -> None:
        if destination.name == "result.json":
            raise OSError("injected primary replace failure")
        real_replace(staged, destination)

    monkeypatch.setattr(result_commit, "_replace_exact", fail_primary)
    with pytest.raises(ResultCommitError, match="replace"):
        write_result_json(tmp_path, {"status": "success"}, command="opt")

    assert (tmp_path / "result.json").read_bytes() == old_result
    assert json.loads((tmp_path / "result.json").read_text())["run_id"] == "old"
    assert json.loads((tmp_path / "summary.json").read_text())["run_id"] == "current"
    status, payload, _hint = _read_current_summary(
        tmp_path / "summary.json", run_id="current", expected_primary=True
    )
    assert status == "summary_run_mismatch"
    assert payload == {}
    _assert_no_staged_siblings(tmp_path)


@pytest.mark.parametrize(
    ("filename", "also_write_summary_json"),
    [("summary.json", True), ("result.json", False)],
)
def test_one_file_replace_failure_preserves_old_file(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    filename: str,
    also_write_summary_json: bool,
) -> None:
    destination = tmp_path / filename
    destination.write_bytes(b"old-one-file")

    def fail_replace(_staged: Path, _destination: Path) -> None:
        raise OSError("injected one-file replace failure")

    monkeypatch.setattr(result_commit, "_replace_exact", fail_replace)
    with pytest.raises(ResultCommitError, match="replace"):
        write_result_json(
            tmp_path,
            {"status": "success"},
            command="test",
            filename=filename,
            also_write_summary_json=also_write_summary_json,
        )

    assert destination.read_bytes() == b"old-one-file"
    _assert_no_staged_siblings(tmp_path)


def test_conflicting_active_run_id_is_rejected_before_publication(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    old_result, old_summary = _seed_pair(tmp_path)
    monkeypatch.setenv(RUN_ID_ENV, "current")

    with pytest.raises(ValueError, match="conflicts"):
        write_result_json(
            tmp_path,
            {"status": "success", "run_id": "different"},
            command="opt",
        )

    assert (tmp_path / "result.json").read_bytes() == old_result
    assert (tmp_path / "summary.json").read_bytes() == old_summary


def test_primary_is_published_last_even_if_listed_as_a_mirror(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    primary = tmp_path / "primary.json"
    mirror = tmp_path / "mirror.json"
    published: list[Path] = []
    real_replace = result_commit._replace_exact

    def record_replace(staged: Path, destination: Path) -> None:
        published.append(destination)
        real_replace(staged, destination)

    monkeypatch.setattr(result_commit, "_replace_exact", record_replace)
    result_commit.commit_exact(
        primary,
        b"payload",
        mirrors=(primary, mirror, mirror),
    )

    assert published == [mirror, primary]


def test_current_aggregate_ignores_unrelated_stale_result_json(tmp_path: Path) -> None:
    (tmp_path / "summary.json").write_text(
        json.dumps({"run_id": "current", "status": "success"}),
        encoding="utf-8",
    )
    (tmp_path / "result.json").write_text(
        json.dumps({"run_id": "old", "status": "success"}),
        encoding="utf-8",
    )

    status, payload, _hint = _read_current_summary(
        tmp_path / "summary.json",
        run_id="current",
        expected_primary=False,
    )

    assert status == "ok"
    assert payload["run_id"] == "current"


def test_leaf_pair_with_same_run_id_but_different_content_is_rejected(
    tmp_path: Path,
) -> None:
    (tmp_path / "summary.json").write_text(
        json.dumps({"run_id": "current", "value": "summary"}),
        encoding="utf-8",
    )
    (tmp_path / "result.json").write_text(
        json.dumps({"run_id": "current", "value": "result"}),
        encoding="utf-8",
    )

    status, payload, _hint = _read_current_summary(
        tmp_path / "summary.json",
        run_id="current",
        expected_primary=True,
    )

    assert status == "summary_run_mismatch"
    assert payload == {}
