"""A failed terminal Hessian in a segment's TS optimization stops the whole ``all`` run.

The stubs are those of test_all_stopped_segments.py. The tsopt child of the
segments in ``HESSIAN_FAILED`` goes through the real ``_run_cli_main`` dispatch
and ``completion_guard`` and leaves the record tsopt writes when the terminal
Hessian raises, so the child exits 1 and ``all`` stops with exit 1 instead of
continuing to the later segments.
"""

from __future__ import annotations

import json
from pathlib import Path

import click
import pytest

import test_all_root_summary_golden as golden
from test_all_stopped_segments import CORRECTIONS, ENERGIES
from pdb2reaction.cli.completion import completion_guard
from pdb2reaction.core.utils import write_result_json
from pdb2reaction.workflows import all as all_workflow

REAL_RUN_CLI_MAIN = all_workflow._run_cli_main
HESSIAN_FAILED: tuple = ()


class HessianFailureStubs(golden.PipelineStubs):
    def child(self, name, cli, args, **kwargs):
        if name != "tsopt":
            return super().child(name, cli, args, **kwargs)
        args = [str(arg) for arg in args]
        frame = golden.frame_of(golden.read_first_coords(Path(golden.option_value(args, "-i"))))
        if (frame + 1) // 2 not in HESSIAN_FAILED:
            return super().child(name, cli, args, **kwargs)
        self.children.append([name, *args])
        out = Path(golden.option_value(args, "--out-dir"))

        @click.command(context_settings={"ignore_unknown_options": True, "allow_extra_args": True})
        def tsopt_with_failed_hessian():
            geometry = golden.write(out / "final_geometry.xyz", golden.xyz_frame(frame))
            payload = golden.ts_result(ENERGIES[frame], files={"final_geometry": str(geometry)})
            # The record tsopt leaves when the terminal Hessian raises.
            payload.update(
                hessian_status="failed", hessian_error="RuntimeError: boom",
                saddle_validation="unavailable", n_imaginary_modes=None,
                n_negative_modes=None, imaginary_frequencies_cm=[],
                reaction_mode_index=None, reaction_mode_frequency_cm=None,
                reaction_mode_overlap=None, reaction_mode_source=None,
            )
            payload.pop("execution_status", None)
            payload.pop("scientific_status", None)
            write_result_json(out, payload, command="tsopt")

        tsopt_with_failed_hessian.callback = completion_guard(tsopt_with_failed_hessian.callback)
        return REAL_RUN_CLI_MAIN(name, tsopt_with_failed_hessian, args, **kwargs)


@pytest.mark.parametrize(
    "failed",
    [
        pytest.param((2,), id="segment-2"),
        pytest.param((3,), id="segment-3"),
        pytest.param((1, 2, 3), id="every-segment"),
    ],
)
def test_a_failed_terminal_hessian_stops_the_whole_run(tmp_path, monkeypatch, failed):
    monkeypatch.setattr(golden, "FRAME_ENERGIES", ENERGIES)
    monkeypatch.setattr(golden, "FRAME_G_CORRECTIONS", CORRECTIONS)
    monkeypatch.setattr(golden, "PipelineStubs", HessianFailureStubs)
    monkeypatch.setitem(globals(), "HESSIAN_FAILED", failed)
    run = golden.invoke(
        tmp_path, monkeypatch, golden.Scenario("three_segments", (0, 2, 4, 6), ["--tsopt"])
    )
    summary = json.loads((run.out / "summary.json").read_text(encoding="utf-8"))

    first = min(failed)
    names = [call[0] for call in run.children]
    # No segment after the first failure is post-processed.
    assert names.count("tsopt") == first
    assert names.count("irc") == first - 1
    assert run.exit_code == 1
    assert "tsopt exit code 1" in run.output
    assert summary["execution_status"] == "failed"
    assert summary["scientific_status"] != "success"


def test_resume_after_a_stop_in_post_processing_keeps_the_finished_segments(tmp_path, monkeypatch):
    monkeypatch.setattr(golden, "FRAME_ENERGIES", ENERGIES)
    monkeypatch.setattr(golden, "FRAME_G_CORRECTIONS", CORRECTIONS)
    monkeypatch.setattr(golden, "PipelineStubs", HessianFailureStubs)
    monkeypatch.setitem(globals(), "HESSIAN_FAILED", (2,))
    first = golden.invoke(
        tmp_path, monkeypatch, golden.Scenario("three_segments", (0, 2, 4, 6), ["--tsopt"])
    )
    assert first.exit_code == 1
    summary = json.loads((first.out / "summary.json").read_text(encoding="utf-8"))
    # Segment 1 finished before the run stopped in segment 2.
    assert [item["index"] for item in summary["post_segments"]] == [1]
    assert summary["scientific_status"] != "success"

    monkeypatch.setitem(globals(), "HESSIAN_FAILED", ())
    resumed = golden.invoke(
        tmp_path, monkeypatch,
        golden.Scenario("resume", (0, 2, 4, 6), ["--resume-segment", "2", "--tsopt"]),
    )
    summary = json.loads((resumed.out / "summary.json").read_text(encoding="utf-8"))

    assert [call[0] for call in resumed.children].count("tsopt") == 2
    assert [item["index"] for item in summary["post_segments"]] == [1, 2, 3]
    assert (resumed.exit_code, summary["execution_status"], summary["scientific_status"]) == (
        0, "completed", "success")
