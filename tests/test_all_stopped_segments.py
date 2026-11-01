"""``all`` keeps post-processing the later segments after one stops before IRC.

The stubs are those of test_all_root_summary_golden.py, extended to three
segments (frame 2k + 1 is the saddle of segment k + 1). The TS optimization of
the segments in ``STOPPED`` does not converge, so they stop before IRC.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

import test_all_root_summary_golden as golden

ENERGIES = (-100.0, -99.95, -100.01, -99.955, -100.02, -99.96, -100.03)
CORRECTIONS = (0.010, 0.008, 0.011, 0.009, 0.012, 0.0095, 0.013)
STOPPED: tuple = ()


class ThreeSegmentStubs(golden.PipelineStubs):
    def _tsopt(self, out: Path, args):
        frame = golden.frame_of(golden.read_first_coords(Path(golden.option_value(args, "-i"))))
        geometry = golden.write(out / "final_geometry.xyz", golden.xyz_frame(frame))
        golden.write(out / "result.json", json.dumps(golden.apply_current_run_id(golden.ts_result(
            ENERGIES[frame], converged=(frame + 1) // 2 not in STOPPED,
            files={"final_geometry": str(geometry)},
        ))))


def run_all(tmp_path, monkeypatch, stopped):
    monkeypatch.setattr(golden, "FRAME_ENERGIES", ENERGIES)
    monkeypatch.setattr(golden, "FRAME_G_CORRECTIONS", CORRECTIONS)
    monkeypatch.setattr(golden, "PipelineStubs", ThreeSegmentStubs)
    monkeypatch.setitem(globals(), "STOPPED", tuple(stopped))
    scenario = golden.Scenario("three_segments", (0, 2, 4, 6), ["--tsopt"])
    run = golden.invoke(tmp_path, monkeypatch, scenario)
    summary = json.loads((run.out / "summary.json").read_text(encoding="utf-8"))
    return run, summary


@pytest.mark.parametrize("position", [1, 2, 3])
def test_a_segment_that_stops_before_irc_does_not_stop_the_others(tmp_path, monkeypatch, position):
    run, summary = run_all(tmp_path, monkeypatch, [position])

    assert [call[0] for call in run.children].count("tsopt") == 3
    assert [call[0] for call in run.children].count("irc") == 2
    assert [item["index"] for item in summary["post_segments"]] == [1, 2, 3]
    assert (run.exit_code, summary["execution_status"], summary["scientific_status"]) == (
        0, "completed", "partial")
    reason = "ts_optimization_not_converged"
    assert summary["stopped_before_irc"] is True
    assert summary["pipeline_stop"]["segment"] == position
    assert summary["stopped_segments"] == [
        {"stage": "before_irc", "segment": position, "reason": reason}]
    log = (run.out / "summary.log").read_text(encoding="utf-8")
    assert f"Stopped before IRC  : segment {position} ({reason})" in log
    # A segment without an IRC would leave a gap in the one concatenated curve.
    assert not (run.out / "irc_plot_all.png").exists()
    assert not (run.out / "energy_diagram_MLIP_all.png").exists()


def test_every_segment_stopped_is_listed_and_pipeline_stop_is_the_first(tmp_path, monkeypatch):
    run, summary = run_all(tmp_path, monkeypatch, [1, 2, 3])

    assert [call[0] for call in run.children].count("irc") == 0
    assert [item["index"] for item in summary["post_segments"]] == [1, 2, 3]
    assert [item["segment"] for item in summary["stopped_segments"]] == [1, 2, 3]
    assert summary["pipeline_stop"]["segment"] == 1
    assert (run.exit_code, summary["execution_status"], summary["scientific_status"]) == (
        0, "completed", "partial")
    assert not (run.out / "irc_plot_all.png").exists()


def test_no_stop_draws_the_irc_plot_over_all_segments(tmp_path, monkeypatch):
    run, summary = run_all(tmp_path, monkeypatch, [])

    assert [call[0] for call in run.children].count("irc") == 3
    assert (run.exit_code, summary["scientific_status"]) == (0, "success")
    assert "pipeline_stop" not in summary and "stopped_segments" not in summary
    assert summary["stopped_before_irc"] is False
    assert (run.out / "irc_plot_all.png").is_file()


def test_resume_keeps_only_the_stops_of_the_current_run(tmp_path, monkeypatch):
    first, summary = run_all(tmp_path, monkeypatch, [3])
    assert summary["stopped_segments"][0]["segment"] == 3

    monkeypatch.setattr(golden, "PipelineStubs", ThreeSegmentStubs)
    monkeypatch.setitem(globals(), "STOPPED", ())
    resumed = golden.invoke(
        tmp_path, monkeypatch,
        golden.Scenario("resume", (0, 2, 4, 6), ["--resume-segment", "2", "--tsopt"]),
    )
    summary = json.loads((resumed.out / "summary.json").read_text(encoding="utf-8"))

    assert [call[0] for call in resumed.children].count("tsopt") == 2
    assert summary["stopped_before_irc"] is False
    assert "pipeline_stop" not in summary and "stopped_segments" not in summary
    assert (resumed.exit_code, summary["scientific_status"]) == (0, "success")
