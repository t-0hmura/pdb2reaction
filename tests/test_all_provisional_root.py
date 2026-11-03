"""The root summary.json that ``all`` keeps while post-processing runs.

The stubs are those of test_all_root_summary_golden.py and
test_all_stopped_segments.py. While a segment is post-processed, the root lists
the MEP outputs at the root and the diagrams of the finished segments; a resumed
run publishes this provisional root with the kept segments before
post-processing starts, and refuses a changed post-processing precision.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

import test_all_root_summary_golden as golden
from test_all_stopped_segments import CORRECTIONS, ENERGIES
from pdb2reaction.workflows import all as all_workflow

ROOT_AT_TSOPT: dict = {}


class RootAtTsoptStubs(golden.PipelineStubs):
    """Keep the root summary.json, and which of its listed files exist, as each TS optimization starts."""

    def child(self, name, cli, args, **kwargs):
        if name == "tsopt":
            args = [str(arg) for arg in args]
            frame = golden.frame_of(golden.read_first_coords(Path(golden.option_value(args, "-i"))))
            out = self.tmp_path / "out"
            root = json.loads((out / "summary.json").read_text(encoding="utf-8"))
            listed = list(root.get("current_output_paths", [])) + [
                str(Path(diagram["image"]).relative_to(out))
                for diagram in root.get("energy_diagrams", [])
            ]
            ROOT_AT_TSOPT[(frame + 1) // 2] = (root, [p for p in listed if (out / p).is_file()])
        return super().child(name, cli, args, **kwargs)


@pytest.fixture
def stubs(monkeypatch):
    monkeypatch.setattr(golden, "FRAME_ENERGIES", ENERGIES)
    monkeypatch.setattr(golden, "FRAME_G_CORRECTIONS", CORRECTIONS)
    monkeypatch.setattr(golden, "PipelineStubs", RootAtTsoptStubs)
    ROOT_AT_TSOPT.clear()
    yield ROOT_AT_TSOPT
    ROOT_AT_TSOPT.clear()


def _summary(run):
    return json.loads((run.out / "summary.json").read_text(encoding="utf-8"))


def test_refine_path_run_root_lists_the_root_mep_outputs(tmp_path, monkeypatch, stubs):
    out = tmp_path / "out"
    at_post_processing = []
    echo_section = all_workflow._echo_section

    def record_at_post_processing(message, *args, **kwargs):
        if "Stage 4/" in str(message):
            at_post_processing.append(json.loads((out / "summary.json").read_text(encoding="utf-8")))
        return echo_section(message, *args, **kwargs)

    monkeypatch.setattr(all_workflow, "_echo_section", record_at_post_processing)
    run = golden.invoke(
        tmp_path, monkeypatch,
        golden.Scenario(
            "refine_two_segments", (0, 4), ["--refine-path", "--tsopt"],
            path_search_segments=[golden.SEGMENT_1, golden.SEGMENT_2],
        ),
    )
    assert run.exit_code == 0, run.output[-3000:]
    assert sorted(stubs) == [1, 2]

    # The root published before post-processing starts.
    (root,) = at_post_processing
    paths = root["current_output_paths"]
    assert {"energy_diagram_MEP.png", "mep_trj.xyz"} <= set(paths)
    assert not any(path.startswith("_work") for path in paths)
    assert all((out / path).is_file() for path in paths)
    assert root["scientific_status"] != "success"
    assert root["post_segments"] == []

    root, existing = stubs[1]
    assert {"energy_diagram_MEP.png", "mep_trj.xyz"} <= set(root["current_output_paths"])
    # Every listed output exists at the root.
    assert set(existing) == set(root["current_output_paths"])

    # Segment 2 starts with the record and the diagram of segment 1.
    root, existing = stubs[2]
    assert [item["index"] for item in root["post_segments"]] == [1]
    images = [Path(diagram["image"]).relative_to(run.out).as_posix()
              for diagram in root["energy_diagrams"]]
    assert "segments/seg_01/energy_diagram_MLIP.png" in images
    assert set(images) <= set(existing)


def _completed_three_segments(tmp_path, monkeypatch, extra=()):
    first = golden.invoke(
        tmp_path, monkeypatch,
        golden.Scenario("three_segments", (0, 2, 4, 6), ["--tsopt", *extra]),
    )
    assert first.exit_code == 0, first.output[-3000:]
    summary = _summary(first)
    assert [item["index"] for item in summary["post_segments"]] == [1, 2, 3]
    return first, summary


def test_resume_that_stops_before_its_first_segment_keeps_a_provisional_root(
    tmp_path, monkeypatch, stubs
):
    first, summary = _completed_three_segments(tmp_path, monkeypatch)
    assert summary["scientific_status"] == "success"
    echo_section = all_workflow._echo_section

    def stop_at_post_processing(message, *args, **kwargs):
        if "Stage 4/" in str(message):
            raise RuntimeError("walltime reached before segment 2")
        return echo_section(message, *args, **kwargs)

    with monkeypatch.context() as patch:
        patch.setattr(all_workflow, "_echo_section", stop_at_post_processing)
        with pytest.raises(RuntimeError, match="walltime"):
            golden.invoke(
                tmp_path, monkeypatch,
                golden.Scenario("resume", (0, 2, 4, 6), ["--resume-segment", "2", "--tsopt"]),
            )
    root = json.loads((first.out / "summary.json").read_text(encoding="utf-8"))
    assert root["scientific_status"] != "success"
    assert root["resumed_from_segment"] == 2
    assert [item["index"] for item in root["post_segments"]] == [1]

    resumed = golden.invoke(
        tmp_path, monkeypatch,
        golden.Scenario("resume", (0, 2, 4, 6), ["--resume-segment", "2", "--tsopt"]),
    )
    summary = _summary(resumed)
    assert [call[0] for call in resumed.children].count("tsopt") == 2
    assert [item["index"] for item in summary["post_segments"]] == [1, 2, 3]
    assert (resumed.exit_code, summary["scientific_status"]) == (0, "success")


def test_resume_after_a_stop_in_segment_2_keeps_the_diagram_of_segment_1(
    tmp_path, monkeypatch, stubs
):
    echo_section = all_workflow._echo_section

    def stop_in_segment_2(message, *args, **kwargs):
        if "seg_02" in str(message):
            raise RuntimeError("walltime reached in segment 2")
        return echo_section(message, *args, **kwargs)

    with monkeypatch.context() as patch:
        patch.setattr(all_workflow, "_echo_section", stop_in_segment_2)
        with pytest.raises(RuntimeError, match="walltime"):
            golden.invoke(
                tmp_path, monkeypatch,
                golden.Scenario("three_segments", (0, 2, 4, 6), ["--tsopt"]),
            )
    resumed = golden.invoke(
        tmp_path, monkeypatch,
        golden.Scenario("resume", (0, 2, 4, 6), ["--resume-segment", "2", "--tsopt"]),
    )
    summary = _summary(resumed)
    assert (resumed.exit_code, summary["scientific_status"]) == (0, "success")
    images = [Path(diagram["image"]).relative_to(resumed.out).as_posix()
              for diagram in summary["energy_diagrams"]]
    for index in (1, 2, 3):
        assert f"segments/seg_{index:02d}/energy_diagram_MLIP.png" in images
    assert "segments/seg_01/energy_diagram_MLIP.png" in (
        resumed.out / "summary.log"
    ).read_text(encoding="utf-8")
    assert "segments/seg_01/ts.pdb" in summary["current_output_paths"]


@pytest.mark.parametrize(
    "resume_segment,precision,refused",
    [("2", "fp32", True), ("2", "fp64", False), ("1", "fp32", False)],
)
def test_resume_refuses_a_changed_post_processing_precision(
    tmp_path, monkeypatch, stubs, resume_segment, precision, refused
):
    first, summary = _completed_three_segments(tmp_path, monkeypatch, ("--precision-post", "fp64"))
    assert summary["mlip_precision_post"] == "fp64"
    before = (first.out / "summary.json").read_bytes()

    resumed = golden.invoke(
        tmp_path, monkeypatch,
        golden.Scenario(
            "resume", (0, 2, 4, 6),
            ["--tsopt", "--resume-segment", resume_segment, "--precision-post", precision],
        ),
    )

    if refused:
        assert resumed.exit_code == 2
        assert resumed.children == []
        assert "post-processing precision changed (fp64 in the saved run, fp32 now)" in (
            resumed.output
        )
        assert (first.out / "summary.json").read_bytes() == before
    else:
        assert resumed.exit_code == 0, resumed.output[-3000:]
        assert [call[0] for call in resumed.children].count("tsopt") == 4 - int(resume_segment)
        assert _summary(resumed)["mlip_precision_post"] == precision
