"""``all`` with the path-opt MEP engine and an endpoint HEI: no barrier, no TS dispatch.

The child fixture follows test_all_tsopt_bond_diagnostic.py. path-search builds no
segment for an interval whose highest image is an endpoint; the path-opt branch of
``all`` has to apply the same rule.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from click.testing import CliRunner

from pdb2reaction.cli import cli as root_cli
from pdb2reaction.core.result_commit import RUN_ID_ENV, apply_current_run_id
from pdb2reaction.workflows import all as all_workflow
import test_all_root_summary_golden as golden
from test_all_tsopt_bond_diagnostic import ReachedTSDispatch, _coords, _pdb, _write


def _frame(index, energy):
    atoms = "".join(
        f"{element} {x:.3f} {y:.3f} {z:.3f}\n"
        for element, (x, y, z) in zip(("C", "O"), _coords(index))
    )
    return f"2\nE={energy:.3f} unit=hartree frame{index}\n{atoms}"


@pytest.mark.parametrize(
    "energies,hei",
    [
        pytest.param((0.0, 0.1, 0.02), 1, id="interior-hei-control"),
        pytest.param((0.1, 0.05, 0.0), 0, id="reactant-end-hei"),
        pytest.param((0.0, 0.05, 0.1), 2, id="product-end-hei"),
    ],
)
def test_path_opt_branch_builds_no_segment_for_an_endpoint_hei(tmp_path, monkeypatch, energies, hei):
    monkeypatch.delenv(RUN_ID_ENV, raising=False)
    out = tmp_path / "out"
    inputs = [_write(tmp_path / f"input{i}.pdb", _pdb(i)) for i in (0, 2)]
    ts_calls = []

    def child(name, _cli, args, **kwargs):
        assert name == "path-opt"
        # The parent reads the child's MEP convergence from this result.json.
        assert "--out-json" in args
        child_out = Path(args[args.index("--out-dir") + 1])
        _write(child_out / "final_geometries_trj.xyz", "".join(_frame(i, e) for i, e in enumerate(energies)))
        _write(child_out / "hei.xyz", _frame(hei, energies[hei]))
        _write(child_out / "hei.pdb", _pdb(hei))
        payload = {
            "optimization_status": "converged", "converged": True,
            "preopt_requested": False, "preopt_converged": None,
            "path_optimizers": [],
            "stage_outcomes": [{"stage": "path-opt", "item_id": "gsm_mep",
                                "executed": True, "converged": True,
                                "usable": True, "required": True}],
        }
        _write(child_out / "result.json", json.dumps(apply_current_run_id(payload)))

    def stop_at_ts(hei_path, *args, **kwargs):
        ts_calls.append(Path(hei_path))
        raise ReachedTSDispatch

    def diagram(prefix, *, labels, energies_au, **kwargs):
        return {"name": prefix.stem, "labels": labels, "energies_au": list(energies_au)}

    monkeypatch.setattr(all_workflow, "_run_cli_main", child)
    monkeypatch.setattr(all_workflow, "run_trj2fig", lambda *args, **kwargs: None)
    monkeypatch.setattr(all_workflow, "close_matplotlib_figures", lambda: None)
    monkeypatch.setattr(all_workflow, "_write_segment_energy_diagram", diagram)
    monkeypatch.setattr(all_workflow._path_search, "has_bond_change", lambda *a, **k: (True, "Bond formed"))
    monkeypatch.setattr(all_workflow, "_run_tsopt_on_hei", stop_at_ts)
    args = ["all", *[arg for path in inputs for arg in ("-i", str(path))]]
    args += ["-q", "0", "-m", "1", "--out-dir", str(out),
             "--no-preopt", "--convert-files", "true", "--no-freeze-links",
             "--no-tsopt-from-mep-tan", "--tsopt", "true"]

    if hei == 1:
        with pytest.raises(ReachedTSDispatch):
            CliRunner().invoke(root_cli, args, catch_exceptions=False)
        assert len(ts_calls) == 1 and ts_calls[0].stem == "hei_seg_01"
        summary = json.loads((out / "_work/path_opt/summary.json").read_text())
        segment, = summary["segments"]
        assert segment["index"] == 1 and "stage_outcomes" not in summary
        diagram, = summary["energy_diagrams"]
        assert diagram["labels"] == ["R", "TS1", "P"]
        return

    # The interval is not a reaction step, so the run ends as a failed result.
    result = CliRunner().invoke(root_cli, args)
    assert result.exit_code == 1, result.output + repr(result.exception)
    assert ts_calls == []
    assert f"[seg_001] WARNING: HEI is at an endpoint (idx={hei}). Returning the raw path." in result.output
    summary = json.loads((out / "summary.json").read_text())
    assert summary["segments"] == [] and summary["n_segments"] == 0 and summary["n_images"] == 3
    assert summary["execution_status"] == "completed"
    assert summary["scientific_status"] == "failed"
    raw, = [leaf for leaf in summary["stage_outcomes"] if leaf["item_id"] == "raw_seg_001"]
    assert raw["stage"] == "path" and raw["usable"] is False and raw["reason"] == "endpoint_hei"
    assert "rate_limiting_step" not in summary
    # The MEP diagram draws no TS for an endpoint HEI.
    diagram, = summary["energy_diagrams"]
    assert diagram["labels"] == ["R", "P"]
    assert diagram["energies_au"] == [energies[0], energies[-1]]


@pytest.mark.parametrize("mep_mode", ["gsm", "dmf"])
@pytest.mark.parametrize(
    "second,labels",
    [
        pytest.param((0.02, 0.08, 0.01), ["R", "TS1", "IM1_1", "IM1_2", "TS2", "P"], id="two-interior-heis"),
        pytest.param((0.02, 0.01, 0.0), ["R", "TS1", "IM1_1", "IM1_2", "P"], id="second-endpoint-hei"),
    ],
)
def test_path_opt_mep_diagram_covers_every_interval(tmp_path, monkeypatch, second, labels, mep_mode):
    monkeypatch.delenv(RUN_ID_ENV, raising=False)
    out = tmp_path / "out"
    inputs = [_write(tmp_path / f"input{i}.pdb", _pdb(i)) for i in (0, 2, 4)]
    profiles = {"seg_001_mep": (0.0, 0.1, 0.02), "seg_002_mep": second}

    def child(name, _cli, args, **kwargs):
        assert name == "path-opt"
        # The parent reads the child's MEP convergence from this result.json.
        assert "--out-json" in args
        child_out = Path(args[args.index("--out-dir") + 1])
        energies = profiles[child_out.name]
        start = 2 * (int(child_out.name[4:7]) - 1)
        frames = "".join(_frame(start + i, e) for i, e in enumerate(energies))
        _write(child_out / "final_geometries_trj.xyz", frames)
        _write(child_out / "hei.xyz", _frame(start + 1, energies[1]))
        _write(child_out / "hei.pdb", _pdb(start + 1))
        payload = {
            "optimization_status": "converged", "converged": True,
            "preopt_requested": False, "preopt_converged": None,
            "path_optimizers": [],
            "stage_outcomes": [{"stage": "path-opt", "item_id": "gsm_mep",
                                "executed": True, "converged": True,
                                "usable": True, "required": True}],
        }
        _write(child_out / "result.json", json.dumps(apply_current_run_id(payload)))

    def stop_at_ts(*args, **kwargs):
        raise ReachedTSDispatch

    def diagram(prefix, *, labels, energies_au, title_note, **kwargs):
        return {"name": prefix.stem, "labels": labels, "energies_au": list(energies_au), "title": title_note}

    monkeypatch.setattr(all_workflow, "_run_cli_main", child)
    monkeypatch.setattr(all_workflow, "run_trj2fig", lambda *args, **kwargs: None)
    monkeypatch.setattr(all_workflow, "close_matplotlib_figures", lambda: None)
    monkeypatch.setattr(all_workflow, "_write_segment_energy_diagram", diagram)
    monkeypatch.setattr(all_workflow._path_search, "has_bond_change", lambda *a, **k: (True, "Bond formed"))
    monkeypatch.setattr(all_workflow, "_run_tsopt_on_hei", stop_at_ts)
    args = ["all", *[arg for path in inputs for arg in ("-i", str(path))]]
    args += ["-q", "0", "-m", "1", "--out-dir", str(out), "--mep-mode", mep_mode,
             "--no-preopt", "--convert-files", "true", "--no-freeze-links",
             "--no-tsopt-from-mep-tan", "--tsopt", "true"]

    with pytest.raises(ReachedTSDispatch):
        CliRunner().invoke(root_cli, args, catch_exceptions=False)
    summary = json.loads((out / "_work/path_opt/summary.json").read_text())
    diagram, = summary["energy_diagrams"]
    assert diagram["labels"] == labels
    # The title names the method that built the path.
    assert diagram["title"] == f"({mep_mode.upper()}; all segments)"
    expected = [0.0, 0.1, 0.02, *(second if len(labels) == 6 else (second[0], second[-1]))]
    assert diagram["energies_au"] == pytest.approx(expected)


def test_the_all_segments_diagram_numbers_the_ts_by_its_segment(tmp_path, monkeypatch):
    """The first interval builds no segment, so the one TS of the run is TS2 in every diagram.

    The post-processing runs through the stubs of test_all_root_summary_golden.py.
    """
    # Interval 1 (frames 0-2) rises to its end; frame 3 is the saddle of interval 2.
    monkeypatch.setattr(golden, "FRAME_ENERGIES", (-100.02, -100.01, -100.0, -99.95, -100.03))
    run = golden.invoke(tmp_path, monkeypatch, golden.Scenario("first_endpoint_hei", (0, 2, 4), ["--tsopt"]))
    summary = json.loads((run.out / "summary.json").read_text(encoding="utf-8"))

    assert [segment["index"] for segment in summary["segments"]] == [2]
    labels = {diagram["name"]: diagram["labels"] for diagram in summary["energy_diagrams"]}
    assert labels["energy_diagram_MLIP"] == ["R", "TS2", "P"]
    assert labels["energy_diagram_MLIP_all"] == ["R", "TS2", "P"]
