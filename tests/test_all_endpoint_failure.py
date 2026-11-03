"""Endpoint hard failures stop real all-branch control flow before consumers.

The TS-only and MEP-segment routes run through the real ``all`` CLI with the
children, DFT and calculator stubbed (tests/test_all_root_summary_golden.py);
only the endpoint optimizer is replaced, so exception handlers, retention,
return/continue and cleanup execute unchanged.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from pdb2reaction.workflows import all as workflow
from test_all_root_summary_golden import Scenario, coords, invoke, read_first_coords


def test_failed_ts_only_summary_does_not_invent_a_barrier(tmp_path):
    stop = {"stage": "endpoint_opt", "reason": "endpoint_execution_failed"}
    summary = {
        "segments": [{"index": 1, "tag": "seg_01", "kind": "tsopt"}],
        "pipeline_stop": stop, "energy_diagrams": [],
    }
    post = [{
        "index": 1, "kind": "tsopt", "pipeline_stop": stop,
        "tsopt": {"continue_irc": True, "saddle_validation": "first_order",
                  "optimization_status": "converged", "hessian_status": "completed",
                  "n_imaginary_modes": 1, "energy_valid": True,
                  "structure_valid": True},
        "irc": {"usable": True, "reason": "ok"}, "irc_traj": "irc.xyz",
        "endpoint_opt": {
            "reactant_converged": None,
            "product_converged": True,
            "failures": {"reactant": {"error_type": "ValueError", "error": "boom"}},
        },
    }]
    workflow._enrich_summary(
        summary, pipeline_mode="tsopt-only", out_dir=tmp_path,
        mlip_backend="orb", mlip_model="orb_v3_conservative_omol",
        charge=0, spin=1, post_segments=post,
        config={"tsopt": True, "thermo": True, "dft": True},
    )
    assert summary["scientific_status"] == "partial"
    assert summary["execution_status"] == "failed"
    assert any(
        "reactant_execution_failed" in reason
        for reason in summary["scientific_status_reasons"]
    )
    assert "rate_limiting_step" not in summary


def test_valid_ts1_with_one_nonconverged_endpoint_is_partial_not_failed(tmp_path):
    summary = {
        "segments": [{"index": 1, "kind": "tsopt"}],
        "energy_diagrams": [],
    }
    post = [{
        "index": 1,
        "kind": "tsopt",
        "tsopt": {
            "optimization_status": "converged",
            "hessian_status": "completed",
            "n_imaginary_modes": 1,
            "energy_valid": True,
            "structure_valid": True,
        },
        "endpoint_opt": {
            "reactant_converged": True,
            "product_converged": False,
        },
    }]

    truth = workflow._pipeline_aggregate_truth(
        summary,
        post_segments=post,
        config={"tsopt": True},
        legacy_status="partial",
    )

    assert truth.scientific_status == "partial"
    assert truth.execution_status == "completed"
    assert "all:segment_1:endpoint_opt:product_not_converged" in truth.status_reasons


class Geometry:
    def __init__(self, coordinate):
        self.cart_coords = np.array([coordinate, 0.0, 0.0])
        self.energy = -1.0
        self.calculator = None
        self.freeze_atoms = []

    def set_calculator(self, calculator):
        self.calculator = calculator

    def as_xyz(self):
        return f"1\nobserved\nH {self.cart_coords[0]} 0 0\n"


STATES = ("reactant", "product")


def _endpoint_optimizer(failed_endpoint, outcome_kind, calls):
    """Move each endpoint, then fail one of them, or return its convergence."""

    def optimize(_stubs, geom, _mode, out_dir, tag, *, dump, outcome):
        index = len(calls)
        calls.append(tag)
        if outcome is not None:
            outcome.update({
                "optimization_status": "converged", "converged": True,
                "n_opt_cycles": index + 3, "max_cycles": 3000,
            })
        child = Path(out_dir) / tag
        child.mkdir(parents=True, exist_ok=True)
        moved = np.array(geom.cart_coords, dtype=float)
        moved[0] = (10 + index) / workflow.BOHR2ANG
        geom.cart_coords = moved
        (child / "iterate.xyz").write_text(geom.as_xyz() + "\n")
        if index == failed_endpoint and outcome_kind in {"hard_error", "nonfinite_error"}:
            if outcome_kind == "nonfinite_error":
                geom.cart_coords = np.full_like(moved, np.nan)
            raise ValueError("endpoint numerical failure")
        return geom, child / "iterate.xyz", outcome_kind != "not_converged"

    return optimize


@pytest.mark.parametrize("branch", ["ts", "seg"])
@pytest.mark.parametrize("dump", [False, True])
# The failed endpoint matters only when an endpoint raises.
@pytest.mark.parametrize(
    ("failed_endpoint", "outcome"),
    [
        (0, "hard_error"), (1, "hard_error"),
        (0, "nonfinite_error"), (1, "nonfinite_error"),
        (0, "converged"), (0, "not_converged"),
    ],
)
def test_endpoint_boundary_retains_provenance_and_stops_consumers(
    monkeypatch, tmp_path, branch, failed_endpoint, dump, outcome,
):
    """Both endpoints are tried; a failure keeps its diagnostics and stops freq/DFT."""
    calls = []
    scenario = Scenario(
        f"endpoint_{branch}", (1,) if branch == "ts" else (0, 2),
        ["--tsopt", "--thermo", "--dft", "--dump" if dump else "--no-dump"],
    )
    run = invoke(tmp_path, monkeypatch, scenario,
                 endpoint_opt=_endpoint_optimizer(failed_endpoint, outcome, calls))
    seg_root = run.out / "segments" / "seg_01"
    structures = seg_root / "structures"
    summary = json.loads((run.out / "summary.json").read_text())
    (post,) = summary["post_segments"]
    prefix = "" if branch == "ts" else "seg_01_"

    assert calls == [f"{prefix}{state}" for state in STATES]
    assert (structures / "ts.xyz").is_file()
    # The IRC snapshots keep the IRC coordinates, not the optimizer's iterates.
    for state, frame in zip(STATES, (0, 2)):
        np.testing.assert_allclose(
            read_first_coords(structures / f"{state}_irc.xyz"), coords(frame), atol=1e-5,
        )
    assert run.closed == ["close"]
    failed = outcome in {"hard_error", "nonfinite_error"}
    freq_calls = [call for call in run.children if call[0] == "freq"]
    if failed:
        assert not any((structures / f"{state}.xyz").exists() for state in STATES)
        assert freq_calls == [] and run.dft_calls == []
        assert not (seg_root / "energy_diagram_MLIP.png").exists()
        assert not (seg_root / "reactant.pdb").exists()
        failure = json.loads((seg_root / "endpoint_opt" / "failure.json").read_text())
        info = failure["failures"][STATES[failed_endpoint]]
        assert set(failure["failures"]) == {STATES[failed_endpoint]}
        assert info["error_type"] == "ValueError"
        assert info["error"] == "endpoint numerical failure"
        assert Path(info["irc_structure"]).is_file()
        assert list((seg_root / "endpoint_opt").rglob("iterate.xyz"))
        assert info["geometry_role"] == "last_observed_acceptance_unknown"
        assert info["coordinates_finite"] is (outcome == "hard_error")
        if outcome == "hard_error":
            observed = read_first_coords(Path(info["last_observed_xyz"]))
            assert observed[0][0] == pytest.approx(10 + failed_endpoint)
        else:
            assert "last_observed_xyz" not in info
        assert post["pipeline_stop"]["stage"] == "endpoint_opt"
        if branch == "ts":
            assert summary["pipeline_stop"]["stage"] == "endpoint_opt"
            assert summary["execution_status"] == "failed"
    else:
        assert all((structures / f"{state}.xyz").is_file() for state in STATES)
        assert len(freq_calls) == 3 and run.dft_calls == ["R", "TS", "P"]
        assert (seg_root / "endpoint_opt").exists() is (dump or outcome == "not_converged")
        assert post["endpoint_opt"]["reactant"]["n_opt_cycles"] == 3
        assert post["endpoint_opt"]["product"]["n_opt_cycles"] == 4


def test_ts_only_endpoint_failure_records_irc_trajectory(monkeypatch, tmp_path):
    scenario = Scenario("ts_only_product_failure", (1,), ["--tsopt"],
                        endpoint_failure="product")
    run = invoke(tmp_path, monkeypatch, scenario)
    (stop_log,) = json.loads((run.out / "summary.json").read_text())["post_segments"]
    assert stop_log["pipeline_stop"]["stage"] == "endpoint_opt"
    assert set(stop_log["endpoint_opt"]["failures"]) == {"product"}
    trajectory = run.out / "segments" / "seg_01" / "irc" / "finished_irc_trj.xyz"
    assert trajectory.is_file()
    assert stop_log["irc_traj"] == str(trajectory)


def test_ts_only_endpoint_failure_summary_keeps_provenance_and_config(monkeypatch, tmp_path):
    """The endpoint-failure exit records the same calculator and config as the other TS-only exits."""
    failed = invoke(tmp_path / "failed", monkeypatch, Scenario(
        "ts_only_product_failure", (1,), ["--tsopt"], endpoint_failure="product",
    ))
    stopped = invoke(tmp_path / "stopped", monkeypatch, Scenario(
        "ts_only_stopped", (1,), ["--tsopt"], tsopt_converged=False,
    ))
    failed_summary = json.loads((failed.out / "summary.json").read_text())
    stopped_summary = json.loads((stopped.out / "summary.json").read_text())
    for key in ("primary_method", "primary_method_label"):
        assert failed_summary[key] == stopped_summary[key]
    assert failed_summary["primary_method"] == "mlip"
    assert failed_summary["config"] == stopped_summary["config"]


@pytest.mark.parametrize(
    "outcome",
    ["converged", "not_converged", "handled_stop", "missing", "stale", "nonfinite_coords", "nonfinite_energy"],
)
def test_endpoint_helper_requires_current_finite_output(monkeypatch, tmp_path, outcome):
    import click
    from pdb2reaction.workflows._run_session import ArtifactClaimError

    final_path = tmp_path / "reactant_lbfgs_opt" / "final.xyz"
    if outcome == "stale":
        final_path.parent.mkdir()
        final_path.write_text("prior invocation")
    geom = Geometry(1)

    class Optimizer:
        def __init__(self, geometry, **_kwargs):
            self.final_fn = final_path
            self.is_converged = outcome == "converged"
            self.is_stalled = False
            self.cur_cycle = 6
            self.stop_reason = ""

        def run(self):
            geom.cart_coords[0] = 2
            if outcome == "handled_stop":
                raise workflow.ZeroStepLength("test stop")
            if outcome not in {"missing", "stale"}:
                final_path.write_text(geom.as_xyz())

    terminal = Geometry(np.nan if outcome == "nonfinite_coords" else 2)
    if outcome == "nonfinite_energy":
        terminal.energy = np.inf
    monkeypatch.setattr(workflow, "LBFGS", Optimizer)
    monkeypatch.setattr(workflow, "geom_loader", lambda *_a, **_kw: terminal)

    if outcome in {"missing", "stale", "nonfinite_coords", "nonfinite_energy"}:
        with pytest.raises((click.ClickException, ArtifactClaimError)):
            workflow._optimize_endpoint_geom(
                geom, "grad", tmp_path, "reactant", dump=False, thresh=None,
            )
    else:
        endpoint_outcome = {}
        result, path, converged = workflow._optimize_endpoint_geom(
            geom, "grad", tmp_path, "reactant", dump=False, thresh=None,
            outcome=endpoint_outcome,
        )
        assert result is terminal and path == final_path
        assert converged is (outcome == "converged")
        assert path.is_file()
        assert endpoint_outcome["optimization_status"] == (
            "converged" if outcome == "converged" else "not_converged"
        )
        assert endpoint_outcome["n_opt_cycles"] == 7
        assert "max_cycles" in endpoint_outcome
