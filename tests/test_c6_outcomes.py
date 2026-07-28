"""C6 — truthful scientific outcomes.

Every test in this file asserts the load-bearing C6 invariant: no artifact
existence and no finite fallback may promote a required nonconverged / missing
scientific leaf to success, while a genuinely converged leaf's public output is
unchanged aside from the additive outcome fields.

The falsifiers are grouped by the false-promotion path each one closes; every one
of them would have reported SUCCESS (or a wrong minimum / a cached Hessian / a
Gibbs diagram) under the pre-C6 fallback and now correctly reports
nonconverged / missing / excluded.
"""

from __future__ import annotations

import json
import tempfile
from pathlib import Path

import pytest

from pdb2reaction.workflows._outcomes import (
    AggregateTruth,
    LeafOutcome,
    ScanPointOutcome,
    aggregate_workflow_truth,
    attach_outcomes,
    eligible_points,
    make_leaf,
    make_scan_point,
    scan_scientific_status,
    seed_eligible_mask,
)


# ---------------------------------------------------------------------------
# 1. Pure outcome types + serializer compatibility (P07)
# ---------------------------------------------------------------------------


def test_leaf_outcome_fail_closed_usability() -> None:
    # Executed + explicitly converged + finite -> usable.
    assert make_leaf("s", "a", executed=True, converged=True).usable is True
    # Explicitly not converged -> unusable, even though it executed.
    nc = make_leaf("s", "a", executed=True, converged=False)
    assert nc.usable is False and nc.reason == "not_converged"
    # Unknown convergence (None) fails closed -> unusable.
    unk = make_leaf("s", "a", executed=True, converged=None)
    assert unk.usable is False and unk.reason == "convergence_unknown"
    # Converged but energy invalid -> unusable.
    bad = make_leaf("s", "a", executed=True, converged=True, energy_valid=False)
    assert bad.usable is False and bad.reason == "energy_invalid"


def test_scan_point_seed_eligibility_fail_closed() -> None:
    ok = make_scan_point("p", executed=True, converged=True, energy=-1.0, artifact_written=True)
    assert ok.seed_eligible is True
    # Nonconverged is never a seed even with a finite energy + artifact.
    nc = make_scan_point("p", executed=True, converged=False, energy=-9.0, artifact_written=True)
    assert nc.seed_eligible is False
    # Converged but the artifact write failed -> not eligible (artifact_missing).
    na = make_scan_point("p", executed=True, converged=True, energy=-1.0, artifact_written=False)
    assert na.seed_eligible is False and na.reason == "artifact_missing"
    # NaN energy -> not eligible.
    nan = make_scan_point("p", executed=True, converged=True, energy=float("nan"), artifact_written=True)
    assert nan.seed_eligible is False and nan.reason == "energy_invalid"


def test_aggregate_success_partial_failed() -> None:
    # All required usable, no missing expected -> success.
    t = aggregate_workflow_truth(
        [make_leaf("p", "seg_1", executed=True, converged=True)], ["seg_1"]
    )
    assert (t.scientific_status, t.execution_status) == ("success", "completed")
    # A usable diagnostic artifact + a missing required leaf -> partial.
    raw = LeafOutcome("p", "raw", required=True, executed=True, converged=True,
                      usable=False, reason="endpoint_hei", artifacts=("mep.pdb",))
    t2 = aggregate_workflow_truth([raw], ["seg_1"])
    assert t2.scientific_status == "partial"
    assert "missing:seg_1" in t2.status_reasons
    # Nothing usable, no artifact, execution failed -> failed.
    t3 = aggregate_workflow_truth(
        [make_leaf("p", "seg_1", required=True, executed=False, converged=None)], ["seg_1"]
    )
    assert (t3.scientific_status, t3.execution_status) == ("failed", "failed")


def test_serializer_roundtrip_and_additive_only(tmp_path: Path) -> None:
    from pdb2reaction.core.utils import write_result_json, RESULT_JSON_SCHEMA_VERSION

    leaf = make_leaf("scan", "stage_1", executed=True, converged=True)
    truth = aggregate_workflow_truth([leaf], ["stage_1"])
    # A leaf-command result: legacy status stays "completed", outcomes are added.
    data = {"status": "completed", "min_energy_hartree": -1.5}
    attach_outcomes(data, truth=truth, stage_outcomes=[leaf])
    write_result_json(tmp_path, dict(data), command="scan")
    r = json.loads((tmp_path / "result.json").read_text())

    # Legacy contract preserved.
    assert r["status"] == "completed"
    assert r["schema_version"] == RESULT_JSON_SCHEMA_VERSION
    assert r["min_energy_hartree"] == -1.5
    # Additive outcome fields present and truthful.
    assert r["scientific_status"] == "success"
    assert r["execution_status"] == "completed"
    assert r["expected_item_ids"] == ["stage_1"]
    assert r["stage_outcomes"][0]["usable"] is True

    # An old reader that only knows `status` obtains the same type/value and is
    # unaffected by the additive fields.
    assert isinstance(r["status"], str) and r["status"] == "completed"


def test_dataclasses_are_json_safe() -> None:
    assert LeafOutcome("s", "a").to_dict()["item_id"] == "a"
    assert ScanPointOutcome("p").to_dict()["seed_eligible"] is False
    assert AggregateTruth("completed", "success").to_dict()["scientific_status"] == "success"


# ---------------------------------------------------------------------------
# 2. FALSIFIER — scan non-convergence (M50 + M48)
#    A failed point with the numerically lowest energy would have become the
#    reported minimum / baseline under the old raw-min; it must be excluded.
# ---------------------------------------------------------------------------


def test_scan_failed_low_energy_point_excluded_from_minimum() -> None:
    # These carry the seed-eligibility fields of the records scan2d/scan3d build
# (tri-state bias_converged).
    records = [
        {"i": 0, "j": 0, "energy_hartree": -1.00, "bias_converged": True, "artifact_written": True},
        {"i": 0, "j": 1, "energy_hartree": -1.20, "bias_converged": True, "artifact_written": True},
        # A FAILED point that happens to have the lowest raw energy.
        {"i": 0, "j": 2, "energy_hartree": -9.99, "bias_converged": False, "artifact_written": True},
    ]
    mask = seed_eligible_mask(records)
    assert mask == [True, True, False]

    raw_min = min(r["energy_hartree"] for r in records)
    eligible = [r for r, ok in zip(records, mask) if ok]
    eligible_min = min(r["energy_hartree"] for r in eligible)

    # OLD behavior would have reported the failed point's energy.
    assert raw_min == -9.99
    # NEW behavior: the reported minimum comes only from converged points.
    assert eligible_min == -1.20

    status, reasons = scan_scientific_status(
        [
            make_scan_point(f"{r['i']}_{r['j']}", executed=True,
                            converged=r["bias_converged"], energy=r["energy_hartree"],
                            artifact_written=r["artifact_written"])
            for r in records
        ]
    )
    assert status == "partial"
    assert reasons == ("unusable_points:1",)


def test_scan_normal_return_is_not_convergence() -> None:
    # M50: an optimizer that returns normally but reports is_converged=False (a
    # cycle-limit stop) must not be recorded as converged. seed_eligible_mask
    # reads the recorded bit; only an explicit True survives.
    normal_return_but_not_converged = {
        "energy_hartree": -1.0, "bias_converged": False,
        "artifact_written": True,
    }
    unknown_convergence = {
        "energy_hartree": -1.0, "bias_converged": None,
        "artifact_written": True,
    }
    converged = {
        "energy_hartree": -1.0, "bias_converged": True,
        "artifact_written": True,
    }
    assert seed_eligible_mask(
        [normal_return_but_not_converged, unknown_convergence, converged]
    ) == [False, False, True]


def test_scan_stage_leaf_partial_when_middle_step_fails() -> None:
    # Binds to the extracted _outcomes.combine_step_convergence helper that
    # scan.py's production stage-leaf fold calls (this test exercises the shared
    # helper, not scan.py's own code path):
    # stage_1 all steps converged; stage_2 has a middle step that failed
    # (OptimizationError) while the final step converged. Legacy per-stage
    # `status` would still read "completed", but the aggregate scientific_status
    # must be partial because a converged final step cannot hide the failure.
    from pdb2reaction.workflows._outcomes import combine_step_convergence as _combine

    assert _combine([True, True, True]) is True
    assert _combine([True, False, True]) is False  # middle-step failure survives
    assert _combine([True, None, True]) is None     # unknown fails closed
    assert _combine([]) is None                      # no steps -> unknown

    stage1 = make_leaf("scan", "stage_1", executed=True,
                       converged=_combine([True, True, True]))
    # middle step failed, final step converged -> combined is False.
    stage2 = make_leaf("scan", "stage_2", executed=True,
                       converged=_combine([True, False, True]))
    assert stage1.usable is True and stage2.usable is False

    truth = aggregate_workflow_truth([stage1, stage2], ["stage_1", "stage_2"])
    assert truth.scientific_status == "partial"
    assert any("stage_2" in r for r in truth.status_reasons)


# ---------------------------------------------------------------------------
# 3. FALSIFIER — path expected-segment miss / endpoint HEI (M09 + M59)
#    An endpoint-HEI raw diagram (segments=[]) would have been promoted to
#    success because a diagram exists; it must be partial.
# ---------------------------------------------------------------------------


def test_path_endpoint_hei_zero_segments_is_partial_not_success() -> None:
    from pdb2reaction.workflows.path_search import _path_leaves_and_expected

    # Endpoint-HEI branch returns segments=[] but a raw R/P diagram can be drawn.
    leaves, expected = _path_leaves_and_expected(
        [], raw_artifacts=["mep.pdb", "energy_diagram_MEP.png"]
    )
    truth = aggregate_workflow_truth(leaves, expected)
    assert truth.scientific_status == "partial"  # would have been "success"
    # The raw path is retained as a reportable-but-unusable diagnostic artifact.
    raw = [leaf for leaf in leaves if leaf.item_id == "raw_path"][0]
    assert raw.usable is False and raw.reason == "endpoint_hei"
    assert "mep.pdb" in raw.artifacts


def test_path_engine_nonconverged_endpoint_reason_retained() -> None:
    from pdb2reaction.workflows.path_search import _path_leaves_and_expected

    leaves, expected = _path_leaves_and_expected(
        [], raw_artifacts=["mep.pdb"], engine_converged=False
    )
    raw = [leaf for leaf in leaves if leaf.item_id == "raw_path"][0]
    assert "endpoint_hei" in raw.reason and "engine_nonconverged" in raw.reason


def test_path_single_reactive_segment_is_success() -> None:
    from pdb2reaction.workflows.path_search import SegmentReport, _path_leaves_and_expected

    seg = SegmentReport(tag="seg_01", barrier_kcal=10.0, delta_kcal=-2.0,
                        summary="", kind="seg", seg_index=1, converged=True)
    leaves, expected = _path_leaves_and_expected([seg])
    truth = aggregate_workflow_truth(leaves, expected)
    assert truth.scientific_status == "success"  # legacy behavior preserved


def test_path_nonconverged_reactive_segment_is_not_success() -> None:
    # M09: a reactive segment whose StringOptimizer hit its cycle limit
    # (is_converged=False) writes its trajectory but must NOT count toward
    # completeness. The pre-C6 _path_leaves_and_expected hardcoded converged=True
    # and would have reported success; the real threading now reads the field.
    from pdb2reaction.workflows.path_search import SegmentReport, _path_leaves_and_expected

    seg = SegmentReport(tag="seg_01", barrier_kcal=10.0, delta_kcal=-2.0,
                        summary="", kind="seg", seg_index=1, converged=False)
    leaves, expected = _path_leaves_and_expected([seg])
    truth = aggregate_workflow_truth(leaves, expected)
    assert truth.scientific_status != "success"
    seg_leaf = [leaf for leaf in leaves if leaf.item_id == "segment_1"][0]
    assert seg_leaf.usable is False and seg_leaf.reason == "not_converged"

    # Convergence-unknown (None: no readable signal) also fails closed.
    seg_unk = SegmentReport(tag="seg_01", barrier_kcal=10.0, delta_kcal=-2.0,
                            summary="", kind="seg", seg_index=1, converged=None)
    leaves_u, expected_u = _path_leaves_and_expected([seg_unk])
    assert aggregate_workflow_truth(leaves_u, expected_u).scientific_status != "success"


def test_path_bridge_only_is_not_success() -> None:
    from pdb2reaction.workflows.path_search import SegmentReport, _path_leaves_and_expected

    bridge = SegmentReport(tag="bridge_01", barrier_kcal=0.0, delta_kcal=0.0,
                           summary="", kind="bridge", seg_index=1)
    leaves, expected = _path_leaves_and_expected([bridge], raw_artifacts=["mep.pdb"])
    truth = aggregate_workflow_truth(leaves, expected)
    # A path made only of a non-reactive bridge has no usable reactive segment.
    assert truth.scientific_status != "success"


# ---------------------------------------------------------------------------
# 4. FALSIFIER — IRC directional non-convergence (M42)
#    A nonconverged direction (whose trajectory + Hessian still exist) must not
#    be promoted; only the converged direction is usable.
# ---------------------------------------------------------------------------


def _irc_direction_leaves(*, forward, forward_conv, backward, backward_conv,
                          n_fwd=10, n_bwd=10):
    """Exercise _outcomes.irc_direction_leaves, the directional-leaf builder that
    irc.py calls in production (bind to the shared helper, don't copy it)."""
    from pdb2reaction.workflows._outcomes import irc_direction_leaves
    return irc_direction_leaves(
        (
            ("forward", bool(forward), forward_conv, n_fwd,
             ["forward_irc.pdb"] if forward else []),
            ("backward", bool(backward), backward_conv, n_bwd,
             ["backward_irc.pdb"] if backward else []),
        )
    )


def test_irc_forward_converged_backward_not_is_partial() -> None:
    # Both directions requested and both trajectories/Hessians exist, but only
    # forward converged.
    leaves, expected = _irc_direction_leaves(
        forward=True, forward_conv=True, backward=True, backward_conv=False
    )
    truth = aggregate_workflow_truth(leaves, expected)
    assert truth.scientific_status == "partial"
    backward = [leaf for leaf in leaves if leaf.item_id == "backward"][0]
    assert backward.usable is False  # not promoted despite artifacts existing


def test_irc_one_sided_request_succeeds() -> None:
    # Backward explicitly disabled: its absence is optional, not a failure.
    leaves, expected = _irc_direction_leaves(
        forward=True, forward_conv=True, backward=False, backward_conv=None
    )
    truth = aggregate_workflow_truth(leaves, expected)
    assert truth.scientific_status == "success"
    assert expected == ["forward"]


def test_irc_convergence_attribute_absent_fails_closed() -> None:
    # A missing convergence attribute (None) must fail closed, not read as True.
    leaves, expected = _irc_direction_leaves(
        forward=True, forward_conv=None, backward=False, backward_conv=None
    )
    truth = aggregate_workflow_truth(leaves, expected)
    assert truth.scientific_status != "success"


def test_irc_hessian_cache_gate_condition() -> None:
    # Binds to the PRODUCTION gate irc.py uses to decide endpoint-Hessian caching
    # (irc_hessian_cache_eligible == `getattr(obj, attr, None) is True`).
    from pdb2reaction.workflows._outcomes import irc_hessian_cache_eligible

    class _FakeEuler:
        forward_is_converged = False
        backward_is_converged = True

    e = _FakeEuler()
    # forward: nonconverged -> the Hessian key must NOT be cached.
    assert irc_hessian_cache_eligible(e, "forward_is_converged") is False
    # backward: converged -> eligible for caching.
    assert irc_hessian_cache_eligible(e, "backward_is_converged") is True
    # An absent attribute also fails closed.
    assert irc_hessian_cache_eligible(e, "never_ran_is_converged") is False
    # A non-boolean truthy value (e.g. the integer 1) is not convergence.
    e.forward_is_converged = 1  # type: ignore[assignment]
    assert irc_hessian_cache_eligible(e, "forward_is_converged") is False


# ---------------------------------------------------------------------------
# 5. FALSIFIER — FREQ/DFT name/number fallback (M28)
#    A nonzero freq exit must not be treated as success just because a
#    thermoanalysis.yaml with finite fields exists on disk.
# ---------------------------------------------------------------------------


def test_freq_nonzero_exit_is_not_usable_despite_finite_yaml(tmp_path: Path, monkeypatch) -> None:
    from pdb2reaction.workflows import all as all_workflow

    fdir = tmp_path / "R"
    fdir.mkdir(parents=True, exist_ok=True)
    # A thermoanalysis.yaml with finite thermochemistry exists (e.g. a prior run
    # or a partial write) — under the old code its finite numbers would be used.
    (fdir / "thermoanalysis.yaml").write_text(
        "sum_EE_and_thermal_free_energy_ha: -123.456\n"
        "thermal_correction_free_energy_ha: 0.05\n"
    )
    pdb = tmp_path / "R.xyz"
    pdb.write_text("1\n\nH 0 0 0\n")

    # Freq child exits nonzero.
    monkeypatch.setattr(all_workflow, "_run_cli_main", lambda *a, **k: 1)
    monkeypatch.setattr(all_workflow, "_echo", lambda *a, **k: None)

    out = all_workflow._run_freq_for_state(
        pdb, 0, 1, fdir, None, False, None, False, overrides={}
    )
    # The finite YAML is NOT promoted; thermochemistry is unusable.
    assert out == {}


def test_freq_zero_exit_returns_parsed_thermo(tmp_path: Path, monkeypatch) -> None:
    from pdb2reaction.workflows import all as all_workflow

    fdir = tmp_path / "R"
    fdir.mkdir(parents=True, exist_ok=True)
    (fdir / "thermoanalysis.yaml").write_text(
        "sum_EE_and_thermal_free_energy_ha: -123.456\n"
    )
    pdb = tmp_path / "R.xyz"
    pdb.write_text("1\n\nH 0 0 0\n")

    # Freq child succeeds (exit 0).
    monkeypatch.setattr(all_workflow, "_run_cli_main", lambda *a, **k: 0)
    monkeypatch.setattr(all_workflow, "_echo", lambda *a, **k: None)

    out = all_workflow._run_freq_for_state(
        pdb, 0, 1, fdir, None, False, None, False, overrides={}
    )
    assert out.get("sum_EE_and_thermal_free_energy_ha") == pytest.approx(-123.456)


def test_freq_no_dump_does_not_consume_stale_thermo(tmp_path: Path, monkeypatch) -> None:
    from pdb2reaction.workflows import all as all_workflow

    fdir = tmp_path / "R"
    fdir.mkdir(parents=True)
    (fdir / "thermoanalysis.yaml").write_text(
        "sum_EE_and_thermal_free_energy_ha: -123.456\n",
        encoding="utf-8",
    )
    structure = tmp_path / "R.xyz"
    structure.write_text("1\n\nH 0 0 0\n", encoding="utf-8")
    monkeypatch.setattr(all_workflow, "_run_cli_main", lambda *a, **k: 0)

    assert all_workflow._run_freq_for_state(
        structure,
        0,
        1,
        fdir,
        None,
        False,
        None,
        False,
        overrides={"dump": False},
    ) == {}


def test_freq_preserves_explicit_default_runtime_values_over_yaml(
    tmp_path: Path, monkeypatch
) -> None:
    from pdb2reaction.workflows import all as all_workflow

    pdb = tmp_path / "R.xyz"
    pdb.write_text("1\n\nH 0 0 0\n", encoding="utf-8")
    config = tmp_path / "runtime.yaml"
    config.write_text(
        "calc:\n  workers: 8\n  workers_per_node: 4\n"
        "  backend: orb\n  solvent: water\n  solvent_model: cpcm\n",
        encoding="utf-8",
    )
    captured: list[str] = []

    def _capture(_name, _command, argv, **_kwargs):
        captured.extend(argv)
        return 1

    monkeypatch.setattr(all_workflow, "_run_cli_main", _capture)
    monkeypatch.setattr(all_workflow, "_echo", lambda *a, **k: None)

    all_workflow._run_freq_for_state(
        pdb,
        0,
        1,
        tmp_path / "freq",
        config,
        False,
        None,
        False,
        overrides={
            "workers": 1,
            "workers_per_node": 1,
            "backend": "uma",
            "solvent": "none",
            "solvent_model": "alpb",
        },
    )

    for flag, value in (
        ("--workers", "1"),
        ("--workers-per-node", "1"),
        ("--backend", "uma"),
        ("--solvent", "none"),
        ("--solvent-model", "alpb"),
    ):
        assert captured[captured.index(flag) + 1] == value
    assert captured[captured.index("--config") + 1] == str(config)


@pytest.mark.parametrize("overrides, expected_count", [({}, 0), ({"symmetry_number": 3}, 1)])
def test_all_freq_forwards_symmetry_number_only_when_overridden(
    tmp_path: Path,
    monkeypatch,
    overrides: dict,
    expected_count: int,
) -> None:
    from pdb2reaction.workflows import all as all_workflow

    structure = tmp_path / "R.xyz"
    structure.write_text("1\n\nH 0 0 0\n", encoding="utf-8")
    captured: list[str] = []

    def _capture(_name, _command, argv, **_kwargs):
        captured.extend(argv)
        return 1

    monkeypatch.setattr(all_workflow, "_run_cli_main", _capture)
    monkeypatch.setattr(all_workflow, "_echo", lambda *a, **k: None)

    all_workflow._run_freq_for_state(
        structure,
        0,
        1,
        tmp_path / "freq",
        None,
        False,
        None,
        False,
        overrides=overrides,
    )

    assert captured.count("--symmetry-number") == expected_count
    if expected_count:
        assert captured[captured.index("--symmetry-number") + 1] == "3"


# ---------------------------------------------------------------------------
# 6. LEGACY-COMPAT — a genuinely converged leaf's public output is unchanged
#    except for the additive new fields.
# ---------------------------------------------------------------------------


def test_legacy_converged_output_is_byte_compatible(tmp_path: Path) -> None:
    from pdb2reaction.core.utils import write_result_json

    # A representative subset of the legacy scan2d result payload for a fully
# converged run.
    legacy = {
        "status": "completed",
        "charge": 0,
        "spin": 1,
        "min_energy_hartree": -1.2345,
        "n_grid_points": 9,
        "files": {"surface_csv": "surface.csv"},
    }

    # Post-C6, the additive outcome fields are appended; every legacy key must be
    # bit-identical.
    points = [
        make_scan_point(f"p{i}", executed=True, converged=True, energy=-1.0, artifact_written=True)
        for i in range(9)
    ]
    sci, reasons = scan_scientific_status(points)
    assert sci == "success"

    enriched = dict(legacy)
    attach_outcomes(enriched, point_outcomes=points, scientific_status=sci,
                    scientific_status_reasons=reasons)

    write_result_json(tmp_path, dict(enriched), command="scan2d")
    r = json.loads((tmp_path / "result.json").read_text())

    for key, value in legacy.items():
        assert r[key] == value, f"legacy key {key} changed"
    # Only additive keys are new; scientific_status reports the truthful success.
    assert r["scientific_status"] == "success"
    assert "scientific_status_reasons" not in r  # no reasons on a clean success
    new_keys = set(r) - set(legacy)
    # The only new keys are additive outcome / envelope fields.
    assert "point_outcomes" in new_keys and "scientific_status" in new_keys


def test_eligible_points_helper() -> None:
    pts = [
        make_scan_point("a", executed=True, converged=True, energy=-1.0, artifact_written=True),
        make_scan_point("b", executed=True, converged=False, energy=-2.0, artifact_written=True),
    ]
    assert [p.point_id for p in eligible_points(pts)] == ["a"]


# ---------------------------------------------------------------------------
# 7. FALSIFIER — DMF (path-opt) non-convergence via IPOPT status (M09)
#    Binds path_opt._run_dmf_mep's status parsing + the DMF LeafOutcome. A
#    nonconverged solve (status 2) retains its trajectory but is unusable; a
#    converged solve (status 0/1) is a usable success.
# ---------------------------------------------------------------------------


def test_dmf_ipopt_status_maps_to_convergence() -> None:
    from pdb2reaction.workflows._outcomes import ipopt_status_to_converged

    assert ipopt_status_to_converged(0) == (True, "ipopt_converged")
    assert ipopt_status_to_converged(1) == (True, "ipopt_converged")
    assert ipopt_status_to_converged(2) == (False, "ipopt_status_2")
    assert ipopt_status_to_converged(-1) == (False, "ipopt_status_-1")
    # A missing/unreadable status fails closed to unknown (never converged).
    assert ipopt_status_to_converged(None) == (None, "convergence_unknown")


def test_dmf_nonconverged_leaf_unusable_artifact_retained() -> None:
    # path_opt.cli builds `_mk_leaf("path-opt", "dmf_mep", converged=_dmf_converged,
    # artifacts=["final_geometries_trj.xyz"], reason=_dmf_reason)`.
    from pdb2reaction.workflows._outcomes import ipopt_status_to_converged

    conv, reason = ipopt_status_to_converged(2)
    leaf = make_leaf("path-opt", "dmf_mep", executed=True, converged=conv,
                     artifacts=["final_geometries_trj.xyz"], reason=reason)
    assert leaf.usable is False                    # not promoted by artifact
    assert "final_geometries_trj.xyz" in leaf.artifacts  # artifact retained
    # Unusable required leaf whose trajectory is retained -> partial (a reportable
    # diagnostic), never success.
    assert aggregate_workflow_truth([leaf], ["dmf_mep"]).scientific_status == "partial"


def test_dmf_converged_leaf_is_success() -> None:
    from pdb2reaction.workflows._outcomes import ipopt_status_to_converged

    conv, reason = ipopt_status_to_converged(0)
    leaf = make_leaf("path-opt", "dmf_mep", executed=True, converged=conv,
                     artifacts=["final_geometries_trj.xyz"], reason=reason)
    assert leaf.usable is True
    assert aggregate_workflow_truth([leaf], ["dmf_mep"]).scientific_status == "success"


def test_dmf_result_status_is_legacy_byte_compatible() -> None:
    # The pre-C6 DMFMepResult exposed no readable is_converged, so the legacy DMF
    # result.json `status` field always read "completed" and its `converged`
    # field always read null. The C6 convergence truth is carried by
    # scientific_status / stage_outcomes; NEITHER legacy field may flip on a
    # genuinely-converged run (a non-additive value change would break byte-compat).
    from pdb2reaction.workflows import path_opt

    src = Path(path_opt.__file__).read_text(encoding="utf-8")
    # The DMF result_data must pin the legacy status literal, not derive it from
    # _dmf_converged (which would flip "completed"->"converged").
    assert '"status": "completed",' in src
    assert '"status": "converged" if _dmf_converged' not in src
    # The DMF `converged` legacy field must stay pinned to null (its base value),
    # never the real IPOPT bit (which would flip null->true/false non-additively).
    assert '"converged": None,' in src
    assert '"converged": _dmf_converged' not in src
    # ...but the IPOPT truth must still feed the additive scientific_status leaf.
    assert 'converged=_dmf_converged' in src


# ---------------------------------------------------------------------------
# 8. FALSIFIER — M50 scan producer records truthful convergence
#    A fake optimizer that returns normally but reports is_converged=False must
#    record bias_converged=False, so a normal (non-raising) run is not mistaken
#    for convergence and the point is excluded from seeds.
# ---------------------------------------------------------------------------


def test_m50_producer_records_nonconvergence_from_optimizer() -> None:
    from pdb2reaction.workflows._outcomes import optimizer_converged_bit

    class _FakeOpt:
        def __init__(self, conv):
            self.is_converged = conv

    # scan2d/scan3d record `bias_converged = optimizer_converged_bit(opt)`; scan
    # records the same bit per step.
    assert optimizer_converged_bit(_FakeOpt(False)) is False
    assert optimizer_converged_bit(_FakeOpt(True)) is True
    # A non-boolean 1 (truthy but not convergence) collapses to unknown.
    assert optimizer_converged_bit(_FakeOpt(1)) is None
    # A missing attribute fails closed.

    class _NoAttr:
        pass

    assert optimizer_converged_bit(_NoAttr()) is None

    # The recorded bit drives seed eligibility: a normal-return-but-nonconverged
    # point is excluded even with a finite (lowest) energy.
    record = {"energy_hartree": -9.9, "bias_converged": optimizer_converged_bit(_FakeOpt(False))}
    assert seed_eligible_mask([record]) == [False]


# ---------------------------------------------------------------------------
# 9. FALSIFIER — the ALL-pipeline aggregate consumes truthful leaves (M42 + item 1)
#    A never_stop / max-cycle IRC (trajectory present, direction nonconverged)
#    must not yield scientific_status=success in the all-pipeline aggregate,
#    while the legacy `status` string is unchanged.
# ---------------------------------------------------------------------------


def test_read_irc_outcome_gates_on_scientific_status(tmp_path: Path) -> None:
    from pdb2reaction.workflows.all import _read_irc_outcome

    irc_dir = tmp_path / "irc"
    irc_dir.mkdir()
    # A converged both-direction IRC child result.
    (irc_dir / "result.json").write_text(json.dumps({
        "status": "completed",
        "scientific_status": "success",
        "forward_converged": True,
        "backward_converged": True,
        "files": {"finished_irc": "finished_irc_trj.xyz"},
    }))
    ok = _read_irc_outcome(irc_dir)
    assert ok["usable"] is True and ok["traj"] == "finished_irc_trj.xyz"

    # A backward direction that hit its cycle limit: trajectory still exists.
    (irc_dir / "result.json").write_text(json.dumps({
        "status": "completed",              # the IRC PROCESS ran (legacy)
        "scientific_status": "partial",     # but a direction did not converge
        "scientific_status_reasons": ["irc:backward:not_converged"],
        "forward_converged": True,
        "backward_converged": False,
        "files": {"finished_irc": "finished_irc_trj.xyz"},
    }))
    bad = _read_irc_outcome(irc_dir)
    assert bad["usable"] is False
    assert "backward" in bad["reason"]

    # A missing result.json fails closed.
    (irc_dir / "result.json").unlink()
    assert _read_irc_outcome(irc_dir)["usable"] is False


def test_all_pipeline_aggregate_excludes_nonconverged_irc() -> None:
    from pdb2reaction.workflows.all import _pipeline_aggregate_truth

    summary = {"segments": [{
        "index": 1, "kind": "seg", "barrier_kcal": 10.0, "converged": True,
    }]}
    config = {"tsopt": True, "thermo": False, "dft": False}

    # Legacy status="success" (a trajectory exists, so _derive_pipeline_status is
    # satisfied). The IRC leaf reports backward nonconverged -> the aggregate must
    # demote scientific_status while leaving the legacy axis unchanged.
    nonconverged = [{
        "index": 1,
        "irc_traj": "finished_irc_trj.xyz",
        "irc": {"usable": False, "reason": "irc:backward:not_converged",
                "traj": "finished_irc_trj.xyz"},
        "endpoint_opt": {"reactant_converged": True, "product_converged": True},
    }]
    truth = _pipeline_aggregate_truth(
        summary, post_segments=nonconverged, config=config, legacy_status="success",
    )
    assert truth.scientific_status != "success"           # would have been success
    assert any("segment_1" in r for r in truth.status_reasons)

    # Every requested IRC direction converged + endpoints converged -> success.
    converged = [{
        "index": 1,
        "irc_traj": "finished_irc_trj.xyz",
        "irc": {"usable": True, "reason": "ok", "traj": "finished_irc_trj.xyz"},
        "endpoint_assignment": {"connectivity_validated": True},
        "endpoint_opt": {"reactant_converged": True, "product_converged": True},
    }]
    truth_ok = _pipeline_aggregate_truth(
        summary, post_segments=converged, config=config, legacy_status="success",
    )
    assert truth_ok.scientific_status == "success"


def test_all_pipeline_tsopt_only_does_not_require_mep_convergence() -> None:
    """A direct-TS segment has no MEP convergence field to gate on."""
    from pdb2reaction.workflows.all import _pipeline_aggregate_truth

    summary = {"segments": [{"index": 1, "kind": "tsopt", "barrier_kcal": 10.0}]}
    post = [{
        "index": 1,
        "irc": {"usable": True, "reason": "ok"},
        "endpoint_assignment": {"connectivity_validated": True},
        "endpoint_opt": {"reactant_converged": True, "product_converged": True},
    }]
    truth = _pipeline_aggregate_truth(
        summary, post_segments=post, config={"tsopt": True},
        legacy_status="success",
    )
    assert truth.scientific_status == "success"
    assert "mep_convergence_unknown" not in truth.status_reasons


def test_all_pipeline_aggregate_gates_on_endpoint_opt() -> None:
    from pdb2reaction.workflows.all import _pipeline_aggregate_truth

    summary = {"segments": [{
        "index": 1, "kind": "seg", "barrier_kcal": 10.0, "converged": True,
    }]}
    config = {"tsopt": True}
    # IRC usable but the product endpoint optimization did not converge.
    post = [{
        "index": 1,
        "irc": {"usable": True, "reason": "ok"},
        "endpoint_opt": {"reactant_converged": True, "product_converged": False},
    }]
    truth = _pipeline_aggregate_truth(
        summary, post_segments=post, config=config, legacy_status="success",
    )
    assert truth.scientific_status != "success"


def test_all_pipeline_aggregate_preserves_legacy_severity() -> None:
    # The composed scientific_status is never LESS severe than the legacy axis:
    # a legacy `partial` (e.g. DFT failed) with a fully-converged IRC stays partial.
    from pdb2reaction.workflows.all import _pipeline_aggregate_truth

    summary = {"segments": [{"index": 1, "kind": "seg", "converged": True}]}
    post = [{
        "index": 1,
        "irc": {"usable": True, "reason": "ok"},
        "endpoint_assignment": {"connectivity_validated": True},
        "endpoint_opt": {"reactant_converged": True, "product_converged": True},
    }]
    truth = _pipeline_aggregate_truth(
        summary, post_segments=post, config={"tsopt": True, "dft": True},
        legacy_status="partial", legacy_reasons=["segment 1: DFT failed (TS)"],
    )
    assert truth.scientific_status == "partial"
    assert "segment 1: DFT failed (TS)" in truth.status_reasons


def test_all_pipeline_requires_validated_mep_irc_connectivity() -> None:
    from pdb2reaction.workflows.all import _pipeline_aggregate_truth

    summary = {
        "segments": [
            {"index": 1, "kind": "seg", "converged": True}
        ]
    }
    base = {
        "index": 1,
        "irc": {"usable": True, "reason": "ok"},
        "endpoint_opt": {
            "reactant_converged": True,
            "product_converged": True,
        },
    }
    unvalidated = {
        **base,
        "endpoint_assignment": {"connectivity_validated": False},
    }
    validated = {
        **base,
        "endpoint_assignment": {"connectivity_validated": True},
    }

    assert _pipeline_aggregate_truth(
        summary,
        post_segments=[unvalidated],
        config={"tsopt": True},
        legacy_status="success",
    ).scientific_status != "success"
    assert _pipeline_aggregate_truth(
        summary,
        post_segments=[validated],
        config={"tsopt": True},
        legacy_status="success",
    ).scientific_status == "success"


def test_all_pipeline_aggregate_post_missing_fails_closed_when_tsopt_requested() -> None:
    # item 1 FALSIFIER (shipped-artifact fail-open): an intermediate MEP summary
    # (post_segments not yet assembled) with tsopt requested must NOT be promoted
    # to success on the MEP trajectory's existence alone. The reactive leaf fails
    # closed (post_missing) until IRC/endpoint post-processing actually runs.
    from pdb2reaction.workflows.all import _pipeline_aggregate_truth

    summary = {"segments": [{"index": 1, "kind": "seg", "barrier_kcal": 10.0}]}
    # post_segments=[] (post ran but produced no record for this segment).
    truth = _pipeline_aggregate_truth(
        summary, post_segments=[], config={"tsopt": True}, legacy_status="success",
    )
    assert truth.scientific_status != "success"        # would have been success (fail-open)
    assert any("segment_1" in r for r in truth.status_reasons)
    # post_segments=None (the very first intermediate write, before post-processing)
    # fails closed too.
    truth_none = _pipeline_aggregate_truth(
        summary, post_segments=None, config={"tsopt": True}, legacy_status="success",
    )
    assert truth_none.scientific_status != "success"


def test_all_pipeline_aggregate_no_tsopt_uses_segment_converged() -> None:
    # item 1 FALSIFIER (no convergence gate on the no-tsopt final summary): a
    # path-only final summary must gate on the segment's OWN reported convergence,
    # never a silent default-True.
    from pdb2reaction.workflows.all import _pipeline_aggregate_truth

    cfg = {"tsopt": False}
    # A nonconverged segment must not be a success even with a barrier band.
    nc = {"segments": [{"index": 1, "kind": "seg", "barrier_kcal": 10.0, "converged": False}]}
    assert _pipeline_aggregate_truth(
        nc, post_segments=None, config=cfg, legacy_status="success",
    ).scientific_status != "success"
    # A missing convergence field fails closed (unknown), not success.
    unk = {"segments": [{"index": 1, "kind": "seg", "barrier_kcal": 10.0}]}
    assert _pipeline_aggregate_truth(
        unk, post_segments=None, config=cfg, legacy_status="success",
    ).scientific_status != "success"
    # A genuinely-converged segment IS a success (legacy behavior preserved).
    ok = {"segments": [{"index": 1, "kind": "seg", "barrier_kcal": 10.0, "converged": True}]}
    assert _pipeline_aggregate_truth(
        ok, post_segments=None, config=cfg, legacy_status="success",
    ).scientific_status == "success"


@pytest.mark.parametrize("mep_converged", [False, None])
def test_all_pipeline_post_success_cannot_promote_bad_mep(
    mep_converged: bool | None,
) -> None:
    """Successful IRC/endpoints cannot overwrite false/unknown MEP truth."""
    from pdb2reaction.workflows.all import _pipeline_aggregate_truth

    summary = {"segments": [{
        "index": 1, "kind": "seg", "barrier_kcal": 10.0,
        "converged": mep_converged,
    }]}
    post = [{
        "index": 1,
        "irc": {"usable": True, "reason": "ok"},
        "endpoint_opt": {"reactant_converged": True, "product_converged": True},
    }]
    truth = _pipeline_aggregate_truth(
        summary, post_segments=post, config={"tsopt": True},
        legacy_status="success",
    )
    assert truth.scientific_status != "success"


def test_read_path_opt_segment_converged_is_tristate(tmp_path: Path) -> None:
    from pdb2reaction.workflows.all import _read_path_opt_segment_converged

    assert _read_path_opt_segment_converged(tmp_path) is None
    result = tmp_path / "result.json"
    result.write_text(json.dumps({"stage_outcomes": [{"converged": True}]}))
    assert _read_path_opt_segment_converged(tmp_path) is True
    result.write_text(json.dumps({"stage_outcomes": [{"converged": False}]}))
    assert _read_path_opt_segment_converged(tmp_path) is False
    result.write_text("{ not json")
    assert _read_path_opt_segment_converged(tmp_path) is None


def test_all_path_opt_child_emits_machine_result() -> None:
    """The no-refine path-opt producer must enable its result.json contract."""
    from pdb2reaction.workflows import all as all_workflow

    source = Path(all_workflow.__file__).read_text(encoding="utf-8")
    branch = source[source.index("if not refine_path:"):source.index("final_trj = path_dir", source.index("if not refine_path:"))]
    assert '"--out-json"' in branch
    assert 'seg_dir / "result.json"' in branch
