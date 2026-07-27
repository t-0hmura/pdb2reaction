"""Unit tests for pdb2reaction.workflows._all_helpers."""

from __future__ import annotations

import tempfile
from pathlib import Path

import click
import pytest
import yaml

from pdb2reaction.workflows._all_helpers import (
    build_energy_level_dict,
    build_pipeline_summary_payload,
    build_thermo_mode_validation,
    build_thermo_symmetry_provenance,
    validated_thermo_triplet,
)


def test_tsopt_result_validation_requires_verified_first_order_saddle() -> None:
    from pdb2reaction.workflows.all import _validate_tsopt_result_payload

    _validate_tsopt_result_payload(
        {"status": "converged", "n_imaginary_modes": 1}
    )

    for payload in (
        {"status": "not_converged", "n_imaginary_modes": 1},
        {"status": "converged", "n_imaginary_modes": 0},
        {"status": "converged", "n_imaginary_modes": 2},
        {"status": "unknown"},
    ):
        with pytest.raises(click.ClickException, match="IRC was not started"):
            _validate_tsopt_result_payload(payload)


def test_all_stops_before_irc_when_tsopt_result_is_invalid(
    tmp_path: Path, monkeypatch,
) -> None:
    import json

    from pdb2reaction.workflows import all as all_workflow

    hei = tmp_path / "hei.xyz"
    hei.write_text("1\nHEI\nH 0.0 0.0 0.0\n", encoding="utf-8")
    captured_args: list[str] = []

    def fake_run_cli_main(_name, _cli, args, **_kwargs) -> None:
        captured_args.extend(args)
        ts_dir = Path(args[args.index("--out-dir") + 1])
        (ts_dir / "result.json").write_text(
            json.dumps(
                {"status": "not_converged", "n_imaginary_modes": 2}
            ),
            encoding="utf-8",
        )

    monkeypatch.setattr(all_workflow, "_run_cli_main", fake_run_cli_main)
    monkeypatch.setattr(all_workflow, "_echo_detail", lambda *_a, **_k: None)

    with pytest.raises(click.ClickException, match="IRC was not started"):
        all_workflow._run_tsopt_on_hei(
            hei,
            charge=0,
            spin=1,
            calc_cfg={"backend": "uma"},
            args_yaml=None,
            out_dir=tmp_path / "segment",
            freeze_links=False,
            opt_mode_default="hess",
            ref_pdb=None,
            convert_files=False,
        )

    assert "--out-json" in captured_args


def test_all_tr_projection_is_injected_into_child_config() -> None:
    from pdb2reaction.workflows.all import _write_args_yaml_with_freeze_atoms
    from pdb2reaction.workflows._run_session import RunSession

    session = RunSession()
    path = _write_args_yaml_with_freeze_atoms(
        None, [], tr_projection="legacy-active", session=session
    )
    assert path is not None
    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    assert payload["geom"]["tr_projection"] == "legacy-active"
    session.close()
    assert not path.exists()


def test_build_energy_level_dict_kcal_projection() -> None:
    d = build_energy_level_dict(
        labels=["R", "TS", "P"],
        energies_au=[-100.0, -99.5, -100.2],
        ref_energy=-100.0,
        au_to_kcal=627.5095,
        diagram_path="/tmp/diag.png",
        structures={"R": "r.pdb", "TS": "ts.pdb", "P": "p.pdb"},
    )
    assert d["labels"] == ["R", "TS", "P"]
    assert d["energies_kcal"][0] == 0.0
    assert abs(d["energies_kcal"][1] - 313.75475) < 1e-3
    assert d["barrier_kcal"] == d["energies_kcal"][1]
    assert d["delta_kcal"] == d["energies_kcal"][-1]


def test_build_energy_level_dict_no_input_mutation() -> None:
    labels = ["R", "TS", "P"]
    structs = {"R": "r.pdb"}
    d = build_energy_level_dict(
        labels=labels,
        energies_au=[-1.0, -0.5, -1.1],
        ref_energy=-1.0,
        au_to_kcal=627.5,
        diagram_path="/tmp/x.png",
        structures=structs,
    )
    d["labels"].append("EXTRA")
    d["structures"]["NEW"] = "n.pdb"
    assert labels == ["R", "TS", "P"]
    assert "NEW" not in structs


def test_energy_diagram_metadata_does_not_claim_failed_png(
    tmp_path: Path, monkeypatch,
) -> None:
    from pdb2reaction.workflows import all as all_workflow

    class FailingFigure:
        def update_layout(self, **kwargs) -> None:
            pass

        def write_image(self, *args, **kwargs) -> None:
            raise RuntimeError("renderer unavailable")

    monkeypatch.setattr(all_workflow, "build_energy_diagram", lambda **kwargs: FailingFigure())
    payload = all_workflow._write_segment_energy_diagram(
        tmp_path / "energy_diagram_MLIP",
        labels=["R", "TS", "P"],
        energies_au=[-1.0, -0.9, -1.1],
        title_note="test",
    )

    assert payload is not None
    assert payload["image"] is None
    assert payload["image_written"] is False
    assert "renderer unavailable" in payload["image_error"]


def test_build_pipeline_summary_payload_shape() -> None:
    with tempfile.TemporaryDirectory() as d:
        out_dir = Path(d) / "out"
        path_dir = Path(d) / "path"
        out_dir.mkdir()
        path_dir.mkdir()
        summary = {
            "n_images": 5,
            "n_segments": 1,
            "segments": [{"kind": "seg", "bond_changes": "A->B"}],
            "energy_diagrams": [{"name": "MEP", "x": [0, 1]}],
        }
        payload = build_pipeline_summary_payload(
            out_dir=out_dir,
            path_dir=path_dir,
            summary=summary,
            refine_path=True,
            do_tsopt=True,
            do_thermo=False,
            do_dft=True,
            dft_func_basis_use="wb97m-v/def2-tzvpd",
            opt_mode="GRAD",
            opt_mode_post="HESS",
            path_opt_mode="GRAD",
            post_opt_mode="HESS",
            ts_opt_mode="HESS",
            endpoint_opt_mode="GRAD",
            mep_mode_kind="dmf",
            dmf_correlated=True,
            mlip_backend="mace",
            mlip_model="mace-off23-small",
            mlip_precision="fp64",
            command_str="pdb2reaction all -i foo.pdb",
            q_int=-1,
            spin=1,
            freeze_atoms=[10, 20, 30],
            post_segment_logs=[{"seg": 1, "status": "ok"}],
        )
    assert payload["pipeline_mode"] == "path-search"
    assert payload["dft_func_basis"] == "wb97m-v/def2-tzvpd"
    assert payload["opt_mode"] == "grad"
    assert payload["opt_mode_post"] == "hess"
    assert payload["path_opt_mode"] == "grad"
    assert payload["post_opt_mode"] == "hess"
    assert payload["ts_opt_mode"] == "hess"
    assert payload["endpoint_opt_mode"] == "grad"
    assert payload["mep_mode"] == "dmf"
    assert payload["dmf_correlated"] is True
    assert payload["mlip_backend"] == "mace"
    assert payload["mlip_model"] == "mace-off23-small"
    assert payload["mlip_precision"] == "fp64"
    assert "uma_model" not in payload
    assert payload["charge"] == -1
    assert payload["spin"] == 1
    assert payload["freeze_atoms"] == [10, 20, 30]
    assert payload["mep"]["diagram"]["name"] == "MEP"


def test_build_pipeline_summary_payload_dft_disabled_drops_basis() -> None:
    payload = build_pipeline_summary_payload(
        out_dir=Path("."),
        path_dir=Path("."),
        summary={},
        refine_path=False,
        do_tsopt=False,
        do_thermo=False,
        do_dft=False,
        dft_func_basis_use="wb97m-v/def2-tzvpd",
        opt_mode=None,
        opt_mode_post=None,
        path_opt_mode=None,
        post_opt_mode=None,
        ts_opt_mode=None,
        endpoint_opt_mode=None,
        mep_mode_kind=None,
        dmf_correlated=False,
        mlip_backend="orb",
        mlip_model=None,
        mlip_precision="fp64",
        command_str="",
        q_int=0,
        spin=1,
        freeze_atoms=[],
        post_segment_logs=[],
    )
    assert payload["pipeline_mode"] == "path-opt"
    assert payload["dft_func_basis"] is None
    assert payload["opt_mode"] is None


def test_thermo_symmetry_provenance_copies_only_complete_valid_states() -> None:
    payload = build_thermo_symmetry_provenance(
        {
            "R": {"symmetry_number": 2, "symmetry_number_source": "config"},
            "TS": {"symmetry_number": 3, "symmetry_number_source": "cli"},
            "P": {"symmetry_number": 1},
            "other": {"symmetry_number": 9, "symmetry_number_source": "test"},
        }
    )

    assert payload == {
        "R": {"symmetry_number": 2, "symmetry_number_source": "config"},
        "TS": {"symmetry_number": 3, "symmetry_number_source": "cli"},
    }


@pytest.mark.parametrize("invalid", [True, 0, -1, 2.0, "2", None])
def test_thermo_symmetry_provenance_rejects_invalid_values(invalid) -> None:
    assert build_thermo_symmetry_provenance(
        {
            "R": {
                "symmetry_number": invalid,
                "symmetry_number_source": "default",
            }
        }
    ) == {}


def test_thermo_mode_validation_requires_minimum_ts_minimum_order() -> None:
    result = build_thermo_mode_validation(
        {
            "R": {"num_imag_freq": 0},
            "TS": {"num_imag_freq": 1},
            "P": {"num_imag_freq": 0},
        }
    )
    assert result["valid"] is True
    assert result["reasons"] == []


def test_thermo_mode_validation_rejects_frozen_legacy_active_projection() -> None:
    payloads = {
        "R": {"num_imag_freq": 0},
        "TS": {
            "num_imag_freq": 1,
            "n_freeze_atoms": 2,
            "rigid_projection": {
                "treatment": "legacy-active",
                "frozen_atom_count": 2,
            },
        },
        "P": {"num_imag_freq": 0},
    }

    result = build_thermo_mode_validation(payloads)

    assert result["valid"] is False
    assert any("cannot certify" in reason for reason in result["reasons"])


@pytest.mark.parametrize(
    "payloads",
    [
        {"R": {"num_imag_freq": 1}, "TS": {"num_imag_freq": 1}, "P": {"num_imag_freq": 0}},
        {"R": {"num_imag_freq": 0}, "TS": {"num_imag_freq": 0}, "P": {"num_imag_freq": 0}},
        {"R": {"num_imag_freq": 0}, "TS": {"num_imag_freq": 1}, "P": {}},
    ],
)
def test_thermo_mode_validation_fails_closed(payloads) -> None:
    result = build_thermo_mode_validation(payloads)
    assert result["valid"] is False
    assert result["reasons"]


def test_validated_thermo_triplet_requires_modes_and_finite_values() -> None:
    payloads = {
        "R": {"num_imag_freq": 0, "g": -10.0},
        "TS": {"num_imag_freq": 1, "g": -9.0},
        "P": {"num_imag_freq": 0, "g": -11.0},
    }
    assert validated_thermo_triplet(payloads, "g") == (-10.0, -9.0, -11.0)
    payloads["P"]["g"] = float("nan")
    assert validated_thermo_triplet(payloads, "g") is None
    payloads["P"] = {"num_imag_freq": 1, "g": -11.0}
    assert validated_thermo_triplet(payloads, "g") is None


def test_enrich_summary_uses_backend_neutral_refined_energy_keys(tmp_path: Path) -> None:
    from pdb2reaction.workflows.all import _enrich_summary

    summary = {
        "segments": [{"index": 1, "kind": "seg", "barrier_kcal": 5.0}],
        "energy_diagrams": [{"name": "energy_diagram_G_MLIP_all"}],
    }
    result = _enrich_summary(
        summary,
        version="",
        pipeline_mode="path-opt",
        mlip_backend="mace",
        mlip_model="MACE-OMOL-0",
        charge=0,
        spin=1,
        post_segments=[
            {
                "index": 1,
                "gibbs_mlip": {"barrier_kcal": 12.5, "delta_kcal": -1.0},
                "endpoint_opt": {"reactant_converged": True},
            }
        ],
        config={
            "tsopt": True,
            "thermo": True,
            "path_opt_mode": "grad",
            "post_opt_mode": "hess",
            "ts_opt_mode": "hess",
            "endpoint_opt_mode": "grad",
            "mep_mode": "dmf",
            "dmf_correlated": True,
        },
        out_dir=tmp_path,
    )

    assert result["rate_limiting_step"]["barrier_kcal"] == 12.5
    assert result["rate_limiting_step"]["method"] == "MLIP_Gibbs"
    assert result["schema_version"]
    assert any(ref["method"] == "pdb2reaction" for ref in result["references"])
    assert any(
        ref["method"] == "Direct Max Flux (DMF)"
        for ref in result["references"]
    )
    assert any(
        ref["method"] == "quasi-RRHO thermochemistry"
        for ref in result["references"]
    )
    assert any(
        ref["method"] == "Correlated FB-ENM (CFB-ENM) initialization for DMF"
        for ref in result["references"]
    )
    assert any(
        ref["method"] == "Limited-memory BFGS (L-BFGS)"
        for ref in result["references"]
    )


def test_aggregate_summary_writer_injects_current_run_id_atomically(
    tmp_path: Path, monkeypatch
) -> None:
    import json

    from pdb2reaction.core.result_commit import RUN_ID_ENV
    from pdb2reaction.workflows.all import _write_summary_json

    monkeypatch.setenv(RUN_ID_ENV, "aggregate-current")
    summary = {"status": "success", "segments": []}

    path = _write_summary_json(tmp_path / "summary.json", summary)

    assert path == tmp_path / "summary.json"
    assert json.loads(path.read_text(encoding="utf-8"))["run_id"] == "aggregate-current"
    assert summary["run_id"] == "aggregate-current"


def test_enrich_summary_labels_overall_reaction_energy_method(tmp_path: Path) -> None:
    from pdb2reaction.workflows.all import _enrich_summary

    result = _enrich_summary(
        {
            "segments": [{"index": 1, "kind": "seg", "barrier_kcal": 5.0}],
            "energy_diagrams": [
                {
                    "name": "energy_diagram_G_DFT_plus_MLIP_all",
                    "energies_kcal": [0.0, 5.0, -2.5],
                }
            ],
        },
        version="",
        pipeline_mode="path-search",
        mlip_backend="orb",
        mlip_model="orb-v3",
        charge=0,
        spin=1,
        post_segments=[{
            "index": 1,
            "gibbs_dft_mlip": {"barrier_kcal": 5.0},
        }],
        out_dir=tmp_path,
    )

    assert result["overall_reaction_energy_kcal"] == -2.5
    assert result["overall_reaction_energy_method"] == "DFT//MLIP_Gibbs"


def test_tsonly_reenrichment_uses_local_diagram_and_post_method(
    tmp_path: Path,
) -> None:
    from pdb2reaction.workflows.all import _enrich_summary

    summary = {
        "segments": [
            {"index": 1, "kind": "tsopt", "barrier_kcal": 8.0},
        ],
        "energy_diagrams": [
            {
                "name": "energy_diagram_DFT",
                "energies_kcal": [0.0, 7.5, -1.25],
            }
        ],
    }
    common = {
        "version": "",
        "pipeline_mode": "tsopt-only",
        "mlip_backend": "mace",
        "mlip_model": "MACE-OMOL-0",
        "charge": 0,
        "spin": 1,
        "out_dir": tmp_path,
    }
    _enrich_summary(summary, **common)
    assert summary["rate_limiting_step"]["method"] == "MLIP"

    _enrich_summary(
        summary,
        post_segments=[
            {
                "index": 1,
                "dft": {"barrier_kcal": 7.5, "delta_kcal": -1.25},
            }
        ],
        **common,
    )

    assert summary["rate_limiting_step"] == {
        "segment": 1,
        "barrier_kcal": 7.5,
        "method": "DFT",
    }
    assert summary["overall_reaction_energy_kcal"] == -1.25
    assert summary["overall_reaction_energy_method"] == "DFT"


def test_enrich_summary_rate_limiting_uses_highest_complete_method(tmp_path: Path) -> None:
    """DFT names containing MLIP must not be downgraded by substring matching."""
    from pdb2reaction.workflows.all import _enrich_summary

    result = _enrich_summary(
        {
            "segments": [{"index": 1, "kind": "seg", "barrier_kcal": 10.0}],
            # Deliberately use a non-priority list order: method selection must
            # be exact and independent of diagram order or PNG export.
            "energy_diagrams": [
                {"name": "energy_diagram_G_DFT_plus_MLIP_all", "energies_kcal": [0, 41, -4]},
                {"name": "energy_diagram_MLIP_all", "energies_kcal": [0, 11, -1]},
                {"name": "energy_diagram_DFT_all", "energies_kcal": [0, 31, -3]},
            ],
        },
        version="",
        pipeline_mode="path-search",
        mlip_backend="orb",
        mlip_model="orb-v3",
        charge=0,
        spin=1,
        post_segments=[{
            "index": 1,
            "mlip": {"barrier_kcal": 11.0},
            "gibbs_mlip": {"barrier_kcal": 21.0},
            "dft": {"barrier_kcal": 31.0},
            "gibbs_dft_mlip": {"barrier_kcal": 41.0},
        }],
        out_dir=tmp_path,
    )

    assert result["rate_limiting_step"]["barrier_kcal"] == 41.0
    assert result["rate_limiting_step"]["method"] == "DFT//MLIP_Gibbs"
    assert result["overall_reaction_energy_method"] == "DFT//MLIP_Gibbs"


def test_enrich_summary_does_not_mix_partial_high_level_methods(tmp_path: Path) -> None:
    from pdb2reaction.workflows.all import _enrich_summary

    result = _enrich_summary(
        {
            "segments": [
                {"index": 1, "kind": "seg", "barrier_kcal": 8.0},
                {"index": 2, "kind": "seg", "barrier_kcal": 9.0},
            ],
            "energy_diagrams": [
                {"name": "energy_diagram_DFT_all", "energies_kcal": [0, 30, -3]},
                {"name": "energy_diagram_MLIP_all", "energies_kcal": [0, 12, -2]},
            ],
        },
        version="",
        pipeline_mode="path-search",
        mlip_backend="mace",
        mlip_model="MACE-OMOL-0",
        charge=0,
        spin=1,
        post_segments=[
            {"index": 1, "mlip": {"barrier_kcal": 11.0}, "dft": {"barrier_kcal": 30.0}},
            {"index": 2, "mlip": {"barrier_kcal": 12.0}},
        ],
        out_dir=tmp_path,
    )

    assert result["rate_limiting_step"]["barrier_kcal"] == 12.0
    assert result["rate_limiting_step"]["method"] == "MLIP"
    assert result["overall_reaction_energy_method"] == "MLIP"


def test_enrich_summary_marks_requested_dft_failure_partial(tmp_path: Path) -> None:
    from pdb2reaction.workflows.all import _enrich_summary

    result = _enrich_summary(
        {
            "segments": [{"index": 1, "kind": "seg", "barrier_kcal": 5.0}],
            "energy_diagrams": [{"name": "energy_diagram_MEP"}],
        },
        version="",
        pipeline_mode="path-search",
        mlip_backend="orb",
        mlip_model="orb-v3",
        charge=0,
        spin=1,
        post_segments=[
            {
                "index": 1,
                "dft": {"status": "failed", "failed_states": ["TS"]},
            }
        ],
        config={"tsopt": False, "thermo": False, "dft": True},
        out_dir=tmp_path,
    )

    assert result["status"] == "partial"
    assert result["status_reasons"] == ["segment 1: DFT failed (TS)"]


def test_enrich_summary_marks_non_saddle_thermochemistry_partial(tmp_path: Path) -> None:
    from pdb2reaction.workflows.all import _enrich_summary

    result = _enrich_summary(
        {
            "segments": [{"index": 1, "kind": "seg", "barrier_kcal": 5.0}],
            "energy_diagrams": [{"name": "energy_diagram_G_MLIP_all"}],
        },
        version="",
        pipeline_mode="path-search",
        mlip_backend="mace",
        mlip_model="MACE-OMOL-0",
        charge=0,
        spin=1,
        post_segments=[
            {
                "index": 1,
                "gibbs_mlip": {"barrier_kcal": 4.0},
                "ts_imag": {"n_imag": 0},
            }
        ],
        config={"tsopt": False, "thermo": True, "dft": False},
        out_dir=tmp_path,
    )

    assert result["status"] == "partial"
    assert "n_imag=0, expected 1" in result["status_reasons"][0]


def test_enrich_summary_success_requires_all_requested_post_results(tmp_path: Path) -> None:
    from pdb2reaction.workflows.all import _enrich_summary

    result = _enrich_summary(
        {
            "segments": [{"index": 1, "kind": "seg", "barrier_kcal": 5.0}],
            "energy_diagrams": [{"name": "energy_diagram_G_DFT_plus_MLIP_all"}],
        },
        version="",
        pipeline_mode="path-search",
        mlip_backend="uma",
        mlip_model="uma-s-1p2",
        charge=0,
        spin=1,
        post_segments=[
            {
                "index": 1,
                "mlip": {"barrier_kcal": 5.0},
                "irc_traj": "finished_irc_trj.xyz",
                "gibbs_mlip": {"barrier_kcal": 4.0},
                "ts_imag": {"n_imag": 1},
                "dft": {"barrier_kcal": 6.0},
                "gibbs_dft_mlip": {"barrier_kcal": 5.5},
            }
        ],
        config={"tsopt": True, "thermo": True, "dft": True},
        out_dir=tmp_path,
    )

    assert result["status"] == "success"
    assert "status_reasons" not in result
def test_irc_endpoint_topology_tie_uses_rmsd_and_records_provenance(
    monkeypatch,
    tmp_path,
):
    from pdb2reaction.workflows import all as workflow

    left = object()
    right = object()
    mep_left = object()
    mep_right = object()
    generated = iter((mep_left, mep_right))
    monkeypatch.setattr(
        workflow,
        "read_xyz_first_last",
        lambda _path: (["C"], [[0.0, 0.0, 0.0]], [[1.0, 0.0, 0.0]]),
    )
    monkeypatch.setattr(
        workflow,
        "_geom_from_angstrom",
        lambda *_args, **_kwargs: next(generated),
    )
    monkeypatch.setattr(
        workflow._path_search,
        "has_bond_change",
        lambda *_args, **_kwargs: (False, ""),
    )
    scores = {
        (left, mep_left): 4.0,
        (left, mep_right): 1.0,
        (right, mep_left): 1.0,
        (right, mep_right): 4.0,
    }
    monkeypatch.setattr(
        workflow._path_search,
        "_rmsd_between",
        lambda first, second: scores[(first, second)],
    )

    (
        oriented_left,
        oriented_right,
        _left_tag,
        _right_tag,
        reversed_irc,
        assignment,
    ) = workflow._orient_irc_endpoints(
        left,
        right,
        endpoint_trajectory=tmp_path / "mep.xyz",
        freeze_atoms=[],
        seg_tag="seg_01",
    )

    assert (oriented_left, oriented_right) == (right, left)
    assert reversed_irc is True
    assert assignment["method"] == "rmsd_topology_tie"
    assert assignment["rmsd_swapped"] < assignment["rmsd_direct"]
    assert assignment["connectivity_validated"] is True
