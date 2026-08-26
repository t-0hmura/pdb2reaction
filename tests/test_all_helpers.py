"""Unit tests for pdb2reaction.workflows._all_helpers and the ``workflows.all``
orchestration helpers."""

from __future__ import annotations

import tempfile
from pathlib import Path
from types import SimpleNamespace

import click
import pytest
import yaml

from pdb2reaction.workflows._all_helpers import (
    build_energy_level_dict,
    build_pipeline_summary_payload,
    build_thermo_symmetry_provenance,
)


@pytest.mark.parametrize(
    ("names", "expected"),
    [
        (("reactant.pdb", "product.pdb"), ".pdb"),
        (("reactant.cif", "product.cif"), ".cif"),
        (("reactant.mmcif", "product.mmcif"), ".cif"),
        (("reactant.pdb", "product.cif"), ".cif"),
    ],
)
def test_public_merged_coordinate_suffix_tracks_original_input_format(
    names: tuple[str, ...], expected: str
) -> None:
    from pdb2reaction.workflows.all import _public_merged_coordinate_suffix

    assert _public_merged_coordinate_suffix(tuple(Path(name) for name in names)) == expected


def test_tsopt_continuation_separates_numerical_status_and_saddle_order() -> None:
    from pdb2reaction.workflows.all import _tsopt_continuation_decision

    first_order = _tsopt_continuation_decision({
        "optimization_status": "converged",
        "hessian_status": "completed",
        "saddle_validation": "first_order",
        "n_imaginary_modes": 1,
        "reaction_mode_index": 0,
        "reaction_mode_frequency_cm": -450.0,
    })
    assert first_order["continue_irc"] is True
    assert first_order["reason"] == "first_order_saddle"

    higher_order = _tsopt_continuation_decision({
        "optimization_status": "converged",
        "hessian_status": "completed",
        "saddle_validation": "higher_order",
        "n_imaginary_modes": 2,
        "reaction_mode_index": 1,
        "reaction_mode_frequency_cm": -120.0,
    })
    assert higher_order["continue_irc"] is True
    assert higher_order["reason"] == "higher_order_saddle"

    invalid_positive_root = _tsopt_continuation_decision({
        "optimization_status": "converged",
        "hessian_status": "completed",
        "n_imaginary_modes": 2,
        "imaginary_frequencies_cm": [-450.0, -120.0],
        "reaction_mode_index": 8,
        "reaction_mode_frequency_cm": 25.0,
    })
    assert invalid_positive_root["continue_irc"] is True
    assert invalid_positive_root["reaction_mode_index"] == 0
    assert invalid_positive_root["reaction_mode_frequency_cm"] == -450.0
    assert invalid_positive_root["reaction_mode_fallback"] is True

    for payload in (
        {"optimization_status": "not_converged", "hessian_status": "completed", "n_imaginary_modes": 1},
        {"optimization_status": "converged", "hessian_status": "completed", "n_imaginary_modes": 0},
        {"optimization_status": "converged", "hessian_status": "failed", "n_imaginary_modes": None},
        {"status": "unknown"},
    ):
        assert _tsopt_continuation_decision(payload)["continue_irc"] is False

def test_all_dft_child_omits_wrapper_defaults_for_yaml_resolution(
    tmp_path: Path, monkeypatch,
) -> None:
    from pdb2reaction.workflows import all as all_workflow
    import subprocess

    captured = {}

    def fake_run(cmd, **_kwargs):
        captured["cmd"] = list(cmd)
        return SimpleNamespace(stdout="", stderr="", returncode=0)

    monkeypatch.setattr(subprocess, "run", fake_run)
    config = tmp_path / "config.yaml"
    config.write_text("dft:\n  engine: cpu\n  func_basis: pbe/def2-svp\n")

    all_workflow._run_dft_for_state(
        tmp_path / "input.pdb",
        0,
        1,
        tmp_path / "dft",
        config,
        func_basis=None,
        engine=None,
    )

    assert "--func-basis" not in captured["cmd"]
    assert "--engine" not in captured["cmd"]
    assert captured["cmd"][-2:] == ["--config", str(config)]


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
        ts_dir.mkdir(parents=True, exist_ok=True)
        (ts_dir / "final_geometry.xyz").write_text(
            "1\nTS\nH 0.0 0.0 0.0\n", encoding="utf-8"
        )
        (ts_dir / "result.json").write_text(
            json.dumps({
                "status": "not_converged",
                "optimization_status": "not_converged",
                "hessian_status": "completed",
                "saddle_validation": "higher_order",
                "n_imaginary_modes": 2,
                "imaginary_frequencies_cm": [-450.0, -120.0],
                "reaction_mode_index": 0,
                "reaction_mode_frequency_cm": -450.0,
            }),
            encoding="utf-8",
        )

    monkeypatch.setattr(all_workflow, "_run_cli_main", fake_run_cli_main)
    monkeypatch.setattr(all_workflow, "_echo_detail", lambda *_a, **_k: None)

    _ts_path, ts_geom = all_workflow._run_tsopt_on_hei(
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

    assert ts_geom._tsopt_continuation["continue_irc"] is False
    assert ts_geom._tsopt_continuation["reason"] == "ts_optimization_not_converged"
    assert "--out-json" in captured_args


def test_all_coord_type_is_injected_into_child_config() -> None:
    from pdb2reaction.workflows.all import _write_args_yaml_with_freeze_atoms
    from pdb2reaction.workflows._run_session import RunSession

    session = RunSession()
    path = _write_args_yaml_with_freeze_atoms(
        None, [], coord_type="dlc", session=session
    )
    assert path is not None
    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    assert payload["geom"]["coord_type"] == "dlc"
    session.close()
    assert not path.exists()


def test_all_cli_freeze_atoms_reach_the_child_config() -> None:
    """``all --freeze-atoms`` must land in every child stage's geom block.

    The CLI list is merged with YAML ``geom.freeze_atoms`` (never replaced) and
    stays 1-based on the way out, exactly like the single-stage subcommands.
    """
    from pdb2reaction.workflows import all as all_mod
    from pdb2reaction.workflows.all import _write_args_yaml_with_freeze_atoms
    from pdb2reaction.workflows._run_session import RunSession

    all_mod._set_yaml_freeze_atoms({"geom": {"freeze_atoms": [3]}}, "14, 7")
    assert all_mod._freeze_atoms_for_log() == [2, 6, 13]  # 0-based internally

    session = RunSession()
    path = _write_args_yaml_with_freeze_atoms(
        None, all_mod._freeze_atoms_for_log(), session=session
    )
    assert path is not None
    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    assert payload["geom"]["freeze_atoms"] == [3, 7, 14]
    session.close()
    all_mod._set_yaml_freeze_atoms(None)


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
            "status": "partial",
            "status_reasons": ["legacy incomplete"],
            "execution_status": "completed",
            "scientific_status": "failed",
            "scientific_status_reasons": ["endpoint optimization failed"],
        }
        payload = build_pipeline_summary_payload(
            out_dir=out_dir,
            path_dir=path_dir,
            summary=summary,
            refine_path=True,
            flatten=True,
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
    assert payload["flatten"] is True
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
    assert payload["status"] == "partial"
    assert payload["status_reasons"] == ["legacy incomplete"]
    assert payload["execution_status"] == "completed"
    assert payload["scientific_status"] == "failed"
    assert payload["scientific_status_reasons"] == [
        "endpoint optimization failed"
    ]
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
        flatten=False,
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
            "R": {
                "point_group": "C2",
                "point_group_source": "auto",
                "symmetry_number": 2,
                "symmetry_number_source": "auto",
            },
            "TS": {"symmetry_number": 3, "symmetry_number_source": "config"},
            "P": {"symmetry_number": 1},
            "other": {"symmetry_number": 9, "symmetry_number_source": "test"},
        }
    )

    assert payload == {
        "R": {
            "point_group": "C2",
            "point_group_source": "auto",
            "symmetry_number": 2,
            "symmetry_number_source": "auto",
        },
        "TS": {
            "symmetry_number": 3,
            "symmetry_number_source": "config",
        },
    }


@pytest.mark.parametrize("invalid", [True, 0, -1, 2.0, "2", None])
def test_thermo_symmetry_provenance_rejects_invalid_values(invalid) -> None:
    assert build_thermo_symmetry_provenance(
        {
            "R": {
                "symmetry_number": invalid,
                "symmetry_number_source": "auto",
            }
        }
    ) == {}


def test_thermo_value_uses_finite_result_independent_of_mode_count() -> None:
    from pdb2reaction.workflows.all import _thermo_value_ha

    assert _thermo_value_ha({"num_imag_freq": 3, "g": -11.0}, "g") == -11.0
    assert _thermo_value_ha({"num_imag_freq": 0, "g": float("nan")}, "g") is None
    assert _thermo_value_ha({}, "g") is None


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
    assert result["mlip_model_label"] == "MACE-OMOL-0"
    assert any(ref["method"] == "pdb2reaction" for ref in result["references"])
    assert any(ref["method"] == "MACE" for ref in result["references"])
    assert any(
        ref["method"] == "OMol25" for ref in result["references"]
    )
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


def test_summary_frequency_reader_uses_the_resolved_zero_cutoff(tmp_path) -> None:
    from pdb2reaction.workflows.all import _read_imaginary_frequency

    freq_dir = tmp_path / "freq"
    freq_dir.mkdir()
    (freq_dir / "frequencies_cm-1.txt").write_text(
        "1 -4.99\n2 -5.00\n3 -5.01\n4 12.0\n",
        encoding="utf-8",
    )

    result = _read_imaginary_frequency(freq_dir, 5.0)

    assert result is not None
    assert result["n_imag"] == 1
    assert result["nu_imag_max_cm"] == -5.01
    assert result["frequency_zero_cutoff_cm"] == 5.0
