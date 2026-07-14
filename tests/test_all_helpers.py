"""Unit tests for pdb2reaction.workflows._all_helpers."""

from __future__ import annotations

import tempfile
from pathlib import Path

import click
import pytest
import yaml

from pdb2reaction.workflows._all_helpers import (
    AllContext,
    build_energy_level_dict,
    build_pipeline_summary_payload,
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


def test_all_context_fields_match_cli_signature() -> None:
    """Regression guard: AllContext field names must mirror cli() params."""
    import dataclasses
    import inspect
    from pdb2reaction.workflows.all import cli
    callback = cli.callback
    assert callback is not None
    sig = inspect.signature(callback)
    cli_param_names = {n for n in sig.parameters if n != "ctx"}
    ctx_field_names = {f.name for f in dataclasses.fields(AllContext)}
    missing = cli_param_names - ctx_field_names
    extra = ctx_field_names - cli_param_names
    assert not missing, f"AllContext missing fields: {sorted(missing)}"
    assert not extra, f"AllContext has extra fields: {sorted(extra)}"
    # Absolute count guard catches symmetric add+remove pairs where the
    # name sets stay equal but the field count drifts.
    assert len(ctx_field_names) == 72, (
        f"AllContext field count drifted from 72 to {len(ctx_field_names)}; "
        f"update this assertion alongside the cli() signature change."
    )


def test_all_tr_projection_is_injected_into_child_config() -> None:
    from pdb2reaction.workflows.all import _write_args_yaml_with_freeze_atoms

    path = _write_args_yaml_with_freeze_atoms(
        None, [], tr_projection="legacy-active"
    )
    assert path is not None
    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    assert payload["geom"]["tr_projection"] == "legacy-active"


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
            mep_mode_kind="dmf",
            mlip_backend="mace",
            mlip_model="mace-off23-small",
            command_str="pdb2reaction all -i foo.pdb",
            q_int=-1,
            spin=1,
            freeze_atoms=[10, 20, 30],
            post_segment_logs=[{"seg": 1, "status": "ok"}],
        )
    assert payload["pipeline_mode"] == "path-search"
    assert payload["dft_func_basis"] == "wb97m-v/def2-tzvpd"
    assert payload["opt_mode"] == "grad"
    assert payload["mep_mode"] == "dmf"
    assert payload["mlip_backend"] == "mace"
    assert payload["mlip_model"] == "mace-off23-small"
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
        mep_mode_kind=None,
        mlip_backend="orb",
        mlip_model=None,
        command_str="",
        q_int=0,
        spin=1,
        freeze_atoms=[],
        post_segment_logs=[],
    )
    assert payload["pipeline_mode"] == "path-opt"
    assert payload["dft_func_basis"] is None
    assert payload["opt_mode"] is None


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
            }
        ],
        out_dir=tmp_path,
    )

    assert result["rate_limiting_step"]["barrier_kcal"] == 12.5
    assert result["rate_limiting_step"]["method"] == "MLIP_Gibbs"
    assert result["schema_version"]


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
