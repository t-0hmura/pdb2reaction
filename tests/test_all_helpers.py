"""Unit tests for pdb2reaction.workflows._all_helpers."""

from __future__ import annotations

import tempfile
from pathlib import Path

from pdb2reaction.workflows._all_helpers import (
    AllContext,
    build_energy_level_dict,
    build_pipeline_summary_payload,
)


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
    assert len(ctx_field_names) == 65, (
        f"AllContext field count drifted from 65 to {len(ctx_field_names)}; "
        f"update this assertion alongside the cli() signature change."
    )


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
            uma_model="uma-s-1p1",
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
    assert payload["uma_model"] == "uma-s-1p1"
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
        uma_model=None,
        command_str="",
        q_int=0,
        spin=1,
        freeze_atoms=[],
        post_segment_logs=[],
    )
    assert payload["pipeline_mode"] == "path-opt"
    assert payload["dft_func_basis"] is None
    assert payload["opt_mode"] is None
