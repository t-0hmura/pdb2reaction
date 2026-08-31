import click
import pytest
from click.testing import CliRunner

from pdb2reaction.io.summary import write_summary_log
from pdb2reaction.workflows.all import (
    _derive_pipeline_status,
    _enrich_summary,
    _pipeline_aggregate_truth,
    _reject_redundant_dft_postprocessing,
    _ts_imag_record,
    _validate_postprocessing_dependencies,
)


def test_all_rejects_post_dft_when_primary_backend_is_dft() -> None:
    with pytest.raises(click.UsageError, match="separate process/job"):
        _reject_redundant_dft_postprocessing(
            effective_backend="dft", do_dft=True
        )
    _reject_redundant_dft_postprocessing(
        effective_backend="uma", do_dft=True
    )


@pytest.mark.parametrize("via_yaml", [False, True])
def test_all_cli_rejects_redundant_dft_before_pipeline(tmp_path, via_yaml) -> None:
    from pdb2reaction.workflows.all import cli as all_cli

    xyz = tmp_path / "h2.xyz"
    xyz.write_text("2\nH2\nH 0 0 0\nH 0 0 0.74\n", encoding="utf-8")
    out_dir = tmp_path / "result"
    args = [
        "-i", str(xyz), "-i", str(xyz), "-q", "0", "--dft", "true", "--dry-run",
        "--out-dir", str(out_dir),
    ]
    if via_yaml:
        config = tmp_path / "dft.yaml"
        config.write_text("calc:\n  backend: dft\n", encoding="utf-8")
        args.extend(["--config", str(config)])
    else:
        args.extend(["-b", "dft"])

    result = CliRunner().invoke(all_cli, args)

    assert result.exit_code == 2
    assert "separate process/job" in result.output
    assert not out_dir.exists()


@pytest.mark.parametrize(
    ("do_thermo", "do_dft"),
    [(True, False), (False, True), (True, True)],
)
def test_all_rejects_ts_labeled_postprocessing_without_tsopt(
    do_thermo, do_dft
):
    with pytest.raises(click.UsageError, match="require `--tsopt`"):
        _validate_postprocessing_dependencies(
            do_tsopt=False,
            do_thermo=do_thermo,
            do_dft=do_dft,
        )


@pytest.mark.parametrize(
    ("do_tsopt", "do_thermo", "do_dft"),
    [(False, False, False), (True, False, False), (True, True, True)],
)
def test_all_accepts_scientifically_defined_postprocessing_combinations(
    do_tsopt, do_thermo, do_dft
):
    _validate_postprocessing_dependencies(
        do_tsopt=do_tsopt,
        do_thermo=do_thermo,
        do_dft=do_dft,
    )


def test_all_show_config_separates_canonical_primary_and_post_dft(tmp_path) -> None:
    from pdb2reaction.workflows.all import cli as all_cli

    xyz = tmp_path / "h2.xyz"
    xyz.write_text("2\nH2\nH 0 0 0\nH 0 0 0.74\n", encoding="utf-8")
    post_out = tmp_path / "post_dft"
    result = CliRunner().invoke(
        all_cli,
        [
            "-i", str(xyz), "-i", str(xyz), "-q", "-1", "-m", "2",
            "--tsopt", "true", "--dft", "true",
            "--dft-func-basis", "hf/sto-3g",
            "--dft-engine", "cpu", "--dft-out-dir", str(post_out),
            "--show-config", "--dry-run", "--out-dir", str(tmp_path / "result"),
        ],
    )

    assert result.exit_code == 0, result.output
    assert "primary_calculator:" in result.output
    assert "primary_method_label: MLIP" in result.output
    assert "post_dft:" in result.output
    assert "functional: hf" in result.output
    assert "basis: sto-3g" in result.output
    assert "charge: -1" in result.output
    assert "multiplicity: 2" in result.output
    assert "engine: cpu" in result.output
    assert f"out_dir_override: {post_out}" in result.output


def test_all_show_config_primary_dft_uses_explicit_state(tmp_path) -> None:
    from pdb2reaction.workflows.all import cli as all_cli

    xyz = tmp_path / "h2.xyz"
    xyz.write_text("2\nH2\nH 0 0 0\nH 0 0 0.74\n", encoding="utf-8")
    result = CliRunner().invoke(
        all_cli,
        [
            "-i", str(xyz), "-i", str(xyz), "-q", "-1", "-m", "2",
            "-b", "dft", "--func-basis", "hf/sto-3g", "--engine", "cpu",
            "--show-config", "--dry-run", "--out-dir", str(tmp_path / "result"),
        ],
    )

    assert result.exit_code == 0, result.output
    assert "primary_method_label: DFT" in result.output
    assert "charge: -1" in result.output
    assert "multiplicity: 2" in result.output
    assert "post_dft: null" in result.output


@pytest.mark.parametrize(
    "primary_backend, expected",
    [("dft", "DFT thermochemistry"), ("uma", "MLIP thermochemistry")],
)
def test_missing_thermochemistry_reason_uses_primary_calculator_label(
    primary_backend, expected
) -> None:
    status, reasons = _derive_pipeline_status(
        {
            "segments": [
                {"index": 1, "kind": "seg", "bond_changes": "C1-O2"}
            ],
            "energy_diagrams": [{"name": "MEP"}],
        },
        post_segments=[{"index": 1}],
        config={"tsopt": False, "thermo": True, "dft": False},
        primary_backend=primary_backend,
    )

    assert status == "partial"
    assert any(expected in reason for reason in reasons)


@pytest.mark.parametrize("bond_changes", ["", "(no covalent changes detected)", "forming 1-2"])
def test_bond_diagnostics_do_not_suppress_requested_postprocessing(bond_changes):
    from pdb2reaction.workflows.all import _is_reactive_segment

    segment = {"index": 1, "kind": "seg", "converged": True,
               "bond_changes": bond_changes}
    assert _is_reactive_segment(segment)
    assert not _is_reactive_segment({**segment, "kind": "bridge"})
    summary = {"segments": [segment], "energy_diagrams": [{"name": "MEP"}]}
    config = {"tsopt": True}
    status, reasons = _derive_pipeline_status(summary, post_segments=[], config=config)
    truth = _pipeline_aggregate_truth(summary, post_segments=[], config=config,
                                      legacy_status=status, legacy_reasons=reasons)
    assert status == "partial"
    assert "segment 1: requested post-processing record is missing" in reasons
    assert truth.scientific_status != "success"
    assert truth.expected_item_ids == ("segment_1",)


def test_tsopt_frequency_counts_do_not_replace_optimizer_completion() -> None:
    summary = {
        "segments": [{"index": 1, "kind": "seg", "converged": True}],
        "energy_diagrams": [{"name": "MEP"}],
    }
    post = [{
        "index": 1,
        "mlip": {},
        "irc_traj": "finished_irc_trj.xyz",
        "ts_imag": {"n_imag": 2},
    }]

    status, reasons = _derive_pipeline_status(
        summary,
        post_segments=post,
        config={"tsopt": True, "thermo": False, "dft": False},
    )

    assert status == "success"
    assert reasons == []


def test_tsopt_imaginary_mode_record_carries_certification_details() -> None:
    assert _ts_imag_record(1, [-512.31], 5.0) == {
        "n_imag": 1,
        "imag_freqs_cm": [-512.31],
        "nu_imag_max_cm": -512.31,
        "min_abs_imag_cm": 512.31,
        "frequency_zero_cutoff_cm": 5.0,
    }


def test_path_postprocessing_uses_terminal_tsopt_imaginary_mode_record(
    tmp_path,
) -> None:
    """A validated path TS must not require the legacy ``ts_imag`` duplicate."""
    summary = {
        "segments": [{
            "index": 1,
            "tag": "seg_001",
            "kind": "seg",
            "converged": True,
            "barrier_kcal": 40.49,
            "delta_kcal": 10.11,
        }],
        "energy_diagrams": [{
            "name": "energy_diagram_MLIP_all",
            "labels": ["R", "TS1", "P"],
            "energies_kcal": [0.0, 40.10, 10.93],
        }],
    }
    post = [{
        "index": 1,
        "tag": "seg_001",
        "kind": "seg",
        "mlip": {"barrier_kcal": 40.10, "delta_kcal": 10.93},
        "irc_traj": "finished_irc_trj.xyz",
        "tsopt": {
            "optimization_status": "converged",
            "continue_irc": True,
            "saddle_validation": "first_order",
            "n_imaginary_modes": 1,
            "reaction_mode_frequency_cm": -512.31,
        },
        "irc": {"usable": True, "reason": "ok"},
        "endpoint_assignment": {"connectivity_validated": True},
        "endpoint_opt": {
            "reactant_converged": True,
            "product_converged": True,
            "connectivity_validated": True,
        },
    }]

    _enrich_summary(
        summary,
        version="",
        pipeline_mode="path-opt",
        out_dir=tmp_path,
        mlip_backend="orb",
        mlip_model="orb_v3_conservative_omol",
        charge=0,
        spin=1,
        post_segments=post,
        config={"tsopt": True, "thermo": False, "dft": False},
    )

    assert summary["status"] == "success"
    assert summary["scientific_status"] == "success"
    assert "status_reasons" not in summary
    assert "scientific_status_reasons" not in summary

    post[0]["ts_imag"] = _ts_imag_record(1, [-512.31], 5.0)
    summary["root_out_dir"] = str(tmp_path)
    summary["post_segments"] = post
    destination = tmp_path / "summary.log"
    write_summary_log(destination, summary)
    text = destination.read_text(encoding="utf-8")
    assert "Scientific status   : success" in text
    assert "TS imaginary-mode validation is missing" not in text
    assert "n_imag       : 1" in text
    assert "ν_imag (max) : -512.3 cm^-1" in text


def test_postprocessing_reports_each_missing_reactive_segment() -> None:
    summary = {
        "segments": [
            {"index": 1, "kind": "seg", "converged": True},
            {"index": 2, "kind": "seg", "converged": True},
        ],
        "energy_diagrams": [{"name": "MEP"}],
    }

    status, reasons = _derive_pipeline_status(
        summary,
        post_segments=[{
            "index": 1,
            "mlip": {},
            "irc_traj": "finished_irc_trj.xyz",
            "ts_imag": {"n_imag": 1},
        }],
        config={"tsopt": True},
    )

    assert status == "partial"
    assert "segment 2: requested post-processing record is missing" in reasons
