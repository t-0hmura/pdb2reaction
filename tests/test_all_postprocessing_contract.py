import click
import pytest

from pdb2reaction.io.summary import write_summary_log
from pdb2reaction.workflows.all import (
    _derive_pipeline_status,
    _enrich_summary,
    _pipeline_aggregate_truth,
    _ts_imag_record,
    _validate_postprocessing_dependencies,
)


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


def test_explicit_no_change_segment_does_not_require_postprocessing() -> None:
    summary = {
        "segments": [
            {
                "index": 1,
                "kind": "seg",
                "converged": True,
                "bond_changes": "(no covalent changes detected)",
            }
        ],
        "energy_diagrams": [{"name": "MEP"}],
    }
    config = {"tsopt": True, "thermo": False, "dft": False}

    status, reasons = _derive_pipeline_status(
        summary, post_segments=[], config=config
    )
    truth = _pipeline_aggregate_truth(
        summary,
        post_segments=[],
        config=config,
        legacy_status=status,
        legacy_reasons=reasons,
    )

    assert status == "success"
    assert truth.scientific_status == "success"
    assert truth.expected_item_ids == ()


def test_tsopt_validates_imaginary_mode_without_thermochemistry() -> None:
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

    assert status == "partial"
    assert any("n_imag=2, expected 1" in reason for reason in reasons)


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
