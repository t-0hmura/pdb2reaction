import click
import pytest

from pdb2reaction.workflows.all import (
    _derive_pipeline_status,
    _pipeline_aggregate_truth,
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
