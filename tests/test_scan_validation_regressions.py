"""Fail-closed validation for scan axes and scientific eligibility."""

import click
import numpy as np
import pytest

from pdb2reaction.core.utils import (
    parse_scan_list_quads,
    parse_scan_list_triples,
    parse_scan_spec_stages,
)
from pdb2reaction.workflows.scan2d import _rbf_support
from pdb2reaction.workflows.scan_common import collect_staged_scan_values


@pytest.mark.parametrize(
    "raw",
    [
        "[]",
        "[(1, 1, 1.5)]",
        "[(1, 2, 1.5), (2, 1, 1.6)]",
    ],
)
def test_scan_stage_rejects_empty_self_and_duplicate_pairs(raw) -> None:
    with pytest.raises(click.BadParameter):
        parse_scan_list_triples(
            raw,
            one_based=False,
            atom_meta=None,
            option_name="--scan-lists",
        )


@pytest.mark.parametrize(
    "raw",
    [
        "[(1, 1, 1.0, 2.0), (2, 3, 1.0, 2.0)]",
        "[(1, 2, 1.0, 2.0), (2, 1, 1.0, 2.0)]",
    ],
)
def test_grid_scan_rejects_self_and_duplicate_axes(raw) -> None:
    with pytest.raises(click.BadParameter):
        parse_scan_list_quads(
            raw,
            expected_len=2,
            one_based=False,
            atom_meta=None,
            option_name="--scan-lists",
        )


def test_scan_spec_expands_bidirectional_stage_with_reset_markers(tmp_path) -> None:
    spec = tmp_path / "scan.yaml"
    spec.write_text(
        "one_based: true\n"
        "stages:\n"
        "  - [[1, 2, 1.2, 1.8]]\n",
        encoding="utf-8",
    )

    stages, one_based, snapshots, resets = parse_scan_spec_stages(
        spec,
        one_based_default=False,
        atom_meta=None,
        return_bidirectional_markers=True,
    )

    assert one_based is True
    assert stages == [[(0, 1, 1.2)], [(0, 1, 1.8)]]
    assert snapshots == frozenset({0})
    assert resets == frozenset({1})


def test_mixed_grouped_and_repeated_scan_options_are_rejected() -> None:
    with pytest.raises(click.BadParameter, match="Do not mix repeated"):
        collect_staged_scan_values(
            ("[(1,2,1.2)]", "[(1,2,1.4)]"),
            ("[(1,2,1.6)]",),
        )


def test_grouped_scan_values_keep_their_order() -> None:
    assert collect_staged_scan_values(
        ("[(1,2,1.2)]",),
        ("[(1,2,1.4)]", "[(1,2,1.6)]"),
    ) == ("[(1,2,1.2)]", "[(1,2,1.4)]", "[(1,2,1.6)]")


def test_all_rejects_bidirectional_scan_tuple_cleanly() -> None:
    from pdb2reaction.workflows.all import _parse_scan_lists_literals

    with pytest.raises(click.BadParameter, match="only .* scan triples"):
        _parse_scan_lists_literals(
            ("[(1,2,1.2,1.6)]",),
            one_based=True,
        )


def test_bidirectional_trajectory_contains_one_center_frame() -> None:
    from pdb2reaction.workflows.scan import _assemble_bidirectional_trajectory

    assert _assemble_bidirectional_trajectory(
        ["center-to-start-1", "center-to-start-2"],
        "center",
        ["center-to-end-1"],
    ) == [
        "center-to-start-2",
        "center-to-start-1",
        "center",
        "center-to-end-1",
    ]


def test_bidirectional_zero_step_second_leg_keeps_first_leg() -> None:
    from pdb2reaction.workflows.scan import _assemble_bidirectional_trajectory

    assert _assemble_bidirectional_trajectory(
        ["center-to-start"],
        "center",
        [],
    ) == ["center-to-start", "center"]


@pytest.mark.parametrize(
    ("x", "y", "expected"),
    [
        ([1.0], [2.0], (1, 0)),
        ([1.0, 2.0, 3.0], [2.0, 2.0, 2.0], (3, 1)),
        ([1.0, 2.0, 1.0], [2.0, 2.0, 3.0], (3, 2)),
    ],
)
def test_scan2d_rbf_support_requires_non_collinear_points(
    x, y, expected,
) -> None:
    assert _rbf_support(np.asarray(x), np.asarray(y)) == expected
