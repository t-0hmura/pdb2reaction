"""Fail-closed validation for scan axes and scientific eligibility."""

import click
import inspect
import numpy as np
import pytest

from pdb2reaction.core.utils import (
    claim_unique_scan_stem,
    parse_scan_list_quads,
    parse_scan_list_triples,
    parse_scan_spec_stages,
)
from pdb2reaction.workflows.scan2d import _rbf_support
from pdb2reaction.workflows.scan3d import _rbf_support_3d
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


@pytest.mark.parametrize(
    ("x", "y", "z", "expected"),
    [
        ([1.0], [2.0], [3.0], (1, 0)),
        ([0.0, 1.0, 0.0, 1.0], [0.0, 0.0, 1.0, 1.0], [2.0] * 4, (4, 2)),
        (
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
            (4, 3),
        ),
    ],
)
def test_scan3d_rbf_support_requires_non_coplanar_points(
    x, y, z, expected,
) -> None:
    assert _rbf_support_3d(
        np.asarray(x),
        np.asarray(y),
        np.asarray(z),
    ) == expected


def test_scan_artifact_stem_disambiguates_rounded_distance_collision() -> None:
    used: set[str] = set()

    first = claim_unique_scan_stem("point_i100_j200", (0, 0), used)
    second = claim_unique_scan_stem("point_i100_j200", (1, 0), used)

    assert first == "point_i100_j200"
    assert second == "point_i100_j200_grid_001_000"


def test_scan_html_plots_are_responsive_and_scan2d_keeps_native_projection() -> None:
    from pdb2reaction.workflows import scan2d, scan3d

    scan2d_source = inspect.getsource(scan2d)
    scan3d_source = inspect.getsource(scan3d)

    assert "plane_proj = go.Surface" in scan2d_source
    assert "plane_z = z_bottom + 0.005" in scan2d_source
    assert "z=np.full_like(ZI, plane_z)" in scan2d_source
    assert "bottom_contours = go.Scatter3d" in scan2d_source
    assert "contour_z = z_bottom + 0.01" in scan2d_source
    assert "go.Figure(data=[surface3d, plane_proj, bottom_contours])" in scan2d_source
    assert 'name="2D Contour Projection (Bottom)"' in scan2d_source
    assert "Computed grid points" not in scan2d_source
    assert '"project": {"z": True}' not in scan2d_source
    for source in (scan2d_source, scan3d_source):
        assert 'config={"responsive": True, "displaylogo": False}' in source
        assert 'default_width="100%"' in source
        assert 'default_height="100%"' in source


def test_scan2d_bottom_contours_are_explicit_line_segments() -> None:
    from pdb2reaction.workflows.scan2d import _contour_line_segments

    x_grid, y_grid = np.meshgrid(
        np.asarray([1.0, 1.5, 2.0]),
        np.asarray([1.0, 1.5, 2.0]),
    )
    z_grid = x_grid + y_grid
    line_x, line_y = _contour_line_segments(
        x_grid, y_grid, z_grid, np.asarray([2.5, 3.0])
    )

    assert line_x
    assert len(line_x) == len(line_y)
    assert any(value is None for value in line_x)
    assert all(value is None or 1.0 <= value <= 2.0 for value in line_x)
    assert all(value is None or 1.0 <= value <= 2.0 for value in line_y)
