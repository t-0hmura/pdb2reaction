import numpy as np
import pytest
from types import SimpleNamespace

from pdb2reaction.workflows.path_opt import _select_hei_index as select_path_opt_hei
from pdb2reaction.workflows.path_search import (
    CombinedPath,
    _frame_ranges_by_segment,
    _select_hei_index as select_path_search_hei,
    _stitch_paths,
)


@pytest.mark.parametrize("selector", [select_path_opt_hei, select_path_search_hei])
@pytest.mark.parametrize(
    ("energies", "expected"),
    [
        ([0.0, 1.0, 2.0], 2),
        ([0.0, 2.0, 1.0, 3.0], 3),
        ([3.0, 2.0, 1.0], 0),
        ([0.0, 3.0, 1.0], 1),
    ],
)
def test_hei_is_the_global_energy_maximum(selector, energies, expected):
    assert selector(energies) == expected


@pytest.mark.parametrize("selector", [select_path_opt_hei, select_path_search_hei])
@pytest.mark.parametrize("energies", [[], [0.0, np.nan], [0.0, np.inf]])
def test_hei_selection_rejects_missing_or_nonfinite_energies(selector, energies):
    with pytest.raises(ValueError):
        selector(energies)


def test_frame_ranges_preserve_bridge_and_disjoint_segment_provenance():
    images = [
        *[SimpleNamespace(mep_seg_index=1, mep_seg_kind="seg") for _ in range(3)],
        *[SimpleNamespace(mep_seg_index=2, mep_seg_kind="bridge") for _ in range(2)],
        *[SimpleNamespace(mep_seg_index=3, mep_seg_kind="seg") for _ in range(3)],
        SimpleNamespace(mep_seg_index=1, mep_seg_kind="seg"),
    ]

    ranges = _frame_ranges_by_segment(images)

    assert ranges[1] == {
        "kind": "seg",
        "frame_ranges": [[0, 3], [8, 9]],
    }
    assert ranges[2] == {
        "kind": "bridge",
        "frame_ranges": [[3, 5]],
        "frame_start": 3,
        "frame_stop": 5,
    }
    assert ranges[3]["frame_ranges"] == [[5, 8]]


def test_cross_pair_recursive_segment_receives_interface_pair_index(
    monkeypatch, tmp_path,
):
    from pdb2reaction.workflows import path_search

    left = SimpleNamespace(coords3d=np.zeros((1, 3)))
    right = SimpleNamespace(coords3d=np.ones((1, 3)))
    captured = []

    monkeypatch.setattr(
        path_search,
        "has_bond_change",
        lambda *_args, **_kwargs: (True, "changed"),
    )

    def builder(tail, head, tag, pair_index):
        captured.append((tag, pair_index))
        return CombinedPath(
            images=[tail, head],
            energies=[0.0, 1.0],
            segments=[],
        )

    _stitch_paths(
        [([left], [0.0]), ([right], [1.0])],
        stitch_rmsd_thresh=1.0e-6,
        bridge_rmsd_thresh=1.0e-6,
        shared_calc=None,
        gs_cfg={},
        stopt_cfg={},
        out_dir=tmp_path,
        tag="pair_01",
        ref_pdb_path=None,
        bond_cfg={},
        segment_builder=builder,
        bridge_pair_index=1,
    )

    assert captured == [("pair_01_mid", 1)]
