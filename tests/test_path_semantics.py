import numpy as np
import pytest
from types import SimpleNamespace

from pdb2reaction.workflows.path_opt import _select_hei_index as select_path_opt_hei
from pdb2reaction.workflows.path_search import (
    CombinedPath,
    SegmentReport,
    _build_mep_diagram_data,
    _echo_mep_diagram_inputs,
    _frame_ranges_by_segment,
    _select_hei_index as select_path_search_hei,
    _stitch_paths,
)
from pysisyphus.constants import AU2KCALPERMOL


@pytest.mark.parametrize("selector", [select_path_opt_hei, select_path_search_hei])
@pytest.mark.parametrize(
    ("energies", "expected"),
    [
        ([0.0, 1.0, 2.0], 2),
        ([0.0, 2.0, 1.0, 3.0], 1),
        ([3.0, 2.0, 1.0], 0),
        ([0.0, 3.0, 1.0], 1),
    ],
)
def test_hei_prefers_a_local_peak_over_higher_endpoints(selector, energies, expected):
    assert selector(energies) == expected


@pytest.mark.parametrize("drop_shared_start", [False, True])
def test_segment_peak_survives_stitching_a_shared_first_frame(drop_shared_start):
    from pdb2reaction.workflows.path_search import SegmentReport, _segment_hei_index

    images = [SimpleNamespace() for _ in range(4)]
    energies = [0.0, 2.0, 1.0, 3.0]
    report = SegmentReport("rxn", 2.0, 3.0, "Bond formed", peak_image=images[1])
    start = int(drop_shared_start)
    assert _segment_hei_index(images[start:], energies[start:], report) == 1 - start
    assert report.barrier_kcal == 2.0


def test_monotone_bridge_does_not_overwrite_a_reactive_segment_barrier(tmp_path, monkeypatch):
    from pdb2reaction.workflows import path_search

    left = SimpleNamespace(coords3d=np.zeros((1, 3)))
    right = SimpleNamespace(coords3d=np.ones((1, 3)))
    bridge_images = [left, SimpleNamespace(coords3d=np.full((1, 3), 0.5)), right]
    bridge = SimpleNamespace(images=bridge_images, energies=[0.0, 1.0, 2.0], hei_idx=2, is_converged=True)
    monkeypatch.setattr(path_search, "_bridge_segments", lambda *_a, **_k: bridge)
    report = path_search.SegmentReport("rxn", 5.0, -1.0, "Bond formed")
    reports = [report]
    path_search._stitch_paths(
        [([left], [0.0]), ([right], [2.0])],
        stitch_rmsd_thresh=1e-6, bridge_rmsd_thresh=1e-6,
        shared_calc=None, gs_cfg={}, stopt_cfg={}, out_dir=tmp_path,
        tag="pair", ref_pdb_path=None, segments_out=reports,
    )
    assert report.barrier_kcal == 5.0
    assert len(reports) == 2
    assert reports[1].kind == "bridge" and reports[1].barrier_kcal is None
    assert reports[1].peak_image is right


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


def test_mep_diagram_anchors_au_energies_to_full_path_endpoints():
    first_energy = -100.0
    final_delta_kcal = 4.5
    combined = CombinedPath(
        images=[SimpleNamespace(), SimpleNamespace()],
        energies=[first_energy, first_energy + final_delta_kcal / AU2KCALPERMOL],
        segments=[
            SegmentReport("lead", 0.0, 2.0, "", kind="kink"),
            SegmentReport("rxn1", 10.0, 5.0, "Bond formed", kind="seg"),
            SegmentReport("middle", 0.0, 3.0, "", kind="bridge"),
            SegmentReport("rxn2", 12.0, -4.0, "Bond broken", kind="seg"),
            SegmentReport("tail", 0.0, -1.5, "", kind="kink"),
        ],
    )

    payload, _ = _build_mep_diagram_data(combined)

    assert payload["energies_au"][0] == pytest.approx(combined.energies[0])
    assert payload["energies_au"][-1] == pytest.approx(combined.energies[-1])
    for energy_au, energy_kcal in zip(
        payload["energies_au"], payload["energies_kcal"]
    ):
        assert (energy_au - first_energy) * AU2KCALPERMOL == pytest.approx(
            energy_kcal
        )


def test_mep_diagram_uses_peak_energy_independently_of_the_local_barrier():
    first = -100.0
    combined = CombinedPath(
        images=[SimpleNamespace(), SimpleNamespace()],
        energies=[first, first + 5.0 / AU2KCALPERMOL],
        segments=[
            SegmentReport(
                "rxn1", 12.0, 3.0, "Bond formed", kind="seg", seg_index=1,
                peak_energy_hartree=first + 8.0 / AU2KCALPERMOL,
            ),
            SegmentReport(
                "middle", 6.0, -2.0, "", kind="bridge", seg_index=2,
                peak_energy_hartree=first + 7.0 / AU2KCALPERMOL,
            ),
            SegmentReport(
                "rxn2", 10.0, 4.0, "Bond broken", kind="seg", seg_index=3,
                peak_energy_hartree=first + 9.0 / AU2KCALPERMOL,
            ),
        ],
    )
    payload, _ = _build_mep_diagram_data(combined)
    by_label = dict(zip(payload["labels"], payload["energies_kcal"]))
    assert by_label["TS1"] == pytest.approx(8.0)
    assert by_label["IM1_TS"] == pytest.approx(7.0)
    assert by_label["TS3"] == pytest.approx(9.0)


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


@pytest.mark.parametrize("level", [0, 1, 2, 3])
def test_mep_diagram_inputs_print_only_at_level_three(monkeypatch, capsys, level):
    from pdb2reaction.core import output

    monkeypatch.setattr(output, "_VERBOSE_LEVEL", level)
    _echo_mep_diagram_inputs(["R", "TS1", "P"], [0.0, 12.5, -3.25])

    out = capsys.readouterr().out
    if level >= 3:
        assert out == (
            '[diagram] build_energy_diagram.labels = ["R", "TS1", "P"]\n'
            "[diagram] build_energy_diagram.energies_kcal = "
            "[0.000000, 12.500000, -3.250000]\n"
        )
    else:
        assert out == ""
