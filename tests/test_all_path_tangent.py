"""Tests for carrying the MEP reaction direction into TS optimization."""

from __future__ import annotations

import numpy as np
from ase import Atoms
from ase.io import write

from pdb2reaction.core.utils import write_xyz_trj_with_energy
from pdb2reaction.workflows.all import _ensure_hei_path_tangent


def test_hei_path_tangent_uses_neighbors_of_matching_image(tmp_path) -> None:
    images = [
        Atoms("H2", positions=positions)
        for positions in (
            [[0.0, 0.0, 0.0], [0.7, -0.2, 0.0]],
            [[0.0, 0.0, 0.0], [0.7, 0.0, 0.0]],
            [[0.0, 0.0, 0.0], [1.0, 0.3, 0.0]],
        )
    ]
    mep = tmp_path / "mep_trj.xyz"
    hei = tmp_path / "hei.xyz"
    mode = tmp_path / "hei_mode.txt"
    write(mep, images)
    write(hei, images[1])
    mode.write_text("1\n0\n0\n", encoding="utf-8")

    result = _ensure_hei_path_tangent(mep, hei, mode)

    assert result == mode
    tangent = np.loadtxt(mode)
    incoming = (images[1].positions - images[0].positions).reshape(-1)
    outgoing = (images[2].positions - images[1].positions).reshape(-1)
    incoming /= np.linalg.norm(incoming)
    outgoing /= np.linalg.norm(outgoing)
    expected = incoming + outgoing
    expected /= np.linalg.norm(expected)
    np.testing.assert_allclose(tangent, expected)


def test_hei_path_tangent_uses_energy_upwinding_when_available(tmp_path) -> None:
    images = [
        Atoms("H", positions=positions)
        for positions in (
            [[-2.0, 0.0, 0.0]],
            [[0.0, 0.0, 0.0]],
            [[0.0, 1.0, 0.0]],
        )
    ]
    energies = [0.0, 2.0, 0.5]
    mep = tmp_path / "mep_trj.xyz"
    hei = tmp_path / "hei.xyz"
    mode = tmp_path / "hei_mode.txt"
    write_xyz_trj_with_energy(images, energies, mep)
    write(hei, images[1])

    result = _ensure_hei_path_tangent(mep, hei, mode)

    assert result == mode
    tangent = np.loadtxt(mode)
    incoming = (images[1].positions - images[0].positions).reshape(-1)
    outgoing = (images[2].positions - images[1].positions).reshape(-1)
    # At this local maximum the next neighbor is higher than the previous
    # one, so the standard improved tangent weights the outgoing secant by
    # delta_max=2.0 and the incoming secant by delta_min=1.5.
    expected = 2.0 * outgoing + 1.5 * incoming
    expected /= np.linalg.norm(expected)
    np.testing.assert_allclose(tangent, expected)
