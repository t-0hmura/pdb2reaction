"""A soft leading imaginary mode must warn without changing the saddle status.

Saddle certification counts imaginary modes and does not weigh them, so a
few-cm^-1 soft mode certifies exactly like a real reaction coordinate. Two
runs of the same input from bit-identical starting geometries were observed to
diverge and land on a real saddle (~-450 cm^-1) and on a soft-mode structure
(~-15 cm^-1) respectively, with the soft one still reported as n_imag=1. The
warning is diagnostic only: the terminal status must stay exactly as it was.
"""

import pytest

from pdb2reaction.core.defaults import TS_IMAG_SOFT_WARN_CM
from pdb2reaction.workflows.tsopt import _warn_if_leading_imaginary_mode_is_soft


def test_soft_leading_mode_warns(capsys) -> None:
    _warn_if_leading_imaginary_mode_is_soft([-15.72])
    out = capsys.readouterr().out
    assert "WARNING" in out and "-15.72 cm^-1" in out
    assert str(int(TS_IMAG_SOFT_WARN_CM)) in out


def test_real_reaction_coordinate_is_silent(capsys) -> None:
    _warn_if_leading_imaginary_mode_is_soft([-447.78])
    assert capsys.readouterr().out == ""


def test_leading_mode_is_the_most_negative_one(capsys) -> None:
    # n_imag=2 with a real coordinate plus a soft companion: the certifying
    # mode is the stiff one, so this must stay silent.
    _warn_if_leading_imaginary_mode_is_soft([-447.78, -5.72])
    assert capsys.readouterr().out == ""


@pytest.mark.parametrize("ims", [None, []])
def test_no_imaginary_modes_is_silent(capsys, ims) -> None:
    _warn_if_leading_imaginary_mode_is_soft(ims)
    assert capsys.readouterr().out == ""


def test_threshold_boundary_is_not_warned(capsys) -> None:
    _warn_if_leading_imaginary_mode_is_soft([-TS_IMAG_SOFT_WARN_CM])
    assert capsys.readouterr().out == ""
