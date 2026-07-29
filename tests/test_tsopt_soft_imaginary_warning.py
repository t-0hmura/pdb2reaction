"""A soft leading imaginary mode must warn without changing the saddle status.

Saddle certification counts imaginary modes without a magnitude cutoff. The
soft-mode warning is diagnostic only and does not alter terminal status.
"""

import pytest

from pdb2reaction.core.defaults import TS_IMAG_SOFT_WARN_CM
from pdb2reaction.workflows.tsopt import _warn_if_leading_imaginary_mode_is_soft


def test_soft_leading_mode_warns(capsys) -> None:
    _warn_if_leading_imaginary_mode_is_soft([-15.72])
    out = capsys.readouterr().out
    assert "WARNING" in out and "-15.72 cm^-1" in out
    assert str(int(TS_IMAG_SOFT_WARN_CM)) in out


def test_large_magnitude_mode_is_silent(capsys) -> None:
    _warn_if_leading_imaginary_mode_is_soft([-447.78])
    assert capsys.readouterr().out == ""


def test_leading_mode_is_the_most_negative_one(capsys) -> None:
    # The warning evaluates the most negative mode, not a soft companion.
    _warn_if_leading_imaginary_mode_is_soft([-447.78, -5.72])
    assert capsys.readouterr().out == ""


@pytest.mark.parametrize("ims", [None, []])
def test_no_imaginary_modes_is_silent(capsys, ims) -> None:
    _warn_if_leading_imaginary_mode_is_soft(ims)
    assert capsys.readouterr().out == ""


def test_threshold_boundary_is_not_warned(capsys) -> None:
    _warn_if_leading_imaginary_mode_is_soft([-TS_IMAG_SOFT_WARN_CM])
    assert capsys.readouterr().out == ""
