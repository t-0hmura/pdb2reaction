"""Regression tests for opt-in IRC energy-stop bypass controls."""

from __future__ import annotations

from pdb2reaction.core.defaults import IRC_KW
from pysisyphus.irc.IRC import IRC


def _irc_stop_probe(*, never_stop: bool, increased: bool, converged: bool) -> IRC:
    irc = object.__new__(IRC)
    irc.never_stop = never_stop
    irc.energy_increased = increased
    irc.energy_converged = converged
    return irc


def test_never_stop_is_opt_in() -> None:
    assert IRC_KW["never_stop"] is False


def test_default_irc_stops_on_energy_increase_or_plateau() -> None:
    assert _irc_stop_probe(
        never_stop=False, increased=True, converged=False
    )._energy_stop_message() == "Energy increased!"
    assert _irc_stop_probe(
        never_stop=False, increased=False, converged=True
    )._energy_stop_message() == "Energy converged!"


def test_never_stop_ignores_energy_only_stops() -> None:
    assert _irc_stop_probe(
        never_stop=True, increased=True, converged=False
    )._energy_stop_message() == ""
    assert _irc_stop_probe(
        never_stop=True, increased=False, converged=True
    )._energy_stop_message() == ""
