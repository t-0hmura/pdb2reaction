"""Flags and settings that do nothing must say so.

Requested settings that cannot be honored must emit a diagnostic instead of
appearing to have been applied. These output-only guards do not change
behavior, thresholds, or defaults.
"""

from pathlib import Path

import click

import pdb2reaction.backends as backends_mod
import pdb2reaction.backends.orb as orb_mod
import pdb2reaction.cli.common_options as common_options_mod
import pdb2reaction.core.utils as utils_mod


def _src(mod) -> str:
    return Path(mod.__file__).read_text(encoding="utf-8")


def test_auto_backend_says_which_one_it_picked() -> None:
    # Provenance records the resolved name, so without this line an automatic
    # choice is indistinguishable from an explicit one.
    src = _src(backends_mod)
    assert "backend='auto' resolved to" in src


def test_dropped_calc_keys_are_reported() -> None:
    # `_filter_kwargs` silently drops calc: keys a backend does not accept.
    src = _src(backends_mod)
    assert "ignored calc setting(s)" in src
    assert "k not in filtered and k not in _seeded" in src


def test_unusable_calc_keys_are_reported_but_seeded_defaults_are_not() -> None:
    """Every workflow seeds `calc_cfg` from the UMA-flavoured block, so a key UMA
    accepts and this backend does not is an untouched default, not a choice.
    Reporting those buried the line under seven or eight settings per run. A key
    outside UMA's own accepted set is genuinely unusable and still reported.
    """
    accepted_map = backends_mod._BACKEND_ACCEPTED_KEYS
    seeded = accepted_map["uma"]

    def dropped(backend: str, kwargs: dict) -> list:
        filtered = backends_mod._filter_kwargs(kwargs, accepted_map[backend])
        return sorted(k for k in kwargs if k not in filtered and k not in seeded)

    # The measured seeded noise: UMA accepts these, the others do not.
    noise = {
        "max_neigh": 30,
        "print_vram": True,
        "r_edges": True,
        "radius": 6.0,
        "task_name": "omol",
        "workers": 1,
        "workers_per_node": 1,
    }
    for backend in ("orb", "mace", "aimnet2"):
        assert set(noise) - accepted_map[backend], f"{backend} accepts all of {noise}"
        assert dropped(backend, dict(noise)) == []
        # A key no backend knows is still surfaced.
        assert dropped(backend, {**noise, "nonexistent_knob": 1}) == ["nonexistent_knob"]


def test_orb_reports_an_unhonoured_compile_request() -> None:
    # The second loader base drops `compile`; reaching it means
    # calc.compile_model=true was not honoured.
    src = _src(orb_mod)
    assert "does not accept" in src and "compile_model=true" in src


def test_freeze_links_reports_having_no_cap_hydrogen() -> None:
    # --freeze-links defaults on; with no LKH/HL the run freezes nothing.
    src = _src(utils_mod)
    assert "no LKH/HL cap hydrogen found" in src
    assert "no cap parent is frozen despite --freeze-links" in src


def test_no_deterministic_cannot_override_the_environment(capsys) -> None:
    """`--no-deterministic` does not switch off what the env var turned on."""
    import os

    from pdb2reaction.cli.common_options import _deterministic_callback

    class _Ctx:
        resilient_parsing = False

        @staticmethod
        def get_parameter_source(_name):
            return click.core.ParameterSource.COMMANDLINE

    class _Param:
        name = "deterministic"

    prev = os.environ.get("PDB2REACTION_STRICT_DETERMINISTIC")
    os.environ["PDB2REACTION_STRICT_DETERMINISTIC"] = "1"
    try:
        _deterministic_callback(_Ctx(), _Param(), False)
        err = capsys.readouterr().err
        assert "stays strictly deterministic despite" in err
    finally:
        if prev is None:
            os.environ.pop("PDB2REACTION_STRICT_DETERMINISTIC", None)
        else:
            os.environ["PDB2REACTION_STRICT_DETERMINISTIC"] = prev


def test_no_deterministic_is_silent_without_the_environment(capsys) -> None:
    import os

    from pdb2reaction.cli.common_options import _deterministic_callback

    class _Ctx:
        resilient_parsing = False

    prev = os.environ.pop("PDB2REACTION_STRICT_DETERMINISTIC", None)
    try:
        _deterministic_callback(_Ctx(), None, False)
        assert capsys.readouterr().err == ""
    finally:
        if prev is not None:
            os.environ["PDB2REACTION_STRICT_DETERMINISTIC"] = prev
