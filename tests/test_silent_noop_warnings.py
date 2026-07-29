"""Flags and settings that do nothing must say so.

Requested settings that cannot be honored must emit a diagnostic instead of
appearing to have been applied. These output-only guards do not change
behavior, thresholds, or defaults.
"""

from pathlib import Path

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
    assert "_dropped = sorted(k for k in kwargs if k not in filtered)" in src


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

    prev = os.environ.get("PDB2REACTION_STRICT_DETERMINISTIC")
    os.environ["PDB2REACTION_STRICT_DETERMINISTIC"] = "1"
    try:
        _deterministic_callback(_Ctx(), None, False)
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
