"""Unit tests for pdb2reaction.backends.apply_precision_to_calc_cfg.

Pins the --precision routing contract: per-backend kwarg translation,
the fp64 ⇒ hessian_double coupling, and the warning emitted when a
config deliberately requests the inconsistent fp64 + hessian_double=False
combination.
"""

from __future__ import annotations

import warnings

import pytest

from pdb2reaction.backends import (
    BackendError,
    apply_effective_precision,
    apply_precision_to_calc_cfg,
)


def test_fp64_routes_uma_precision_kwarg() -> None:
    cfg = {"backend": "uma"}
    apply_precision_to_calc_cfg(cfg, "fp64")
    assert cfg["precision"] == "fp64"


def test_fp64_default_forces_hessian_double_silently() -> None:
    cfg = {"backend": "uma"}
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        apply_precision_to_calc_cfg(cfg, "fp64")
    assert cfg["hessian_double"] is True
    assert caught == []


def test_fp64_overrides_explicit_false_with_warning() -> None:
    cfg = {"backend": "uma", "hessian_double": False}
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        apply_precision_to_calc_cfg(cfg, "fp64")
    assert cfg["hessian_double"] is True
    assert len(caught) == 1
    assert "fp64" in str(caught[0].message)


def test_fp32_leaves_hessian_double_untouched() -> None:
    cfg = {"backend": "uma", "hessian_double": False}
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        apply_precision_to_calc_cfg(cfg, "fp32")
    assert cfg["hessian_double"] is False
    assert caught == []


def test_fp64_rejected_for_aimnet2() -> None:
    with pytest.raises(BackendError):
        apply_precision_to_calc_cfg({"backend": "aimnet2"}, "fp64")


def test_invalid_precision_rejected() -> None:
    with pytest.raises(BackendError):
        apply_precision_to_calc_cfg({"backend": "uma"}, "fp16")


def test_orb_fp32_maps_to_reduced_tf32_precision() -> None:
    """Root cause of the ORB freq HOSP inflation: ORB's fp32 is the reduced
    'float32-high' (TF32) matmul mode, not true single precision — so an ORB
    finite-difference Hessian at fp32 carries TF32 force noise (spurious imag)."""
    cfg = {"backend": "orb"}
    apply_precision_to_calc_cfg(cfg, "fp32")
    assert cfg["precision"] == "float32-high"


def test_orb_fp64_maps_to_true_float64() -> None:
    """The fix's target: dispatching fp64 to ORB selects true float64, which
    suppresses the TF32 finite-difference noise (13 imaginary -> 3 at a c04 TS)."""
    cfg = {"backend": "orb"}
    apply_precision_to_calc_cfg(cfg, "fp64")
    assert cfg["precision"] == "float64"


def test_config_borne_precision_dispatched_when_no_cli_flag() -> None:
    """freq.py cli(): when the --precision *flag* is absent, the config's unified
    ``calc.precision`` (written by ``all`` via _write_args_yaml_with_freeze_atoms)
    must still be dispatched, so ``all --precision fp64`` reaches the ORB freq
    Hessian instead of being dropped on the way to the child stage. Mirrors the
    cli() resolution: flag wins, else the config value if it is a unified token."""
    calc_cfg = {"backend": "orb", "precision": "fp64"}
    apply_effective_precision(calc_cfg, None)  # no --precision flag → honor config
    assert calc_cfg["precision"] == "float64"


def test_already_dispatched_config_precision_is_passed_through() -> None:
    """The ('fp32','fp64') guard: a hand-edited config carrying an already
    backend-dispatched value ('float64') must NOT be re-dispatched — re-dispatch
    would raise BackendError since 'float64' is not a unified token."""
    calc_cfg = {"backend": "orb", "precision": "float64"}
    apply_effective_precision(calc_cfg, None)  # 'float64' not a unified token → untouched
    assert calc_cfg["precision"] == "float64"


def test_cli_precision_flag_wins_over_config() -> None:
    """The ``--precision`` flag overrides the config's ``calc.precision``."""
    calc_cfg = {"backend": "orb", "precision": "fp32"}
    apply_effective_precision(calc_cfg, "fp64")  # explicit flag wins
    assert calc_cfg["precision"] == "float64"


@pytest.mark.parametrize("token", [None, "auto"])
def test_unspecified_precision_defaults_to_fp64_for_orb(token) -> None:
    """No --precision flag and no config token → ORB gets true float64, not the
    reduced 'float32-high' (TF32) matmul mode."""
    calc_cfg = {"backend": "orb", "precision": token}
    apply_effective_precision(calc_cfg, None)
    assert calc_cfg["precision"] == "float64"


@pytest.mark.parametrize("token", [None, "auto"])
def test_unspecified_precision_defaults_to_fp64_for_mace(token) -> None:
    """MACE ships fp64 upstream; the unspecified default must reach
    ``default_dtype`` (a literal 'fp32' default used to dispatch float32 here and
    shadow MACE_BACKEND_DEFAULTS, breaking imaginary-mode counts)."""
    calc_cfg = {"backend": "mace", "precision": token}
    apply_effective_precision(calc_cfg, None)
    assert calc_cfg["default_dtype"] == "float64"
    assert calc_cfg["hessian_double"] is True


@pytest.mark.parametrize("token", [None, "auto"])
def test_unspecified_precision_stays_fp32_for_uma(token) -> None:
    """UMA's upstream (fairchem) baseline is fp32 and the paper's numbers rest on
    it: the orb/mace fp64 default must not leak to UMA."""
    calc_cfg = {"backend": "uma", "precision": token}
    apply_effective_precision(calc_cfg, None)
    assert calc_cfg["precision"] == "fp32"


def test_explicit_fp32_still_reaches_orb_and_mace() -> None:
    """The sentinel exists to keep an explicit ``--precision fp32`` distinguishable
    from 'unspecified'; fp32 must still downgrade both backends."""
    orb_cfg = {"backend": "orb", "precision": "auto"}
    apply_effective_precision(orb_cfg, "fp32")
    assert orb_cfg["precision"] == "float32-high"

    mace_cfg = {"backend": "mace", "precision": "auto"}
    apply_effective_precision(mace_cfg, "fp32")
    assert mace_cfg["default_dtype"] == "float32"


def test_auto_sentinel_never_survives_dispatch() -> None:
    """``"auto"`` is internal: every backend must overwrite it with a concrete
    token, since ``calc.precision`` is echoed in the run summary and propagated
    into the ``all`` pipeline's child configs."""
    for backend in ("uma", "orb", "mace", "aimnet2"):
        calc_cfg = {"backend": backend, "precision": "auto"}
        apply_effective_precision(calc_cfg, None)
        assert calc_cfg["precision"] != "auto", backend
