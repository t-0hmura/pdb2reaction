"""Pins the default models and the routes for the newly supported models.

Defaults follow each backend's recommendation (UMA uma-s-1p2p1, ORB orbmol-v2,
MACE MACE-OMOL-0, AIMNet2 aimnet2). OrbMol-v2's name carries no
"conservative", so ORB checks the loaded regressor's class; MACE-POLAR-1 loads
through mace_polar. None of these tests needs a GPU, a network or a backend
package.
"""

from __future__ import annotations

import inspect
import sys
from types import ModuleType, SimpleNamespace

import pytest

from pdb2reaction.backends.base import BackendError
from pdb2reaction.core import defaults
from pdb2reaction.core.output import mlip_model_label


def test_default_models() -> None:
    assert defaults.DEFAULT_UMA_MODEL == "uma-s-1p2p1"
    assert defaults.CALC_KW_DEFAULT["model"] == "uma-s-1p2p1"
    assert defaults.ORB_BACKEND_DEFAULTS["model"] == "orbmol-v2"
    assert defaults.MACE_BACKEND_DEFAULTS["model"] == "MACE-OMOL-0"
    from pdb2reaction.backends.orb import OrbASECalculator, OrbCalculator

    for factory in (OrbCalculator.__init__, OrbASECalculator.__new__):
        assert inspect.signature(factory).parameters["model"].default == "orbmol-v2"


@pytest.mark.parametrize(
    ("backend", "model", "label"),
    [
        ("uma", "uma-s-1p2p1", "UMA-S-1.2.1 (OMol)"),
        ("uma", "uma-s-1p2p1-ft-d4rev", "UMA-S-1.2.1-FT-D4rev (OMol)"),
        ("uma", "uma-s-1p2", "UMA-S-1.2 (OMol)"),
        ("orb", "orbmol-v2", "OrbMol-v2"),
        ("orb", "orb_v3_conservative_omol", "ORB-v3-conservative-OMol"),
        ("mace", "MACE-POLAR-1-M", "MACE-POLAR-1-M"),
    ],
)
def test_model_labels(backend, model, label) -> None:
    assert mlip_model_label(backend, model, "omol") == label


def test_orb_name_precheck_and_regressor_check() -> None:
    from pdb2reaction.backends.orb import (
        _is_conservative_orb_model,
        _require_conservative_regressor,
    )

    assert _is_conservative_orb_model("orbmol-v2")
    assert _is_conservative_orb_model("orb_v3_conservative_omol")
    assert not _is_conservative_orb_model("orb_v3_direct_omol")

    class ConservativeForcefieldRegressor:
        pass

    class DirectForcefieldRegressor:
        pass

    _require_conservative_regressor(ConservativeForcefieldRegressor(), "orbmol-v2")
    with pytest.raises(BackendError, match="conservative"):
        _require_conservative_regressor(DirectForcefieldRegressor(), "orb-v2")


def test_orb_resolves_orbmol_v2_and_explains_old_orb_models() -> None:
    from pdb2reaction.backends.orb import OrbCalculator

    loader = object()
    calc = object.__new__(OrbCalculator)
    calc._pretrained = SimpleNamespace(ORB_PRETRAINED_MODELS={"orbmol-v2": loader})
    assert calc._resolve_loader("orbmol-v2") is loader

    calc._pretrained = SimpleNamespace(ORB_PRETRAINED_MODELS={"orb-v3-conservative-omol": loader})
    with pytest.raises(BackendError, match="orb-models >= 0.7.0"):
        calc._resolve_loader("orbmol-v2")


@pytest.mark.parametrize(
    ("spec", "key"),
    [("MACE-POLAR-1-S", "polar-1-s"), ("mace-polar-1-m", "polar-1-m"), ("polar-1-L", "polar-1-l")],
)
def test_mace_polar_names_route_to_mace_polar(monkeypatch, spec, key) -> None:
    from pdb2reaction.backends.mace import MACECalculator

    captured = {}
    calculators = ModuleType("mace.calculators")
    for name in ("mace_anicc", "mace_mp", "mace_off", "mace_omol"):
        setattr(calculators, name, lambda **_kw: pytest.fail("wrong MACE factory"))
    foundations = ModuleType("mace.calculators.foundations_models")
    foundations.mace_mp_urls = {}
    foundations.mace_polar = lambda **kw: captured.update(kw) or "polar-calc"
    monkeypatch.setitem(sys.modules, "mace", ModuleType("mace"))
    monkeypatch.setitem(sys.modules, "mace.calculators", calculators)
    monkeypatch.setitem(sys.modules, "mace.calculators.foundations_models", foundations)

    calc = object.__new__(MACECalculator)
    calc.device_str, calc.default_dtype = "cpu", "float64"
    assert calc._build_calc(spec) == "polar-calc"
    assert captured == {"model": key, "device": "cpu", "default_dtype": "float64"}


@pytest.mark.parametrize(
    ("backend", "model", "expected"),
    [
        ("orb", "orbmol-v2", {"orb_v3", "omol25"}),
        ("orb", "orb_v3_conservative_omol", {"orb_v3", "omol25"}),
        ("mace", "MACE-POLAR-1-M", {"mace", "mace_design", "omol25"}),
        ("uma", "uma-s-1p2p1-ft-d4rev", {"uma", "omol25"}),
    ],
)
def test_new_models_cite_their_model_and_dataset(backend, model, expected) -> None:
    from pdb2reaction.io.summary import _method_citation_record_keys

    keys = set(_method_citation_record_keys({"mlip_backend": backend, "mlip_model": model}))
    assert expected <= keys
