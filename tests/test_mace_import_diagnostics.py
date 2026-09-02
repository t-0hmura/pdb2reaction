"""MACE import diagnostics preserve installed-runtime failures."""

from __future__ import annotations

import builtins
import importlib.metadata

import pytest

from pdb2reaction.backends.base import BackendError
from pdb2reaction.backends.mace import MACECalculator


def test_installed_mace_import_failure_reports_root_cause(monkeypatch) -> None:
    real_import = builtins.__import__

    def _import(name, *args, **kwargs):
        if name == "mace.calculators":
            raise RuntimeError("lmdb.h is unavailable")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", _import)
    monkeypatch.setattr(
        importlib.metadata,
        "version",
        lambda name: "0.3.16" if name == "mace-torch" else "unknown",
    )

    calculator = object.__new__(MACECalculator)
    with pytest.raises(BackendError) as exc_info:
        calculator._build_calc("MACE-OMOL-0")

    message = str(exc_info.value)
    assert "Installed mace-torch 0.3.16" in message
    assert "RuntimeError: lmdb.h is unavailable" in message
    assert "pip install mace-torch" not in message
