"""Pins the announcement that brackets the first load of each MLIP model.

Weight downloads happen inside the backend constructor and print nothing of
their own, so a first run looked indistinguishable from a hang. Both factories
bracket the construction with ``Preparing MLIP model (...)`` / ``Done.`` once per
distinct ``(backend, model)``: a cached model shows both lines back to back, a
download stops between them, a second construction of the same model stays
silent, and a construction that raises is not remembered as announced.
"""

from __future__ import annotations

import pytest

from pdb2reaction import backends


class _StubBackend:
    """Stands in for a real calculator class; records the kwargs it received."""

    def __init__(self, **kwargs):
        self.kwargs = kwargs


@pytest.fixture(autouse=True)
def _fresh_announcement_state(monkeypatch):
    """Each test starts with nothing announced (the set is process-global)."""
    monkeypatch.setattr(backends, "_ANNOUNCED_MODEL_LOADS", set())


@pytest.fixture()
def _stub_cls(monkeypatch):
    monkeypatch.setattr(
        backends, "_import_cls", lambda backend, cls_key: _StubBackend
    )


@pytest.mark.parametrize(
    "factory",
    ["create_calculator", "create_ase_calculator"],
)
def test_first_load_is_bracketed_then_silent(factory, capsys, _stub_cls) -> None:
    calc = getattr(backends, factory)("uma", model="uma-s-1p2", device="cpu")
    assert isinstance(calc, _StubBackend)
    out = capsys.readouterr().out
    assert "[backend] Preparing MLIP model (uma / uma-s-1p2)..." in out
    assert "[backend] Done." in out
    # The notice must precede the load it brackets, not follow it.
    assert out.index("Preparing MLIP model") < out.index("[backend] Done.")

    # Same (backend, model) again: already loaded once, so no second notice.
    getattr(backends, factory)("uma", model="uma-s-1p2", device="cpu")
    assert "Preparing MLIP model" not in capsys.readouterr().out


def test_each_model_is_announced_once(capsys, _stub_cls) -> None:
    backends.create_ase_calculator("uma", model="uma-s-1p2", device="cpu")
    backends.create_ase_calculator("uma", model="uma-m-1p1", device="cpu")
    out = capsys.readouterr().out
    assert out.count("Preparing MLIP model") == 2
    assert "(uma / uma-s-1p2)" in out and "(uma / uma-m-1p1)" in out


def test_failed_load_is_announced_again(capsys, monkeypatch) -> None:
    class _Boom:
        def __init__(self, **kwargs):
            raise RuntimeError("checkpoint download failed")

    monkeypatch.setattr(backends, "_import_cls", lambda backend, cls_key: _Boom)
    with pytest.raises(RuntimeError):
        backends.create_ase_calculator("uma", model="uma-s-1p2", device="cpu")
    out = capsys.readouterr().out
    assert "Preparing MLIP model" in out
    assert "[backend] Done." not in out  # nothing was loaded

    # The retry must announce again: the first attempt never completed.
    monkeypatch.setattr(backends, "_import_cls", lambda backend, cls_key: _StubBackend)
    backends.create_ase_calculator("uma", model="uma-s-1p2", device="cpu")
    retry = capsys.readouterr().out
    assert "Preparing MLIP model" in retry and "[backend] Done." in retry
