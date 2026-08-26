"""Pins the announcement that brackets the first load of each MLIP model."""

from __future__ import annotations

import pytest

from pdb2reaction import backends
from pdb2reaction.core.output import mlip_model_label


class _StubBackend:
    def __init__(self, **kwargs):
        self.kwargs = kwargs


@pytest.fixture(autouse=True)
def _fresh_announcement_state(monkeypatch):
    monkeypatch.setattr(backends, "_ANNOUNCED_MODEL_LOADS", set())


@pytest.fixture()
def _stub_cls(monkeypatch):
    monkeypatch.setattr(
        backends, "_import_cls", lambda backend, cls_key: _StubBackend
    )


@pytest.mark.parametrize("factory", ["create_calculator", "create_ase_calculator"])
def test_first_load_is_bracketed_then_silent(factory, capsys, _stub_cls) -> None:
    calc = getattr(backends, factory)("uma", model="uma-s-1p2", device="cpu")
    assert isinstance(calc, _StubBackend)
    out = capsys.readouterr().out
    assert "[backend] Preparing MLIP model (UMA / UMA-S-1.2 (OMol))..." in out
    assert "[backend] Done." in out
    assert out.index("Preparing MLIP model") < out.index("[backend] Done.")
    lines = out.splitlines()
    preparing = next(i for i, line in enumerate(lines) if "Preparing MLIP model" in line)
    done = lines.index("[backend] Done.")
    assert preparing == 0 or lines[preparing - 1] == ""
    assert done + 1 < len(lines) and lines[done + 1] == ""

    getattr(backends, factory)("uma", model="uma-s-1p2", device="cpu")
    assert "Preparing MLIP model" not in capsys.readouterr().out


def test_each_model_is_announced_once(capsys, _stub_cls) -> None:
    backends.create_ase_calculator("uma", model="uma-s-1p2", device="cpu")
    backends.create_ase_calculator("uma", model="uma-m-1p1", device="cpu")
    out = capsys.readouterr().out
    assert out.count("Preparing MLIP model") == 2
    assert "(UMA / UMA-S-1.2 (OMol))" in out
    assert "(UMA / UMA-M-1.1 (OMol))" in out


def test_uma_task_is_part_of_the_model_announcement(capsys, _stub_cls) -> None:
    backends.create_ase_calculator(
        "uma", model="uma-s-1p2", task_name="omol", device="cpu"
    )
    backends.create_ase_calculator(
        "uma", model="uma-s-1p2", task_name="omat", device="cpu"
    )
    out = capsys.readouterr().out
    assert "UMA-S-1.2 (OMol)" in out
    assert "UMA-S-1.2 (OMat)" in out


def test_uma_default_and_explicit_omol_share_one_announcement(
    capsys, _stub_cls
) -> None:
    backends.create_ase_calculator("uma", model="uma-s-1p2", device="cpu")
    backends.create_ase_calculator(
        "uma", model="uma-s-1p2", task_name="omol", device="cpu"
    )
    assert capsys.readouterr().out.count("Preparing MLIP model") == 1


def test_public_model_labels_canonicalize_supported_aliases() -> None:
    assert mlip_model_label("mace", "off:small") == "MACE-OFF23-small"
    assert (
        mlip_model_label("orb", "orb-v3-conservative-omol")
        == "ORB-v3-conservative-OMol"
    )


def test_failed_load_is_announced_again(capsys, monkeypatch) -> None:
    class _Boom:
        def __init__(self, **kwargs):
            raise RuntimeError("checkpoint download failed")

    monkeypatch.setattr(backends, "_import_cls", lambda backend, cls_key: _Boom)
    with pytest.raises(RuntimeError):
        backends.create_ase_calculator("uma", model="uma-s-1p2", device="cpu")
    out = capsys.readouterr().out
    assert "Preparing MLIP model" in out
    assert "[backend] Done." not in out

    monkeypatch.setattr(backends, "_import_cls", lambda backend, cls_key: _StubBackend)
    backends.create_ase_calculator("uma", model="uma-s-1p2", device="cpu")
    retry = capsys.readouterr().out
    assert "Preparing MLIP model" in retry and "[backend] Done." in retry
