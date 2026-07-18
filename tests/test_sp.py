"""Single-point workflow regression tests."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest
from click.testing import CliRunner


class _FakeGeometry:
    atoms = ["C"]
    cart_coords = np.zeros(3, dtype=float)

    def __init__(self) -> None:
        self._calc = None
        self.cached_results = None

    def set_calculator(self, calc) -> None:
        self._calc = calc

    @property
    def energy(self) -> float:
        return 1.25

    @property
    def gradient(self) -> np.ndarray:
        return np.zeros(3, dtype=float)

    def set_results(self, results) -> None:
        self.cached_results = results


class _FakeCalculator:
    def get_hessian(self, atoms, coords):
        assert atoms == ["C"]
        assert np.array_equal(coords, np.zeros(3))
        return {
            "energy": 1.25,
            "forces": np.zeros(3),
            "hessian": np.eye(3),
        }


@pytest.mark.parametrize(
    ("backend", "mode_args", "expected_mode"),
    [
        ("orb", ["--hessian-calc-mode", "Analytical"], "Analytical"),
        ("mace", ["--hessian-calc-mode", "Analytical"], "Analytical"),
        ("aimnet2", ["--hessian-calc-mode", "Analytical"], "Analytical"),
        ("orb", [], "FiniteDifference"),
        ("uma", [], "Analytical"),
    ],
)
def test_sp_propagates_resolved_hessian_mode_to_backend(
    monkeypatch,
    tmp_path: Path,
    backend: str,
    mode_args: list[str],
    expected_mode: str,
) -> None:
    from pdb2reaction.cli import cli as root_cli
    from pdb2reaction.workflows import sp

    inp = tmp_path / "geom.xyz"
    inp.write_text("1\nframe\nC 0.0 0.0 0.0\n")
    out = tmp_path / "out"
    created: list[dict] = []

    monkeypatch.setattr(sp, "geom_loader", lambda *_args, **_kwargs: _FakeGeometry())

    def fake_create_calculator(**kwargs):
        created.append(dict(kwargs))
        return _FakeCalculator()

    monkeypatch.setattr(sp, "create_calculator", fake_create_calculator)

    result = CliRunner().invoke(
        root_cli,
        [
            "sp", "-i", str(inp), "-q", "0", "-m", "1", "-b", backend,
            "--hess", "--out-json", "-o", str(out), *mode_args,
        ],
    )

    assert result.exit_code == 0, result.output
    assert len(created) == 1
    assert created[0]["hessian_calc_mode"] == expected_mode
    assert np.array_equal(np.load(out / "hessian.npy"), np.eye(3))
    payload = json.loads((out / "result.json").read_text())
    assert payload["status"] == "ok"


def test_sp_converts_yaml_freeze_atoms_to_internal_indices(
    monkeypatch, tmp_path: Path,
) -> None:
    from pdb2reaction.cli import cli as root_cli
    from pdb2reaction.workflows import sp

    inp = tmp_path / "geom.xyz"
    inp.write_text("2\nframe\nC 0.0 0.0 0.0\nC 0.0 0.0 1.0\n")
    cfg = tmp_path / "config.yaml"
    cfg.write_text("geom:\n  freeze_atoms: [1]\n")
    loaded: list[dict] = []
    created: list[dict] = []

    def fake_loader(*_args, **kwargs):
        loaded.append(dict(kwargs))
        return _FakeGeometry()

    def fake_create_calculator(**kwargs):
        created.append(dict(kwargs))
        return _FakeCalculator()

    monkeypatch.setattr(sp, "geom_loader", fake_loader)
    monkeypatch.setattr(sp, "create_calculator", fake_create_calculator)

    result = CliRunner().invoke(
        root_cli,
        [
            "sp", "-i", str(inp), "-q", "0", "-m", "1",
            "--config", str(cfg), "--hess", "-o", str(tmp_path / "out"),
        ],
    )

    assert result.exit_code == 0, result.output
    assert loaded[0]["freeze_atoms"] == [0]
    assert created[0]["freeze_atoms"] == [0]
    assert created[0]["return_partial_hessian"] is True
