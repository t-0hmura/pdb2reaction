"""Strict Hessian-method dispatch at the backend boundary."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import torch
from click.testing import CliRunner

from pdb2reaction.backends.base import BackendError, MLIPCalculator


class _QuadraticCalculator(MLIPCalculator):
    def __init__(self, **kwargs):
        self.energy_force_calls = 0
        self.fd_calls = 0
        self.analytical_calls = 0
        super().__init__(print_timing=False, out_hess_torch=False, **kwargs)

    def _compute_energy_forces_ev(self, elem, coord_ang):
        self.energy_force_calls += 1
        coords = np.asarray(coord_ang, dtype=float).reshape(-1, 3)
        return 0.5 * float(np.sum(coords * coords)), -coords

    def _build_fd_hessian_cpu(self, elem, coord_ang, *, eps_ang=1.0e-3):
        self.fd_calls += 1
        return super()._build_fd_hessian_cpu(
            elem, coord_ang, eps_ang=eps_ang
        )


class _AnalyticalQuadraticCalculator(_QuadraticCalculator):
    def _supports_analytical_hessian(self) -> bool:
        return True

    def _compute_analytical_hessian_ev(self, elem, coord_ang):
        self.analytical_calls += 1
        return np.eye(3 * len(elem), dtype=float)


def test_unsupported_analytical_request_fails_before_fd_or_energy_force() -> None:
    calc = _QuadraticCalculator(hessian_calc_mode="Analytical")

    with pytest.raises(BackendError, match="Analytical Hessian is not available"):
        calc.get_hessian(["H"], np.array([0.1, 0.2, 0.3]))

    assert calc.fd_calls == 0
    assert calc.energy_force_calls == 0


def test_finite_difference_request_uses_fd() -> None:
    calc = _QuadraticCalculator(hessian_calc_mode="finitedifference")

    result = calc.get_hessian(["H"], np.array([0.1, 0.2, 0.3]))

    assert calc.hessian_calc_mode == "FiniteDifference"
    assert calc.fd_calls == 1
    assert result["hessian"].shape == (3, 3)


def test_supported_analytical_request_uses_analytical_branch() -> None:
    calc = _AnalyticalQuadraticCalculator(hessian_calc_mode="analytical")

    result = calc.get_hessian(["H"], np.array([0.1, 0.2, 0.3]))

    assert calc.hessian_calc_mode == "Analytical"
    assert calc.analytical_calls == 1
    assert calc.fd_calls == 0
    assert result["hessian"].shape == (3, 3)


def test_torch_hessian_trim_and_symmetrization_match_dense_reference() -> None:
    calc = _AnalyticalQuadraticCalculator(
        freeze_atoms=[1],
        return_partial_hessian=True,
        hessian_double=True,
    )
    source = torch.arange(36, dtype=torch.float32).reshape(6, 6)
    expected = source[:3, :3].to(torch.float64)
    expected = 0.5 * (expected + expected.T)
    from pdb2reaction.backends.base import H_EVAA_2_AU

    expected *= H_EVAA_2_AU

    trimmed = calc._apply_active_trim_torch(source, n_atoms=2)
    actual = calc._au_hessian_torch(trimmed)

    torch.testing.assert_close(actual, expected)
    assert actual.dtype == torch.float64


@pytest.mark.parametrize("mode", ["", "FD", "analytic", "typo"])
def test_unknown_direct_api_method_is_rejected(mode: str) -> None:
    with pytest.raises(BackendError, match="Choose from: FiniteDifference, Analytical"):
        _QuadraticCalculator(hessian_calc_mode=mode)


def test_factory_rejects_unknown_method_before_backend_import(monkeypatch) -> None:
    import pdb2reaction.backends as backends

    imported = False

    def forbidden_import(*args, **kwargs):
        nonlocal imported
        imported = True
        raise AssertionError("invalid method must fail before backend import")

    monkeypatch.setattr(backends, "_import_cls", forbidden_import)

    with pytest.raises(BackendError, match="Choose from: FiniteDifference, Analytical"):
        backends.create_calculator(
            backend="orb", hessian_calc_mode="typo"
        )
    assert imported is False


def test_custom_sp_analytical_failure_creates_no_hessian(
    tmp_path: Path,
) -> None:
    from pdb2reaction.cli import cli as root_cli

    inp = tmp_path / "geom.xyz"
    inp.write_text("1\nframe\nHe 0.1 0.0 0.0\n")
    calc_file = tmp_path / "harmonic_calc.py"
    calc_file.write_text(
        "from ase.calculators.calculator import Calculator, all_changes\n"
        "import numpy as np\n"
        "class Harmonic(Calculator):\n"
        "    implemented_properties = ['energy', 'forces']\n"
        "    def calculate(self, atoms=None, properties=None, "
        "system_changes=all_changes):\n"
        "        super().calculate(atoms, properties, system_changes)\n"
        "        x = np.asarray(atoms.positions, dtype=float)\n"
        "        self.results = {'energy': 0.5 * float((x*x).sum()), "
        "'forces': -x}\n"
        "def get_calculator(**kwargs):\n"
        "    return Harmonic()\n"
    )
    out_dir = tmp_path / "out"

    result = CliRunner().invoke(
        root_cli,
        [
            "sp",
            "-i",
            str(inp),
            "-q",
            "0",
            "-m",
            "1",
            "--calc-file",
            str(calc_file),
            "--hess",
            "--hessian-calc-mode",
            "Analytical",
            "-o",
            str(out_dir),
        ],
    )

    assert result.exit_code != 0
    assert "Analytical Hessian is not available" in result.output
    assert not (out_dir / "hessian.npy").exists()
