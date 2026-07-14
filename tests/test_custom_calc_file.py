"""Smoke + unit tests for the ``--calc-file`` custom ASE Calculator backend.

Exercises loading an arbitrary ASE Calculator from a user Python file and using
it as the energy/gradient source (the R1 reviewer point: couple GFN-xTB / DFTB+
/ any ASE engine). Uses a dependency-free, element-agnostic toy harmonic
calculator so the test needs no MLIP weights or GPU.
"""

from __future__ import annotations

import textwrap
from pathlib import Path

import numpy as np
from click.testing import CliRunner

from pdb2reaction.cli import cli as root_cli
from pdb2reaction.backends import apply_calc_file_to_calc_cfg

# Toy ASE calculator: V = 0.5 * sum(pos**2) eV, F = -pos eV/Ang. No external
# deps beyond ASE+numpy, no element restrictions.
TOY_CALC = textwrap.dedent(
    '''
    import numpy as np
    from ase.calculators.calculator import Calculator, all_changes


    class ToyHarmonic(Calculator):
        implemented_properties = ["energy", "forces"]

        def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
            super().calculate(atoms, properties, system_changes)
            pos = atoms.get_positions()
            self.results["energy"] = 0.5 * float(np.sum(pos ** 2))
            self.results["forces"] = -pos


    def get_calculator(charge=0, spin=1, device="auto", **kwargs):
        return ToyHarmonic()
    '''
)

WATER_XYZ = "3\n\nO 0.0 0.0 0.0\nH 0.0 0.757 0.587\nH 0.0 -0.757 0.587\n"


def _write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


def test_load_ase_calculator_and_compute(tmp_path: Path) -> None:
    from pdb2reaction.backends.custom import CustomCalculator, load_ase_calculator

    calc_file = _write(tmp_path / "toy.py", TOY_CALC)

    ase_calc = load_ase_calculator(str(calc_file))
    assert hasattr(ase_calc, "get_potential_energy")
    assert hasattr(ase_calc, "get_forces")

    cc = CustomCalculator(calc_file=str(calc_file), charge=0, spin=1)
    elem = ["O", "H", "H"]
    coord = np.array(
        [[0.0, 0.0, 0.0], [0.0, 0.757, 0.587], [0.0, -0.757, 0.587]]
    )
    e_ev, f_ev = cc._compute_energy_forces_ev(elem, coord)
    assert abs(e_ev - 0.5 * float(np.sum(coord ** 2))) < 1e-9
    assert np.allclose(f_ev, -coord)


def test_load_ase_calculator_errors(tmp_path: Path) -> None:
    from pdb2reaction.backends.base import BackendError
    from pdb2reaction.backends.custom import load_ase_calculator

    import pytest

    missing_factory = _write(tmp_path / "no_factory.py", "x = 1\n")
    with pytest.raises(BackendError):
        load_ase_calculator(str(missing_factory))

    not_a_calc = _write(
        tmp_path / "bad.py", "def get_calculator(**kw):\n    return 42\n"
    )
    with pytest.raises(BackendError):
        load_ase_calculator(str(not_a_calc))


def test_apply_calc_file_switches_backend() -> None:
    cfg = {"backend": "uma", "model": "uma-s-1p1"}
    apply_calc_file_to_calc_cfg(cfg, "/path/to/toy.py", "get_calculator")
    assert cfg["backend"] == "custom"
    assert cfg["calc_file"] == "/path/to/toy.py"
    assert cfg["calc_factory"] == "get_calculator"
    assert "model" not in cfg  # inherited UMA default dropped

    # No calc-file -> the --backend selection is untouched.
    cfg2 = {"backend": "uma"}
    apply_calc_file_to_calc_cfg(cfg2, None, None)
    assert cfg2["backend"] == "uma"


def test_sp_cli_with_calc_file(tmp_path: Path) -> None:
    calc_file = _write(tmp_path / "toy.py", TOY_CALC)
    xyz = _write(tmp_path / "water.xyz", WATER_XYZ)
    out_dir = tmp_path / "sp_out"

    runner = CliRunner()
    result = runner.invoke(
        root_cli,
        [
            "sp", "-i", str(xyz),
            "--calc-file", str(calc_file),
            "-q", "0", "-m", "1",
            "-o", str(out_dir),
        ],
        catch_exceptions=False,
    )
    assert result.exit_code == 0, result.output
    assert "custom" in result.output


def test_sp_yaml_custom_factory_is_not_overwritten_by_cli_default(
    tmp_path: Path,
) -> None:
    calc_file = _write(
        tmp_path / "toy_build.py", TOY_CALC.replace("def get_calculator", "def build")
    )
    xyz = _write(tmp_path / "water.xyz", WATER_XYZ)
    config = _write(
        tmp_path / "config.yaml",
        f"calc:\n  calc_file: {calc_file}\n  calc_factory: build\n",
    )

    result = CliRunner().invoke(
        root_cli,
        [
            "sp",
            "-i",
            str(xyz),
            "--config",
            str(config),
            "-q",
            "0",
            "-m",
            "1",
            "--show-config",
            "--dry-run",
        ],
        catch_exceptions=False,
    )
    assert result.exit_code == 0, result.output
    assert "backend: custom" in result.output
    assert "calc_factory: build" in result.output
