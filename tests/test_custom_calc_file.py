"""Smoke + unit tests for the ``--calc-file`` custom ASE Calculator backend.

Exercises loading an arbitrary ASE Calculator from a user Python file and using
it as the energy/gradient source (the R1 reviewer point: couple GFN-xTB / DFTB+
/ any ASE engine). Uses a dependency-free, element-agnostic toy harmonic
calculator so the test needs no MLIP weights or GPU.
"""

from __future__ import annotations

import json
import textwrap
from pathlib import Path

import numpy as np
import pytest
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
    from pysisyphus.constants import ANG2BOHR

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

    # Default finite-difference Hessian emits a tagged timing line.  Direct
    # library use must not require CLI bootstrap for that output path.
    hessian_result = cc.get_hessian(elem, coord.reshape(-1) * ANG2BOHR)
    assert np.asarray(hessian_result["hessian"]).shape == (9, 9)


def test_calculator_class_export_is_instantiated(tmp_path: Path) -> None:
    from ase.calculators.calculator import Calculator

    from pdb2reaction.backends.custom import load_ase_calculator

    calc_file = _write(
        tmp_path / "class_calc.py",
        "from ase.calculators.emt import EMT\nget_calculator = EMT\n",
    )

    assert isinstance(load_ase_calculator(str(calc_file)), Calculator)


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


def test_custom_calculator_provenance_uses_loader_default_factory() -> None:
    from pdb2reaction.core.utils import calculator_provenance

    provenance = calculator_provenance({
        "backend": "custom",
        "calc_file": "/tmp/toy.py",
    })

    assert provenance == {
        "mlip_backend": "custom",
        "mlip_model": "toy.py:get_calculator",
        "mlip_precision": None,
    }


@pytest.mark.parametrize(
    ("backend", "cfg", "expected"),
    [
        ("uma", {"precision": "fp64"}, "fp64"),
        ("orb", {"precision": "float32-high"}, "fp32"),
        ("orb", {"precision": "float32-highest"}, "fp32"),
        ("mace", {"default_dtype": "float64"}, "fp64"),
        ("aimnet2", {}, "fp32"),
    ],
)
def test_calculator_provenance_uses_canonical_precision(
    backend: str, cfg: dict, expected: str
) -> None:
    from pdb2reaction.core.utils import calculator_provenance

    assert calculator_provenance({"backend": backend, **cfg})["mlip_precision"] == expected


def test_calculator_provenance_uses_backend_defaults() -> None:
    from pdb2reaction.core.utils import calculator_provenance

    assert calculator_provenance({"backend": "orb"}) == {
        "mlip_backend": "orb",
        "mlip_model": "orb_v3_conservative_omol",
        "mlip_precision": "fp64",
    }


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
            "--out-json",
            "-o", str(out_dir),
        ],
        catch_exceptions=False,
    )
    assert result.exit_code == 0, result.output
    assert "custom" in result.output
    payload = json.loads((out_dir / "result.json").read_text(encoding="utf-8"))
    assert payload["mlip_backend"] == "custom"
    assert payload["mlip_model"] == "toy.py:get_calculator"
    assert payload["mlip_precision"] is None
    assert payload["custom_calculator"] == "toy.py:get_calculator"
    assert payload["n_atoms"] == 3


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


def test_opt_json_uses_yaml_resolved_custom_provenance(tmp_path: Path) -> None:
    calc_file = _write(tmp_path / "toy.py", TOY_CALC)
    xyz = _write(tmp_path / "water.xyz", WATER_XYZ)
    config = _write(
        tmp_path / "config.yaml",
        f"calc:\n  calc_file: {calc_file}\n  calc_factory: get_calculator\n",
    )
    out_dir = tmp_path / "opt_out"

    result = CliRunner().invoke(
        root_cli,
        [
            "opt",
            "-i",
            str(xyz),
            "--config",
            str(config),
            "-q",
            "0",
            "-m",
            "1",
            "--max-cycles",
            "1",
            "--out-json",
            "-o",
            str(out_dir),
        ],
        catch_exceptions=False,
    )

    assert result.exit_code == 0, result.output
    payload = json.loads((out_dir / "result.json").read_text(encoding="utf-8"))
    assert payload["backend"] == payload["mlip_backend"] == "custom"
    assert payload["model"] == payload["mlip_model"] == "toy.py:get_calculator"
    assert payload["mlip_precision"] is None
