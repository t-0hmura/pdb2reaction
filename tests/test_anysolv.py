"""AnySolv uses one additive energy surface in ordinary and ASE calculations."""
from types import SimpleNamespace
import sys

import numpy as np
import pytest
import torch
from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes

from pdb2reaction.backends.anysolv import AnySolvCorrection, wrap_ase_calculator
from pdb2reaction.backends.solvent import SolventCorrectedCalculator
from pysisyphus.constants import AU2EV, BOHR2ANG


@pytest.fixture
def predictor(monkeypatch):
    calls = []
    def predict(structure, **kwargs):
        numbers, x = structure
        calls.append((numbers.copy(), dict(kwargs)))
        return float(0.5 * np.sum(x*x) + 2), -x
    monkeypatch.setitem(sys.modules, "anysolv", SimpleNamespace(predict_solvation_energy=predict))
    return calls


@pytest.mark.parametrize("checkpoint", ["model_compact", "model_moe", "weights.pt"])
def test_model_charge_spin_precision_and_force_hessian(predictor, checkpoint):
    correction = AnySolvCorrection(solvent="acetonitrile", charge=-1, spin=2,
                                    checkpoint=checkpoint, dtype="fp64", device="cpu")
    x = np.array([[0.1, 0.2, -0.3], [0.9, 0.1, 0.4]])
    e, f, h = correction.compute(["O", "H"], x, need_hessian=True)
    assert e == pytest.approx(2 + .5 * np.sum(x*x))
    np.testing.assert_allclose(f, -x)
    np.testing.assert_allclose(h, np.eye(6), atol=1e-12)
    numbers, settings = predictor[0]
    assert numbers.tolist() == [8, 1]
    assert settings == dict(charge=-1, spin=2, solvent="acetonitrile", checkpoint=checkpoint,
                            device="cpu", dtype=torch.float64, inference_settings="default")


def test_compact_default_is_explicit(predictor):
    AnySolvCorrection(solvent="water").compute(["He"], np.zeros((1, 3)))
    assert predictor[0][1]["checkpoint"] == "model_compact"


def test_optional_dependency_error_is_actionable(monkeypatch):
    monkeypatch.setitem(sys.modules, "anysolv", None)
    correction = AnySolvCorrection(solvent="water")
    with pytest.raises(ImportError, match=r"pdb2reaction\[anysolv\]"):
        correction.compute(["He"], np.zeros((1, 3)))


@pytest.mark.parametrize("step", [0, -1, float("nan"), float("inf")])
def test_invalid_hessian_step_rejected(step):
    with pytest.raises(ValueError, match="positive and finite"):
        AnySolvCorrection(solvent="water", hessian_step=step)


class ZeroBase:
    charge, mult, device_str, freeze_atoms = 0, 1, "cpu", []
    hessian_calc_mode, return_partial_hessian = "Analytical", False
    hessian_double, out_hess_torch, print_timing = True, False, False
    def get_energy(self, elem, coords):
        return {"energy": 0.}
    def get_forces(self, elem, coords):
        return {"energy": 0., "forces": np.zeros(len(coords))}
    def get_hessian(self, elem, coords):
        return {"energy": 0., "forces": np.zeros(len(coords)), "hessian": np.zeros((len(coords),len(coords)))}


class ZeroASE(Calculator):
    implemented_properties = ["energy", "forces"]
    def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        self.results = {"energy": 0., "forces": np.zeros((len(atoms),3))}


def test_regular_and_ase_surfaces_units_and_constraints(predictor):
    symbols = ["O", "H"]
    x = np.array([[.1,.2,.3], [.9,.2,.1]])
    ordinary = SolventCorrectedCalculator(ZeroBase(), solvent="water", solvent_model="anysolv")
    atoms = Atoms(symbols, positions=x)
    atoms.calc = wrap_ase_calculator(ZeroASE(), solvent="water")
    result = ordinary.get_hessian(symbols, (x/BOHR2ANG).reshape(-1))
    assert result["energy"] * AU2EV == pytest.approx(atoms.get_potential_energy())
    np.testing.assert_allclose(result["forces"].reshape(-1,3)*AU2EV/BOHR2ANG,
                               atoms.get_forces(), atol=1e-12)
    np.testing.assert_allclose(result["hessian"]*AU2EV/BOHR2ANG**2, np.eye(6), atol=1e-12)
    ordinary.freeze_atoms = [0]
    constrained = ordinary.get_hessian(symbols, (x/BOHR2ANG).reshape(-1))
    np.testing.assert_array_equal(constrained["forces"][:3], np.zeros(3))
    np.testing.assert_array_equal(constrained["hessian"][:3], np.zeros((3,6)))


def test_disabled_correction_does_not_import_anysolv(monkeypatch):
    monkeypatch.setitem(sys.modules, "anysolv", None)
    calc = SolventCorrectedCalculator(ZeroBase(), solvent="none", solvent_model="anysolv")
    assert calc.get_energy(["He"], np.zeros(3))["energy"] == 0


@pytest.mark.parametrize("backend,expected",[("uma",torch.float32),("orb",torch.float64),
                                             ("mace",torch.float64),("aimnet2",torch.float32)])
def test_factory_default_dtype_matches_base_in_both_routes(monkeypatch,predictor,backend,expected):
    import pdb2reaction.backends as backends
    monkeypatch.setattr(backends,"_import_cls",lambda name,kind:
                        (lambda **kwargs: ZeroBase()) if kind=="pysis_cls" else ZeroASE)
    settings = dict(solvent="water",solvent_model="anysolv",device="cpu")
    ordinary = backends.create_calculator(backend,**settings)
    ordinary.get_energy(["He"],np.zeros(3))
    assert predictor[-1][1]["dtype"] is expected
    atoms = Atoms("He",positions=np.zeros((1,3)))
    atoms.calc = backends.create_ase_calculator(backend,**settings)
    atoms.get_potential_energy()
    assert predictor[-1][1]["dtype"] is expected


def test_checkpoint_and_hessian_settings_are_recorded_for_resume():
    from pdb2reaction.core.utils import calculator_provenance
    config = dict(backend="uma",solvent="water",solvent_model="anysolv",precision="fp64",
                  anysolv_checkpoint="model_moe",anysolv_hessian_step=.002)
    result = calculator_provenance(config)["solvation_correction"]
    assert result == dict(method="anysolv",solvent="water",checkpoint="model_moe",dtype="float64",
                          hessian_method="force-finite-difference",hessian_step_angstrom=.002)


@pytest.mark.parametrize("factory",["create_calculator","create_ase_calculator"])
def test_unsupported_dft_combination_fails_before_model_loading(monkeypatch,factory):
    import pdb2reaction.backends as backends
    def forbidden(*args):
        raise AssertionError("DFT must not be loaded for unsupported AnySolv")
    monkeypatch.setattr(backends,"_import_cls",forbidden)
    with pytest.raises(backends.BackendError,match="use PCM/SMD"):
        getattr(backends,factory)("dft",solvent="water",solvent_model="anysolv")
