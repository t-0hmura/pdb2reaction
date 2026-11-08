"""Optional additive AnySolv correction, in eV and Angstrom units."""

from __future__ import annotations

import numpy as np


class AnySolvCorrection:
    """Reuse the public inference API for energies, forces and force-difference Hessians."""

    def __init__(self, *, solvent, charge=0, spin=1, checkpoint="model_compact",
                 device="cpu", dtype="float64", hessian_step=1e-3):
        self.solvent = solvent
        self.charge, self.spin = int(charge), int(spin)
        self.checkpoint = str(checkpoint)
        self.device = str(device)
        token = str(dtype).lower()
        if token in {"float64", "fp64", "double", "torch.float64"}:
            self.dtype = "float64"
        elif token in {"float32", "fp32", "single", "torch.float32", "tf32", "float32-high", "float32-highest"}:
            self.dtype = "float32"
        else:
            raise ValueError(f"Unsupported AnySolv precision: {dtype!r}")
        self.hessian_step = float(hessian_step)
        if checkpoint is None or not self.checkpoint.strip():
            raise ValueError("AnySolv checkpoint must not be empty")
        if not np.isfinite(self.hessian_step) or self.hessian_step <= 0:
            raise ValueError("AnySolv Hessian step must be positive and finite")

    def _predict(self, symbols, coordinates):
        try:
            from anysolv import predict_solvation_energy
        except ImportError as exc:
            raise ImportError("AnySolv requires pip install 'pdb2reaction[anysolv]'.") from exc
        import torch
        from ase.data import atomic_numbers

        dtype = torch.float64 if self.dtype in ("float64", "fp64") else torch.float32
        device = self.device
        if device in ("auto", "None"):
            device = "cuda" if torch.cuda.is_available() else "cpu"
        numbers = np.array([atomic_numbers[str(s).capitalize()] for s in symbols])
        energy, forces = predict_solvation_energy(
            (numbers, np.asarray(coordinates)), charge=self.charge, spin=self.spin,
            solvent=self.solvent, checkpoint=self.checkpoint, device=device,
            dtype=dtype, inference_settings="default",
        )
        forces = np.asarray(forces, dtype=np.float64)
        if forces.shape != (len(symbols), 3) or not np.isfinite(forces).all() or not np.isfinite(energy):
            raise ValueError("AnySolv returned non-finite energy/forces or an invalid force shape")
        return float(energy), forces

    def compute(self, symbols, coordinates, *, need_forces=False, need_hessian=False):
        coordinates = np.asarray(coordinates, dtype=np.float64).reshape(len(symbols), 3)
        if not np.isfinite(coordinates).all():
            raise ValueError("AnySolv coordinates must be finite")
        energy, forces = self._predict(symbols, coordinates)
        hessian = None
        if need_hessian:
            n = coordinates.size
            hessian = np.empty((n, n), dtype=np.float64)
            step = self.hessian_step
            for column in range(n):
                plus, minus = coordinates.copy(), coordinates.copy()
                plus.flat[column] += step
                minus.flat[column] -= step
                _, fp = self._predict(symbols, plus)
                _, fm = self._predict(symbols, minus)
                hessian[:, column] = -(fp - fm).reshape(-1) / (2 * step)
            hessian = 0.5 * (hessian + hessian.T)
        return energy, forces if need_forces or need_hessian else None, hessian


def wrap_ase_calculator(base, **settings):
    """Use the same correction on the ASE/DMF energy surface."""
    from ase.calculators.calculator import Calculator, all_changes
    from ase.calculators.mixing import SumCalculator

    correction = AnySolvCorrection(**settings)

    class DeltaCalculator(Calculator):
        implemented_properties = ["energy", "forces"]

        def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
            super().calculate(atoms, properties, system_changes)
            e, f, _ = correction.compute(
                self.atoms.get_chemical_symbols(), self.atoms.get_positions(), need_forces=True,
            )
            self.results = {"energy": e, "forces": f}

    class SolvatedCalculator(SumCalculator):
        def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
            source = atoms if atoms is not None else self.atoms
            prepared = source.copy()
            prepared.info.update(charge=correction.charge, spin=correction.spin)
            super().calculate(prepared, properties, system_changes)

    return SolvatedCalculator([base, DeltaCalculator()])
