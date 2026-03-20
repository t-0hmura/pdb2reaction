# pdb2reaction/backends/aimnet2.py

"""
AIMNet2 backend for pdb2reaction.

Requires: ``pip install pdb2reaction[aimnet2]`` (aimnet).
"""

from __future__ import annotations

from typing import Optional, Sequence

import numpy as np

from .base import MLIPCalculator, BackendError


def _unique_ordered(items):
    seen = set()
    out = []
    for item in items:
        if item not in seen:
            seen.add(item)
            out.append(item)
    return out


class AIMNet2Calculator(MLIPCalculator):
    """AIMNet2 backend via aimnet."""

    def __init__(
        self,
        *,
        model: str = "aimnet2",
        # Base class parameters
        charge: int = 0,
        spin: int = 1,
        device: str = "auto",
        freeze_atoms: Optional[Sequence[int]] = None,
        hessian_calc_mode: str = "FiniteDifference",
        return_partial_hessian: bool = False,
        hessian_double: bool = True,
        print_timing: bool = True,
        **kwargs,
    ):
        try:
            import torch
        except Exception as exc:
            raise BackendError(
                "AIMNet2 backend requires torch and aimnet. "
                "Install with: pip install 'pdb2reaction[aimnet2]'"
            ) from exc

        super().__init__(
            charge=charge,
            spin=spin,
            device=device,
            freeze_atoms=freeze_atoms,
            hessian_calc_mode=hessian_calc_mode,
            return_partial_hessian=return_partial_hessian,
            hessian_double=hessian_double,
            print_timing=print_timing,
            **kwargs,
        )

        self._torch = torch
        if str(device).lower() == "auto":
            device = "cuda" if torch.cuda.is_available() else "cpu"
        self.device_str = str(device)
        self.model_name = str(model)

        self._calculator = self._load_calculator(self.model_name)

    def _load_calculator(self, model_name: str):
        try:
            from aimnet.calculators import AIMNet2Calculator as _AIMNet2Calc
        except Exception as exc:
            raise BackendError(
                "AIMNet2 backend requires aimnet. "
                "Install with: pip install 'pdb2reaction[aimnet2]'"
            ) from exc

        kwargs_variants = _unique_ordered(
            tuple(sorted(kw.items())) for kw in [
                {"device": self.device_str},
                {},
            ]
        )

        attempts = []
        for kw_tuple in kwargs_variants:
            kw = dict(kw_tuple)
            attempts.append(((model_name,), kw))
            attempts.append(((), {"model": str(model_name), **kw}))
            if kw:
                attempts.append(((str(model_name),), {}))

        last_exc = None
        for args, kwargs in attempts:
            try:
                return _AIMNet2Calc(*args, **kwargs)
            except Exception as exc:
                last_exc = exc

        raise BackendError(
            f"Failed to initialize AIMNet2 model '{model_name}' via aimnet."
        ) from last_exc

    @staticmethod
    def _to_scalar(value):
        if type(value).__module__.startswith("torch"):
            value = value.detach().cpu().numpy()
        if hasattr(value, "item"):
            return float(value.item())
        return float(np.asarray(value).reshape(-1)[0])

    @staticmethod
    def _extract_array(value, force_2d: bool):
        if type(value).__module__.startswith("torch"):
            value = value.detach().cpu().numpy()
        arr = np.asarray(value, dtype=np.float64)
        if force_2d:
            if arr.ndim == 3 and arr.shape[0] == 1:
                arr = arr[0]
            return arr.reshape(-1, 3)
        return arr

    @staticmethod
    def _pick_first_available(mapping, names):
        for name in names:
            if name in mapping:
                return mapping[name]
        lower = {str(k).lower(): v for k, v in mapping.items()}
        for name in names:
            val = lower.get(str(name).lower())
            if val is not None:
                return val
        return None

    def _call(self, symbols, coords_ang, with_hessian: bool):
        from ase import Atoms

        atoms = Atoms(symbols=symbols, positions=np.asarray(coords_ang, dtype=np.float64))
        numbers = np.asarray(atoms.get_atomic_numbers(), dtype=np.int64)
        coord_np = np.asarray(coords_ang, dtype=np.float32).reshape(-1, 3)

        data = {
            "coord": coord_np,
            "numbers": numbers.reshape(-1),
            "charge": np.asarray([float(self.charge)], dtype=np.float32),
            "mult": np.asarray([float(self.mult)], dtype=np.float32),
        }

        out = self._calculator(data, forces=True, hessian=bool(with_hessian))

        if isinstance(out, tuple):
            out = list(out)

        if isinstance(out, (list, tuple)):
            if len(out) < 2:
                raise BackendError(f"Unexpected AIMNet2 output tuple length {len(out)}")
            energy = self._to_scalar(out[0])
            forces = self._extract_array(out[1], force_2d=True)
            hess = self._extract_array(out[2], force_2d=False) if with_hessian and len(out) > 2 else None
            return energy, forces, hess

        if not isinstance(out, dict):
            raise BackendError(f"Unexpected AIMNet2 output type: {type(out)}")

        energy = self._pick_first_available(out, ("energy", "E", "e", "total_energy", "free_energy"))
        if energy is None:
            raise BackendError(f"AIMNet2 output missing energy key. Keys: {sorted(out.keys())}")
        energy = self._to_scalar(energy)

        forces = self._pick_first_available(out, ("forces", "force"))
        if forces is None:
            raise BackendError(f"AIMNet2 output missing forces key. Keys: {sorted(out.keys())}")
        forces = self._extract_array(forces, force_2d=True)

        hess = self._pick_first_available(out, ("hessian", "Hessian", "hess", "hessians"))
        if hess is not None:
            hess = self._extract_array(hess, force_2d=False)
        if with_hessian and hess is None:
            raise BackendError(f"AIMNet2 output missing hessian key. Keys: {sorted(out.keys())}")

        return energy, forces, hess

    def _compute_energy_forces_ev(self, elem, coord_ang):
        energy, forces, _ = self._call(list(elem), coord_ang, with_hessian=False)
        return energy, np.asarray(forces, dtype=np.float64)


class AIMNet2ASECalculator:
    """Factory that returns an AIMNet2-backed ASE calculator for DMF.

    Since aimnet provides its own ASE interface, this creates a
    lightweight ASE calculator that delegates to the aimnet calculator.
    """

    def __new__(cls, *, model: str = "aimnet2", device: str = "auto"):
        try:
            import torch
            from aimnet.calculators import AIMNet2Calculator as _AIMNet2Calc
        except Exception as exc:
            raise BackendError(
                "AIMNet2 backend requires aimnet. "
                "Install with: pip install 'pdb2reaction[aimnet2]'"
            ) from exc

        if str(device).lower() == "auto":
            device = "cuda" if torch.cuda.is_available() else "cpu"

        # Try different constructor signatures
        for args, kwargs in [
            ((model,), {"device": device}),
            ((), {"model": model, "device": device}),
            ((model,), {}),
        ]:
            try:
                return _AIMNet2Calc(*args, **kwargs)
            except Exception:
                continue

        raise BackendError(f"Failed to initialize AIMNet2 ASE calculator for model '{model}'.")
