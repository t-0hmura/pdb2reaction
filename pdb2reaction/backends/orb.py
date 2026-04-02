# pdb2reaction/backends/orb.py

"""
ORB (orb-models) backend for pdb2reaction.

Requires: ``pip install "pdb2reaction[orb]"`` (orb-models).
"""

from __future__ import annotations

from typing import Optional, Sequence

import numpy as np

from .base import MLIPCalculator, BackendError

# Deprecated model aliases
ORB_DEPRECATED_MODEL_ALIASES = {
    "orb-v1": "orb-v2",
    "orb-d3-v1": "orb-d3-v2",
    "orb-d3-sm-v1": "orb-d3-sm-v2",
    "orb-d3-xs-v1": "orb-d3-xs-v2",
    "orb-v1-mptraj-only": "orb-mptraj-only-v2",
    "orb-mptraj-only-v1": "orb-mptraj-only-v2",
}


def _is_conservative_orb_model(model_name: str) -> bool:
    norm = str(model_name).replace("_", "-").lower()
    return ("conservative" in norm) and ("direct" not in norm)


def _unique_ordered(items):
    seen = set()
    out = []
    for item in items:
        if item not in seen:
            seen.add(item)
            out.append(item)
    return out


class OrbCalculator(MLIPCalculator):
    """ORB backend via orb-models."""

    def __init__(
        self,
        *,
        model: str = "orb_v3_conservative_omol",
        precision: str = "float32",
        compile_model: bool = False,
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
            from orb_models.forcefield import pretrained as orb_pretrained
        except Exception as exc:
            raise BackendError(
                "ORB backend requires orb-models and torch. "
                "Install with: pip install \"pdb2reaction[orb]\""
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
        self._pretrained = orb_pretrained

        if str(device).lower() == "auto":
            device = "cuda" if torch.cuda.is_available() else "cpu"
        self.device_str = device

        self.model_name = str(model)
        self.precision = str(precision)
        self.compile_model = bool(compile_model)

        if not _is_conservative_orb_model(self.model_name):
            raise BackendError(
                f"Only conservative Orb models are supported. Requested '{self.model_name}'."
            )

        self._loader = self._resolve_loader(self.model_name)
        self._model_obj, self._adapter = self._load_model()
        self._ase_calc = self._build_ase_calculator()

    def _resolve_loader(self, model_name: str):
        norm_dash = str(model_name).replace("_", "-").lower()
        if norm_dash in ORB_DEPRECATED_MODEL_ALIASES:
            model_name = ORB_DEPRECATED_MODEL_ALIASES[norm_dash]

        if not _is_conservative_orb_model(model_name):
            raise BackendError(f"Only conservative Orb models are supported. Requested '{model_name}'.")

        if hasattr(self._pretrained, "ORB_PRETRAINED_MODELS"):
            model_map = getattr(self._pretrained, "ORB_PRETRAINED_MODELS")
            cands = [model_name, model_name.replace("_", "-"), model_name.replace("-", "_")]
            cands.extend([x.lower() for x in cands])
            for cand in cands:
                if cand in model_map:
                    return model_map[cand]
            lower_map = {str(k).lower(): v for k, v in model_map.items()}
            for cand in cands:
                if str(cand).lower() in lower_map:
                    return lower_map[str(cand).lower()]

        for cand in (model_name, model_name.replace("-", "_"), model_name.lower().replace("-", "_")):
            if hasattr(self._pretrained, cand):
                return getattr(self._pretrained, cand)

        raise BackendError(f"Unknown Orb model '{model_name}'.")

    def _load_model(self):
        bases = [
            {"device": self.device_str, "precision": self.precision, "compile": self.compile_model},
            {"device": self.device_str, "precision": self.precision},
            {"device": self.device_str},
            {},
        ]
        last_exc = None
        for kwargs in bases:
            try:
                out = self._loader(**kwargs)
                if isinstance(out, tuple) and len(out) >= 2:
                    return out[0], out[1]
                return out, None
            except Exception as exc:
                last_exc = exc
                continue
        raise BackendError(f"Failed to load Orb model '{self.model_name}': {last_exc}")

    def _build_ase_calculator(self):
        def _try_construct(cls):
            arg_combos = []
            if self._adapter is not None:
                arg_combos.append(((self._model_obj, self._adapter), {"device": self.device_str}))
                arg_combos.append(((self._model_obj, self._adapter), {}))
            arg_combos.append(((self._model_obj,), {"device": self.device_str}))
            arg_combos.append(((self._model_obj,), {}))
            for args, kw in arg_combos:
                try:
                    return cls(*args, **kw)
                except TypeError:
                    continue
            return None

        try:
            from orb_models.forcefield.inference.calculator import ORBCalculator
            calc = _try_construct(ORBCalculator)
            if calc is not None:
                return calc
        except ImportError:
            pass

        try:
            from orb_models.forcefield.calculator import ORBCalculator
            calc = _try_construct(ORBCalculator)
            if calc is not None:
                return calc
        except ImportError:
            pass

        raise BackendError("Failed to build ORBCalculator.")

    def _compute_energy_forces_ev(self, elem, coord_ang):
        from ase import Atoms

        atoms = Atoms(symbols=list(elem), positions=np.asarray(coord_ang, dtype=np.float64))
        atoms.info["charge"] = int(self.charge)
        atoms.info["spin"] = int(self.mult)
        atoms.calc = self._ase_calc

        energy = float(atoms.get_potential_energy())
        forces = np.asarray(atoms.get_forces(), dtype=np.float64)
        return energy, forces


class OrbASECalculator:
    """Factory that returns an ORB ASE calculator for DMF."""

    def __new__(
        cls,
        *,
        model: str = "orb_v3_conservative_omol",
        device: str = "auto",
        precision: str = "float32",
        compile_model: bool = False,
    ):
        # Build and return the ORB ASE calculator directly
        calc = OrbCalculator(
            model=model, device=device, precision=precision,
            compile_model=compile_model, charge=0, spin=1,
        )
        return calc._ase_calc
