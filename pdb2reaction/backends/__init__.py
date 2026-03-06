# pdb2reaction/backends/__init__.py

"""
MLIP backend factory and registry.

Usage::

    from pdb2reaction.backends import create_calculator, create_ase_calculator

    calc = create_calculator(backend="uma", charge=0, spin=1, ...)
    ase_calc = create_ase_calculator(backend="uma", model="uma-s-1p2", ...)
"""

from __future__ import annotations

import importlib
import warnings
from typing import Any, Dict

from .base import BackendError, MLIPCalculator

# Lazy-import registry.  Modules are imported only when the backend is used.
BACKEND_REGISTRY: Dict[str, Dict[str, str]] = {
    "uma": {
        "module": "pdb2reaction.backends.uma",
        "pysis_cls": "UMACalculator",
        "ase_cls": "UMAASECalculator",
    },
    "orb": {
        "module": "pdb2reaction.backends.orb",
        "pysis_cls": "OrbCalculator",
        "ase_cls": "OrbASECalculator",
    },
    "mace": {
        "module": "pdb2reaction.backends.mace",
        "pysis_cls": "MACECalculator",
        "ase_cls": "MACEASECalculator",
    },
    "aimnet2": {
        "module": "pdb2reaction.backends.aimnet2",
        "pysis_cls": "AIMNet2Calculator",
        "ase_cls": "AIMNet2ASECalculator",
    },
}

# Keys accepted by each backend (used to filter calc_cfg kwargs)
_BACKEND_ACCEPTED_KEYS: Dict[str, set] = {
    "uma": {
        "charge", "spin", "device", "freeze_atoms", "hessian_calc_mode",
        "return_partial_hessian", "hessian_double", "print_timing",
        "model", "task_name", "workers", "workers_per_node",
        "max_neigh", "radius", "r_edges", "out_hess_torch", "print_vram",
    },
    "orb": {
        "charge", "spin", "device", "freeze_atoms", "hessian_calc_mode",
        "return_partial_hessian", "hessian_double", "print_timing",
        "model", "precision", "compile_model",
    },
    "mace": {
        "charge", "spin", "device", "freeze_atoms", "hessian_calc_mode",
        "return_partial_hessian", "hessian_double", "print_timing",
        "model", "default_dtype",
    },
    "aimnet2": {
        "charge", "spin", "device", "freeze_atoms", "hessian_calc_mode",
        "return_partial_hessian", "hessian_double", "print_timing",
        "model",
    },
}

# Keys accepted by ASE calculator factories
_ASE_ACCEPTED_KEYS: Dict[str, set] = {
    "uma": {"model", "device", "task_name", "workers", "workers_per_node"},
    "orb": {"model", "device", "precision", "compile_model"},
    "mace": {"model", "device", "default_dtype"},
    "aimnet2": {"model", "device"},
}

VALID_BACKENDS = tuple(BACKEND_REGISTRY.keys())


def _import_cls(backend: str, cls_key: str):
    """Import a class from the backend registry."""
    reg = BACKEND_REGISTRY.get(backend)
    if reg is None:
        raise BackendError(
            f"Unknown backend '{backend}'. Choose from: {', '.join(VALID_BACKENDS)}"
        )
    mod = importlib.import_module(reg["module"])
    return getattr(mod, reg[cls_key])


def _filter_kwargs(kwargs: Dict[str, Any], accepted: set) -> Dict[str, Any]:
    """Return only keys present in *accepted*."""
    return {k: v for k, v in kwargs.items() if k in accepted}


def resolve_backend(backend: str) -> str:
    """Resolve ``'auto'`` to the best available backend (UMA-first)."""
    if backend != "auto":
        return backend
    for key in ("uma", "orb", "mace", "aimnet2"):
        try:
            _import_cls(key, "pysis_cls")
            return key
        except Exception:
            continue
    raise BackendError(
        "No MLIP backend available. Install one of: "
        "fairchem-core (UMA), orb-models (ORB), mace-torch (MACE), aimnet (AIMNet2)."
    )


def create_calculator(backend: str = "uma", **kwargs) -> MLIPCalculator:
    """Create a PySisyphus-compatible MLIP calculator.

    Parameters
    ----------
    backend : str
        One of ``'uma'``, ``'orb'``, ``'mace'``, ``'aimnet2'``, ``'auto'``.
    **kwargs
        Backend-specific and common parameters.  Unknown keys for the
        selected backend are silently ignored.

        Solvent correction keys (``solvent``, ``solvent_model``, ``xtb_cmd``,
        ``xtb_acc``) are extracted and used to wrap the base calculator with
        ``SolventCorrectedCalculator`` when ``solvent != 'none'``.
    """
    # Extract solvent keys before backend filtering
    solvent = str(kwargs.pop("solvent", "none"))
    solvent_model = str(kwargs.pop("solvent_model", "alpb"))
    xtb_cmd = str(kwargs.pop("xtb_cmd", "xtb"))
    xtb_acc = float(kwargs.pop("xtb_acc", 0.2))

    backend = resolve_backend(backend)
    accepted = _BACKEND_ACCEPTED_KEYS.get(backend, set())
    filtered = _filter_kwargs(kwargs, accepted)
    cls = _import_cls(backend, "pysis_cls")
    calc = cls(**filtered)

    # Wrap with solvent correction if enabled
    from .solvent import solvent_correction_enabled
    if solvent_correction_enabled(solvent):
        from .solvent import SolventCorrectedCalculator
        calc = SolventCorrectedCalculator(
            calc,
            solvent=solvent,
            solvent_model=solvent_model,
            xtb_cmd=xtb_cmd,
            xtb_acc=xtb_acc,
        )

    return calc


def create_ase_calculator(backend: str = "uma", **kwargs):
    """Create an ASE-compatible MLIP calculator (used for DMF).

    Parameters
    ----------
    backend : str
        One of ``'uma'``, ``'orb'``, ``'mace'``, ``'aimnet2'``, ``'auto'``.
    **kwargs
        Backend-specific parameters.
    """
    backend = resolve_backend(backend)
    accepted = _ASE_ACCEPTED_KEYS.get(backend, set())
    filtered = _filter_kwargs(kwargs, accepted)
    cls = _import_cls(backend, "ase_cls")
    return cls(**filtered)


__all__ = [
    "BACKEND_REGISTRY",
    "VALID_BACKENDS",
    "BackendError",
    "MLIPCalculator",
    "create_ase_calculator",
    "create_calculator",
    "resolve_backend",
]
