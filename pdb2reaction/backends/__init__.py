"""
MLIP backend factory and registry.

Usage::

    from pdb2reaction.backends import create_calculator, create_ase_calculator

    calc = create_calculator(backend="uma", charge=0, spin=1, ...)
    ase_calc = create_ase_calculator(backend="uma", model="uma-s-1p1", ...)
"""

from __future__ import annotations

import importlib
import warnings
from typing import Any, Dict, Optional

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
        "precision",  # fp32 (default) | fp64 (full-precision base inference)
    },
    "orb": {
        "charge", "spin", "device", "freeze_atoms", "hessian_calc_mode",
        "return_partial_hessian", "hessian_double", "out_hess_torch",
        "print_timing",
        "model", "precision", "compile_model",
    },
    "mace": {
        "charge", "spin", "device", "freeze_atoms", "hessian_calc_mode",
        "return_partial_hessian", "hessian_double", "out_hess_torch",
        "print_timing",
        "model", "default_dtype",
    },
    "aimnet2": {
        "charge", "spin", "device", "freeze_atoms", "hessian_calc_mode",
        "return_partial_hessian", "hessian_double", "out_hess_torch",
        "print_timing",
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


# Backend-specific value domain for the unified ``--precision`` CLI flag.
# Keyed by canonical user value ('fp32' | 'fp64') -> { backend: (kw_name, kw_value) }.
# UMA accepts 'fp32'/'fp64' directly (its own enum). ORB expects orb_models'
# matmul-precision strings ('float32-high' | 'float64' — 'float32' alone is
# silently demoted to 'highest', see backends/orb.py:55-62). MACE wants
# numpy-style dtype strings.
_PRECISION_DISPATCH: Dict[str, Dict[str, tuple]] = {
    "fp32": {
        "uma":  ("precision", "fp32"),
        "orb":  ("precision", "float32-high"),
        "mace": ("default_dtype", "float32"),
    },
    "fp64": {
        "uma":  ("precision", "fp64"),
        "orb":  ("precision", "float64"),
        "mace": ("default_dtype", "float64"),
    },
}


def apply_precision_to_calc_cfg(calc_cfg: Dict[str, Any], precision: str) -> None:
    """Route the unified ``--precision`` CLI value into backend-specific kwargs.

    Mutates ``calc_cfg`` in place. For aimnet2 (no precision knob) ``fp32`` is a
    no-op and ``fp64`` is rejected (its model inputs are cast to float32
    upstream). Raises ``BackendError`` for invalid values or aimnet2 + fp64.
    """
    val = str(precision or "").lower()
    if val not in _PRECISION_DISPATCH:
        raise BackendError(
            f"--precision must be 'fp32' or 'fp64', got {precision!r}"
        )
    backend = resolve_backend(calc_cfg.get("backend") or "uma")
    mapping = _PRECISION_DISPATCH[val]
    if backend not in mapping:
        # AIMNet2 (or any future backend with no precision knob) cannot honour
        # `--precision fp64`: upstream `aimnet` casts model inputs to float32
        # (`aimnet/calculators/calculator.py keys_in: torch.float`), so silently
        # accepting fp64 would lie about the actual numeric precision of the
        # run. Reject the combination loudly; fp32 stays a no-op so users can
        # swap `--backend` without changing scripts.
        if val == "fp64":
            raise BackendError(
                f"--precision fp64 is not supported by backend {backend!r}: "
                f"its model inputs are cast to float32 upstream, so the run "
                f"would not actually be fp64."
            )
        return
    kw_name, kw_val = mapping[backend]
    calc_cfg[kw_name] = kw_val
    if val == "fp64":
        # fp64 model precision implies a fp64 Hessian. Leaving ``hessian_double``
        # off while the forward pass is fp64 is internally inconsistent: the
        # optimizer / eigen linear algebra would discard exactly the precision
        # the user paid for. ``hessian_double`` is not a CLI flag (always True
        # by default), so the only way to reach the mismatch is a hand-edited
        # config. Warn that the deliberate inconsistency is being overridden,
        # then force fp64 on so the run stays self-consistent.
        if calc_cfg.get("hessian_double", True) is False:
            warnings.warn(
                "--precision fp64 forces a fp64 Hessian: overriding the "
                "hessian_double=False set in the config. A fp32 Hessian under "
                "fp64 model precision is internally inconsistent (the optimizer "
                "and eigen linear algebra would discard the fp64 precision). "
                "Use --precision fp32 if a fp32 Hessian is intended.",
                stacklevel=2,
            )
        calc_cfg["hessian_double"] = True


def apply_backend_model_to_calc_cfg(calc_cfg: Dict[str, Any], backend_model: Optional[str] = None) -> None:
    """Route the unified ``--backend-model`` CLI value into the ``model`` kwarg
    (the active backend's model variant). Mutates ``calc_cfg`` in place; a no-op
    when unset, so the backend keeps its built-in default model.

    Also consumes (pops) any raw ``backend_model`` token already in ``calc_cfg``
    (e.g. propagated through a ``--config`` YAML), since it is not a Calculator
    kwarg; that token is used when no explicit argument is passed.
    """
    raw = calc_cfg.pop("backend_model", None)
    val = backend_model if backend_model is not None else raw
    if val is None or str(val).strip() == "":
        return  # keep the backend default model
    calc_cfg["model"] = str(val).strip()


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
    # Env-var / direct-API entry point for strict determinism (the CLI uses the
    # --deterministic flag callback). Idempotent; no-op unless requested.
    from ._determinism import (
        is_deterministic_active,
        is_deterministic_requested,
        setup_deterministic,
    )
    if is_deterministic_requested():
        setup_deterministic()
    if backend == "aimnet2" and is_deterministic_active():
        raise BackendError(
            "AIMNet2 forces are not bit-reproducible under strict determinism: "
            "they are produced by a custom CUDA kernel (torch.ops.aimnet, "
            "conv_sv_2d_sp_bwd) outside torch.use_deterministic_algorithms "
            "control (~1e-9 residual across runs; the energy is reproducible)."
        )
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
    from ._determinism import (
        is_deterministic_active,
        is_deterministic_requested,
        setup_deterministic,
    )
    if is_deterministic_requested():
        setup_deterministic()
    if backend == "aimnet2" and is_deterministic_active():
        raise BackendError(
            "AIMNet2 forces are not bit-reproducible under strict determinism: "
            "they are produced by a custom CUDA kernel (torch.ops.aimnet, "
            "conv_sv_2d_sp_bwd) outside torch.use_deterministic_algorithms "
            "control (~1e-9 residual across runs; the energy is reproducible)."
        )
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
