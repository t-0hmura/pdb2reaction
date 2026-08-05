"""
MLIP backend factory and registry.

Usage::

    from pdb2reaction.backends import create_calculator, create_ase_calculator
import click

    calc = create_calculator(backend="uma", charge=0, spin=1, ...)
    ase_calc = create_ase_calculator(backend="uma", model="uma-s-1p2", ...)
"""

from __future__ import annotations

import importlib
import warnings
from typing import Any, Dict, Optional

from .base import BackendError, MLIPCalculator, normalize_hessian_calc_mode
import click

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
    # User-supplied ASE Calculator loaded from a Python file (``--calc-file``).
    # Not auto-resolvable (needs a file); selected by passing ``--calc-file``.
    "custom": {
        "module": "pdb2reaction.backends.custom",
        "pysis_cls": "CustomCalculator",
        "ase_cls": "make_custom_ase_calculator",
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
    "custom": {
        "charge", "spin", "device", "freeze_atoms", "hessian_calc_mode",
        "return_partial_hessian", "hessian_double", "out_hess_torch",
        "print_timing",
        "calc_file", "calc_factory",
    },
}

# Keys accepted by ASE calculator factories
_ASE_ACCEPTED_KEYS: Dict[str, set] = {
    # `precision` is accepted so the ASE/DMF path can match the pysisyphus eval
    # PES (e.g. under --precision fp64); UMAASECalculator honours it. Without it
    # the DMF path optimizer was pinned to fp32 while the HEI energy ranker ran
    # fp64 → frames optimized on one PES, ranked on another.
    "uma": {"model", "device", "task_name", "workers", "workers_per_node", "precision"},
    "orb": {"model", "device", "precision", "compile_model"},
    "mace": {"model", "device", "default_dtype"},
    "aimnet2": {"model", "device", "charge", "spin"},
    "custom": {"calc_file", "calc_factory", "charge", "spin", "device"},
}

VALID_BACKENDS = tuple(BACKEND_REGISTRY.keys())

_BACKEND_AVAILABILITY_MODULES = {
    "uma": "fairchem.core",
    "orb": "orb_models.forcefield.pretrained",
    "mace": "mace.calculators",
    "aimnet2": "aimnet.calculators",
}


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

# Precision used when the user names none (`--precision` absent and no
# ``calc.precision`` in the config, i.e. the ``"auto"`` token of
# ``CALC_KW_DEFAULT``). Defaults follow each backend's established numerical
# mode: UMA uses its upstream fp32 baseline, while ORB and MACE use fp64.
_BACKEND_DEFAULT_PRECISION: Dict[str, str] = {
    "uma": "fp32",
    "orb": "fp64",
    "mace": "fp64",
    "aimnet2": "fp32",  # no precision knob; fp32 is a no-op
    "custom": "fp32",   # the user's own Calculator owns its dtype
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
        calc_cfg["precision"] = val
        return
    kw_name, kw_val = mapping[backend]
    calc_cfg[kw_name] = kw_val
    if kw_name != "precision":
        # MACE carries its choice in ``default_dtype``; keep ``precision`` as the
        # canonical token so the run-summary line and the ``all`` pipeline's
        # child-config propagation never surface the raw ``"auto"`` sentinel.
        # Filtered out of the MACE kwargs by _BACKEND_ACCEPTED_KEYS.
        calc_cfg["precision"] = val
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


def apply_effective_precision(
    calc_cfg: Dict[str, Any], cli_precision: Optional[str]
) -> None:
    """Apply CLI, config, or backend-default precision in that order.

    The unified ``calc.precision`` value is dispatched for pipeline children as
    well as standalone commands. An absent or ``auto`` value resolves through
    ``_BACKEND_DEFAULT_PRECISION``; already backend-specific values are left
    untouched.
    """
    if resolve_backend(calc_cfg.get("backend") or "uma") == "custom":
        calc_cfg.pop("precision", None)
        return
    eff = cli_precision if cli_precision is not None else calc_cfg.get("precision")
    if eff is None or str(eff).lower() == "auto":
        backend = resolve_backend(calc_cfg.get("backend") or "uma")
        eff = _BACKEND_DEFAULT_PRECISION.get(backend, "fp32")
    if str(eff).lower() in ("fp32", "fp64"):
        apply_precision_to_calc_cfg(calc_cfg, str(eff).lower())


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


def apply_calc_file_to_calc_cfg(
    calc_cfg: Dict[str, Any],
    calc_file: Optional[str] = None,
    calc_factory: Optional[str] = None,
) -> None:
    """Route the ``--calc-file`` CLI value into the ``custom`` backend.

    When a calc-file is given, switch the backend to ``custom`` (loading an
    arbitrary ASE Calculator from the user Python file) and store the file path
    and factory name. Mutates ``calc_cfg`` in place; a no-op when no calc-file
    is set, so the normal ``--backend`` selection is untouched.

    Also consumes (pops) any raw ``calc_file`` / ``calc_factory`` already in
    ``calc_cfg`` (e.g. propagated through a ``--config`` YAML); an explicit CLI
    argument takes precedence over the YAML token.
    """
    raw_file = calc_cfg.pop("calc_file", None)
    raw_factory = calc_cfg.pop("calc_factory", None)
    chosen = calc_file if calc_file is not None else raw_file
    factory = calc_factory if calc_factory is not None else raw_factory
    has_factory = factory is not None and str(factory).strip() != ""
    if chosen is None or str(chosen).strip() == "":
        # An effective factory name without an effective calculator file would
        # be dropped here and a built-in backend would run instead.
        if has_factory:
            raise click.BadParameter(
                f"calc_factory {str(factory).strip()!r} names a factory inside a "
                "custom calculator file, but no calculator file is effective. "
                "Supply --calc-file or calc.calc_file, or drop calc_factory."
            )
        return  # no calc-file: keep the --backend selection
    calc_cfg["backend"] = "custom"
    calc_cfg["calc_file"] = str(chosen)
    if has_factory:
        calc_cfg["calc_factory"] = str(factory).strip()
    # A user ASE Calculator has no MLIP model variant: drop the inherited UMA
    # default so run headers don't mislabel it as 'uma-s-1p2'.
    calc_cfg.pop("model", None)


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
            importlib.import_module(_BACKEND_AVAILABILITY_MODULES[key])
            _import_cls(key, "pysis_cls")
        except Exception:
            continue
        # Say which one won. Provenance records the resolved name, so without
        # this the user cannot tell an explicit choice from an automatic one.
        click.echo(
            f"[backend] NOTE: backend='auto' resolved to {key!r} "
            "(first importable of uma, orb, mace, aimnet2).",
            err=True,
        )
        return key
    raise BackendError(
        "No MLIP backend available. Install one of: "
        "fairchem-core (UMA), orb-models (ORB), mace-torch (MACE), aimnet (AIMNet2)."
    )


def create_calculator(backend: str = "uma", **kwargs) -> MLIPCalculator:
    """Create a PySisyphus-compatible MLIP calculator.

    Parameters
    ----------
    backend : str
        One of ``'uma'``, ``'orb'``, ``'mace'``, ``'aimnet2'``, ``'custom'``,
        or ``'auto'``. ``'custom'`` requires ``calc_file``.
    **kwargs
        Backend-specific and common parameters. Unknown keys for the selected
        backend are ignored with a warning.

        Solvent correction keys (``solvent``, ``solvent_model``, ``xtb_cmd``,
        ``xtb_acc``) are extracted and used to wrap the base calculator with
        ``SolventCorrectedCalculator`` when ``solvent != 'none'``.
    """
    # Validate the cross-backend numerical method before importing or loading a
    # backend. Direct API callers receive the same strict enum as CLI callers.
    if "hessian_calc_mode" in kwargs:
        kwargs["hessian_calc_mode"] = normalize_hessian_calc_mode(
            kwargs["hessian_calc_mode"]
        )

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
    # A `calc:` key this backend does not accept is dropped here. Silence reads
    # as "applied", which is how a documented-but-ignored setting survives.
    _dropped = sorted(k for k in kwargs if k not in filtered)
    if _dropped:
        click.echo(
            f"[backend] WARNING: {backend} ignored calc setting(s) "
            f"{', '.join(_dropped)}; they are not accepted by this backend.",
            err=True,
        )
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
        One of ``'uma'``, ``'orb'``, ``'mace'``, ``'aimnet2'``, ``'custom'``,
        or ``'auto'``. ``'custom'`` requires ``calc_file``.
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
    "normalize_hessian_calc_mode",
    "create_ase_calculator",
    "create_calculator",
    "resolve_backend",
    "apply_precision_to_calc_cfg",
    "apply_backend_model_to_calc_cfg",
    "apply_calc_file_to_calc_cfg",
]
