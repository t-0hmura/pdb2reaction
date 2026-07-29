# pdb2reaction/backends/custom.py

"""Custom user-supplied calculator backend.

Load an arbitrary `ASE <https://wiki.fysik.dtu.dk/ase/>`_ Calculator from a user
Python file and use it as the energy/gradient source for any ``pdb2reaction``
subcommand. This lets a user couple the pipeline with GFN-xTB (via tblite /
``xtb-python``), DFTB+, ORCA, or any ASE-compatible engine without modifying
``pdb2reaction`` — the boundary is the standard ASE Calculator interface.

The user file must expose a factory callable (default name ``get_calculator``)
that returns an :class:`ase.calculators.calculator.Calculator`::

    # my_calc.py
    from ase.calculators.emt import EMT

    def get_calculator(charge=0, spin=1, device="auto", **kwargs):
        return EMT()

Then::

    pdb2reaction sp -i model.xyz --calc-file my_calc.py

The factory is passed ``charge``, ``spin`` (multiplicity; also offered as
``mult`` / ``multiplicity``), and ``device`` when its signature accepts them
(or unconditionally if it declares ``**kwargs``), so backends that need the
total charge (e.g. xTB) can be configured. A module-level Calculator instance
(or a Calculator object bound to the factory name) is also accepted directly.

The pysisyphus-facing :class:`CustomCalculator` only implements
``_compute_energy_forces_ev``; freeze-atom masking, finite-difference Hessian
assembly, and unit conversion are inherited from :class:`MLIPCalculator`.
"""

from __future__ import annotations

import importlib.util
import inspect
import sys
from pathlib import Path
from typing import Any, Dict, Optional, Sequence

import numpy as np

from .base import BackendError, MLIPCalculator

_MODULE_COUNTER = 0


def _is_ase_calculator(obj: Any) -> bool:
    """Duck-typed check for an ASE-compatible Calculator."""
    if inspect.isclass(obj):
        return False
    try:
        from ase.calculators.calculator import Calculator as _ASECalc
    except ImportError:
        # ASE not importable in this environment; fall back to the duck type
        # (energy + forces accessors) so a non-ASE-subclass calculator still works.
        _ASECalc = None
    if _ASECalc is not None and isinstance(obj, _ASECalc):
        return True
    return hasattr(obj, "get_potential_energy") and hasattr(obj, "get_forces")


def _load_user_module(calc_file: str):
    """Import the user's Python file as an isolated module by path."""
    path = Path(calc_file).expanduser()
    if not path.is_file():
        raise BackendError(f"--calc-file not found: {path}")
    global _MODULE_COUNTER
    _MODULE_COUNTER += 1
    mod_name = f"pdb2reaction_calc_file_{_MODULE_COUNTER}"
    spec = importlib.util.spec_from_file_location(mod_name, str(path))
    if spec is None or spec.loader is None:
        raise BackendError(f"Could not load --calc-file as a Python module: {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[mod_name] = module
    try:
        spec.loader.exec_module(module)
    except Exception as exc:  # surface the user's import-time error clearly
        sys.modules.pop(mod_name, None)
        raise BackendError(
            f"Error while importing --calc-file {path}: {exc}"
        ) from exc
    return module, path


def _call_factory(factory, *, charge: int, spin: int, device: str, extra: Dict[str, Any]):
    """Call the user factory, passing only the kwargs its signature accepts."""
    candidate: Dict[str, Any] = {
        "charge": charge,
        "spin": spin,
        "mult": spin,
        "multiplicity": spin,
        "device": device,
    }
    candidate.update(extra or {})
    try:
        sig = inspect.signature(factory)
    except (TypeError, ValueError):
        # Builtins / C callables expose no signature; pass the full candidate set.
        return factory(**candidate)
    params = sig.parameters.values()
    if any(p.kind == p.VAR_KEYWORD for p in params):
        return factory(**candidate)
    names = {p.name for p in params}
    return factory(**{k: v for k, v in candidate.items() if k in names})


def load_ase_calculator(
    calc_file: str,
    calc_factory: str = "get_calculator",
    *,
    charge: int = 0,
    spin: int = 1,
    device: str = "auto",
    **extra: Any,
):
    """Load and return an ASE Calculator from a user Python file.

    Parameters
    ----------
    calc_file
        Path to a Python file exposing ``calc_factory``.
    calc_factory
        Name of a callable returning an ASE Calculator, or of a module-level
        Calculator instance. Defaults to ``"get_calculator"``.
    charge, spin, device
        Forwarded to the factory when its signature accepts them.
    """
    module, path = _load_user_module(calc_file)
    if not hasattr(module, calc_factory):
        exported = [n for n in vars(module) if not n.startswith("_")]
        raise BackendError(
            f"--calc-file {path} has no attribute '{calc_factory}'. "
            f"Define `def {calc_factory}(charge=0, spin=1, **kwargs)` returning "
            f"an ASE Calculator (or use --calc-factory NAME). "
            f"Found top-level names: {', '.join(exported) or '(none)'}."
        )
    attr = getattr(module, calc_factory)

    # A Calculator instance bound directly to the factory name is accepted as-is.
    if _is_ase_calculator(attr):
        return attr
    if not callable(attr):
        raise BackendError(
            f"--calc-file {path}: '{calc_factory}' is neither a callable factory "
            f"nor an ASE Calculator instance (got {type(attr).__name__})."
        )
    try:
        calc = _call_factory(attr, charge=charge, spin=spin, device=device, extra=extra)
    except Exception as exc:
        raise BackendError(
            f"--calc-file {path}: calling {calc_factory}(...) failed: {exc}"
        ) from exc
    if not _is_ase_calculator(calc):
        raise BackendError(
            f"--calc-file {path}: {calc_factory}(...) returned {type(calc).__name__}, "
            f"not an ASE Calculator (needs get_potential_energy / get_forces)."
        )
    return calc


def make_custom_ase_calculator(
    *,
    calc_file: str,
    calc_factory: str = "get_calculator",
    charge: int = 0,
    spin: int = 1,
    device: str = "auto",
    **kwargs: Any,
):
    """ASE-calculator factory entry point (registry ``ase_cls``).

    Returns the user's ASE Calculator directly (used by the DMF ``path-opt`` /
    ``path-search`` ASE path).
    """
    return load_ase_calculator(
        calc_file, calc_factory, charge=charge, spin=spin, device=device, **kwargs
    )


class CustomCalculator(MLIPCalculator):
    """PySisyphus-compatible wrapper around a user-supplied ASE Calculator.

    Energy and forces are obtained from the ASE Calculator (eV and eV/Å), which
    is exactly the unit contract expected by :class:`MLIPCalculator`. Hessians
    fall back to the inherited finite-difference path.
    """

    def __init__(
        self,
        *,
        calc_file: str,
        calc_factory: str = "get_calculator",
        charge: int = 0,
        spin: int = 1,
        device: str = "auto",
        freeze_atoms: Optional[Sequence[int]] = None,
        hessian_calc_mode: str = "FiniteDifference",
        return_partial_hessian: bool = True,
        hessian_double: bool = True,
        out_hess_torch: bool = True,
        print_timing: bool = True,
        **kwargs: Any,
    ):
        super().__init__(
            charge=charge,
            spin=spin,
            device=device,
            freeze_atoms=freeze_atoms,
            hessian_calc_mode=hessian_calc_mode,
            return_partial_hessian=return_partial_hessian,
            hessian_double=hessian_double,
            out_hess_torch=out_hess_torch,
            print_timing=print_timing,
        )
        self.calc_file = str(calc_file)
        self.calc_factory = str(calc_factory)
        self._ase_calc = load_ase_calculator(
            calc_file, calc_factory, charge=charge, spin=spin, device=device
        )

    def _compute_energy_forces_ev(self, elem: Sequence[str], coord_ang: np.ndarray):
        from ase import Atoms

        atoms = Atoms(
            symbols=list(elem),
            positions=np.asarray(coord_ang, dtype=float).reshape(-1, 3),
        )
        # Expose charge/spin per frame the same way the ORB/MACE backends do
        # (orb.py:211-212): a user Calculator that reads them (any OMol-style
        # wrapper) would otherwise run neutral/singlet. The construction-time
        # kwargs reach the factory only once, so a factory that ignores its
        # charge= argument has no other channel.
        atoms.info["charge"] = int(self.charge)
        atoms.info["spin"] = int(self.mult)
        atoms.calc = self._ase_calc
        energy_ev = float(atoms.get_potential_energy())
        forces_ev_ang = np.asarray(atoms.get_forces(), dtype=float).reshape(-1, 3)
        return energy_ev, forces_ev_ang
