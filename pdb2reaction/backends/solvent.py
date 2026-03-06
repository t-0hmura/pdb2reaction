# pdb2reaction/backends/solvent.py

"""
xTB implicit-solvent delta-correction wrapper.

Wraps any ``MLIPCalculator`` and adds an xTB-based implicit-solvent
correction::

    ΔE = E_xTB(solvent) − E_xTB(vacuum)
    ΔF = F_xTB(solvent) − F_xTB(vacuum)
    ΔH = H_xTB(solvent) − H_xTB(vacuum)

The correction is added to the MLIP energy/forces/Hessian.

Usage::

    from pdb2reaction.backends.solvent import SolventCorrectedCalculator

    base = create_calculator(backend="uma", charge=0, spin=1)
    calc = SolventCorrectedCalculator(base, solvent="water", solvent_model="alpb")
"""

from __future__ import annotations

from typing import Any, Dict, Optional, Sequence

import numpy as np

from pysisyphus.constants import BOHR2ANG, AU2EV

from .base import MLIPCalculator, BackendError

# Unit conversion (atomic units ↔ eV/Å)
_EV2AU = 1.0 / AU2EV
_ANG2BOHR = 1.0 / BOHR2ANG
_F_EVAA_2_AU = _EV2AU / _ANG2BOHR
_H_EVAA_2_AU = _EV2AU / _ANG2BOHR / _ANG2BOHR

# xTB correction helpers (lazy-imported from bundled module)
_xtb_mod = None


def _get_xtb_mod():
    global _xtb_mod
    if _xtb_mod is None:
        from . import xtb_alpb_correction as _mod
        _xtb_mod = _mod
    return _xtb_mod


def solvent_correction_enabled(solvent: str) -> bool:
    """Return True if solvent is not 'none'/empty."""
    mod = _get_xtb_mod()
    return mod.solvent_correction_enabled(solvent)


class SolventCorrectedCalculator(MLIPCalculator):
    """Decorator that adds xTB implicit-solvent correction to any MLIP calculator.

    If ``solvent`` is ``'none'`` (default), all calls are passed through
    to the base calculator unchanged.
    """

    def __init__(
        self,
        base: MLIPCalculator,
        *,
        solvent: str = "none",
        solvent_model: str = "alpb",
        xtb_cmd: str = "xtb",
        xtb_acc: float = 0.2,
    ):
        # Inherit charge/spin/device from base
        super().__init__(
            charge=base.charge,
            spin=base.mult,
            device=getattr(base, "device_str", "auto"),
            freeze_atoms=base.freeze_atoms,
            hessian_calc_mode=base.hessian_calc_mode,
            return_partial_hessian=base.return_partial_hessian,
            hessian_double=base.hessian_double,
            print_timing=base.print_timing,
        )

        self.base = base
        self.solvent = str(solvent)
        self.solvent_model = str(solvent_model)
        self.xtb_cmd = str(xtb_cmd)
        self.xtb_acc = float(xtb_acc)

        self._enabled = solvent_correction_enabled(self.solvent)

        if self._enabled:
            mod = _get_xtb_mod()
            # Validate solvent model
            mod.normalize_solvent_model(self.solvent_model)

    def _solvent_delta(
        self,
        symbols: Sequence[str],
        coords_ang: np.ndarray,
        need_forces: bool = False,
        need_hessian: bool = False,
    ):
        """Compute xTB solvent delta (eV/Å units)."""
        mod = _get_xtb_mod()
        de_ev, df_ev_ang, dh_ev_ang2 = mod.delta_alpb_minus_vac(
            symbols=list(symbols),
            coords_ang=np.asarray(coords_ang, dtype=np.float64),
            charge=self.base.charge,
            multiplicity=self.base.mult,
            solvent=self.solvent,
            solvent_model=self.solvent_model,
            need_forces=need_forces,
            need_hessian=need_hessian,
            xtb_cmd=self.xtb_cmd,
            xtb_acc=self.xtb_acc,
            xtb_workdir="tmp",
            xtb_keep_files=False,
            ncores=mod.resolve_xtb_ncores(),
        )
        return float(de_ev), df_ev_ang, dh_ev_ang2

    # ------------------------------------------------------------------
    # PySisyphus API — delegate + add solvent correction
    # ------------------------------------------------------------------

    def get_energy(self, elem, coords):
        result = self.base.get_energy(elem, coords)
        if not self._enabled:
            return result
        coord_ang = np.asarray(coords, dtype=np.float64).reshape(-1, 3) * BOHR2ANG
        de_ev, _, _ = self._solvent_delta(elem, coord_ang, need_forces=False)
        result["energy"] += de_ev * _EV2AU
        return result

    def get_forces(self, elem, coords):
        result = self.base.get_forces(elem, coords)
        if not self._enabled:
            return result
        coord_ang = np.asarray(coords, dtype=np.float64).reshape(-1, 3) * BOHR2ANG
        de_ev, df_ev_ang, _ = self._solvent_delta(elem, coord_ang, need_forces=True)
        result["energy"] += de_ev * _EV2AU
        if df_ev_ang is not None:
            df_au = np.asarray(df_ev_ang, dtype=np.float64).reshape(-1) * _F_EVAA_2_AU
            result["forces"] = np.asarray(result["forces"], dtype=np.float64) + df_au
        return result

    def get_hessian(self, elem, coords):
        result = self.base.get_hessian(elem, coords)
        if not self._enabled:
            return result
        coord_ang = np.asarray(coords, dtype=np.float64).reshape(-1, 3) * BOHR2ANG
        de_ev, df_ev_ang, dh_ev_ang2 = self._solvent_delta(
            elem, coord_ang, need_forces=True, need_hessian=True
        )
        result["energy"] += de_ev * _EV2AU
        if df_ev_ang is not None:
            df_au = np.asarray(df_ev_ang, dtype=np.float64).reshape(-1) * _F_EVAA_2_AU
            result["forces"] = np.asarray(result["forces"], dtype=np.float64) + df_au
        if dh_ev_ang2 is not None:
            dh_au = np.asarray(dh_ev_ang2, dtype=np.float64) * _H_EVAA_2_AU
            n_atoms = len(elem)
            dh_au = self._apply_active_trim_np(dh_au, n_atoms)
            result["hessian"] = np.asarray(result["hessian"], dtype=np.float64) + dh_au
        return result

    # Pass through subclass hooks (unused since we override get_* directly)
    def _compute_energy_forces_ev(self, elem, coord_ang):
        return self.base._compute_energy_forces_ev(elem, coord_ang)
