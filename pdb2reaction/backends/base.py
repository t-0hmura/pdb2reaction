# pdb2reaction/backends/base.py

"""
Base MLIP calculator class for PySisyphus integration.

All MLIP backends (UMA, ORB, MACE, AIMNet2) extend ``MLIPCalculator``,
which inherits from ``pysisyphus.calculators.Calculator``.  Common
functionality — freeze-atom handling, finite-difference Hessian assembly,
unit conversion, and timing/VRAM logging — lives here.
"""

from __future__ import annotations

import time
from typing import Any, Dict, List, Optional, Sequence

import click
import numpy as np

from pysisyphus.calculators.Calculator import Calculator
from pysisyphus.constants import BOHR2ANG, ANG2BOHR, AU2EV

# ---- Unit conversion constants ----
EV2AU = 1.0 / AU2EV                       # eV → Hartree
F_EVAA_2_AU = EV2AU / ANG2BOHR            # eV Å⁻¹ → Hartree Bohr⁻¹
H_EVAA_2_AU = EV2AU / ANG2BOHR / ANG2BOHR # eV Å⁻² → Hartree Bohr⁻²


class BackendError(RuntimeError):
    """Raised for backend-specific runtime failures."""


def _prepare_model_for_autograd_hessian(model_obj, torch_mod):
    """Prepare a torch module for ``torch.autograd.functional.hessian``.

    Saves and mutates model state so autograd can traverse the forward pass
    without tracking parameter gradients or running stochastic dropout:

    * every parameter's ``requires_grad`` is saved and set to ``False``
    * ``model_obj.train(True)`` is called (autograd backward through modules
      that short-circuit in eval-mode needs training mode active)
    * every Dropout / AlphaDropout / FeatureAlphaDropout module has its
      ``p`` saved and forced to ``0.0``, and is switched to eval mode
      ``train(False)``, so dropout is effectively a no-op during the
      Hessian pass

    Returns a ``state`` dict which must be passed to
    ``_restore_model_after_autograd_hessian`` in a ``finally`` block.
    """
    state = {
        "was_training": bool(getattr(model_obj, "training", False)),
        "param_flags": [],
        "dropout_states": [],
    }

    if hasattr(model_obj, "parameters"):
        for param in model_obj.parameters():
            state["param_flags"].append((param, bool(param.requires_grad)))
            param.requires_grad_(False)

    if hasattr(model_obj, "train"):
        model_obj.train(True)

    dropout_types = []
    nn_mod = getattr(torch_mod, "nn", None)
    if nn_mod is not None:
        for name in (
            "Dropout",
            "Dropout1d",
            "Dropout2d",
            "Dropout3d",
            "AlphaDropout",
            "FeatureAlphaDropout",
        ):
            cls = getattr(nn_mod, name, None)
            if cls is not None:
                dropout_types.append(cls)

    if dropout_types and hasattr(model_obj, "modules"):
        dtypes = tuple(dropout_types)
        for module in model_obj.modules():
            if not isinstance(module, dtypes):
                continue
            old_p = getattr(module, "p", None)
            state["dropout_states"].append(
                (module, bool(getattr(module, "training", False)), old_p)
            )
            if old_p is not None:
                try:
                    module.p = 0.0
                except Exception:
                    pass
            module.train(False)

    return state


def _restore_model_after_autograd_hessian(model_obj, state):
    """Undo the mutations performed by :func:`_prepare_model_for_autograd_hessian`.

    Safe to call from a ``finally`` block: operates on whatever partial
    state was recorded and gracefully handles the empty-state case.
    """
    for module, was_training, old_p in state.get("dropout_states", []):
        if old_p is not None:
            try:
                module.p = old_p
            except Exception:
                pass
        module.train(was_training)

    if hasattr(model_obj, "train"):
        model_obj.train(state.get("was_training", False))

    for param, req_grad in state.get("param_flags", []):
        param.requires_grad_(req_grad)


class MLIPCalculator(Calculator):
    """PySisyphus-compatible MLIP calculator base class.

    Subclasses must implement:

    * ``_compute_energy_forces_ev(elem, coord_ang)``
      → ``(energy_eV, forces_eV_Ang)``
    * optionally ``_compute_analytical_hessian_ev(elem, coord_ang)``
      → ``hessian_eV_Ang2`` (3N×3N ndarray or torch Tensor)

    The base class provides ``get_energy``, ``get_forces``, ``get_hessian``
    methods that handle freeze-atom masking, unit conversion, and FD Hessian
    assembly automatically.
    """

    implemented_properties = ["energy", "forces", "hessian"]

    def __init__(
        self,
        *,
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
        super().__init__(charge=charge, mult=spin, **kwargs)

        self.device_str = device
        self.hessian_calc_mode = hessian_calc_mode
        self.return_partial_hessian = bool(return_partial_hessian)
        self.hessian_double = bool(hessian_double)
        self.print_timing = bool(print_timing)

        # Normalise freeze_atoms
        if freeze_atoms is None:
            freeze_iter: List[int] = []
        else:
            try:
                freeze_iter = [int(i) for i in list(freeze_atoms)]
            except Exception:
                freeze_iter = []
        self.freeze_atoms = sorted(set(freeze_iter))

    # ------------------------------------------------------------------
    # Subclass hooks (override in backend implementations)
    # ------------------------------------------------------------------

    def _compute_energy_forces_ev(
        self, elem: Sequence[str], coord_ang: np.ndarray
    ) -> tuple:
        """Return ``(energy_eV: float, forces_eV_Ang: ndarray(N,3))``."""
        raise NotImplementedError

    def _compute_analytical_hessian_ev(
        self, elem: Sequence[str], coord_ang: np.ndarray
    ) -> Any:
        """Return Hessian in eV/Å² as (3N,3N) ndarray or torch Tensor.

        Raise ``NotImplementedError`` if analytical Hessian is not available
        for this backend.  The base class will fall back to FD.
        """
        raise NotImplementedError(
            f"Analytical Hessian is not available for {self.__class__.__name__}."
        )

    def _supports_analytical_hessian(self) -> bool:
        """Return True if this backend can compute analytical Hessians."""
        return False

    # ------------------------------------------------------------------
    # Freeze-atom helpers
    # ------------------------------------------------------------------

    def _active_and_frozen_dof_idx(self, n_atoms: int):
        frozen_set = set(self.freeze_atoms)
        active_atoms = [i for i in range(n_atoms) if i not in frozen_set]
        active_dof_idx = [3 * i + j for i in active_atoms for j in range(3)]
        frozen_dof_idx = [3 * i + j for i in self.freeze_atoms for j in range(3)]
        return active_atoms, active_dof_idx, frozen_dof_idx

    def _build_partial_hessian_meta(self, n_atoms: int) -> Optional[dict]:
        """Return the partial-Hessian metadata dict (to be stored as
        ``out['within_partial_hessian']``) when partial Hessian is active,
        otherwise ``None``."""
        if not self.return_partial_hessian or len(self.freeze_atoms) == 0:
            return None
        active_atoms, active_dof_idx, _ = self._active_and_frozen_dof_idx(n_atoms)
        return {
            "active_atoms": np.array(active_atoms, dtype=int),
            "active_dofs": np.array(active_dof_idx, dtype=int),
            "active_n_dof": len(active_dof_idx),
            "full_n_dof": n_atoms * 3,
        }

    def _zero_frozen_forces_ev(self, F: np.ndarray) -> np.ndarray:
        """Zero forces (eV/Å) on frozen atoms."""
        if (F is None) or (len(self.freeze_atoms) == 0):
            return F
        Fz = F.copy()
        Fz[np.asarray(self.freeze_atoms, dtype=int)] = 0.0
        return Fz

    # ------------------------------------------------------------------
    # FD Hessian (CPU, used by non-UMA backends)
    # ------------------------------------------------------------------

    def _build_fd_hessian_cpu(
        self,
        elem: Sequence[str],
        coord_ang: np.ndarray,
        *,
        eps_ang: float = 1.0e-3,
    ) -> Dict[str, Any]:
        """Central-difference Hessian assembled on CPU.

        Frozen-atom DOFs are never displaced, so their rows and columns remain
        zero in the full (3N×3N) output. When ``return_partial_hessian=True``
        the active-DOF sub-block is extracted and returned instead.
        """
        n_atoms = len(elem)
        dof = n_atoms * 3

        frozen_set = set(self.freeze_atoms)
        active_atoms = [i for i in range(n_atoms) if i not in frozen_set]
        active_dof_idx = [3 * i + j for i in active_atoms for j in range(3)]

        # Base point
        e0, F0 = self._compute_energy_forces_ev(elem, coord_ang)
        F0 = np.asarray(F0, dtype=np.float64).reshape(n_atoms, 3)

        H = np.zeros((dof, dof), dtype=np.float64)
        coord_plus = coord_ang.copy()
        coord_minus = coord_ang.copy()

        for k in active_dof_idx:
            a = k // 3
            c = k % 3

            coord_plus[a, c] = coord_ang[a, c] + eps_ang
            _, Fp = self._compute_energy_forces_ev(elem, coord_plus)
            Fp = np.asarray(Fp, dtype=np.float64).reshape(-1)

            coord_minus[a, c] = coord_ang[a, c] - eps_ang
            _, Fm = self._compute_energy_forces_ev(elem, coord_minus)
            Fm = np.asarray(Fm, dtype=np.float64).reshape(-1)

            H[:, k] = -(Fp - Fm) / (2.0 * eps_ang)

            coord_plus[a, c] = coord_ang[a, c]
            coord_minus[a, c] = coord_ang[a, c]

        # Symmetrise
        H = 0.5 * (H + H.T)

        # Reduce to active-DOF block if requested
        if self.return_partial_hessian:
            idx = np.array(active_dof_idx, dtype=int)
            H = H[np.ix_(idx, idx)]

        return {"energy": e0, "forces": F0, "hessian": H}

    # ------------------------------------------------------------------
    # Hessian postprocessing
    # ------------------------------------------------------------------

    def _hessian_ev_to_au(self, H_ev_ang2: np.ndarray) -> np.ndarray:
        """Convert (3N,3N) Hessian from eV/Å² to Hartree/Bohr² and symmetrise."""
        H = np.asarray(H_ev_ang2, dtype=np.float64)
        if H.ndim == 4:
            n = H.shape[0]
            H = H.reshape(n * 3, n * 3)
        H = 0.5 * (H + H.T)
        H *= H_EVAA_2_AU
        return H

    def _apply_active_trim_np(self, H: np.ndarray, n_atoms: int) -> np.ndarray:
        """Extract the active-DOF sub-block when ``return_partial_hessian=True``,
        or zero frozen rows/columns when ``False`` (numpy version)."""
        if len(self.freeze_atoms) == 0:
            return H
        _, active_dof_idx, frozen_dof_idx = self._active_and_frozen_dof_idx(n_atoms)
        dof = n_atoms * 3
        H = H.reshape(dof, dof)
        if self.return_partial_hessian:
            idx = np.array(active_dof_idx, dtype=int)
            return H[np.ix_(idx, idx)]
        else:
            if frozen_dof_idx:
                H[:, frozen_dof_idx] = 0.0
                H[frozen_dof_idx, :] = 0.0
            return H

    # ------------------------------------------------------------------
    # PySisyphus API
    # ------------------------------------------------------------------

    def get_energy(self, elem, coords):
        coord_ang = np.asarray(coords, dtype=np.float64).reshape(-1, 3) * BOHR2ANG
        e_ev, _ = self._compute_energy_forces_ev(elem, coord_ang)
        return {"energy": e_ev * EV2AU}

    def get_forces(self, elem, coords):
        coord_ang = np.asarray(coords, dtype=np.float64).reshape(-1, 3) * BOHR2ANG
        e_ev, F_ev = self._compute_energy_forces_ev(elem, coord_ang)
        F_ev = self._zero_frozen_forces_ev(np.asarray(F_ev, dtype=np.float64).reshape(-1, 3))
        return {
            "energy": e_ev * EV2AU,
            "forces": (F_ev.reshape(-1) * F_EVAA_2_AU),
        }

    def get_hessian(self, elem, coords):
        coord_ang = np.asarray(coords, dtype=np.float64).reshape(-1, 3) * BOHR2ANG
        n_atoms = coord_ang.shape[0]

        hess_total_start = time.perf_counter()
        mode_label = "FiniteDifference"

        mode = (self.hessian_calc_mode or "FiniteDifference").strip().lower()
        use_analytical = mode in ("analytical", "analytic")

        if use_analytical and self._supports_analytical_hessian():
            mode_label = "Analytical"
            t0 = time.perf_counter()
            e_ev, F_ev = self._compute_energy_forces_ev(elem, coord_ang)
            H_ev = self._compute_analytical_hessian_ev(elem, coord_ang)
            H_ev = np.asarray(H_ev, dtype=np.float64)
            if H_ev.ndim == 4:
                H_ev = H_ev.reshape(n_atoms * 3, n_atoms * 3)
            mode_elapsed = time.perf_counter() - t0

            F_ev = self._zero_frozen_forces_ev(
                np.asarray(F_ev, dtype=np.float64).reshape(-1, 3)
            )
            H_ev = self._apply_active_trim_np(H_ev, n_atoms)
            H_au = self._hessian_ev_to_au(H_ev)

            out = {
                "energy": e_ev * EV2AU,
                "forces": (F_ev.reshape(-1) * F_EVAA_2_AU),
                "hessian": H_au,
            }
        else:
            t0 = time.perf_counter()
            res = self._build_fd_hessian_cpu(elem, coord_ang)
            mode_elapsed = time.perf_counter() - t0

            F_ev = self._zero_frozen_forces_ev(
                np.asarray(res["forces"], dtype=np.float64).reshape(-1, 3)
            )
            H_au = self._hessian_ev_to_au(res["hessian"])

            out = {
                "energy": res["energy"] * EV2AU,
                "forces": (F_ev.reshape(-1) * F_EVAA_2_AU),
                "hessian": H_au,
            }

        total_elapsed = time.perf_counter() - hess_total_start
        if self.print_timing:
            click.echo(f"[HessianTiming] mode: {mode_label} | elapsed: {mode_elapsed:.2f} s")

        partial_meta = self._build_partial_hessian_meta(n_atoms)
        if partial_meta is not None:
            out["within_partial_hessian"] = partial_meta

        return out
