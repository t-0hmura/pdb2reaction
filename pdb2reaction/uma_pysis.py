# pdb2reaction/uma_pysis.py

"""
uma_pysis — UMA calculator wrapper for PySisyphus
====================================================================

Usage (API)
-----
    from pdb2reaction.uma_pysis import uma_pysis

Examples::
    >>> from pdb2reaction.uma_pysis import uma_pysis
    >>> calc = uma_pysis(charge=0, spin=1, model="uma-s-1p1")
    >>> calc.get_energy(["C", "O"], coords)
    {'energy': -228.123456}

Description
-----
- Provides energy, forces, and Hessian for molecular systems using FAIR‑Chem UMA pretrained ML potentials via ASE/AtomicData.
- Two Hessian modes:
  - "Analytical": second‑order autograd of the UMA energy on the GPU.
    Requires more VRAM/time than finite differences. During evaluation,
    model parameter gradients are disabled; dropout layers are force‑disabled;
    the model is temporarily toggled to train mode only to build the autograd
    graph and then restored to eval.
      * if `return_partial_hessian=True`, the Hessian is reduced to the Active‑DOF
        block (non‑frozen atoms);
      * otherwise the full matrix is returned and columns for frozen DOF are zeroed.
  - "FiniteDifference": central differences of forces assembled on the GPU.
    Columns for frozen DOF are skipped; optionally return only the active‑DOF block.
    Default step: ε = 1.0e‑3 Å.
- Device handling: device="auto" selects CUDA if available, else CPU; tensors use torch.float32 by default.
- **Precision option**: set `double_precision=True` to run energy/forces/Hessian in float64 (double) end‑to‑end.
  By default (`False`), float32 (single) is used. This toggles both model/graph tensors and working arrays.
- Neighborhood/graph: optional overrides for `max_neigh`, `radius`, `r_edges`.
  On‑the‑fly graphs are built (`otf_graph=True`), and `task_name`, `charge`, `spin`
  are attached to the batch.
- Robustness: analytical path catches CUDA out‑of‑memory and advises switching to
  finite differences.
- Default Hessian mode at construction is "Analytical". If a falsy mode reaches
  `get_hessian`, a fallback "FiniteDifference" is used (implementation detail).
- Units: UMA runs in Å and eV; PySisyphus public API converts to Hartree/Bohr:
  energy eV→Hartree, forces eV·Å⁻¹→Hartree·Bohr⁻¹, Hessian eV·Å⁻²→Hartree·Bohr⁻².

Outputs
-----
PySisyphus calculator interface (`implemented_properties = ["energy", "forces", "hessian"]`):

- get_energy(elem, coords)
  Returns: {"energy": E}
  - E: float, Hartree.
  - coords: Bohr in, internally converted to Å.

- get_forces(elem, coords)
  Returns: {"energy": E, "forces": F}
  - F: shape (3N,), Hartree/Bohr. Components for `freeze_atoms` are zeroed.
  - E: float, Hartree.
  - coords: Bohr in, internally converted to Å.

- get_hessian(elem, coords)
  Returns: {"energy": E, "forces": F, "hessian": H}
  - F: shape (3N,), Hartree/Bohr. Components for `freeze_atoms` are zeroed.
  - H: shape (3N, 3N), Hartree/Bohr² (or (3N_active, 3N_active) if returning the
    active‑DOF submatrix in either mode).
  - If `return_partial_hessian=True`: H contains only the Active‑DOF block
    (non‑frozen atoms). Otherwise H is full sized and columns corresponding to
    frozen DOF are zeroed.

Notes:
-----
- freeze_atoms: list of 0‑based atom indices; **applies to both Analytical and
  FiniteDifference**. Forces on frozen atoms are returned as 0. In Hessians,
  either the matrix is reduced to the Active‑DOF block (`return_partial_hessian=True`)
  or (for full size) columns for frozen DOF are zeroed.
- return_partial_hessian: if True, return only the Active‑DOF submatrix in both modes.
- UMA loader: pretrained_mlip.get_predict_unit(model). The predictor is moved to
  the selected device, set to eval, and all nn.Dropout layers are disabled (p=0).
- During analytical Hessian evaluation, model parameters have requires_grad=False;
  the model is briefly set to train() to enable autograd and then restored to eval().
  CUDA caches are cleared if needed.
- Neighborhood defaults come from the model backbone (e.g., max_neighbors, cutoff)
  unless explicitly overridden.
- CLI entry point: run_pysis() registers the calculator, enabling: `uma_pysis input.yaml`.
"""

from __future__ import annotations
from typing import Sequence, Optional, Dict, Any, List

import numpy as np
import torch
import torch.nn as nn
from ase import Atoms

from fairchem.core import pretrained_mlip
from fairchem.core.datasets.atomic_data import AtomicData
from fairchem.core.datasets import data_list_collater

from pysisyphus.calculators.Calculator import Calculator
from pysisyphus.constants import BOHR2ANG, ANG2BOHR, AU2EV
from pysisyphus import run

# ------------ unit conversion constants ----------------------------
EV2AU          = 1.0 / AU2EV                     # eV → Hartree
F_EVAA_2_AU    = EV2AU / ANG2BOHR                # eV Å⁻¹ → Hartree Bohr⁻¹
H_EVAA_2_AU    = EV2AU / ANG2BOHR / ANG2BOHR     # eV Å⁻² → Hartree Bohr⁻²

# Unified host/GPU dtypes for this implementation
# NOTE: Will be overridden to float64 when `double_precision=True` in uma_pysis.__init__
_ndtype = np.float32
_tdtype = torch.float32

# Default for geom_loader
GEOM_KW_DEFAULT: Dict[str, Any] = {
    "coord_type": "cart",       # str, coordinate representation for geom_loader (Cartesian recommended)
    "freeze_atoms": [],         # list[int], 0-based atom indices to freeze
}

# UMA calculator defaults
CALC_KW: Dict[str, Any] = {
    # Charge and multiplicity
    "charge": 0,              # int, total charge
    "spin": 1,                # int, multiplicity (2S+1)

    # Model selection
    "model": "uma-s-1p1",     # str, UMA pretrained model ID
    "task_name": "omol",      # str, dataset/task tag carried into UMA's AtomicData

    # Device & graph construction
    "device": "auto",         # str, "cuda" | "cpu" | "auto"
    "max_neigh": None,        # Optional[int], override model's neighbor cap
    "radius": None,           # Optional[float], cutoff radius (Å)
    "r_edges": False,         # bool, store edge vectors in graph (UMA option)
    "out_hess_torch": True,   # bool, return Hessian as torch.Tensor (RSIRFO expects GPU Hessian when available)

    # Freeze atoms
    "freeze_atoms": None,     # Optional[Sequence[int]], list of freeze atoms. Corresponding forces are zeroed. Non-DOF Hessian columns are set to 0 (or trimmed).

    # Hessian interfaces to UMA
    "hessian_calc_mode": "Analytical",        # str, "Analytical" (default) | "FiniteDifference"
    "return_partial_hessian": True,           # bool, receive only the active-DOF Hessian block

    # Precision
    "double_precision": False,                # bool, if True use float64 for energy/forces/Hessian
}

# ===================================================================
#                         UMA core wrapper
# ===================================================================
class UMAcore:
    """
    Thin wrapper around fairchem-UMA predict_unit.
    """

    def __init__(
        self,
        elem: Sequence[str],
        *,
        charge: int = 0,
        spin: int = 1,
        model: str = "uma-s-1p1",
        task_name: str = "omol",
        device: str = "auto",
        max_neigh: Optional[int] = None,
        radius:    Optional[float] = None,
        r_edges:   bool = False,
    ):
        # Select device ------------------------------------------------
        if device == "auto":
            device = "cuda" if torch.cuda.is_available() else "cpu"
        self.device_str = device
        self.device = torch.device(device)

        self._AtomicData = AtomicData
        self._collater   = data_list_collater

        # Load predictor and place it on target device & dtype --------
        self.predict = pretrained_mlip.get_predict_unit(model, device=self.device_str)
        # <<< precision-aware >>>  (uses module-level _tdtype configured by uma_pysis)
        self.predict.model.to(self.device, dtype=_tdtype)
        self.predict.model.eval()
        # Hard-disable dropout
        for m in self.predict.model.modules():
            if isinstance(m, nn.Dropout):
                m.p = 0.0

        self.elem      = [e.capitalize() for e in elem]
        self.charge    = charge
        self.spin      = spin
        self.task_name = task_name

        self._max_neigh_user = max_neigh
        self._radius_user    = radius
        self._r_edges_user   = r_edges

    # ----------------------------------------------------------------
    def _model_backbone(self):
        """
        Handle both plain and (D)DP-wrapped models.
        """
        mdl = self.predict.model
        mod = getattr(mdl, "module", mdl)
        return mod.backbone

    # ----------------------------------------------------------------
    def _ase_to_batch(self, atoms: Atoms):
        """
        Convert ASE Atoms → UMA AtomicData (batched).
        """
        backbone = self._model_backbone()
        default_max_neigh = getattr(backbone, "max_neighbors", None)
        default_radius    = getattr(backbone, "cutoff", None)

        max_neigh = self._max_neigh_user if self._max_neigh_user is not None else default_max_neigh
        radius    = self._radius_user    if self._radius_user    is not None else default_radius
        r_edges   = self._r_edges_user

        atoms.info.update({"charge": self.charge, "spin": self.spin})
        data = self._AtomicData.from_ase(
            atoms,
            max_neigh=max_neigh,
            radius   =radius,
            r_edges  =r_edges,
        ).to(self.device, dtype=_tdtype)
        data.dataset = self.task_name
        batch = self._collater([data], otf_graph=True).to(self.device, dtype=_tdtype)
        return batch

    # ----------------------------------------------------------------
    def compute(
        self,
        coord_ang: np.ndarray,
        *,
        forces: bool = False,
        hessian: bool = False,
    ) -> Dict[str, Any]:
        """
        Evaluate energy and (optionally) forces/Hessian at given coordinates.

        Parameters
        ----------
        coord_ang : (N, 3) array-like in Å
        forces : bool
            If True, return forces (and enable grads).
        hessian : bool
            If True, keep autograd path so the analytical Hessian can be built.

        Returns
        -------
        dict with keys:
          - "energy"  : float, eV
          - "forces"  : np.ndarray, shape (N, 3), eV/Å or None
          - "hessian" : torch.Tensor, shape (N,3,N,3), eV/Å² or None
        """
        atoms = Atoms(self.elem, positions=coord_ang)
        batch = self._ase_to_batch(atoms)

        need_grad = bool(forces or hessian)
        # <<< precision-aware >>>
        pos = batch.pos.detach().clone().to(self.device, dtype=_tdtype)
        pos.requires_grad_(need_grad)
        batch.pos = pos

        if need_grad:
            res = self.predict.predict(batch)
        else:
            with torch.no_grad():
                res = self.predict.predict(batch)

        energy = float(res["energy"].squeeze().detach().item())
        forces_np = (
            res["forces"].detach().cpu().numpy().astype(_ndtype, copy=False)
            if (forces or hessian)
            else None
        )

        if hessian:
            # Disable parameter grads during Hessian evaluation
            p_flags = [p.requires_grad for p in self.predict.model.parameters()]
            for p in self.predict.model.parameters():
                p.requires_grad_(False)

            self.predict.model.train()
            try:
                def e_fn(flat):
                    batch.pos = flat.view(-1, 3)
                    out = self.predict.predict(batch)["energy"].squeeze()
                    return out

                H = torch.autograd.functional.hessian(e_fn, batch.pos.view(-1), vectorize=False)
                H = H.view(len(atoms), 3, len(atoms), 3).detach()
            finally:
                self.predict.model.eval()
                for p, flag in zip(self.predict.model.parameters(), p_flags):
                    p.requires_grad_(flag)
                if self.device.type == "cuda":
                    torch.cuda.empty_cache()
        else:
            H = None

        return {"energy": energy, "forces": forces_np, "hessian": H}


# ===================================================================
#                    PySisyphus calculator class
# ===================================================================
class uma_pysis(Calculator):
    """
    PySisyphus-compatible UMA calculator.
    """

    implemented_properties = ["energy", "forces", "hessian"]

    def __init__(
        self,
        *,
        charge: int = CALC_KW["charge"],
        spin: int = CALC_KW["spin"],
        model: str = CALC_KW["model"],
        task_name: str = CALC_KW["task_name"],
        device: str = CALC_KW["device"],
        out_hess_torch: bool = CALC_KW["out_hess_torch"],
        max_neigh: Optional[int] = CALC_KW["max_neigh"],
        radius:    Optional[float] = CALC_KW["radius"],
        r_edges:   bool = CALC_KW["r_edges"],
        freeze_atoms: Optional[Sequence[int]] = CALC_KW["freeze_atoms"],
        hessian_calc_mode: str = CALC_KW["hessian_calc_mode"],
        return_partial_hessian: bool = CALC_KW["return_partial_hessian"],
        double_precision: bool = CALC_KW["double_precision"],
        # -------------------------------------------------------------------
        **kwargs,
    ):
        """
        Parameters
        ----------
        hessian_calc_mode : {"Analytical", "FiniteDifference"}, default "Analytical"
            Select how to compute the Hessian in `get_hessian()`.
            - "Analytical": autograd-based analytical Hessian (GPU).
            - "FiniteDifference": central-difference Hessian; the matrix is
              assembled on the GPU.
        freeze_atoms : list[int], optional
            Atom indices (0-based). In both modes, DOFs of these atoms are
            treated as frozen.
        return_partial_hessian : bool, default True
            If True, return only the Active-DOF Hessian (submatrix for non-frozen atoms).
            If False, return a full (3N×3N) matrix where frozen-DOF columns are 0.
        double_precision : bool, default False
            If True, run energy/forces/Hessian in float64 end-to-end (model, graph tensors,
            working arrays). If False (default), use float32 for performance.
        """
        global _ndtype, _tdtype
        if bool(double_precision):
            _ndtype = np.float64
            _tdtype = torch.float64
        else:
            _ndtype = np.float32
            _tdtype = torch.float32

        super().__init__(charge=charge, mult=spin, **kwargs)
        self._core: Optional[UMAcore] = None
        self._core_kw = dict(
            charge=charge,
            spin=spin,
            model=model,
            task_name=task_name,
            device=device,
            max_neigh=max_neigh,
            radius=radius,
            r_edges=r_edges,
        )
        self.out_hess_torch = out_hess_torch
        self.hessian_calc_mode = hessian_calc_mode
        self.freeze_atoms: List[int] = sorted(set(int(i) for i in (freeze_atoms or [])))
        self.return_partial_hessian = bool(return_partial_hessian)
        self.double_precision = bool(double_precision)

    # ---------- helpers ---------------------------------------------
    def _ensure_core(self, elem: Sequence[str]):
        if self._core is None:
            self._core = UMAcore(elem, **self._core_kw)

    @staticmethod
    def _au_energy(E: float) -> float:
        return E * EV2AU

    @staticmethod
    def _au_forces(F: np.ndarray) -> np.ndarray:
        return (F * F_EVAA_2_AU).reshape(-1)

    def _au_hessian(self, H: torch.Tensor):
        """
        Convert Hessian from eV/Å² to Hartree/Bohr² and format the output.

        Parameters
        ----------
        H : torch.Tensor
            Shape (N,3,N,3), on device.

        Returns
        -------
        numpy.ndarray or torch.Tensor
            If `out_hess_torch` is False, returns a NumPy array on CPU.
            Otherwise returns a torch.Tensor on the current device.
        """
        n = H.size(0)
        H = H.view(n * 3, n * 3)
        # Use a same-device/dtype scalar to avoid unintended FP64 promotion.
        scale = torch.tensor(H_EVAA_2_AU, device=H.device, dtype=H.dtype)
        H = H * scale
        return H.detach().cpu().numpy() if not self.out_hess_torch else H.detach()

    # ---- Common utilities for freeze/active DOF ----------------
    def _active_and_frozen_dof_idx(self, n_atoms: int):
        frozen_set = set(self.freeze_atoms)
        active_atoms = [i for i in range(n_atoms) if i not in frozen_set]
        active_dof_idx = [3 * i + j for i in active_atoms for j in range(3)]
        frozen_dof_idx = [3 * i + j for i in self.freeze_atoms for j in range(3)]
        return active_atoms, active_dof_idx, frozen_dof_idx

    def _zero_frozen_forces_ev(self, F: np.ndarray) -> np.ndarray:
        """Zero forces (eV/Å) on frozen atoms."""
        if (F is None) or (len(self.freeze_atoms) == 0):
            return F
        Fz = F.copy()
        Fz[np.asarray(self.freeze_atoms, dtype=int)] = 0.0
        return Fz

    def _apply_analytical_active_trim(self, H: torch.Tensor) -> torch.Tensor:
        """
        Apply Active‑DOF trimming/column-zeroing to an Analytical Hessian.

        If `return_partial_hessian=True`, reduce to Active‑DOF block.
        Else keep full size and zero columns corresponding to frozen DOF.
        """
        n_atoms = H.size(0)
        if len(self.freeze_atoms) == 0:
            return H  # nothing to do

        active_atoms, active_dof_idx, frozen_dof_idx = self._active_and_frozen_dof_idx(n_atoms)
        H2d = H.view(n_atoms * 3, n_atoms * 3)

        if self.return_partial_hessian:
            idx = torch.tensor(active_dof_idx, device=H.device, dtype=torch.long)
            H_sub = H2d.index_select(0, idx).index_select(1, idx)
            n_act = len(active_atoms)
            return H_sub.view(n_act, 3, n_act, 3)
        else:
            if frozen_dof_idx:
                cols = torch.tensor(frozen_dof_idx, device=H.device, dtype=torch.long)
                H2d.index_fill_(1, cols, 0.0)  # zero columns of frozen DOF
            return H2d.view(n_atoms, 3, n_atoms, 3)

    # ---------- Finite-Difference Hessian (GPU assembly) -------------
    def _build_fd_hessian_gpu(
        self,
        elem: Sequence[str],
        coord_ang: np.ndarray,
        *,
        eps_ang: float = 1.0e-3,  # central-difference step (Å)
    ) -> Dict[str, Any]:
        """
        Assemble the Hessian by central differences of forces, with the matrix
        stored and operated on the GPU:

            H_[:, k] = -(F(x + h e_k) - F(x - h e_k)) / (2 h)

        where F is the flattened force vector (eV/Å). Columns corresponding to
        `freeze_atoms` DOF are skipped. If `return_partial_hessian` is True, the
        returned matrix contains only the Active-DOF block (non-frozen atoms);
        otherwise the full (3N×3N) matrix is returned.

        Returns
        -------
        dict with keys:
          - "energy"  : float, eV
          - "forces"  : np.ndarray, shape (N, 3), eV/Å
          - "hessian" : torch.Tensor, shape (N_out,3,N_out,3), eV/Å² on device
        """
        dev = torch.device("cuda" if torch.cuda.is_available() else "cpu") if self._core is None else self._core.device
        dtype = _tdtype

        n_atoms = len(elem)
        dof = n_atoms * 3

        # Active atoms / DOF
        frozen_set = set(self.freeze_atoms)
        active_atoms = [i for i in range(n_atoms) if i not in frozen_set]
        active_dof_idx = [3 * i + j for i in active_atoms for j in range(3)]

        # Base point evaluation
        res0 = self._core.compute(coord_ang, forces=True, hessian=False)
        energy0_eV = res0["energy"]
        F0 = res0["forces"].astype(_ndtype, copy=False)  # (N,3) eV/Å

        # GPU-side Hessian storage (2D for easy column insertion)
        H = torch.zeros((dof, dof), device=dev, dtype=dtype)

        # Host-side work arrays for coordinate perturbations
        coord_plus = coord_ang.copy()
        coord_minus = coord_ang.copy()

        # Compute columns only for active DOF
        for k in active_dof_idx:
            a = k // 3
            c = k % 3

            # x + h
            coord_plus[a, c] = coord_ang[a, c] + eps_ang
            res_p = self._core.compute(coord_plus, forces=True, hessian=False)
            Fp = res_p["forces"].reshape(-1).astype(_ndtype, copy=False)

            # x - h
            coord_minus[a, c] = coord_ang[a, c] - eps_ang
            res_m = self._core.compute(coord_minus, forces=True, hessian=False)
            Fm = res_m["forces"].reshape(-1).astype(_ndtype, copy=False)

            # Column on GPU
            Fp_t = torch.from_numpy(Fp).to(dev, dtype=dtype)
            Fm_t = torch.from_numpy(Fm).to(dev, dtype=dtype)
            col = -(Fp_t - Fm_t) / (2.0 * eps_ang)  # (3N,) eV/Å²

            H[:, k] = col

            # Restore
            coord_plus[a, c]  = coord_ang[a, c]
            coord_minus[a, c] = coord_ang[a, c]

        # Reduce to Active-DOF block if requested
        if self.return_partial_hessian:
            idx = torch.tensor(active_dof_idx, device=dev, dtype=torch.long)
            H = H.index_select(0, idx).index_select(1, idx)
            n_active_atoms = len(active_atoms)
            H = H.view(n_active_atoms, 3, n_active_atoms, 3)
        else:
            H = H.view(n_atoms, 3, n_atoms, 3)

        return {"energy": energy0_eV, "forces": F0, "hessian": H}

    # ---------- PySisyphus API --------------------------------------
    def get_energy(self, elem, coords):
        self._ensure_core(elem)
        coord_ang = np.asarray(coords, dtype=_ndtype).reshape(-1, 3) * BOHR2ANG
        res = self._core.compute(coord_ang, forces=False, hessian=False)
        return {"energy": self._au_energy(res["energy"])}

    def get_forces(self, elem, coords):
        self._ensure_core(elem)
        coord_ang = np.asarray(coords, dtype=_ndtype).reshape(-1, 3) * BOHR2ANG
        res = self._core.compute(coord_ang, forces=True, hessian=False)

        # Zero forces on frozen atoms (in eV/Å before conversion)
        F_ev = self._zero_frozen_forces_ev(res["forces"])

        return {
            "energy": self._au_energy(res["energy"]),
            "forces": self._au_forces(F_ev),
        }

    def get_hessian(self, elem, coords):
        """
        Compute the Hessian according to `hessian_calc_mode`.

        Modes
        -----
        - "Analytical": autograd-based analytical Hessian (GPU).
          If a CUDA out-of-memory error occurs, a clear message is raised
          instructing you to switch to finite differences.
          Active‑DOF trimming/column‑zeroing is applied to match FD semantics.
        - "FiniteDifference": central-difference Hessian assembled on the GPU.
          * Columns for `freeze_atoms` DOF are skipped.
          * If `return_partial_hessian` is True, return the Active-DOF submatrix;
            otherwise return a full-size matrix (frozen columns remain zero).
        """
        self._ensure_core(elem)
        coord_ang = np.asarray(coords, dtype=_ndtype).reshape(-1, 3) * BOHR2ANG

        mode = (self.hessian_calc_mode or "FiniteDifference").strip().lower()
        if mode in ("analytical", "analytic"):
            try:
                res = self._core.compute(coord_ang, forces=True, hessian=True)
            except (torch.cuda.OutOfMemoryError, RuntimeError) as e:
                # Detect CUDA OOM patterns
                msg = str(e).lower()
                if "out of memory" in msg and "cuda" in msg:
                    raise RuntimeError(
                        "Analytical Hessian computation failed due to CUDA out-of-memory. "
                        "Your GPU memory appears to be limited. Please switch to the finite-"
                        "difference Hessian by specifying `--hessian-calc-mode FiniteDifference` "
                        "in the external CLI, or `hessian_calc_mode=\"FiniteDifference\"` in this calculator."
                    ) from e
                raise

            # Zero forces for frozen atoms (eV/Å)
            res_forces_ev = self._zero_frozen_forces_ev(res["forces"])

            # Apply Active‑DOF trimming/column‑zeroing for Analytical
            H = self._apply_analytical_active_trim(res["hessian"])

            return {
                "energy": self._au_energy(res["energy"]),
                "forces": self._au_forces(res_forces_ev),
                "hessian": self._au_hessian(H),
            }

        elif mode in ("finitedifference", "finite-difference", "fd"):
            res = self._build_fd_hessian_gpu(elem, coord_ang)

            # Zero forces for frozen atoms (eV/Å)
            res_forces_ev = self._zero_frozen_forces_ev(res["forces"])

            return {
                "energy": self._au_energy(res["energy"]),
                "forces": self._au_forces(res_forces_ev),
                "hessian": self._au_hessian(res["hessian"]),
            }
        else:
            # Fallback: treat unknown mode as FiniteDifference
            res = self._build_fd_hessian_gpu(elem, coord_ang)

            # Zero forces for frozen atoms (eV/Å)
            res_forces_ev = self._zero_frozen_forces_ev(res["forces"])

            return {
                "energy": self._au_energy(res["energy"]),
                "forces": self._au_forces(res_forces_ev),
                "hessian": self._au_hessian(res["hessian"]),
            }


# ---------- CLI ----------------------------------------
def run_pysis():
    """
    Enable `uma_pysis input.yaml`.
    """
    run.CALC_DICT["uma_pysis"] = uma_pysis
    run.run()


if __name__ == "__main__":
    run_pysis()
