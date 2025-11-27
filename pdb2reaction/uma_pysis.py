"""
uma_pysis — UMA calculator wrapper for PySisyphus
====================================================================

Usage (API)
-----------
    from pdb2reaction.uma_pysis import uma_pysis

Examples
--------
    >>> from pdb2reaction.uma_pysis import uma_pysis
    >>> calc = uma_pysis(charge=0, spin=1, model="uma-s-1p1")
    >>> calc.get_energy(["C", "O"], coords)
    {'energy': -228.123456}

Description
-----------
- Provides energy, forces, and Hessian for molecular systems using FAIR‑Chem UMA
  pretrained ML potentials via ASE/AtomicData.
- Two Hessian modes (on the selected device; GPU if available):
  - **Analytical**: second‑order autograd of the UMA energy.
    Requires more VRAM/time than finite differences. During evaluation,
    model parameter gradients are disabled; dropout layers are force‑disabled;
    the model is temporarily toggled to `train()` only to build the autograd
    graph and then restored to `eval()`.
      * if `return_partial_hessian=True`, the Hessian is reduced to the Active‑DOF
        block (non‑frozen atoms);
      * otherwise the full matrix is returned and columns for frozen DOF are zeroed.
  - **FiniteDifference**: central differences of forces, assembled on the selected
    device. Columns for frozen DOF are skipped; optionally return only the
    active‑DOF block. Default step: ε = 1.0e‑3 Å.
- Device handling: `device="auto"` selects CUDA if available, else CPU.
- **Precision / dtype**:
    * UMA models run in their native precision (currently float32). We no longer
      upcast the full model and graph to float64.
    * Energies and forces returned by the PySisyphus API are always float64
      (Hartree / Hartree·Bohr⁻¹).
    * The Hessian dtype is controlled by `hessian_double`:
        - `hessian_double=True` (default): assemble and return the Hessian in
          float64 (double precision) on the host/device, without changing the
          internal model precision.
        - `hessian_double=False`: assemble and return the Hessian in the model
          dtype (typically float32).
- Neighborhood/graph: optional overrides for `max_neigh`, `radius`, `r_edges`.
  On‑the‑fly graphs are built (`otf_graph=True`), and `task_name` is attached to
  the batch; charge/spin are stored in `Atoms.info`.
- Robustness: analytical path catches CUDA out‑of‑memory and advises switching to
  finite differences.
- Default Hessian mode at construction is `"Analytical"`. If `hessian_calc_mode`
  is falsy *or unrecognized* in `get_hessian`, `"FiniteDifference"` is used.
- Units: UMA runs in Å and eV; PySisyphus public API converts to Hartree/Bohr:
  energy eV→Hartree, forces eV·Å⁻¹→Hartree·Bohr⁻¹, Hessian eV·Å⁻²→Hartree·Bohr⁻².

Outputs (& Directory Layout)
----------------------------
In-memory calculator interface (`implemented_properties = ["energy", "forces", "hessian"]`)
  ├─ get_energy(elem, coords)  → ``{"energy": E}`` (E in Hartree; coords supplied in Bohr and converted to Å internally)
  ├─ get_forces(elem, coords)  → ``{"energy": E, "forces": F}`` (F has 3N components in Hartree/Bohr; frozen DOF set to 0)
  └─ get_hessian(elem, coords) → ``{"energy": E, "forces": F, "hessian": H}``
        • ``H`` has shape (3N, 3N) in Hartree/Bohr² or (3N_active, 3N_active) when ``return_partial_hessian=True``
        • Columns for frozen DOF are zeroed in the full-size form


Notes
-----
- `freeze_atoms`: list of 0‑based atom indices; **applies to both Analytical and
  FiniteDifference**. Forces on frozen atoms are returned as 0. In Hessians,
  either the matrix is reduced to the Active‑DOF block (`return_partial_hessian=True`)
  or (for full size) columns for frozen DOF are zeroed.
- `return_partial_hessian`: if True, return only the Active‑DOF submatrix in both modes.
- UMA loader: `pretrained_mlip.get_predict_unit(model)`. The predictor is moved to
  the selected device and set to eval; all `nn.Dropout` layers are disabled (`p=0`).
- During analytical Hessian evaluation, model parameters have `requires_grad=False`;
  the model is briefly set to `train()` to enable autograd and then restored to `eval()`.
  CUDA caches are cleared if needed.
- Neighborhood defaults come from the model backbone (e.g., `max_neighbors`, `cutoff`)
  unless explicitly overridden.
- CLI entry point: `run_pysis()` registers the calculator, enabling:
  `uma_pysis input.yaml`.
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

# Default for potential geometry loader integrations (kept for completeness)
GEOM_KW_DEFAULT: Dict[str, Any] = {
    "coord_type": "cart",       # coordinate representation (Cartesian recommended)
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
    "out_hess_torch": True,   # bool, return Hessian as torch.Tensor

    # Freeze atoms
    "freeze_atoms": None,     # Optional[Sequence[int]], list of freeze atoms

    # Hessian interfaces to UMA
    "hessian_calc_mode": "Analytical",        # "Analytical" (default) | "FiniteDifference"
    "return_partial_hessian": True,           # receive only the active-DOF Hessian block

    # Hessian precision (energy/forces are always returned as float64)
    "hessian_double": True,                   # if True, assemble/return Hessian in float64
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

        # Load predictor and place it on target device ----------------
        self.predict = pretrained_mlip.get_predict_unit(model, device=self.device_str)
        self.predict.model.eval()

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
        ).to(self.device)
        data.dataset = self.task_name
        batch = self._collater([data], otf_graph=True).to(self.device)
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
        pos = batch.pos.detach().clone()
        pos.requires_grad_(need_grad)
        batch.pos = pos

        if need_grad:
            res = self.predict.predict(batch)
        else:
            with torch.no_grad():
                res = self.predict.predict(batch)

        energy = float(res["energy"].squeeze().detach().item())
        forces_np = (
            res["forces"].detach().cpu().numpy()
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
                    return self.predict.predict(batch)["energy"].squeeze()

                H = torch.autograd.functional.hessian(
                    e_fn, batch.pos.view(-1), vectorize=False
                )
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
        hessian_double: bool = CALC_KW["hessian_double"],
        # -------------------------------------------------------------------
        **kwargs,
    ):
        """
        Parameters
        ----------
        hessian_calc_mode : {"Analytical", "FiniteDifference"}, default "Analytical"
            Select how to compute the Hessian in `get_hessian()`.
            - "Analytical": autograd-based analytical Hessian.
            - "FiniteDifference": central-difference Hessian; the matrix is
              assembled on the selected device.
        freeze_atoms : list[int], optional
            Atom indices (0-based). In both modes, DOFs of these atoms are
            treated as frozen.
        return_partial_hessian : bool, default True
            If True, return only the Active-DOF Hessian (submatrix for non-frozen atoms).
            If False, return a full (3N×3N) matrix where frozen-DOF columns are 0.
        hessian_double : bool, default True
            If True, assemble/return the Hessian in float64 (double precision).
            Energy and forces are always returned as float64 regardless of this flag.
        """
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
        self.hessian_double = bool(hessian_double)

    # ---------- helpers ---------------------------------------------
    def _ensure_core(self, elem: Sequence[str]):
        if self._core is None:
            self._core = UMAcore(elem, **self._core_kw)

    @staticmethod
    def _au_energy(E: float) -> float:
        return E * EV2AU

    @staticmethod
    def _au_forces(F: np.ndarray) -> np.ndarray:
        F64 = np.asarray(F, dtype=np.float64)
        return (F64 * F_EVAA_2_AU).reshape(-1)

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

        Notes
        -----
        The dtype of the returned Hessian is controlled solely by
        `self.hessian_double`:

        - If True, the result is float64 (double precision).
        - If False, the result stays in the model/native dtype.
        """
        n = H.size(0)
        H = H.view(n * 3, n * 3)

        # Unit conversion
        H = H * H_EVAA_2_AU

        if self.hessian_double:
            H = H.to(dtype=torch.float64)

        if self.out_hess_torch:
            return H.detach()
        else:
            return H.detach().cpu().numpy()

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
            return H

        _, active_dof_idx, frozen_dof_idx = self._active_and_frozen_dof_idx(n_atoms)
        H = H.view(n_atoms * 3, n_atoms * 3)

        if self.return_partial_hessian:
            idx = torch.tensor(active_dof_idx, device=H.device, dtype=torch.long)
            H_sub = H.index_select(0, idx).index_select(1, idx)
            n_act = len(active_dof_idx) // 3
            return H_sub.view(n_act, 3, n_act, 3)
        else:
            if frozen_dof_idx:
                cols = torch.tensor(frozen_dof_idx, device=H.device, dtype=torch.long)
                H.index_fill_(1, cols, 0.0)  # zero columns of frozen DOF
            return H.view(n_atoms, 3, n_atoms, 3)

    # ---------- Finite-Difference Hessian (GPU/CPU assembly) -------------
    def _build_fd_hessian_gpu(
        self,
        elem: Sequence[str],
        coord_ang: np.ndarray,
        *,
        eps_ang: float = 1.0e-3,  # central-difference step (Å)
    ) -> Dict[str, Any]:
        """
        Assemble the Hessian by central differences of forces, with the matrix
        stored and operated on the selected device:

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
        self._ensure_core(elem)
        core = self._core
        assert core is not None

        dev = core.device

        n_atoms = len(elem)
        dof = n_atoms * 3

        # Active atoms / DOF
        frozen_set = set(self.freeze_atoms)
        active_atoms = [i for i in range(n_atoms) if i not in frozen_set]
        active_dof_idx = [3 * i + j for i in active_atoms for j in range(3)]

        # Base point evaluation
        res0 = core.compute(
            coord_ang, forces=True, hessian=False
        )
        energy0_eV = res0["energy"]
        F0 = res0["forces"]  # (N,3) eV/Å, native numpy dtype

        # Assemble Hessian in the same dtype as the model forces.
        force_dtype = torch.from_numpy(F0).dtype

        # Device-side Hessian storage (2D for easy column insertion)
        H = torch.zeros((dof, dof), device=dev, dtype=force_dtype)

        # Host-side work arrays for coordinate perturbations
        coord_plus = coord_ang.copy()
        coord_minus = coord_ang.copy()

        # Compute columns only for active DOF
        for k in active_dof_idx:
            a = k // 3
            c = k % 3

            # x + h
            coord_plus[a, c] = coord_ang[a, c] + eps_ang
            res_p = core.compute(
                coord_plus, forces=True, hessian=False
            )
            Fp = res_p["forces"].reshape(-1)  # (3N,) eV/Å

            # x - h
            coord_minus[a, c] = coord_ang[a, c] - eps_ang
            res_m = core.compute(
                coord_minus, forces=True, hessian=False
            )
            Fm = res_m["forces"].reshape(-1)  # (3N,) eV/Å

            # Column on device
            Fp_t = torch.from_numpy(Fp).to(dev, dtype=force_dtype)
            Fm_t = torch.from_numpy(Fm).to(dev, dtype=force_dtype)
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
        coord_ang = np.asarray(coords, dtype=np.float64).reshape(-1, 3) * BOHR2ANG
        res = self._core.compute(
            coord_ang, forces=False, hessian=False
        )
        return {"energy": self._au_energy(res["energy"])}

    def get_forces(self, elem, coords):
        self._ensure_core(elem)
        coord_ang = np.asarray(coords, dtype=np.float64).reshape(-1, 3) * BOHR2ANG
        res = self._core.compute(
            coord_ang, forces=True, hessian=False
        )

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
        - "Analytical": autograd-based analytical Hessian.
          If a CUDA out-of-memory error occurs, a clear message is raised
          instructing you to switch to finite differences.
          Active‑DOF trimming/column‑zeroing is applied to match FD semantics.
        - "FiniteDifference": central-difference Hessian assembled on the selected device.
          * Columns for `freeze_atoms` DOF are skipped.
          * If `return_partial_hessian` is True, return the Active-DOF submatrix;
            otherwise return a full-size matrix (frozen columns remain zero).

        If `hessian_calc_mode` is falsy or unrecognized, FiniteDifference is used.
        """
        self._ensure_core(elem)
        coord_ang = np.asarray(coords, dtype=np.float64).reshape(-1, 3) * BOHR2ANG

        mode = (self.hessian_calc_mode or "FiniteDifference").strip().lower()
        if mode in ("analytical", "analytic"):
            try:
                res = self._core.compute(
                    coord_ang,
                    forces=True,
                    hessian=True,
                )
            except (torch.cuda.OutOfMemoryError, RuntimeError) as e:
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

            # Apply Active‑DOF trimming/column-zeroing for Analytical
            H = self._apply_analytical_active_trim(res["hessian"])

            return {
                "energy": self._au_energy(res["energy"]),
                "forces": self._au_forces(res_forces_ev),
                "hessian": self._au_hessian(H),
            }

        # FiniteDifference (default/fallback)
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
